!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA.TestAnalytics.F90
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring 
!                     Warning and Forecasting (GBA-MWF) Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Yuanfu Xie on 2025-02-08 for testing a set of analytic solution numerically:
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This test program tests a set of analytic solution.
!! @copyright (C) 2020 GBA-MWF, All rights reserved.
! @note
! @warning
! @attention
PROGRAM testAnalytics
  USE kinds_m, ONLY: i_kind, r_kind
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE geometry_m, ONLY: geometry_t
  USE gzm_m, ONLY: gzm_t
  USE State_m, ONLY: State_t
  USE State2NC_m
  USE YAMLRead_m

  IMPLICIT NONE

  INCLUDE "mpif.h"

  TYPE(mpddGlob_t) :: mpddGlob
  TYPE(geometry_t) :: geo
  TYPE(gzm_t), ALLOCATABLE :: gzm(:)

  ! For plot:
  TYPE(State_t), ALLOCATABLE :: X(:)

  CHARACTER(LEN=1024) :: yamlFile, ncOutputFile
  CHARACTER(LEN=20) :: task
  INTEGER(i_kind) :: mgStart,mgEnd,i,j,istatus,ic,ie,iq
  INTEGER(i_kind), ALLOCATABLE :: imx(:),viex(:),vimx(:)
  LOGICAL :: pass
  REAL(r_kind) :: a
  REAL(r_kind), ALLOCATABLE :: psi(:,:),chi(:,:),vor(:,:),div(:,:),vor_num(:,:),div_num(:,:), &
    u(:,:,:,:),v(:,:,:,:),u_num(:,:,:,:),v_num(:,:,:,:),amx(:),allmx(:),fun(:),vamx(:),vfun(:)

  ! Get the configuration file
  CALL getarg(1, yamlFile)

  ! Set up MPDD and Geometry:
  CALL mpddGlob%initialize()
  CALL geo%initialize(yamlFile,mpddGlob)
  mgStart = geo%mg%mg_coarsest
  mgEnd = geo%mg%mg_finest

  ! Check the analytic solutions:
  ALLOCATE(gzm(mgStart:mgEnd),X(mgStart:mgEnd),imx(mgStart:mgEnd), &
    amx(mgStart:mgEnd),allmx(mgStart:mgEnd),fun(mgStart:mgEnd))
  ALLOCATE(viex(mgStart:mgEnd),vimx(mgStart:mgEnd),vamx(mgStart:mgEnd),vfun(mgStart:mgEnd))
  
  ! Check all multigrid levels:
  DO i=mgStart,mgEnd
    gzm(i) = gzm_t(geo%mg%sg(i))
    ASSOCIATE(sg => geo%mg%sg(i))
      ALLOCATE(psi(sg%vLevel,sg%num_cell), &
              chi(sg%vLevel,sg%num_cell), &
              vor(sg%vLevel,sg%num_cell), &
              div(sg%vLevel,sg%num_cell), &
              vor_num(sg%vLevel,sg%num_cell), &
              div_num(sg%vLevel,sg%num_cell), &
              u(sg%vLevel,sg%numQuadPerEdge,sg%num_edge,sg%num_cell), &
              v(sg%vLevel,sg%numQuadPerEdge,sg%num_edge,sg%num_cell), &
              u_num(sg%vLevel,sg%numQuadPerEdge,sg%num_edge,sg%num_cell), &
              v_num(sg%vLevel,sg%numQuadPerEdge,sg%num_edge,sg%num_cell))
      u = 0.0D0; v = 0.0D0
      CALL LinearLatCase(sg%num_cell,sg%vLevel,psi,chi,vor,div,u,v,sg%cell_cntr)

      PRINT*,'Got analytics...',i,mgEnd,MAXVAL(v(1,1,1,:)),sg%num_edge, &
        sg%vLevel,sg%gLevel

      CALL gzm(i)%Laplacia(psi,vor_num)

      ! Set the boundary conditions to the analytic function: no check on BC
      DO j=1,sg%num_cell
        IF (sg%cell_type(j) .NE. 0) vor_num(:,j) = vor(:,j)
      END DO

      ! Check velocity:
      CALL psiChi2velocity(psi,chi,u_num,v_num,sg)

      ! Errors:
      imx(i) = 0
      amx(i) = 0.0D0
      DO j=1,sg%num_cell
        IF (ABS(vor_num(1,j)-vor(1,j)) .GT. amx(i)) THEN
          imx(i) = j; amx(i) = ABS(vor(1,j)-vor_num(1,j))
        END IF
      END DO
      fun(i) = MAXVAL(ABS(vor(1,:)))

      CALL MPI_REDUCE(amx(i), allmx(i), 1, mpi_double_precision, mpi_max, &
        mpddGlob%rankBase, mpddGlob%comm, mpddGlob%ierr)

      IF (mpddGlob%isBaseProc()) WRITE(*,1) i,imx(i),amx(i),allmx(i)
  1   FORMAT('Maximum error at G:',I2,' at cell:',I6,' max error: ',E14.6,'  max all: ',E14.6)

      ! Note: velocity check needs to modify the velocity function in analytics as
      !       analytics gives u and v at cell centers but psiChi2velocity gives u and v at edges:
      viex(i) = 0
      vimx(i) = 0
      vamx(i) = 0.0D0
      DO ic=1,sg%num_icell
        IF (sg%cell_type(ic) .EQ. 0) THEN
          DO ie=1,sg%num_edge
            DO iq=1,sg%numQuadPerEdge
              a = DSQRT((u_num(1,iq,ie,ic)-u(1,iq,ie,ic))**2 + &
                        (v_num(1,iq,ie,ic)-v(1,iq,ie,ic))**2)
              IF (a .GT. vamx(i)) THEN
                vimx(i) = ic; viex(i) = ie; vamx(i) = a
              END IF
            END DO
          END DO
        END IF
      END DO
      vfun(i) = MAXVAL(DSQRT((u(1,1,:,:))**2 + &
                             (v(1,1,:,:))**2))

      ! WRITE(*,2) i,sg%mpddInfo_sg%myrank,vimx(i),viex(i),vamx(i),vfun(i)
  2   FORMAT('Velocity error of  at G:',I2,' pc:',I1,' at cell:',I6, &
        ' at edge: ',I2,' max error: ',E14.6,' max fun: ',E14.6)

      ! Plot:
      istatus = yaml_get_var(TRIM(yamlFile), 'IO', 'output_dir', ncOutputFile)
      istatus = yaml_get_var(TRIM(yamlFile), 'RunMode', 'Task', task)
      CALL X(i)%initialize(yamlFile,sg)
      X(i)%fields(1)%DATA(:,:,1) = vor_num-vor
      X(i)%fields(2)%DATA(:,:,1) = vor_num
      X(i)%fields(3)%DATA(:,:,1) = vor
      X(i)%fields(4)%DATA(:,:,1) = v_num(:,1,2,:)
      X(i)%fields(5)%DATA(:,:,1) = v(:,1,2,:)
      CALL Output_NC_State_AV(X(i), TRIM(ncOutputFile), TRIM(task)//"_Analytics", .TRUE., .TRUE.)

      ! Deallocate:
      DEALLOCATE(psi,chi,vor,div,vor_num,div_num,u,v,u_num,v_num)
    END ASSOCIATE
  END DO

  ! Output relative errors:
  pass = .TRUE.

  IF (mpddGlob%isBaseProc()) THEN
    DO i=mgStart+1,mgEnd
      IF (allmx(i-1)/allmx(i) .LT. 3.9D0) pass = .FALSE.
      WRITE(*,3) amx(i-1)/amx(i),i,amx(i)/fun(i)
3     FORMAT('Error reduction ratio: ',E12.4,' at G',I1,' relative error: ',E14.6)
    END DO
  END IF
  ! The following is deactivated until above velocity check is fixed:
!   DO i=mgStart+1,mgEnd
!     IF (vamx(i) .GT. 1.0D-10) THEN
!       IF (vamx(i-1)/vamx(i) .LT. 3.9D0) pass = .FALSE.
!     END IF
!     WRITE(*,4) vamx(i-1)/vamx(i),i,vamx(i)/vfun(i)
! 4   FORMAT('Error reduction of velocity ratio: ',E12.4,' at G',I1,' relative error: ',E14.6)
!   END DO

  ! Check whether test passes:
  IF (pass) THEN
    WRITE(*,5)
5   FORMAT('Test boundary interpolation - Test passed!')
  ELSE
    WRITE(*,6)
6   FORMAT('Test boundary interpolation - Test failed!')
  END IF

  DEALLOCATE(amx,imx,vamx,vimx)

END PROGRAM testAnalytics