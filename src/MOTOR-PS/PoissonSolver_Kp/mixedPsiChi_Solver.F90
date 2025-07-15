!!--------------------------------------------------------------------------------------------------
! PROJECT         : MOTOR-DA.possionSolver_Kp.psMatrix
! AFFILIATION     : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring 
!                   Warning and Forecasting (GBA-MWF) Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Yuanfu Xie, on 2024-12-05
!!--------------------------------------------------------------------------------------------------
!> @brief
!! This module contains a solver of streamfunction and velocity potential from vorticity and divergence
!! with u and v as boundary conditions.
!! @copyright (C) 2024 GBA-MWF, All rights reserved.
!! @note
!! @warning
!! @attention
MODULE mixedPsiChiSolver_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE geometry_m, ONLY: geometry_t
  USE gzm_m, ONLY: gzm_t
  USE State_m, ONLY: State_t
  USE MGOpts_m

  TYPE :: mixedPsiChiSolver_t
    INTEGER(i_kind) :: mgStart,mgEnd
    INTEGER(i_kind) :: numVLevels,numCells
    TYPE(gzm_t), ALLOCATABLE :: gzm(:)
    TYPE(state_t), ALLOCATABLE :: states(:)

    ! Plotting:
    CHARACTER(LEN=1024) :: yamlFile, ncOutputFile
    CHARACTER(LEN=20) :: task
    CONTAINS
      PROCEDURE :: initialize_s
      PROCEDURE :: Solve_s
      PROCEDURE :: rangeSpaceProjection_s
  END TYPE mixedPsiChiSolver_t

  CONTAINS
    SUBROUTINE initialize_s(this,mgStart,mgEnd,geo)
      CLASS(mixedPsiChiSolver_t) :: this
      INTEGER(i_kind), INTENT(IN) :: mgStart,mgEnd
      TYPE(geometry_t), INTENT(IN) :: geo

      ! For plot: to delete when the debugging is done
      !========================================================
      INTEGER(i_kind) :: istatus
      !========================================================

      CALL getarg(1, this%yamlFile)

      ! Plot the solutions: Debugging purpose
      !========================================================
      istatus = yaml_get_var(TRIM(this%yamlFile), 'IO', 'output_dir', this%ncOutputFile)
      istatus = yaml_get_var(TRIM(this%yamlFile), 'RunMode', 'Task', this%task)
      !========================================================

      ! Allocate memory and gzms:
      IF (mgStart .LT. geo%mg%mg_coarsest .OR. mgEnd .GT. geo%mg%mg_finest) THEN
        WRITE(*,1) mgStart,mgEnd,geo%mg%mg_coarsest,geo%mg%mg_finest
1       FORMAT('Multigrid levels: ',2I3,' requested are out of the geometry setting: ',2I3)
        STOP
      END IF

      ALLOCATE(this%gzm(mgStart:mgEnd),this%states(mgStart:mgEnd))
      DO i=mgStart,mgEnd
        this%gzm(i) = gzm_t(geo%mg%sg(i))
        CALL this%states(i)%initialize(this%yamlFile,geo%mg%sg(i))
      END DO
      this%numVLevels = geo%mg%sg(mgEnd)%vLevel
      this%numCells = geo%mg%sg(mgEnd)%num_cell

      this%mgStart = mgStart
      this%mgEnd = mgEnd

    END SUBROUTINE initialize_s

    ! Apply a multigrid method solving a mixed Poisson equations
    ! for Psi and Chi with u and v boundary conditions:
    SUBROUTINE Solve_s(this, vor, div, u, v, psi, chi, geo)
      USE mpddGlob_m, ONLY: mpddGlob_t
      USE geometry_m, ONLY: geometry_t
      USE State_m, ONLY: State_t
      USE State2NC_m
      USE YAMLRead_m
      USE parameters_m, ONLY: EarthRadius

      IMPLICIT NONE
      CLASS(mixedPsiChiSolver_t) :: this
      TYPE(geometry_t), INTENT(IN) :: geo   ! Geometry is used to implement a multigrid solver
      REAL(r_kind), INTENT(IN) :: &
        vor(this%numVLevels,this%numCells), &
        div(this%numVLevels,this%numCells), &
        u(this%numVLevels,this%numCells), &
        v(this%numVLevels,this%numCells)
      
      REAL(r_kind), INTENT(INOUT) :: &
        psi(this%numVLevels,this%numCells), &
        chi(this%numVLevels,this%numCells)

      ! Local variables:
      EXTERNAL :: ddot,dotprd,psolve,AdotX4PsiChi_s
      INTEGER(i_kind) :: glevel,ivlvl,ngmres,mgmres,its,info,i,m
      REAL(r_kind) :: res,del,scaling(2)  ! Scaling for interior and BC
      REAL(r_kind), ALLOCATABLE :: rhs(:),ua(:),va(:)   ! Right hand sides for each multigrid levels
      REAL(r_kind), ALLOCATABLE :: xgmres(:),hgmres(:,:),vgmres(:,:)

      its = 80000
      res = 1.0D-18 ! on input: |Ax-b|/|b|<res;
                  ! on exit:  residual reached
      del = 1.0D1 ! on input: if(del>0) then the x returned is the hookstep
                  ! on exit:  norm of next b predicted by hook
                  ! on exit:  number of its taken
      info = 1    ! on input: if(info==1) print* residuals
                  ! if(info==2) recalc hookstep with new del
                  ! on exit:  0 sucessful, 1 method breakdown, 2 max its

      ! Pass the input states into analysis state:
      this%states(this%mgEnd)%fields(1)%DATA(:,:,1) = psi
      this%states(this%mgEnd)%fields(2)%DATA(:,:,1) = chi
      this%states(this%mgEnd)%fields(3)%DATA(:,:,1) = u
      this%states(this%mgEnd)%fields(4)%DATA(:,:,1) = v
      this%states(this%mgEnd)%fields(5)%DATA(:,:,1) = vor
      this%states(this%mgEnd)%fields(6)%DATA(:,:,1) = div
      ! Interpolate the finest state to coarser grids:
      IF (this%mgEnd .GT. this%mgStart) THEN
        DO glevel=this%mgEnd-1,this%mgStart,-1
          CALL restrictionMG(this%states(glevel), this%states(glevel+1),geo%mg)
        END DO
      END IF

      ! Apply a multigrid solver:
      DO glevel=this%mgEnd,this%mgEnd
        ASSOCIATE(vLevel   => geo%mg%sg(glevel)%vLevel, &
                  num_cell => geo%mg%sg(glevel)%num_cell, &
                  sg       => geo%mg%sg(glevel), &
                  states   => this%states(glevel))

          ! This solver solves the system of equations at each vertical level for better efficiency and accuracy:
          DO ivlvl=1,vLevel

            ngmres = 2*num_cell ! Solve for psi or chi separately
            mgmres = 100      ! GMRES dimension，maximum number of iterations for a restart

            ALLOCATE(rhs(ngmres), xgmres(ngmres), hgmres(mgmres+1,mgmres), vgmres(ngmres,mgmres+1))
            rhs = 0.0D0; xgmres = 0.0D0; hgmres = 0.0D0; vgmres = 0.0D0

            ALLOCATE(ua(vLevel), va(vLevel))

            ! Initial guess:
            xgmres = 1.0D0 !/DSQRT(DBLE(ngmres))
            ! Ideal guess:
            xgmres = 0.0D0

            ! Assign vorticity and divergence to the rhs for interior:
            DO i=1,num_cell
              IF (sg%cell_type(i) .EQ. 0) THEN
                rhs(i         ) = states%fields(5)%DATA(ivlvl,i,1) ! vorticity
                rhs(i+num_cell) = states%fields(6)%DATA(ivlvl,i,1) ! Divergence
              END IF
            END DO

            ! Scale of all vertical levels:
            scaling = DOT_PRODUCT(rhs,rhs)
            scaling = DSQRT(scaling)
            scaling(2) = 1.0D0  ! Scaling parameter for boundary conditions
            WRITE(*,111) scaling,glevel
        111 FORMAT('Scaling params: ',2E14.6,' at level: ',I1)

            ! Scaling the vorticity and divergence:
            rhs = rhs/scaling(1)

            ! Add mixed boundary conditions at the physical boundaries: u and v
            ! DO i=1,sg%num_cell
            !   ! Apply BC to boundary only, find an interior point with a neighbor(s)
            !   IF (sg%cell_type(i) .EQ. 1) THEN
            !     DO j=1,UBOUND(sg%cell_adnb,1)
            !       ! Skip those neighbor cells outside the domain:
            !       IF (sg%cell_adnb(j,i) .NE. 0) THEN
            !         ! Apply the BC to the edge sharing with an interior point:
            !         IF (sg%cell_type(sg%cell_adnb(j,i)) .EQ. 0) THEN
            !           rhs(i         ) = states%fields(3)%DATA(ivlvl,i,1)*scaling(2)
            !           rhs(i+num_cell) = states%fields(4)%DATA(ivlvl,i,1)*scaling(2)
            !         END IF
            !       END IF
            !     END DO
            !   END IF
            ! END DO

            ! Yuanfu Xie modified gmres to gmres_gzm to access gzm for AdotX4PsiChi_s:
            ! Solve for psi:
            psi = 1.0D0/DSQRT(dot_product(psi(ivlvl,:),psi(ivlvl,:)))
            xgmres(1:num_cell) = psi(ivlvl,:) ! Initial guess
            CALL gmresm_gzm(mgmres,ngmres,xgmres,rhs,AdotX4PsiChi_s,psolve,dotprd, &
              hgmres,vgmres,res,del,its,info, this%gzm(glevel),scaling)

            ! Check the output from gmres_gzm:
            WRITE(*,1)
    1       FORMAT('+----------------------------------------+')
            IF (info .EQ. 0) THEN
              WRITE(*,2)
    2         FORMAT('|   GMRES successes:                    |')
            ELSE IF (info .EQ. 1) THEN
              WRITE(*,3) info
    3         FORMAT('| GMRES breakdown, no solution output! ',I1,' |')
            ELSE
              WRITE(*,4)
    4         FORMAT('| GMRES reaches max iterations!          |')
            END IF
            WRITE(*,1)

            ! Reshape the stream function solution:
            states%fields(1)%DATA(ivlvl,:,1) = xgmres(1:num_cell)/scaling(1)
            states%fields(2)%DATA(ivlvl,:,1) = xgmres(num_cell+1:2*num_cell)/scaling(1)
          END DO

          !===================================================================
          ! Check the solutions:
          !===================================================================
          ! Calculate the numerical vorticity and divergence for verification:
          CALL this%gzm(glevel)%Laplacia(states%fields(1)%DATA(:,:,1), &
            this%states(glevel)%fields(5)%DATA(:,:,1))
          CALL this%gzm(glevel)%Laplacia(states%fields(2)%DATA(:,:,1), &
            this%states(glevel)%fields(6)%DATA(:,:,1))

          ! Determining the general solution:
          CALL addLinearSolutions(states%fields(1)%DATA(:,:,1), &
                                  states%fields(2)%DATA(:,:,1), &
                                  u, v, sg)
          psi = states%fields(1)%DATA(:,:,1)
          chi = states%fields(2)%DATA(:,:,1)

          ! Recalculate the boundary values after adding the general solution:
          CALL uvVelocityOnInterior(states%fields(1)%DATA(:,:,1), &
                                    states%fields(2)%DATA(:,:,1), &
                                    states%fields(3)%DATA(:,:,1), &
                                    states%fields(4)%DATA(:,:,1),sg)

          ! Plot the residual of vorticity field:
          BLOCK
            INTEGER(i_kind) :: imx,inx
            REAL(r_kind) :: amx,anx,alow,high

            imx = 0; inx = 0
            amx = 0.0D0; anx = 0.0D0
            DO i=1,num_cell
              IF (sg%cell_type(i) .EQ. 0) THEN
                IF (ABS(states%fields(6)%DATA(1,i,1)) .GT. amx) THEN
                  imx = i; amx = ABS(states%fields(6)%DATA(1,i,1))
                END IF
                IF (ABS(states%fields(5)%DATA(1,i,1)-vor(1,i)) .GT. anx) THEN
                  inx = i; anx = ABS(states%fields(5)%DATA(1,i,1)-vor(1,i))
                END IF
              END IF
            END DO
            IF (inx .GT. 0) THEN
              WRITE(*,33) inx,anx,vor(1,inx),states%fields(5)%DATA(1,inx,1),num_cell,sg%mpddInfo_sg%myrank
            ELSE
              PRINT*,'No vor difference in the interior! '
            END IF
        33  FORMAT('Max residual in vor: ',I6,' res: ',E14.6,' ana/num: ',2E14.6,' num_cell: ',I6,' pc:',I2)
            IF (imx .GT. 0) THEN
              WRITE(*,44) imx,amx,div(1,imx),states%fields(6)%DATA(1,imx,1),num_cell
            ELSE
              PRINT*,'No div difference in the interior! ',sg%mpddInfo_sg%myrank
            END IF
        44  FORMAT('Max residual in div: ',I6,' res: ',E14.6,' ana/num: ',2E14.6,' num_cell: ',I6)

            ! For interior points, plot the residual:
            states%fields(5)%DATA(:,:,2) = states%fields(5)%DATA(:,:,1) - vor

            ! Clean the boundary values:
            DO i=1,num_cell
              IF (sg%cell_type(i) .GT. 2) THEN
                states%fields(5)%DATA(:,i,2) = 0.0D0
              END IF
            END DO

            ! For the boundary points, plot the residual:
            ASSOCIATE(sg => geo%mg%sg(glevel))
              imx = 0; amx = 0.0D0; inx = 0
              alow = 1.0D10; high = -1.0D10
              
              DO i=1,sg%num_cell
                IF (sg%cell_type(i) .LT. 2) THEN
                  states%fields(1)%DATA(:,i,3) = (states%fields(4)%DATA(:,i,1) - v(:,i))
                  states%fields(3)%DATA(:,i,2) = states%fields(3)%DATA(:,i,1)
                  states%fields(3)%DATA(:,i,3) = u(:,i)
                  states%fields(4)%DATA(:,i,2) = states%fields(4)%DATA(:,i,1)
                  states%fields(4)%DATA(:,i,3) = v(:,i)
                  IF (ABS(states%fields(1)%DATA(1,i,3)) .GT. amx) THEN
                    imx = i
                    amx = ABS(states%fields(1)%DATA(1,i,3))
                  END IF
                END IF
              END DO
              IF (imx .NE. 0) THEN
                WRITE(*,333) imx,amx,states%fields(4)%DATA(:,imx,1),v(1,imx),sg%mpddInfo_sg%myrank
  333           FORMAT('Max diff in v: ',I6,' val: ',E14.6,' V: ',2E14.6,' pc:',I2)
              ELSE
                PRINT*,'No difference at the boundaries! ',scaling(2),sg%mpddInfo_sg%myrank
              END IF
            END ASSOCIATE
          END BLOCK
          CALL Output_NC_State_AV(this%states(glevel), TRIM(this%ncOutputFile), TRIM(this%task)//"_rhs", .TRUE., .TRUE.)
        END ASSOCIATE

        ! Deallocate memory:
        DEALLOCATE(rhs,xgmres,hgmres,vgmres,ua,va)
      END DO

    END SUBROUTINE Solve_s

    ! Project the residual in the null space out from the solution:
    SUBROUTINE rangeSpaceProjection_s(this,n,x,y)
      CLASS(mixedPsiChiSolver_t) :: this
      INTEGER(i_kind), INTENT(IN) :: n
      REAL(r_kind), INTENT(IN) :: x(n)
      REAL(r_kind), INTENT(OUT) :: y(n)

      ! Local variables:
      INTEGER(i_kind) :: i,j,k,m,nn,ioff(2)
      REAL(r_kind) :: s(n),z(n),d(n),xx,yy

      m = DSQRT(DBLE(n/2))
      ioff(1) = 0
      ioff(2) = m

      ! Project out the constant:
      d(1:n/2) = x(1:n/2)-SUM(x(1:n/2))/DBLE(n/2)
      d(n/2+1:n) = x(n/2+1:n)-SUM(x(n/2+1:n))/DBLE(n/2)
      y = d
      PRINT*,'MaxMin values: ',maxval(d),minval(d)
      PRINT*,'MaxMin values: ',maxval(x),minval(x)

      ! Project out the linear function:
      k = 0
      s = 0.0D0
      z = 0.0D0
      DO j=1,m
          DO i=1,m
              k = k+1
              s(k+ioff(1)) =  DBLE(i-1)/DBLE(m-1)-0.5D0
              z(k+ioff(1)) = -(DBLE(j-1)/DBLE(m-1)-0.5D0)
              s(k+ioff(2)) =  DBLE(j-1)/DBLE(m-1)-0.5D0
              z(k+ioff(2)) =  DBLE(i-1)/DBLE(m-1)-0.5D0
          END DO
      END DO
      d(1:n/2) = y(1:n/2) - DOT_PRODUCT(y(1:n/2),s(1:n/2))/DOT_PRODUCT(s(1:n/2),s(1:n/2))*s(1:n/2)
      y(1:n/2) = d(1:n/2) - DOT_PRODUCT(d(1:n/2),z(1:n/2))/DOT_PRODUCT(z(1:n/2),z(1:n/2))*z(1:n/2)
      d(n/2+1:n) = y(n/2+1:n) - DOT_PRODUCT(y(n/2+1:n),s(n/2+1:n))/DOT_PRODUCT(s(n/2+1:n),s(n/2+1:n))*s(n/2+1:n)
      y(n/2+1:n) = d(n/2+1:n) - DOT_PRODUCT(d(n/2+1:n),z(n/2+1:n))/DOT_PRODUCT(z(n/2+1:n),z(n/2+1:n))*z(n/2+1:n)

    END SUBROUTINE rangeSpaceProjection_s
END MODULE mixedPsiChiSolver_m