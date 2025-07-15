!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA.poissonSolver_Kp
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Yuanfu Xie, Zilong Qin
! VERSION           : V 0.0
! HISTORY           :
!   Created by Yuanfu Xie, in test_rcsv_mltgrd.F90
!   Reforged by Zilong Qin (zilong.qin@gmail.com), 2020/12/15, @GBA-MWF, Shenzhen, add parallelization
!   Yuanfu Xie modified this test for automatically checking error reductions: 2025-04-16
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module contains the data type for model_states
!! method.
!! @copyright (C) 2020 GBA-MWF, All rights reserved.
! @note
! @warning
! @attention
PROGRAM Test_PSAdjoint
  USE kinds_m, ONLY: i_kind, r_kind
  ! USE RossbyHaurwitzSphere_m, ONLY: RossbyHaurwitzSphere_t ! rossbyhaurwitz function from utility
  USE UnitTestDyn_m, ONLY: UnitTestDyn_t
  USE poissonSolver_m, ONLY: poissonSolver_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE YAMLRead_m
  USE geometry_m, ONLY: geometry_t

!  USE mpi_f08
  INCLUDE "mpif.h"

  TYPE(mpddGlob_t), TARGET :: mpddGlob
  TYPE(geometry_t), TARGET :: geometry
  TYPE(poissonSolver_t) :: ps
  TYPE(UnitTestDyn_t) :: rossby

  REAL(r_kind), ALLOCATABLE :: ValuVort(:, :)         !> Container of Divergence & Vorticity as RHS
  REAL(r_kind), ALLOCATABLE :: tempSoluStrm(:, :)     !> Container of Stream function and Velocity potential as initial value
  REAL(r_kind), ALLOCATABLE :: trueSoluStrm(:, :)     !> Container of Stream function and Velocity potential as true solution
  REAL(r_kind), ALLOCATABLE :: reltSoluStrm(:, :)     !> Container for poissonSolver output
  REAL(r_kind), ALLOCATABLE :: rhs(:, :)              !> Container for poissonSolver output
  REAL(r_kind), ALLOCATABLE :: sigma_3d(:, :)

  INTEGER(i_kind) :: vLevel
  INTEGER(i_kind) :: i, ig
  LOGICAL :: PASS
  REAL(r_kind) :: tempMax, tempMin, mx, mn, t1, t2, amx(10)
  CHARACTER(LEN=1024) :: configFile
  REAL(r_kind) :: YHX_sum, XHTY_sum
  
  ! Get the configFile
  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", configFile)
  configFile = TRIM(configFile)//"/UnitTest/Test_PoissonSolver.yaml"

  IF (yaml_get_var(configFile, 'modelState', 'vLevel', vLevel) /= 0) STOP

  ! Initialize libraries
  CALL mpddGlob%initialize()
  CALL geometry%initialize(configFile, mpddGlob)    ! Initialize the geometry
  ! GenContainers = GenContainers_t(TRIM(configFile))
  ! CALL GenContainers%GenGeometry(geometry)
  ps = poissonSolver_t(configFile, geometry)   ! Initialize the poissonSolver

  ! Yuanfu Xie modified this test for automatically checking error reductions:
  DO ig = 6, 6
    ASSOCIATE (sgFinest => geometry%mg%sg(ig))

      CALL sgFinest%allocateMat(vLevel, ValuVort)
      CALL sgFinest%allocateMat(vLevel, tempSoluStrm)
      CALL sgFinest%allocateMat(vLevel, trueSoluStrm)
      CALL sgFinest%allocateMat(vLevel, reltSoluStrm)
      CALL sgFinest%allocateMat(vLevel, rhs)
      CALL sgFinest%allocateMat(vLevel, sigma_3d)

      ! Initialize the streamfunction and vorticity
      CALL Rossby%initial(sgFinest) ! Beaware that this initialization determine the cache size for generating the corresponding values

      DO i = 1, sgFinest%num_cell
        sigma_3d(:, i) = sgFinest%sigma
      END DO
      ! Undate the inner vars
      CALL Rossby%RH_psi_func(-600.0D0, sgFinest%cell_cntr, sigma_3d, tempSoluStrm, ValuVort)
      CALL Rossby%RH_psi_func(0.0D0, sgFinest%cell_cntr, sigma_3d, trueSoluStrm, ValuVort)

      ! Update the vars on each level
      ! DO i = 2, vLevel
      !    tempSoluStrm(i, :) = tempSoluStrm(1, :)
      !    trueSoluStrm(i, :) = trueSoluStrm(1, :)
      !    ValuVort(i, :) = ValuVort(1, :)
      ! END DO

      ! Define the right hand side, include the boundary infos.
      DO i = 1, sgFinest%num_icell
        IF (sgFinest%cell_type(i) .EQ. 0) THEN              ! On the inner cells, rhs is the vorticity
          rhs(:, i) = ValuVort(:, i)
        ELSE IF ((sgFinest%cell_type(i) .EQ. 1) .OR. (sgFinest%cell_type(i) .EQ. 2)) THEN         ! On the boundary cells, rhs is the stream function
          rhs(:, i) = trueSoluStrm(:, i)
        END IF
      END DO

      ! Update the vars on the boundary using mpi
      CALL sgFinest%ExchangeMatOnHalo2D(vLevel, tempSoluStrm)
      CALL sgFinest%ExchangeMatOnHalo2D(vLevel, trueSoluStrm)
      CALL sgFinest%ExchangeMatOnHalo2D(vLevel, ValuVort)
      CALL sgFinest%ExchangeMatOnHalo2D(vLevel, rhs)

      BLOCK
        REAL(r_kind), ALLOCATABLE :: YY(:, :), XX(:, :), HX(:, :), HTY(:, :)
        REAL(r_kind) :: tempYHX_sum, tempXHTY_sum
        
        CALL sgFinest%allocateMat(vLevel, HX)
        CALL sgFinest%allocateMat(vLevel, HTY)
        CALL sgFinest%allocateMat(vLevel, XX)
        CALL sgFinest%allocateMat(vLevel, YY)

        YY = trueSoluStrm
        XX = rhs
        HX = 0.0_R_KIND
        HTY = 0.0_R_KIND

        CALL sgFinest%ExchangeMatOnHalo2D(vLevel, YY)
        CALL sgFinest%ExchangeMatOnHalo2D(vLevel, XX)

        DO i = 1, sgFinest%num_icell
          IF (sgFinest%cell_type(i) .EQ. 0) THEN              ! On the inner cells, rhs is the vorticity
            rhs(:, i) = ValuVort(:, i)
          ELSE IF ((sgFinest%cell_type(i) .EQ. 1) .OR. (sgFinest%cell_type(i) .EQ. 2)) THEN         ! On the boundary cells, rhs is the stream function
            rhs(:, i) = trueSoluStrm(:, i)
          END IF
        END DO

        DO i = 1, sgFinest%num_icell
          IF ((sgFinest%cell_type(i) .EQ. 1) .OR. (sgFinest%cell_type(i) .EQ. 2)) THEN         ! On the boundary cells, rhs is the stream function
            HX(:, i) = XX(:, i)
          END IF
        END DO
        CALL sgFinest%ExchangeMatOnHalo2D(vLevel, HX)
        CALL ps%PoissonSol(sgFinest%gLevel, vLevel, XX, HX, 'forward')

        DO i = 1, sgFinest%num_icell
          IF ((sgFinest%cell_type(i) .EQ. 1) .OR. (sgFinest%cell_type(i) .EQ. 2)) THEN         ! On the boundary cells, rhs is the stream function
            HTY(:, i) = YY(:, i)
          END IF
        END DO
        CALL sgFinest%ExchangeMatOnHalo2D(vLevel, HTY)
        CALL ps%PoissonSol(sgFinest%gLevel, vLevel, YY, HTY, 'adjoint')

        tempYHX_sum = SUM(YY(:, 1:sgFinest%num_icell) * HX(:, 1:sgFinest%num_icell))
        tempXHTY_sum = SUM(XX(:, 1:sgFinest%num_icell) * HTY(:, 1:sgFinest%num_icell))

        CALL MPI_REDUCE(tempYHX_sum, YHX_sum, 1, mpi_double_precision, MPI_SUM, mpddGlob%rankBase, mpddGlob%comm, mpddGlob%ierr)
        CALL MPI_REDUCE(tempXHTY_sum, XHTY_sum, 1, mpi_double_precision, MPI_SUM, mpddGlob%rankBase, mpddGlob%comm, mpddGlob%ierr)

        IF (mpddGlob%isBaseProc()) PRINT *, 'Test AD YHX:', YHX_sum
        IF (mpddGlob%isBaseProc()) PRINT *, 'Test AD HTYX:', XHTY_sum
      END BLOCK

    END ASSOCIATE

    ! Destroy the structures
    CALL rossby%destroy()
  END DO

  ! Release the ALLOCATED memory
  DEALLOCATE (ValuVort, &
              trueSoluStrm, &
              tempSoluStrm, &
              reltSoluStrm, &
              rhs &
              )

  ! For ctest
  IF (mpddGlob%isBaseProc()) THEN
    print*, 'YHX/XHTY = ', YHX_sum/XHTY_sum
    IF (ABS(YHX_sum/XHTY_sum-1)<1e-7) THEN
      PRINT *, 'Test passed!'
    ELSE
      PRINT *, 'Test failed!'
    END IF
  END IF
  CALL mpddGlob%finalize

END PROGRAM Test_PSAdjoint
