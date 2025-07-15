!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Yuanfu Xie (yuanfu_xie@yahoo.com), 2024/10/11, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------
!> @brief
!! This module defines GZM dynamic core
!!
MODULE dyCoreGZM_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE State_m, ONLY: State_t
  USE TimeIntegrationBase_m, ONLY: TimeIntegrationBase_t
  USE rhsBase_m, ONLY: rhsBase_t
  USE dyCoresBase_m, ONLY: dyCoresBase_t
  USE poissonSolver_m, ONLY: poissonSolver_t
  USE SingleGrid_m, ONLY: SingleGrid_t

  TYPE, EXTENDS(dyCoresBase_t) :: dyCoreGZM_t
    TYPE(poissonSolver_t) :: poissonSolver
  CONTAINS
    PROCEDURE, PUBLIC :: initialize_s => initialGZM_s
    PROCEDURE, PUBLIC :: dyCore_s => dyCoreGZM_s
    PROCEDURE, PUBLIC :: dyAdjoint_s => dyAdjGZM_s
  END TYPE dyCoreGZM_t

CONTAINS
  SUBROUTINE initialGZM_s(this, yamlFile, tim, rhs, sg, dt)
    CLASS(dyCoreGZM_t) :: this
    CHARACTER(LEN=1024), INTENT(IN) :: yamlFile
    CLASS(TimeIntegrationBase_t), TARGET :: tim
    CLASS(rhsBase_t), TARGET :: rhs
    CLASS(SingleGrid_t), TARGET :: sg
    REAL(r_kind), INTENT(IN) :: dt

    ! Pass the model parameters to the dyCoreGZM_t
    this%yamlFile = yamlFile
    this%tim => tim
    this%rhs => rhs
    this%sg => sg
    this%dt = dt
    ! this%poissonSolver = poissonSolver_t(this%yamlFile, geo)
  END SUBROUTINE initialGZM_s

  SUBROUTINE dyCoreGZM_s(this, iC, X)
    IMPLICIT NONE
    CLASS(dyCoreGZM_t) :: this
    TYPE(State_t), INTENT(INOUT) :: iC
    TYPE(State_t), INTENT(INOUT) :: X

    ! Local variables:
    INTEGER(i_kind) :: i, it, intvl, msteps
    TYPE(State_t) :: yk

    ! Time integration step length based on the model frames at this mg level:
    IF (MOD((this%sg%tt(2) - this%sg%tt(1)), this%dt) .NE. 0.0D0) THEN
      WRITE (*, 2) this%dt, this%sg%tt(2) - this%sg%tt(1)     ! dt and a time interval of X
2     FORMAT('Time step is not exactly matching time window of X: ', 2E14.6, /, &
             ' Please choose proper dt and rerun')
      STOP
    END IF
    intvl = (this%sg%tt(2) - this%sg%tt(1)) / this%dt
    this%intval = intvl                           ! pass the interval info

    msteps = UBOUND(iC%fields(1)%DATA, 3) 
    CALL yk%initialize(this%yamlFile, this%sg, msteps)

    ! Model integration:
    WRITE (*, 3)  X%sg%gLevel
3   FORMAT('dyCoreGZM_s - Model is integrating at G ',I2,' ...')
    DO it = 1, X%sg%tSlots
      DO i = 1, intvl
        ! WRITE(*,4) i + intvl * (it-1), intvl * X%sg%tSlots, X%sg%gLevel
4       FORMAT('dyCoreGZM_s - Integrating time step: ', I8, ' of ',I8,' at Glevel: ',I3)
        CALL this%tim%TimeIntegrationBase_fwd(this%dt, i, it, iC, X, yk)
        PRINT*,'dyCoreGZM_s - values: ',(iC .DOT. iC), (X .DOT. X), (yk .DOT. yk), &
          ' at time step: ', it, ' and interval: ', i, ' of ', intvl, ' at Glevel: ', X%sg%gLevel
        ! Calculate the streamfuncion and velocity potential:
        ! CALL solve_poisson(this%poissonSolver, X%sg%vLevel, X%sg%gLevel, &
        !   vorticity(:, :, it), divergence(:, :, it), streamfunc, velocity_pot, X%sg)
        ! CALL solve_poisson(this%poissonSolver,it,X)
      END DO
    END DO

  END SUBROUTINE dyCoreGZM_s

  SUBROUTINE dyAdjGZM_s(this, it, iC, Xadj)
    IMPLICIT NONE
    CLASS(dyCoreGZM_t) :: this
    INTEGER(i_kind), INTENT(IN) :: it
    TYPE(State_t), INTENT(INOUT) :: iC
    TYPE(State_t), INTENT(INOUT) :: Xadj

    ! Local variables:
    INTEGER(i_kind) :: i

    DO i = this%intval, 1, -1
      CALL this%tim%TimeIntegrationBase_adj(this%dt, it, iC, Xadj)
      ! Calculate the adjoint for streamfuncion and velocity potential:
      ! CALL solve_poisson(this%poissonSolver, X%sg%vLevel, X%sg%gLevel, &
      !   vorticity(:, :, it), divergence(:, :, it), streamfunc, velocity_pot, X%sg)
      ! CALL solve_poisson(this%poissonSolver,it,X)
    END DO
  END SUBROUTINE dyAdjGZM_s
END MODULE dyCoreGZM_m
