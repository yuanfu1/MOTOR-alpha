!!---------------------------------------------------------------------------------------
! PROJECT           : PRIVATE-PROJECT: OO Based NWP model
! AFFILIATION       : Self-employed
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : V 1
! HISTORY           :
!   Created  by Yuanfu Xie, 2021 Broomfield, CO, USA
!   Modified by
!!---------------------------------------------------------------------------------------

!> @brief
!! Object of a 4th order Runge-Kutta scheme.
!! @author Yuanfu Xie
!! @copyright (C) 2021 Yuanfu Xie, All rights reserved.
!! @note This allows users to switch from different model grids.
!! @warning
!! @attention
MODULE rk4_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE timeIntegral_m, ONLY: timeIntegral_t
  USE nwp_m, ONLY: nwp_t
  USE nwp_zgrid_m, ONLY: nwp_zgrid_t
  USE modelState_m, ONLY: modelState_t
  USE Helmholtz_m, ONLY: Helmholtz_t
  USE parameters_m, ONLY: Omega

  TYPE, EXTENDS(timeIntegral_t) :: rk4_t
    TYPE(modelState_t) :: yn
    TYPE(modelState_t), ALLOCATABLE :: steps(:)
    ! a Poisson solver:
    TYPE(Helmholtz_t) :: Poisson
  CONTAINS
    FINAL :: destructor
    PROCEDURE :: timeIntegrals => timeIntegrals_rk4
    PROCEDURE :: outputStates => outputStates_rk4
  END TYPE rk4_t

  INTERFACE rk4_t
    PROCEDURE :: constructor
  END INTERFACE rk4_t

CONTAINS

  FUNCTION constructor(nwp) RESULT(this)
    IMPLICIT NONE

    TYPE(rk4_t) :: this
    CLASS(nwp_t), TARGET :: nwp

    ! Local variables:
    LOGICAL :: load
    INTEGER(i_kind) :: i

    this%nwp => nwp

    ! Allocate memory for the step in a RK4:
    ! I intentionally leave this structure here
    ! for possible improvement of removing modelState_t structure
    SELECT TYPE (p => this%nwp)
    TYPE IS (nwp_zgrid_t)
      ALLOCATE (this%steps(4))
      CALL this%yn%copy(p%gzm%states)
      DO i = 1, 4
        CALL this%steps(i)%copy(p%gzm%states)
      END DO
      load = .TRUE.
      this%Poisson = Helmholtz_t(p%geo%mg, load)
    CLASS DEFAULT
      WRITE (*, 1)
1     FORMAT('Currently rk4_t support nwp_zgrid_t only!')
      STOP
    END SELECT
  END FUNCTION constructor

  SUBROUTINE timeIntegrals_rk4(this, dt, numSteps, intOutpt)
    IMPLICIT NONE

    CLASS(rk4_t) :: this
    ! Total time integration steps and output intervals:
    INTEGER(i_kind), INTENT(IN) :: numSteps, intOutpt
    REAL(r_kind), INTENT(IN) :: dt

    ! Local variables:
    INTEGER(i_kind) :: i, it, k
    REAL(r_kind) :: dk(3)
    REAL :: start_time, end_time

    this%stepLength = dt

    ! RK4: y := y + (k1 + 2 k2 + 2 k3 + k4) * (dt/6)
    dk = (/0.5D0 * dt, 0.5D0 * dt, dt/)

    SELECT TYPE (p => this%nwp)
    TYPE IS (nwp_zgrid_t)
      DO it = 1, numSteps
        ! K1 = f(t_n,y_n):
        CALL p%rightHandSide(p%gzm%states, this%steps(1))

        ! RK4 for K2-K4
        DO k = 1, 3
          ! yn + K(k) * dt/2:
          this%yn%fields(p%idx_vort)%DATA = &
            p%gzm%states%fields(p%idx_vort)%DATA + &
            (this%steps(k)%fields(p%idx_vort)%DATA) * dk(k)
          ! Relative vorticity:
          DO i = 1, p%gzm%grid%vLevel
            p%rlv(i, :) = &
              this%yn%fields(p%idx_vort)%DATA(i, :) - &
              2.0D0 * Omega * DSIN(p%gzm%grid%cell_cntr(1, :))
          END DO
          ! Convert the vorticity to stream function:
          CALL CPU_TIME(start_time)
          !CALL this%Poisson%eqnsSolver(1,p%geo%mg%mg_finest, &
          !    p%rlv,this%yn%fields(p%idx_strm)%data)
          ! Fixed streamfunction:
          this%yn%fields(p%idx_strm)%DATA = p%gzm%states%fields(p%idx_strm)%DATA
          CALL CPU_TIME(end_time)
          !WRITE(*,2) end_time-start_time
2         FORMAT('Poisson solver time used: ', D12.4)

          ! K(k+1) = f(t_n+dk(k), yn + k1 * dk(k))
          CALL p%rightHandSide(this%yn, this%steps(k + 1))
        END DO

        ! Relative vorticity:
        p%gzm%states%fields(p%idx_vort)%DATA = &
          p%gzm%states%fields(p%idx_vort)%DATA + &
          (this%steps(1)%fields(p%idx_vort)%DATA + &
           this%steps(2)%fields(p%idx_vort)%DATA * 2.0D0 + &
           this%steps(3)%fields(p%idx_vort)%DATA * 2.0D0 + &
           this%steps(4)%fields(p%idx_vort)%DATA) * (dt / 6.0D0)
        DO k = 1, p%gzm%grid%vLevel
          p%rlv(k, :) = p%gzm%states%fields(p%idx_vort)%DATA(k, :) - &
                        2.0D0 * Omega * DSIN(p%gzm%grid%cell_cntr(1, :))
        END DO
        !CALL this%Poisson%eqnsSolver(1,p%geo%mg%mg_finest, &
        !    p%rlv,p%gzm%states%fields(p%idx_strm)%data)

        ! Output solutions:
        IF (MOD(INT(it * dt), intOutpt) .EQ. 0) THEN
          PRINT *, 'Output?: ', INT(it * dt / intOutpt), &
            MAXVAL(ABS(p%gzm%states%fields(p%idx_vort)%DATA)), &
            SQRT(DOT_PRODUCT(p%gzm%states%fields(p%idx_vort)%DATA(1, :), &
                             p%gzm%states%fields(p%idx_vort)%DATA(1, :)))
          CALL this%outputStates(INT(it * dt / intOutpt))
        END IF
      END DO
    CLASS DEFAULT
      WRITE (*, 1)
1     FORMAT('Currently rk4_t support nwp_zgrid_t only!')
      STOP
    END SELECT
  END SUBROUTINE timeIntegrals_rk4

  SUBROUTINE outputStates_rk4(this, intOutpt)
    IMPLICIT NONE
    CLASS(rk4_t) :: this
    INTEGER(i_kind), INTENT(IN) :: intOutpt  ! Output interval in seconds

    ! Local variables:
    CHARACTER(LEN=256) :: file2plot
    INTEGER(i_kind) :: i

    ! RK4: y := y + (k1 + 2 k2 + 2 k3 + k4) * (dt/6)
    SELECT TYPE (p => this%nwp)
    TYPE IS (nwp_zgrid_t)
      DO i = 1, p%numVars
        IF (p%gzm%states%fields(i)%TYPE .EQ. 'strm') THEN
          file2plot = "model_outpsi"
        ELSE IF (p%gzm%states%fields(i)%TYPE .EQ. 'vort') THEN
          file2plot = "model_outvrt"
        ELSE IF (p%gzm%states%fields(i)%TYPE .EQ. 'hght') THEN
          file2plot = "model_outthk"
        END IF
        WRITE (file2plot(13:14), '(I2.2)') p%gzm%grid%gLevel
        file2plot(15:15) = "_"
        WRITE (file2plot(16:18), '(I3.3)') intOutpt
        OPEN (10, file=file2plot(1:18), form='unformatted')
        WRITE (10) p%gzm%grid%num_cell
        WRITE (10) p%gzm%grid%cell_cntr
        WRITE (10) p%gzm%states%fields(i)%DATA(1, :)
        CLOSE (10)
      END DO
    CLASS DEFAULT
      WRITE (*, 1)
1     FORMAT('Currently rk4_t support nwp_zgrid_t only!')
      STOP
    END SELECT

  END SUBROUTINE outputStates_rk4

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    TYPE(rk4_t), INTENT(INOUT) :: this

    IF (ALLOCATED(this%steps)) DEALLOCATE (this%steps)
  END SUBROUTINE destructor
END MODULE rk4_m
