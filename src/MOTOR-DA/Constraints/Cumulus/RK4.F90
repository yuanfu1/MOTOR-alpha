!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather
! Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Jilong CHEN
! VERSION           : V 0.0
! HISTORY           :
!   Created by Jilong CHEN (jchen@link.cuhk.edu.hk), 2021/1/26, @GBA-MWF,
!   Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!!
SUBROUTINE RK4(this, X)

  USE CumFwd_m, ONLY: CumFwd
  IMPLICIT NONE

  CLASS(PhyCon_t), INTENT(INOUT) :: this
  TYPE(State_t), INTENT(IN) :: X
  INTEGER(i_kind) :: it1, it2, size_q_dim1, size_q_dim2, size_q_dim3
  ! REAL(r_kind) :: cu_dt
  REAL(r_kind), ALLOCATABLE :: dthdt(:, :, :), dqdt(:, :, :), &
                               theta_t(:, :), qvapor_t(:, :), &
                               dprecipdt(:, :)

  size_q_dim1 = SIZE(this%qvapor_out, dim=1)
  size_q_dim2 = SIZE(this%qvapor_out, dim=2)
  size_q_dim3 = SIZE(this%qvapor_out, dim=3)

  ALLOCATE (theta_t(size_q_dim1, size_q_dim2), &
            qvapor_t(size_q_dim1, size_q_dim2), &
            dprecipdt(size_q_dim2, 4))

  ! dthdt = 0.0D0
  ! dqdt = 0.0D0
  dprecipdt = 0.0D0
  DO it1 = 1, size_q_dim3 - 1

    DO it2 = 1, 4
      IF (it2 == 1) THEN
        theta_t = this%theta_out(:, :, it1)
        qvapor_t = this%qvapor_out(:, :, it1)
      ELSEIF (it2 == 4) THEN
        theta_t = this%theta_out(:, :, it1) + this%cu_dt * this%dthdt(:, :, it2 - 1, it1)
        qvapor_t = this%qvapor_out(:, :, it1) + this%cu_dt * this%dqdt(:, :, it2 - 1, it1)
      ELSE
        theta_t = this%theta_out(:, :, it1) + this%cu_dt * this%dthdt(:, :, it2 - 1, it1) / 2.0D0
        qvapor_t = this%qvapor_out(:, :, it1) + this%cu_dt * this%dqdt(:, :, it2 - 1, it1) / 2.0D0
      END IF

      CALL CumFwd(this%uwnd_out(:, :, it1), this%vwnd_out(:, :, it1), &
                  this%wwnd_out(:, :, it1), this%pres_out(:, :, it1), &
                  theta_t, qvapor_t, X, this%ztop, &
                  this%dthdt(:, :, it2, it1), this%dqdt(:, :, it2, it1), dprecipdt(:, it2))

    END DO
    this%theta_out(:, :, it1 + 1) = this%theta_out(:, :, it1) + &
                                    this%cu_dt * (this%dthdt(:, :, 1, it1) + &
                                                  2.0D0 * this%dthdt(:, :, 2, it1) &
                                                  + 2.0D0 * this%dthdt(:, :, 3, it1) &
                                                  + this%dthdt(:, :, 4, it1)) / 6.0D0
    this%qvapor_out(:, :, it1 + 1) = this%qvapor_out(:, :, it1) + &
                                     this%cu_dt * (this%dqdt(:, :, 1, it1) + &
                                                   2.0D0 * this%dqdt(:, :, 2, it1) &
                                                   + 2.0D0 * this%dqdt(:, :, 3, it1) &
                                                   + this%dqdt(:, :, 4, it1)) / 6.0D0

    this%precip(1, :, it1 + 1) = this%precip(1, :, it1) + &
                                 this%cu_dt * (dprecipdt(:, 1) + &
                                               2.0D0 * dprecipdt(:, 2) &
                                               + 2.0D0 * dprecipdt(:, 3) &
                                               + dprecipdt(:, 4)) / 6.0D0

  END DO

  DEALLOCATE (theta_t, &
              qvapor_t, &
              dprecipdt)

END SUBROUTINE
