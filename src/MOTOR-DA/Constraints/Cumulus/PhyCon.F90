!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA-Cumulus
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather
! Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Jilong CHEN, Zilong Qin
! VERSION           : V 0.0
! HISTORY           :
!   Created by Jilong CHEN (jchen@link.cuhk.edu.hk), 2021/9/2, @GBA-MWF,
!   Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!!
MODULE PhyCon_m
  USE State_m, ONLY: State_t
  USE parameters_m, ONLY: g, cp, Lv, k_d, r_k_d, epsw_r, gamma_d, r_d
  USE kinds_m, ONLY: i_kind, r_kind

  TYPE PhyCon_t
    CHARACTER(len=1024) :: configFile
    REAL(r_kind), ALLOCATABLE :: uwnd_out(:, :, :), vwnd_out(:, :, :), wwnd_out(:, :, :), &
                                 pres_out(:, :, :), theta_out(:, :, :), qvapor_out(:, :, :), &
                                 precip(:, :, :), uwnd_outb(:, :, :), vwnd_outb(:, :, :), &
                                 wwnd_outb(:, :, :), pres_outb(:, :, :), theta_outb(:, :, :), &
                                 qvapor_outb(:, :, :), precipb(:, :, :), tt_out(:), dthdt(:, :, :, :), &
                                 dqdt(:, :, :, :)
    REAL(r_kind) :: cu_dt = 3.0D02
    REAL(r_kind) :: ztop = 3.0D04
    INTEGER(i_kind) :: size_3m = 10000
    INTEGER(i_kind) :: kt

  CONTAINS

    PROCEDURE, PUBLIC, PASS(this) :: FwdDataPre
    PROCEDURE, PUBLIC, PASS(this) :: AdjDataPre
    PROCEDURE, PUBLIC, PASS(this) :: RK4
    PROCEDURE, PUBLIC, PASS(this) :: RK4_AD
    FINAL :: destroy
  END TYPE PhyCon_t

  INTERFACE PhyCon_t
    PROCEDURE :: constructor
  END INTERFACE PhyCon_t

CONTAINS

  INCLUDE 'FwdDataPre.F90'
  INCLUDE 'AdjDataPre.F90'
  INCLUDE 'RK4.F90'
  INCLUDE 'RK4_AD.F90'

  IMPURE ELEMENTAL SUBROUTINE destroy(this)
    IMPLICIT NONE
    TYPE(PhyCon_t), INTENT(INOUT) :: this

    IF (ALLOCATED(this%uwnd_out)) DEALLOCATE (this%uwnd_out, &
                                              this%vwnd_out, this%wwnd_out, this%theta_out, &
                                              this%qvapor_out, this%pres_out, this%precip, &
                                              this%uwnd_outb, &
                                              this%vwnd_outb, this%wwnd_outb, this%theta_outb, &
                                              this%qvapor_outb, this%pres_outb, this%precipb, &
                                              this%tt_out, this%dthdt, this%dqdt)
  END SUBROUTINE

  FUNCTION constructor(configFile) RESULT(this)
    IMPLICIT NONE
    TYPE(PhyCon_t) :: this
    CHARACTER(len=1024) :: configFile

    this%configFile = configFile

  END FUNCTION constructor

END MODULE PhyCon_m
