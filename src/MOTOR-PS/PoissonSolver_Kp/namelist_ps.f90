!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   2020-4 created by Yuanpfu Xie
!   This file is reforged from fsw_namelist_m by Yuanfu Xie
!   Reforged by Zilong Qin (zilong.qin@gmail.com), 2020/12/18, @GBA-MWF, Shenzhen
!   Modified by Jiongming Pang (pang.j.m@hotmail.com), 2022/05/12, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!!  This is a module reading in all needed namelist files.
!! @author Zilong Qin, Yuanfu Xie
!! @copyright (C) 2020 GBA-MWF, All rights reserved.
! @note
! @warning
! @attention
MODULE namelist_ps_m

  USE kinds_m, ONLY: i_kind, r_kind
  USE YAMLRead_m

CONTAINS
  SUBROUTINE namelist_ps(configFile, solver, nCycle, nIterPre, nIterPost, nRelax, omegas, max_band)
    ! Model namelist variables:
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    CHARACTER(LEN=3), INTENT(OUT) :: solver
    INTEGER(i_kind), INTENT(OUT) :: &
      nCycle, &
      nIterPre, &
      nIterPost, &
      nRelax, &
      max_band
    REAL(r_kind), ALLOCATABLE, INTENT(OUT) :: omegas(:)
    INTEGER(i_kind)                        :: rc

    INTEGER(i_kind) :: istatus

    istatus = yaml_get_var(TRIM(configFile), 'poissonSovler', 'solver', solver)
    istatus = yaml_get_var(TRIM(configFile), 'poissonSovler', 'nCycle', nCycle)
    istatus = yaml_get_var(TRIM(configFile), 'poissonSovler', 'nRelax', nRelax)
    istatus = yaml_get_var(TRIM(configFile), 'poissonSovler', 'omegas', omegas)
    istatus = yaml_get_var(TRIM(configFile), 'poissonSovler', 'max_band', max_band)
    istatus = yaml_get_var(TRIM(configFile), 'poissonSovler', 'nIterPre', nIterPre)
    istatus = yaml_get_var(TRIM(configFile), 'poissonSovler', 'nIterPost', nIterPost)

  END SUBROUTINE namelist_ps

END MODULE namelist_ps_m
