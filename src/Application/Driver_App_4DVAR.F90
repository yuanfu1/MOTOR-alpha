!!--------------------------------------------------------------------------------------------------
! PROJECT           : Application of MOTOR-DA 4DVAR
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : V 1.0
! HISTORY           : 2024-01-17, created by Yuanfu Xie.
!
!   Created by Yuanfu Xie (yuanfu_xie@yahoo.com), 2024/01/17, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------
!!
!> @brief
!! This is a MOTOR-DA 4DVAR with analysis
PROGRAM App_4DVAR
  USE App_4DVAR_m
  USE bkgMG_m
  USE ObsMG_m
  USE cDefault_m, ONLY: cDefault_t
  USE c4DVar_m, ONLY: c4DVar_t
  USE dyCoreGZM_m, ONLY: dyCoreGZM_t
  USE TimeIntegrationAB3_m, ONLY: TimeIntegrationAB3_t
  USE rhsShallowWater_m, ONLY: rhsShallowWater_t

  IMPLICIT NONE

  TYPE(App_4DVAR_t) :: app
  TYPE(bkgMG_t) :: bkg
  TYPE(ObsMG_t) :: obs

  INTEGER(i_kind) :: istatus, mgStart, mgEnd
  REAL(r_kind) :: dt

  ! Time integral step at the finest grid:
  dt = 30.0D0

  CALL app%initialize_s(istatus)

  ! Initializing bkgd fields:
  WRITE(*,1)
1 FORMAT('Driver_App_4DVAR -- Initializing background fields....')
  CALL bkg%initialize_s(app%geometry, app%yamlVars%yamlFile)
  mgStart = app%geometry%mg%mg_coarsest
  mgEnd = app%geometry%mg%mg_finest

  ! Reading in model background fields:
  WRITE(*,2)
2 FORMAT('Driver_App_4DVAR -- Reading in model background fields...')
  CALL bkg%getBackgrd_s(app%yamlVars%yamlFile)
  CALL bkg%writeBkgd_s(app%yamlVars%yamlFile)

  ! Initializing obs data and reading in data:
  WRITE(*,3)
3 FORMAT('Driver_App_4DVAR -- Initializing and reading obs...')
  CALL obs%initialize_s(app%geometry, app%yamlVars%yamlFile)
  CALL obs%ingestions_s(bkg%XbMG(UBOUND(bkg%XbMG, 1)))

  ! Solve for the analysis:
  WRITE(*,8) mgStart, mgEnd
8 FORMAT('Driver_App_4DVAR -- analyzing...MG from ',I2,' to ',I2)
  CALL app%analysis_s(obs, bkg, mgStart, mgEnd, istatus, dt)

  WRITE(*,4) istatus
4 FORMAT('Driver_App_4DVAR -- Analysis completed with status: ',I2)

  CALL bkg%destroy_s
  CALL app%destroy_s

END PROGRAM App_4DVAR
