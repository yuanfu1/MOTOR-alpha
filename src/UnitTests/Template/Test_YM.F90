!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2021/1/26, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!!
PROGRAM Test_Temp
  USE Template_m, ONLY: Template_t
  !USE NMLRead_m, ONLY: namelist_read
  USE YAMLREAD_m

  CHARACTER(LEN=1024) :: configFile
  CHARACTER(LEN=1024) :: filename_obs
  INTEGER :: istatus

  ! Get the configFile
  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", configFile)
  configFile = TRIM(configFile)//"/template.yaml"

  !CALL namelist_read(configFile, 'filename_obs', filename_obs)
  istatus = yaml_get_var(TRIM(configFile), 'IO', 'ModelFileName', filename_obs)

  PRINT *, TRIM(filename_obs)
END PROGRAM Test_Temp
