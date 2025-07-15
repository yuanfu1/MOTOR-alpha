!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA.EnLoc.Test_EnLoc
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Jiongming Pang
! VERSION           : V 0.0   Calling EnLoc_m for BEC EnLoc calculation.
! HISTORY           :
!   Created by Jiongming Pang (pang.j.m@hotmail.com), 2022/3/30, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!!
PROGRAM Test_EnLoc_model

  USE kinds_m
  USE State_m, ONLY: State_t
  USE EnLoc_m, ONLY: EnLoc_t

  IMPLICIT NONE

  TYPE(EnLoc_t)       :: EnLoc
  TYPE(State_t)       :: Xm(1)
  CHARACTER(LEN=1024) :: configFile

  ! Get the configFile
  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", configFile)
  configFile = TRIM(configFile)//"/UnitTest/testEnLoc_model.yaml"

  CALL EnLoc%b_getEnsdata(configFile, Xm)
  PRINT *, 'Done EnLoc%b_getEnsdata ...'

  CALL EnLoc%b_calEigen()
  PRINT *, 'Done EnLoc%b_calEigen ...'

  CALL EnLoc%b_jacoEigenT()
  PRINT *, 'Done EnLoc%b_jacoEigenT ...'

  CALL EnLoc%b_destroy()

END PROGRAM Test_EnLoc_model
