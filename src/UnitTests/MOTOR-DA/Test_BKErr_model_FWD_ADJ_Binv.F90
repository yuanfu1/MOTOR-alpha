!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA.EnLoc.Test_EnLoc
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Jiongming Pang
! VERSION           : V 0.0   Calling EnLoc_m for BEC EnLoc calculation.
! HISTORY           :
!   Created by Jiongming Pang (pang.j.m@hotmail.com), 2023/1/19, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!!
PROGRAM Test_BKErr_model_FWD_ADJ_Binv

  USE kinds_m
  USE State_m, ONLY: State_t
  USE EnLoc_m, ONLY: EnLoc_t

  IMPLICIT NONE

  TYPE(EnLoc_t)       :: EnLoc
  TYPE(State_t)       :: Xm(1)
  CHARACTER(LEN=1024) :: configFile

  INTEGER(i_kind) :: i, j

  ! Get the configFile
  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", configFile)
  configFile = TRIM(configFile)//"/UnitTest/testEnLoc_model.yaml"

  ! CALL EnLoc%b_getEnsdata(configFile, Xm)
  CALL EnLoc%b_getEnsdata(configFile)
  PRINT *, 'Done EnLoc%b_getEnsdata ...'

  CALL EnLoc%b_qrDecomp()
  PRINT *, 'Done EnLoc%b_qrDecomp ...'

  CALL EnLoc%b_qrSave()
  PRINT *, 'Done EnLoc%b_qrSave ...'

  CALL EnLoc%b_qrLoad()
  PRINT *, 'Done EnLoc%b_qrLoad ...'

  CALL EnLoc%b_qrRinv()
  PRINT *, 'Done EnLoc%b_qrRinv ...'

  CALL EnLoc%b_FWD()
  PRINT *, 'Done EnLoc%b_FWD ...'

  CALL EnLoc%b_ADJ()
  PRINT *, 'Done EnLoc%b_ADJ ...'

  ! PRINT *, "nz:",EnLoc%BKErr(1)%nz
  ! PRINT *, "nc:",EnLoc%BKErr(1)%nc

  ! OPEN(813,FILE="./RinvQinvI.txt",STATUS='REPLACE')
  ! DO i=1,EnLoc%BKErr(1)%nc
  !    WRITE(UNIT=813, FMT="(100F18.8)") EnLoc%BKErr(1)%RinvQinvI(:,i)
  ! END DO
  ! CLOSE(813)

  ! OPEN(813,FILE="./QRinvTI.txt",STATUS='REPLACE')
  ! DO j=1,EnLoc%BKErr(1)%nc
  !    WRITE(UNIT=813, FMT="(100F18.8)") EnLoc%BKErr(1)%QRinvTI(j,:)
  ! END DO
  ! CLOSE(813)

  CALL EnLoc%b_destroy()

END PROGRAM Test_BKErr_model_FWD_ADJ_Binv
