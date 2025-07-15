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
PROGRAM Test_EnLoc_state

  USE kinds_m
  USE geometry_m, ONLY: geometry_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE State_m, ONLY: State_t
  USE IOWRF_m, ONLY: IOWRF_t
  USE EnLoc_m, ONLY: EnLoc_t
  USE State2NC_m
  !USE NMLRead_m
  USE YAMLRead_m

  IMPLICIT NONE

  ! Define types
  TYPE(mpddGlob_t), TARGET :: mpddGlob
  TYPE(geometry_t), TARGET :: geometry

  TYPE(IOWRF_t), ALLOCATABLE, TARGET :: ioWRF
  ! TYPE(IOWRF_t),  target :: ioWRF
  TYPE(State_t), ALLOCATABLE         :: Xm(:)
  TYPE(EnLoc_t)       :: EnLoc
  CHARACTER(LEN=1024) :: configFile, bakFilesDir, outputDir
  CHARACTER(LEN=20), ALLOCATABLE :: varList(:)
  CHARACTER(LEN=20)   :: varName
  CHARACTER(LEN=3)    :: memid
  INTEGER(i_kind)     :: varNum, iv, ensNum, iens, ifile

  ! Get the configFile
  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", configFile)
  configFile = TRIM(configFile)//"/UnitTest/testEnLoc_state.yaml"

  IF (yaml_get_var(configFile, 'IO', 'output_dir', outputDir) /= 0) STOP

  !CALL namelist_read(TRIM(configFile), "varList", varList)
  ifile = yaml_get_var(TRIM(configFile), 'modelState', 'varList', varList)
  varNum = UBOUND(varList, 1)
  PRINT *, "num of analysis vars: ", varNum

  !CALL namelist_read(TRIM(configFile), "ensNum", ensNum)
  ifile = yaml_get_var(TRIM(configFile), 'BMat', 'ensNum', ensNum)
  PRINT *, "num of ensembles: ", ensNum

  ! Initialize
  CALL mpddGlob%initialize()   ! Initialize the mpdd
  CALL geometry%initialize(configFile, mpddGlob)   ! Initialize the geometry

  ALLOCATE (Xm(ensNum))

  ASSOCIATE (sg => geometry%mg%sg(geometry%mg%mg_finest))

    ! Initialize the Xm
    CALL Xm(1)%Initialize(configFile, sg)
    IF (ensNum > 1) THEN
      DO iens = 2, ensNum
        Xm(iens) = Xm(1)
      END DO
    END IF

    ! Initialize the ioWRF
    ALLOCATE (ioWRF)
    CALL ioWRF%initialize(configFile, geometry)

    DO iens = 1, ensNum
      ! Read the wrfout into the Xm
      CALL ioWRF%m_read_bcg_into_Xm_Ens(Xm(iens), sg, iens)

      ! Output Xm to NC files
      WRITE (memid, '(A1,I2.2)') "m", iens
      DO iv = 1, varNum
        varName = varList(iv)
        CALL Output_NC_State_SV(Xm(iens), TRIM(outputDir), "testEnLoc_Ens_bak_"//TRIM(varName)//"_"//memid, TRIM(varName))
      END DO
    END DO
    DEALLOCATE (ioWRF)

  END ASSOCIATE

  IF (mpddGlob%isBaseProc()) THEN
    CALL EnLoc%b_getEnsdata(configFile, Xm)
    PRINT *, 'Done EnLoc%b_getEnsdata ...'

    CALL EnLoc%b_calEigen()
    PRINT *, 'Done EnLoc%b_calEigen ...'

    CALL EnLoc%b_jacoEigenT()
    PRINT *, 'Done EnLoc%b_jacoEigenT ...'

    CALL EnLoc%b_destroy()
  END IF

  DEALLOCATE (Xm)

  ! Destroy the structures
  CALL mpddGlob%finalize

END PROGRAM Test_EnLoc_state
