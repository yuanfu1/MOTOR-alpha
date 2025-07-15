!----------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Centers
!                     for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Jiongming Pang
! VERSION           : V 0.0
! HISTORY           : 2021-09-26, created by Jiongming Pang for testing IOWRF.
!                     2022-01-27, modified by JMPang for testing IOWRF with MultiDA.
!
!   Created by Jiongming Pang (pang.j.m@hotmail.com), 2021-9-26, @GBA-MWF, Shenzhen
!   Modified by Jiongming Pang, 2022-01-27, @GBA-MWF, Shenzhen
!
!----------------------------------------------------------------------------------------

PROGRAM testIOWRF
  USE kinds_m
  USE geometry_m, ONLY: geometry_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE State_m, ONLY: State_t
  USE IOWRF_m, ONLY: IOWRF_t
  USE InterpolateTime_m
  USE State2NC_m
  !USE NMLRead_m
  USE YAMLRead_m

  ! Define types
  TYPE(mpddGlob_t), TARGET :: mpddGlob
  TYPE(geometry_t), TARGET :: geometry
  TYPE(IOWRF_t), TARGET :: ioWRF
  TYPE(State_t)  :: Xm

  CHARACTER(LEN=1024) :: configFile, bakFilesDir, outputDir

  CHARACTER(LEN=20), ALLOCATABLE :: varList(:)
  CHARACTER(LEN=20) :: varName
  CHARACTER(LEN=3)  :: memid
  INTEGER(i_kind)   :: varNum, iv, ifile, istatus

  ! Get the configFile
  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", configFile)
  configFile = TRIM(configFile)//"/UnitTest/testIOWRF.yaml"

  istatus = yaml_get_var(configFile, 'IO', 'output_dir', outputDir)

  !CALL namelist_read(TRIM(configFile), "varList", varList)
  ifile = yaml_get_var(TRIM(configFile), 'modelState', 'varList', varList)
  varNum = UBOUND(varList, 1)
  PRINT *, "num of analysis vars: ", varNum

  ! Initialize
  CALL mpddGlob%initialize()   ! Initialize the mpdd
  CALL geometry%initialize(configFile, mpddGlob)   ! Initialize the geometry

  ! Initialize the ioWRF
  CALL ioWRF%initialize(configFile, geometry)

  ASSOCIATE (sg => geometry%mg%sg(geometry%mg%mg_finest))

    ! Initialize the Xm
    CALL Xm%initialize(configFile, sg)

    ! Read the wrfout into the Xm
    ! CALL ioWRF%m_read_bcg_into_Xm(Xm, sg)
    CALL ioWRF%m_read_bcg_into_Xm_Ens(Xm, sg)

    ! Output Xm to NC files
    ! CALL Output_NC_State_SV(Xm, outputDir, 'testIOWRF_t20220117', "temp")
    DO iv = 1, varNum
      varName = varList(iv)
      CALL Output_NC_State_SV(Xm, TRIM(outputDir), "testIOWRF_MultiDA_bak_"//TRIM(varName), TRIM(varName))
    END DO

  END ASSOCIATE

  ! Destroy the structures
  CALL mpddGlob%finalize
  ! DEALLOCATE(bakFilesNames, varList)

END PROGRAM testIOWRF
