!----------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Centers
!                     for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Jiongming Pang
! VERSION           : V 0.0
! HISTORY           : 2023-08-10, created by Jiongming Pang for testing IOGrapesPostvar.
!
!   Created by Jiongming Pang (pang.j.m@hotmail.com), 2023-08-10, @GBA-MWF, Shenzhen
!
!----------------------------------------------------------------------------------------

PROGRAM testIOGrapesPostvar
  USE kinds_m
  USE geometry_m, ONLY: geometry_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE State_m, ONLY: State_t
  USE IOGrapesPostvar_m, ONLY: IOGrapesPostvar_t
  USE InterpolateTime_m
  USE State2NC_m
  !USE NMLRead_m
  USE YAMLRead_m

  ! Define types
  TYPE(mpddGlob_t), TARGET :: mpddGlob
  TYPE(geometry_t), TARGET :: geometry
  TYPE(IOGrapesPostvar_t), TARGET :: ioGrapesPostvar
  TYPE(State_t)  :: Xm

  CHARACTER(LEN=1024) :: configFile, bakFilesDir, outputDir

  CHARACTER(LEN=20), ALLOCATABLE :: varList(:)
  CHARACTER(LEN=20) :: varName
  CHARACTER(LEN=3)  :: memid
  INTEGER(i_kind)   :: varNum, iv, ifile

  ! Get the configFile
  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", configFile)
  configFile = TRIM(configFile)//"/UnitTest/testIOGrapesPostvar.yaml"
  PRINT *, 'configFile: ', TRIM(configFile)

  ifile = yaml_get_var(TRIM(configFile), 'IO', 'output_dir', outputDir)
  PRINT *, 'outputDir: ', TRIM(outputDir)

  !CALL namelist_read(TRIM(configFile), "varList", varList)
  ifile = yaml_get_var(TRIM(configFile), 'modelState', 'varList', varList)
  varNum = UBOUND(varList, 1)
  PRINT *, "num of analysis vars: ", varNum

  ! Initialize
  CALL mpddGlob%initialize()   ! Initialize the mpdd
  CALL geometry%initialize(configFile, mpddGlob)   ! Initialize the geometry

  ! Initialize the ioGrapesPostvar
  CALL ioGrapesPostvar%initialize(configFile, geometry)

  ASSOCIATE (sg => geometry%mg%sg(geometry%mg%mg_finest))

    ! Initialize the Xm
    CALL Xm%initialize(configFile, sg)

    ! Read the postvar into the Xm
    CALL ioGrapesPostvar%m_read_bcg_into_Xm(Xm, sg)

    ! Output Xm to NC files
    DO iv = 1, varNum
      varName = varList(iv)
      CALL Output_NC_State_SV(Xm, TRIM(outputDir), "testIOGrapesPostvar_MultiDA_bak_"//TRIM(varName), TRIM(varName))
    END DO

  END ASSOCIATE

  ! Destroy the structures
  CALL mpddGlob%finalize
  ! DEALLOCATE(bakFilesNames, varList)

END PROGRAM testIOGrapesPostvar
