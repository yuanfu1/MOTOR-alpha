!----------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Centers
!                     for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Jiongming Pang
! VERSION           : V 0.0
! HISTORY           : 2023-07-12, created by Jiongming Pang for testing IOGrapesModelvar_Ens.
!
!   Created by Jiongming Pang (pang.j.m@hotmail.com), 2023-07-12, @GBA-MWF, Shenzhen
!
!----------------------------------------------------------------------------------------

PROGRAM testIOGrapesModelvar_Ens
  USE kinds_m
  USE geometry_m, ONLY: geometry_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE State_m, ONLY: State_t
  USE IOGrapesModelvar_m, ONLY: IOGrapesModelvar_t
  USE InterpolateTime_m
  USE State2NC_m
  !USE NMLRead_m
  USE YAMLRead_m

  ! Define types
  TYPE(mpddGlob_t), TARGET :: mpddGlob
  TYPE(geometry_t), TARGET :: geometry
  TYPE(IOGrapesModelvar_t), ALLOCATABLE, TARGET :: ioGrapesModelvar
  TYPE(State_t), ALLOCATABLE  :: Xm(:)

  CHARACTER(LEN=1024) :: configFile, bakFilesDir, outputDir

  CHARACTER(LEN=20), ALLOCATABLE :: varList(:)
  CHARACTER(LEN=20) :: varName
  CHARACTER(LEN=3)  :: memid
  INTEGER(i_kind)   :: varNum, iv, ifile, ensNum, iens

  ! Get the configFile
  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", configFile)
  configFile = TRIM(configFile)//"/UnitTest/testIOGrapesModelvar_Ens.yaml"

  ifile = yaml_get_var(configFile, 'IO', 'output_dir', outputDir)

  !CALL namelist_read(TRIM(configFile), "varList", varList)
  ifile = yaml_get_var(TRIM(configFile), 'modelState', 'varList', varList)
  varNum = UBOUND(varList, 1)
  PRINT *, "num of analysis vars: ", varNum

  ifile = yaml_get_var(TRIM(configFile), 'BMat', 'ensNum', ensNum)
  PRINT *, "num of ensembles: ", ensNum

  ! Initialize
  CALL mpddGlob%initialize()   ! Initialize the mpdd
  CALL geometry%initialize(configFile, mpddGlob)   ! Initialize the geometry

  ASSOCIATE (sg => geometry%mg%sg(geometry%mg%mg_finest))

    ALLOCATE (Xm(ensNum))

    DO iens = 1, ensNum

      ALLOCATE (ioGrapesModelvar)
      ! Initialize the ioGrapesModelvar
      CALL ioGrapesModelvar%initialize(configFile, geometry)

      ! Initialize the Xm
      CALL Xm(iens)%initialize(configFile, sg)

      ! Read the modelvar into the Xm
      CALL ioGrapesModelvar%m_read_bcg_into_Xm_Ens(Xm(iens), sg, iens)

      ! Output Xm to NC files
      WRITE (memid, '(A1,I2.2)') "m", iens
      DO iv = 1, varNum
        varName = TRIM(varList(iv))
        CALL Output_NC_State_SV(Xm(iens), TRIM(outputDir), "testIOGrapesModelvar_MultiDA_bak_"//memid, TRIM(varName))
      END DO

      DEALLOCATE (ioGrapesModelvar)

    END DO

  END ASSOCIATE

  ! Destroy the structures
  CALL mpddGlob%finalize
  ! DEALLOCATE(bakFilesNames, varList)

END PROGRAM testIOGrapesModelvar_Ens
