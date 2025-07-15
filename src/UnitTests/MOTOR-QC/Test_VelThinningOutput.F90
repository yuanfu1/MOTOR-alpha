!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2022/3/18, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!!
PROGRAM Test_VelThinningOutput
  USE kinds_m, ONLY: i_kind, r_kind, r_double
  USE RadarRAW_m, ONLY: RadarRAW_t
  USE ObsRadarVel_m, ONLY: ObsRadarVel_t
  USE geometry_m, ONLY: geometry_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE MPObs_m, ONLY: MPObs_t
  USE State_m, ONLY: State_t
  USE ObsSet_m, ONLY: ObsSet_t
  USE C2O_m, ONLY: C2O_t
  USE IOGrapes_m, ONLY: IOGrapes_t
  USE IOERA5_m, ONLY: IOERA5_t
  USE YAMLRead_m
  USE State2NC_m
  USE Obs2State_m

  IMPLICIT NONE

  ! Define types
  TYPE(mpddGlob_t), TARGET :: mpddGlob
  TYPE(geometry_t), TARGET :: geometry

  CHARACTER(LEN=1024) :: inputDir
  CHARACTER(LEN=1024) :: outputDir
  CHARACTER(len=1024) :: filenameInput, filenameOutput
  CHARACTER(len=1024) :: configFile
  CHARACTER(LEN=1024) :: staticDir
  TYPE(ObsRadarVel_t) ::  ObsRadarVel

  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", staticDir)
  configFile = TRIM(staticDir)//"/UnitTest"//"/Test_VelThinningOutputMore.yaml"
  configFile = "/Users/qzl/sources/MOTOR/input/2024110300_T10p1/App_3DVarVerification.yaml"

  IF (yaml_get_var(configFile, 'IO', 'output_dir', outputDir) .NE. 0) THEN
    PRINT *, 'Error: output_dir not found in the configuration file'
    STOP
  END IF

  ! Initializer
  ! Auxtypes
  CALL mpddGlob%initialize()                                   ! Initialize the mpdd

  ! Initialize geometry
  CALL geometry%initialize(configFile, mpddGlob)               ! Initialize the geometry

  ASSOCIATE (sg => geometry%mg%sg(8))
    BLOCK

      TYPE(MPObs_t), TARGET :: mpObs
      TYPE(State_t) :: X, DD
      TYPE(ObsSet_t) :: Y
      TYPE(C2O_t) :: H
      INTEGER(i_kind) :: i, j, k
      TYPE(IOGrapes_t), TARGET :: ioGrapes
      TYPE(IOERA5_t), TARGET :: ioERA5

      CALL mpObs%initializeMPObs(sg)
      CALL X%initialize(configFile, sg)
      Y = ObsSet_t(configFile, mpObs)

      X%Fields(1)%DATA = 0
      ! CALL IOGrapes%initialize(configFile, geometry)
      ! CALL IOGrapes%m_read_bcg_into_Xm(X, sg)
      CALL ioERA5%initialize(configFile, geometry)
      CALL ioERA5%m_read_bcg_into_Xm(X, sg)
      ! CALL Output_NC_State_AVST(X, outputDir, "TestRadarThinning_X_ERA5", sg%tSlots)
      ! STOP

      ! Initial
      CALL ObsRadarVel%ObsInitial(configFile)
      CALL ObsRadarVel%ObsIngest(X)
      CALL ObsRadarVel%ObsPrepareForSg(X)
      CALL ObsRadarVel%ObsThinning(X, Y, mpObs, .TRUE., .FALSE.)
      PRINT *, 'Here'

      CALL H%initialize(configFile, X, Y, X)

      ! PRINT *, 'Y%ObsFields(1)%Get_Name(): ', Y%ObsFields(1)%Get_Name()
      DD = Obs2State_BaseTypeName(sg, Y)

      ! PRINT *, 'DD%Fields(1)%Get_Name(): ', DD%Fields(1)%Get_Name()
      CALL Output_NC_State_AVST(DD, outputDir, "TestRadarVelThinning_rwnd", sg%tSlots)

      ! DD = Obs2State_BaseTypeName(sg, Y - H%transFwdNonLinear_opr(X))

      ! CALL Output_NC_State_AVST(DD, outputDir, "TestRadarVel_Inno_CMAGFS", sg%tSlots)

    END BLOCK
  END ASSOCIATE

  ! Destory
  CALL ObsRadarVel%ObsDestroy()
  PRINT *, 'It is done'

  CALL mpddGlob%finalize
END PROGRAM Test_VelThinningOutput
