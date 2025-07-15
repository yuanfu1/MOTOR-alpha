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
PROGRAM Test_RADPreProcessing
  USE kinds_m, ONLY: i_kind, r_kind, r_double
  USE RadarRAW_m, ONLY: RadarRAW_t
  USE ObsRadarVel_m, ONLY: ObsRadarVel_t
  USE geometry_m, ONLY: geometry_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE MPObs_m, ONLY: MPObs_t
  USE State_m, ONLY: State_t
  USE ObsSet_m, ONLY: ObsSet_t
  USE C2O_m, ONLY: C2O_t
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
  INTEGER(i_kind) :: istatus

  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", staticDir)
  configFile = TRIM(staticDir)//"/UnitTest"//"/Test_VelThinningOutputMore.yaml"
  configFile = "/Users/qzl/sources/MOTOR/input/220527_0000_3km/App_3DVarVerification.yaml"

  !   IF(yaml_get_var(configFile, 'IO', 'output_dir', outputDir) /= 0) STOP
  istatus = yaml_get_var(configFile, 'IO', 'output_dir', outputDir)

  ! Initializer
  ! Auxtypes
  CALL mpddGlob%initialize()                                   ! Initialize the mpdd

  ! Initialize geometry
  CALL geometry%initialize(configFile, mpddGlob)               ! Initialize the geometry

  ASSOCIATE (sg => geometry%mg%sg(7))
    BLOCK

      TYPE(MPObs_t), TARGET :: mpObs
      TYPE(State_t) :: X, DD
      TYPE(ObsSet_t) :: Y
      INTEGER(i_kind) :: i, j, k

      CALL mpObs%initializeMPObs(sg)
      CALL X%initialize(configFile, sg)
      Y = ObsSet_t(configFile, mpObs)

      X%Fields(1)%DATA = 0

      ! Initial
      CALL ObsRadarVel%ObsInitial(configFile)
      CALL ObsRadarVel%ObsIngest(X)
      ! CALL ObsRadarVel%ObsPrepareForSg(X)
      ! CALL ObsRadarVel%ObsThinning(X, Y, mpObs, .TRUE., .FALSE.)

      ! ! PRINT *, 'Y%ObsFields(1)%Get_Name(): ', Y%ObsFields(1)%Get_Name()
      ! DD = Obs2State_BaseTypeName(sg, Y)

      ! ! PRINT *, 'DD%Fields(1)%Get_Name(): ', DD%Fields(1)%Get_Name()
      ! CALL Output_NC_State_AV(DD, outputDir, "TestRadarVelThinning_rwnd")
      DO i = 1, UBOUND(ObsRadarVel%radars, 1)
        CALL ObsRadarVel%radars(i)%dumpVADToFile()
      END DO
      ! CALL ObsRadarVel%radars(1)%dumpVADToFile()
    END BLOCK
  END ASSOCIATE

  ! Destory
  CALL ObsRadarVel%ObsDestroy()
  PRINT *, 'It is done'

  CALL mpddGlob%finalize
END PROGRAM Test_RADPreProcessing
