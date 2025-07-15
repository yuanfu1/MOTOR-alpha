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
PROGRAM Test_RefThinningOutput
  USE kinds_m, ONLY: i_kind, r_kind, r_double
  USE RadarRAW_m, ONLY: RadarRAW_t
  USE ObsRadarRef_m, ONLY: ObsRadarRef_t
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
  TYPE(ObsRadarRef_t) :: ObsRadarRef
  INTEGER(i_kind) :: istatus

  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", staticDir)
  ! configFile = TRIM(staticDir)//"/UnitTest"//"/Test_RefThinningOutput.yaml"
  configFile = '/public/home/simi/optest/3DAnalysis.alpha/task/221103/221103_2200/App_3DAnalysis.yaml'

  istatus = yaml_get_var(configFile, 'IO', 'output_dir', outputDir)
  filenameOutput = TRIM(outputDir)//"/testRadar.nc"

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

      ! Initial
      CALL ObsRadarRef%ObsInitial(configFile)
      CALL ObsRadarRef%ObsIngest(X)

      X%Fields(1)%DATA = 0

      CALL ObsRadarRef%ObsPrepareForSg(X)
      !CALL ObsRadarRef%ObsSuper(X, Y, mpObs)
      CALL ObsRadarRef%ObsThinning(X, Y, mpObs, .TRUE., .FALSE.)

      DD = Obs2State_BaseTypeName(sg, Y)

      CALL Output_NC_State_SV(DD, outputDir, "Test_RefThinningOutput_STMAS_"//TRIM(DD%Fields(1)%Get_Name()), TRIM(DD%Fields(1)%Get_Name()))
    END BLOCK
  END ASSOCIATE

  ! Destory
  CALL ObsRadarRef%ObsDestroy()
  PRINT *, 'It is done'

  CALL mpddGlob%finalize
END PROGRAM Test_RefThinningOutput
