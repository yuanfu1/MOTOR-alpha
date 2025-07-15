PROGRAM Test_RawGNSS

  USE kinds_m, ONLY: i_kind, r_kind
  USE ObsGNSS_m, ONLY: ObsGNSS_t
  USE State_m, ONLY: State_t
  USE geometry_m, ONLY: geometry_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE MPObs_m, ONLY: MPObs_t
  USE ObsSet_m, ONLY: ObsSet_t
  USE ObsField_m, ONLY: ObsField_t
  !USE NMLRead_m
  USE YAMLRead_m
  USE parameters_m

  IMPLICIT NONE

  CHARACTER(LEN=1024) :: configFile
  TYPE(ObsGNSS_t) :: OBS_GNSS
  TYPE(State_t) :: X
  TYPE(ObsSet_t) :: HX
  !TYPE(rttov_nl_t) :: rttov_nl
  TYPE(mpddGlob_t), TARGET :: mpddGlob
  TYPE(geometry_t), TARGET :: geometry
  TYPE(MPObs_t), TARGET :: mpObs

  configFile = "App_3DVarVerification_gnss.yaml"

  PRINT *, '------ RUN TEST RAW GNSS CASE ------ with config file:', configFile

  ! Initialize the mpdd
  CALL mpddGlob%initialize()

  ! Initialize geometry:
  CALL geometry%initialize(configFile, mpddGlob)

  ASSOCIATE (sg => geometry%mg%sg(5))
    CALL X%initialize(configFile, sg)
    CALL mpObs%initializeMPObs(sg)

    CALL OBS_GNSS%ObsInitial(configFile)

    CALL OBS_GNSS%ObsIngest(X)

  END ASSOCIATE

  CALL mpddGlob%barrier

  ! Finalize
  CALL mpddGlob%finalize

  PRINT *, 'JUST FOR TEST!'
END PROGRAM Test_RawGNSS
