PROGRAM Test_Raw_fy4_giirs

  USE kinds_m, ONLY: i_kind, r_kind
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
  TYPE(State_t) :: X
  TYPE(ObsSet_t) :: HX
  !TYPE(rttov_nl_t) :: rttov_nl
  TYPE(mpddGlob_t), TARGET :: mpddGlob
  TYPE(geometry_t), TARGET :: geometry
  TYPE(MPObs_t), TARGET :: mpObs

  ! ============ Call RTTOV for H(x) ============
  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", configFile)
  configFile = TRIM(configFile)//"/UnitTest/Test_Raw_fy4_giirs.yaml"
  PRINT *, 'Test_Sat_fy4_giirs Config: ', TRIM(configFile)

  ! Initialize the mpdd
  CALL mpddGlob%initialize()

  ! Initialize geometry:
  CALL geometry%initialize(configFile, mpddGlob)

  ASSOCIATE (sg => geometry%mg%sg(5))
    CALL X%initialize(configFile, sg)
    CALL mpObs%initializeMPObs(sg)

    PRINT *, 'ObsGIIRS%ObsInitial'
    CALL OBS_GIIRS%ObsInitial(configFile)
    PRINT *, '------------ObsGIIRS%ObsInitial successfully run------------'

    CALL mpddGlob%barrier
    PRINT *, 'ObsGIIRS%ObsIngest'
    CALL OBS_GIIRS%ObsIngest(X)
    PRINT *, '------------ObsGIIRS%ObsIngest successfully run------------'

  END ASSOCIATE

  CALL mpddGlob%barrier

  ! Finalize
  CALL mpddGlob%finalize

  ! For ctest

  IF (X%sg%isBaseProc()) THEN

    PRINT *, MAXVAL(OBS_GIIRS%RadBase(6)%obsData(:, 1))
    PRINT *, MINVAL(OBS_GIIRS%RadBase(6)%obsData(:, 1))

    IF (MAXVAL((OBS_GIIRS%RadBase(6)%obsData(:, 1))) - MINVAL((OBS_GIIRS%RadBase(6)%obsData(:, 1))) > 5) THEN
      PRINT *, 'Test passed!'
    ELSE
      PRINT *, 'Test failed!'
    END IF
  END IF

END PROGRAM Test_Raw_fy4_giirs
