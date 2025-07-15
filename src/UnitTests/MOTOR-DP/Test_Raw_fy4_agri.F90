PROGRAM Test_Raw_fy4_agri

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
  configFile = TRIM(configFile)//"/UnitTest/Test_Raw_fy4_agri.yaml"
  PRINT *, 'Test_Sat_fy4_agri Config: ', TRIM(configFile)

  ! Initialize the mpdd
  CALL mpddGlob%initialize()

  ! Initialize geometry:
  CALL geometry%initialize(configFile, mpddGlob)

  ASSOCIATE (sg => geometry%mg%sg(5))
    CALL X%initialize(configFile, sg)
    CALL mpObs%initializeMPObs(sg)

    PRINT *, 'ObsAGRI%ObsInitial'
    CALL OBS_AGRI%ObsInitial(configFile)
    PRINT *, '------------ObsAGRI%ObsInitial successfully run------------'

    CALL mpddGlob%barrier
    PRINT *, 'ObsAGRI%ObsIngest'
    CALL OBS_AGRI%ObsIngest(X)
    PRINT *, '------------ObsAGRI%ObsIngest successfully run------------'

  END ASSOCIATE

  CALL mpddGlob%barrier

  ! Finalize
  CALL mpddGlob%finalize

  ! For ctest

  IF (X%sg%isBaseProc()) THEN

    PRINT *, MAXVAL(OBS_AGRI%RadBase(6)%obsData(:, 1))
    PRINT *, MINVAL(OBS_AGRI%RadBase(6)%obsData(:, 1))

    IF (MAXVAL((OBS_AGRI%RadBase(6)%obsData(:, 1))) - MINVAL((OBS_AGRI%RadBase(6)%obsData(:, 1))) > 5) THEN
      PRINT *, 'Test passed!'
    ELSE
      PRINT *, 'Test failed!'
    END IF
  END IF

END PROGRAM Test_Raw_fy4_agri
