PROGRAM Test_Raw_fy4_mwts2

  USE kinds_m, ONLY: i_kind, r_kind
  USE ObsMWTS2_m, ONLY: ObsMWTS2_t
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
  TYPE(ObsMWTS2_t) :: OBS_MWTS2
  TYPE(State_t) :: X
  TYPE(ObsSet_t) :: HX
  !TYPE(rttov_nl_t) :: rttov_nl
  TYPE(mpddGlob_t), TARGET :: mpddGlob
  TYPE(geometry_t), TARGET :: geometry
  TYPE(MPObs_t), TARGET :: mpObs

  ! ============ Call RTTOV for H(x) ============
  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", configFile)
  
  ! Get the configuration file
  CALL getarg(1, configFile)
  IF (TRIM(configFile) .EQ. '') THEN
    CALL getarg(0, configFile)
    configFile = TRIM(configFile)//"/UnitTest/Test_Raw_fy3_mwts2_mwhs2.yaml"
  END IF

  PRINT *, 'Test_Sat_fy4_mwts2 Config: ', TRIM(configFile)

  ! Initialize the mpdd
  CALL mpddGlob%initialize()

  ! Initialize geometry:
  CALL geometry%initialize(configFile, mpddGlob)

  ASSOCIATE (sg => geometry%mg%sg(5))
    CALL X%initialize(configFile, sg)
    CALL mpObs%initializeMPObs(sg)

    PRINT *, 'ObsMWTS2%ObsInitial'
    CALL OBS_MWTS2%ObsInitial(configFile)
    PRINT *, '------------ObsMWTS2%ObsInitial successfully run------------'

    CALL mpddGlob%barrier
    PRINT *, 'ObsMWTS2%ObsIngest'
    CALL OBS_MWTS2%ObsIngest(X)
    PRINT *, '------------ObsMWTS2%ObsIngest successfully run------------'

  END ASSOCIATE

  CALL mpddGlob%barrier

  ! Finalize
  CALL mpddGlob%finalize

  ! For ctest

  IF (X%sg%isBaseProc()) THEN

    PRINT *, MAXVAL(OBS_MWTS2%RadBase(6)%obsData(:, 1))
    PRINT *, MINVAL(OBS_MWTS2%RadBase(6)%obsData(:, 1))

    IF (MAXVAL((OBS_MWTS2%RadBase(6)%obsData(:, 1))) - MINVAL((OBS_MWTS2%RadBase(6)%obsData(:, 1))) > 5) THEN
      PRINT *, 'Test passed!'
    ELSE
      PRINT *, 'Test failed!'
    END IF
  END IF

END PROGRAM Test_Raw_fy4_mwts2
