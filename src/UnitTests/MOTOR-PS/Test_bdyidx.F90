PROGRAM test_bdyidx
  USE kinds_m, ONLY: i_kind, r_kind
  USE State_m, ONLY: State_t
  USE GenContainers_m, ONLY: GenContainers_t
  USE geometry_m, ONLY: Geometry_t
  USE YAMLRead_m
  USE State2NC_m
  USE GenBdy_m, ONLY: GenBdy_t

  IMPLICIT NONE

  TYPE(State_t) :: X1, X2
  TYPE(GenContainers_t) :: GenContainers
  TYPE(geometry_t) :: geometry
  TYPE(GenBdy_t) :: GenBdy
  CHARACTER(LEN=1024) :: configFile, outputDir
  INTEGER(i_kind) :: ifile

  PRINT *, 'test1'

  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", configFile)
  configFile = TRIM(configFile)//"/gzm_and_laplace_terrain.yaml"

  ifile = yaml_get_var(TRIM(configFile), 'IO', 'output_dir', outputDir)

  CALL GET_ENVIRONMENT_VARIABLE("OUTPUT_DIR", outputDir)

  PRINT *, 'test2'

  GenContainers = GenContainers_t(TRIM(configFile))
  PRINT *, 'test3'
  CALL GenContainers%GenGeometry(geometry)
  PRINT *, 'test4'
  ASSOCIATE (sg => geometry%mg%sg(5))
    GenBdy = GenBdy_t(TRIM(configFile))
    CALL GenBdy%GenBdyIdx(sg)
    CALL GenBdy%GenBdyBufferWeight(sg)
    CALL GenContainers%GenStates(sg, X1)
    CALL X1%addVar('bdyidx')
    CALL X1%addVar('alpha')
    CALL X1%addVar('beta')

    X1%Fields(X1%getVarIdx('bdyidx'))%DATA(1, :, 1) = sg%bdy_type
    X1%Fields(X1%getVarIdx('alpha'))%DATA(1, :, 1) = sg%BdyAlpha
    X1%Fields(X1%getVarIdx('beta'))%DATA(1, :, 1) = sg%BdyBeta
    CALL Output_NC_State_AV(X1, outputDir, "Test_bdyidx")
  END ASSOCIATE

  CALL GenContainers%destroy()

END PROGRAM
