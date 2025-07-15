PROGRAM test_WarmBubble
  USE kinds_m, ONLY: i_kind, r_kind
  USE rk4_m, ONLY: rk4_t
  USE State_m, ONLY: State_t
  USE WarmBubble_m
  USE GenContainers_m, ONLY: GenContainers_t
  USE geometry_m, ONLY: Geometry_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE InterpValue_m, ONLY: InterpValue_t
  USE MGOpts_m
  USE YAMLRead_m
  USE parameters_m, ONLY: Omega
  USE State2NC_m

  IMPLICIT NONE
  TYPE(rk4_t) :: rk4
  TYPE(State_t) :: XHalf, XFull
  TYPE(GenContainers_t) :: GenContainers
  TYPE(geometry_t) :: geometryFull, geometryHalf
  TYPE(InterpValue_t) :: InterpValue

  INTEGER(i_kind) :: i, j, k, kt, vLevel1, vLevel2, ifile, test_glevel
  REAL(r_kind), ALLOCATABLE :: sigma(:)

  CHARACTER(LEN=1024) :: configFile, outputDir

  CHARACTER(7) :: varNamesHalf(7) = ['psi    ', 'chi    ', 'vor    ', &
                                     'div    ', 'rho    ', 'qvapor ', 'pres   '], &
                  varNamesHalft(7) = ['psit   ', 'chit   ', 'vort   ', &
                                      'divt   ', 'rhot   ', 'qvaport', 'prest  '], &
                  varNamesFull(2) = ['w      ', 'theta  '], &
                  varNamesFullt(2) = ['wt     ', 'thetat ']
  configFile = "/home/jilongchen/research/MOTOR-DA/ZGridModel_new/MOTOR/static/UnitTest/warm_bubble_test.yaml"

  ifile = yaml_get_var(TRIM(configFile), 'IO', 'output_dir', outputDir)

  test_glevel = 7

  GenContainers = GenContainers_t(TRIM(configFile))
  CALL GenContainers%GenGeometry(geometryFull)

  vLevel1 = geometryFull%mg%sg(test_glevel)%vLevel - 1

  ALLOCATE (sigma(vLevel1))

  sigma(1:vLevel1) = (geometryFull%mg%sg(test_glevel)%sigma(1:vLevel1) &
                      + geometryFull%mg%sg(test_glevel)%sigma(2:vLevel1 + 1)) / 2.0D0

  CALL GenContainers%GenGeometry(geometryHalf, sigma)
  ! CALL GenContainers%GenGeometry(geometry1Half)
  DEALLOCATE (sigma)

  ASSOCIATE (sgHalf => geometryHalf%mg%sg(test_glevel), &
             sgFull => geometryFull%mg%sg(test_glevel))
    sgFull%topo = 0.0D0
    sgHalf%topo = 0.0D0
    CALL GenContainers%updateSg(sgFull, sgHalf)
    CALL GenContainers%updateSg(sgHalf, sgFull)
    rk4 = rk4_t(configfile, test_glevel, geometryHalf, geometryFull)

    CALL GenContainers%GenStates(sgHalf, XHalf)
    CALL GenContainers%GenStates(sgFull, XFull)
    InterpValue = InterpValue_t(sgHalf)
    DO kt = 1, SIZE(XFull%Fields(XFull%getVarIdx('theta'))%DATA, dim=3)
      CALL GenWarmBubble(sgFull%cell_cntr(1, :), sgFull%cell_cntr(2, :), sgFull%zHght, &
                         XFull%Fields(XFull%getVarIdx('theta'))%DATA(:, :, kt), &
                         XFull%Fields(XFull%getVarIdx('rho'))%DATA(:, :, kt), &
                         XFull%Fields(XFull%getVarIdx('pres'))%DATA(:, :, kt))
      PRINT *, 'max of pres is: ', MAXVAL(XFull%Fields(XFull%getVarIdx('pres'))%DATA(:, :, kt))

      CALL InterpValue%InterpValue(XFull%Fields(XFull%getVarIdx('rho'))%DATA(:, :, kt), &
                                   XHalf%Fields(XHalf%getVarIdx('rho'))%DATA(:, :, kt))

      CALL InterpValue%InterpValue(XFull%Fields(XFull%getVarIdx('pres'))%DATA(:, :, kt), &
                                   XHalf%Fields(XHalf%getVarIdx('pres'))%DATA(:, :, kt))
    END DO

    CALL GenContainers%GenBdy%GenBdyTendency(XFull, XHalf)

    CALL Output_NC_State_AV(XFull, outputDir, "Test_warmbubble_Full_initial")
    CALL Output_NC_State_AV(XHalf, outputDir, "Test_warmbubble_Half_initial")

    CALL rk4%timeIntegrals(XHalf, XFull)

    CALL Output_NC_State_AV(XFull, outputDir, "Test_warmbubble_Full_after_720s")
    CALL Output_NC_State_AV(XHalf, outputDir, "Test_warmbubble_Half_after_720s")

  END ASSOCIATE

END PROGRAM
