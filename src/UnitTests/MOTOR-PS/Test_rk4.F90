PROGRAM test_rk4
  USE kinds_m, ONLY: i_kind, r_kind
  USE rk4test_m, ONLY: rk4_t
  USE State_m, ONLY: State_t
  USE UnitTestDyn_m, ONLY: UnitTestDyn_t
  USE GenContainers_m, ONLY: GenContainers_t
  USE geometry_m, ONLY: Geometry_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE MGOpts_m
  USE YAMLRead_m
  USE parameters_m, ONLY: Omega
  USE State2NC_m

  IMPLICIT NONE

  TYPE(rk4_t) :: rk4
  TYPE(UnitTestDyn_t) :: Rossby
  TYPE(State_t) :: X1half, X2half, X3half, &
                   X1full, X2full, X3full, &
                   X1thalf, X1tFull
  TYPE(GenContainers_t) :: GenContainers
  TYPE(mpddGlob_t), TARGET :: mpddGlob
  TYPE(SingleGrid_t), TARGET :: sg_t, sg_t2, sg_t3
  TYPE(geometry_t) :: geometry1Full, geometry1Half, geometry2Full, geometry2Half, geometry3Full, geometry3Half

  INTEGER(i_kind) :: i, j, k, vLevel1, vLevel2, ifile, test_glevel
  REAL(r_kind) :: error_1, error_2, error_3
  REAL(r_kind), ALLOCATABLE :: sigma(:), z(:, :), z_s(:), Hz(:, :), dtt(:)

  CHARACTER(LEN=1024) :: configFile, outputDir

  CHARACTER(7) :: varNamesHalf(7) = ['psi    ', 'chi    ', 'vor    ', &
                                     'div    ', 'rho    ', 'qvapor ', 'pres   '], &
                  varNamesHalft(7) = ['psit   ', 'chit   ', 'vort   ', &
                                      'divt   ', 'rhot   ', 'qvaport', 'prest  '], &
                  varNamesFull(2) = ['w      ', 'theta  '], &
                  varNamesFullt(2) = ['wt     ', 'thetat ']

  ! ! Fields to test:
  ! REAL(r_kind), ALLOCATABLE :: tendiv1(:, :), tendiv2(:, :), tendiv3(:, :), &
  !                              tendivt1(:, :), tendivt2(:, :), tendivt3(:, :), &
  !                              vor(:, :), psi(:, :), &
  !                              div(:, :), chi(:, :), &
  !                              wwnd(:, :), rho(:, :), P(:, :), &
  !                              sigma_3d(:, :), sigma(:), sigma2(:), &
  !                              zhght(:, :), z_s(:), Hz(:, :), &
  !                              DM(:, :)

  ! CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", configFile)
  configFile = "/Users/jilongchen/research/MOTOR-DA/MOTOR/static/gzm_and_laplace_terrain.yaml"

  ifile = yaml_get_var(TRIM(configFile), 'IO', 'output_dir', outputDir)
  ! CALL GET_ENVIRONMENT_VARIABLE("OUTPUT_DIR", outputDir)

  test_glevel = 5

  GenContainers = GenContainers_t(TRIM(configFile))
  CALL GenContainers%GenGeometry(geometry1Full)

  vLevel1 = geometry1Full%mg%sg(test_glevel)%vLevel - 1

  ALLOCATE (sigma(vLevel1))

  sigma(1:vLevel1) = (geometry1Full%mg%sg(test_glevel)%sigma(1:vLevel1) &
                      + geometry1Full%mg%sg(test_glevel)%sigma(2:vLevel1 + 1)) / 2.0D0

  CALL GenContainers%GenGeometry(geometry1Half, sigma)
  ! CALL GenContainers%GenGeometry(geometry1Half)
  DEALLOCATE (sigma)

  ASSOCIATE (sg => geometry1Full%mg%sg(test_glevel))
    ALLOCATE (z(sg%vLevel, sg%num_cell), &
              z_s(sg%num_cell), &
              Hz(sg%vLevel, sg%num_cell))
    PRINT *, 'latlon is: ', sg%cell_cntr(:, 1:10)
    CALL Rossby%initial(sg)
    CALL Rossby%Hght_func(geometry1Full%mg%sg(test_glevel)%cell_cntr, z, z_s, Hz)
  END ASSOCIATE

  geometry1Full%mg%sg(test_glevel)%topo = z_s
  geometry1Half%mg%sg(test_glevel)%topo = z_s
  DEALLOCATE (z, z_s, Hz)
  CALL Rossby%destroy()
  PRINT *, 'Debug 3'
  CALL GenContainers%updateSg(geometry1Full%mg%sg(test_glevel), geometry1Half%mg%sg(test_glevel))
  CALL GenContainers%updateSg(geometry1Half%mg%sg(test_glevel), geometry1Full%mg%sg(test_glevel))
  PRINT *, 'Debug 4'
  rk4 = rk4_t(configfile, test_glevel, geometry1Half, geometry1Full)
  PRINT *, 'Debug 5'
  PRINT *, 'size of f in geometry1Half is: ', SIZE(geometry1Half%mg%sg(test_glevel)%f)
  PRINT *, 'size of f in rk4%geometryHalf is: ', SIZE(rk4%geometryHalf%mg%sg(test_glevel)%f)
  CALL GenContainers%GenStates(rk4%geometryHalf%mg%sg(test_glevel), X1Half)
  CALL GenContainers%GenStates(rk4%geometryFull%mg%sg(test_glevel), X1Full)
  CALL GenContainers%GenStates(rk4%geometryHalf%mg%sg(test_glevel), X1tHalf)
  CALL GenContainers%GenStates(rk4%geometryFull%mg%sg(test_glevel), X1tFull)
  DO i = LBOUND(X1Half%fields, 1), UBOUND(X1Half%fields, 1)
    PRINT *, 'varname in X1Half is: ', TRIM(X1Half%fields(i)%Get_Name())
  END DO

  ASSOCIATE (t_in => geometry1Full%mg%sg(test_glevel)%tt)
    ALLOCATE (dtt(SIZE(t_in)))
    PRINT *, 'dtt equals (in Test_rk4): ', SIZE(dtt)
    dtt = t_in - t_in(1)
    PRINT *, 'dtt equals (in Test_rk4): ', dtt
    DO i = 1, SIZE(dtt, 1)
      CALL rk4%RossbyFull%righthandsFull(0.0D0, X1tfull)
      CALL rk4%RossbyHalf%righthandsHalf(0.0D0, X1tHalf)

      DO j = 1, SIZE(varNamesHalf)
        X1half%Fields(X1half%getVarIdx(TRIM(varNamesHalft(j))))%DATA(:, :, i) = &
          X1thalf%Fields(X1thalf%getVarIdx(TRIM(varNamesHalft(j))))%DATA(:, :, 1)
      END DO

      DO j = 1, SIZE(varNamesFull)
        X1Full%Fields(X1Full%getVarIdx(TRIM(varNamesFullt(j))))%DATA(:, :, i) = &
          X1tfull%Fields(X1tFull%getVarIdx(TRIM(varNamesFullt(j))))%DATA(:, :, 1)
      END DO

    END DO
  END ASSOCIATE

  DO i = 1, SIZE(varNamesHalf)
    X1half%Fields(X1half%getVarIdx(TRIM(varNamesHalf(i))))%DATA = X1half%Fields(X1half%getVarIdx(TRIM(varNamesHalft(i))))%DATA
  END DO

  DO i = 1, SIZE(varNamesFull)
    X1Full%Fields(X1Full%getVarIdx(TRIM(varNamesFull(i))))%DATA = X1full%Fields(X1Full%getVarIdx(TRIM(varNamesFullt(i))))%DATA
  END DO
  CALL GenContainers%GenBdy%GenBdyTendency(X1Full, X1Half)

  CALL rk4%timeIntegrals(X1Half, X1Full)

  CALL Output_NC_State_AV(X1Full, outputDir, "Test_rk4_Full_d1")
  CALL Output_NC_State_AV(X1Half, outputDir, "Test_rk4_Half_d1")

END PROGRAM
