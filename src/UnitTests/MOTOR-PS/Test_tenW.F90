PROGRAM test_tenW
  USE kinds_m, ONLY: i_kind, r_kind
  USE TenW_m, ONLY: TenW_t
  USE PreCal_m, ONLY: PreCal_t
  USE State_m, ONLY: State_t
  USE UnitTestDyn_m, ONLY: UnitTestDyn_t
  USE GenContainers_m, ONLY: GenContainers_t
  USE geometry_m, ONLY: Geometry_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE State_m, ONLY: State_t
  !  USE SingleGrid_m, ONLY: SingleGrid_t
  USE MGOpts_m
  USE YAMLRead_m
  USE parameters_m, ONLY: Omega, g
  USE State2NC_m
  USE gzm_m, ONLY: gzm_t

  TYPE(TenW_t) :: TenW
  TYPE(State_t) :: X1, X2, X3
  TYPE(PreCal_t) :: PreCal
  TYPE(UnitTestDyn_t) :: Rossby
  TYPE(GenContainers_t) :: GenContainers
  TYPE(mpddGlob_t), TARGET :: mpddGlob
  !  TYPE(SingleGrid_t) :: sg
  TYPE(geometry_t) :: geometry1Full, geometry1Half, geometry2Full, geometry2Half, geometry3Full, geometry3Half
  TYPE(gzm_t) :: gzm

  INTEGER(i_kind) :: i, k, vLevel1, vLevel2, ifile
  REAL(r_kind) :: errorw_1, errorw_2, errorw_3

  CHARACTER(LEN=1024) :: configFile, outputDir

  ! Fields to test:
  REAL(r_kind), ALLOCATABLE :: tenW1(:, :), tenW2(:, :), tenW3(:, :), &
                               tenWt1(:, :), tenWt2(:, :), tenWt3(:, :), &
                               vor(:, :), psi(:, :), &
                               div(:, :), chi(:, :), &
                               wwnd(:, :), rhohalf(:, :), &
                               rho(:, :), Phalf(:, :), P(:, :), &
                               sigma_3d(:, :), sigma(:), sigma2(:), &
                               zhght(:, :), z_s(:), Hz(:, :), &
                               DM(:, :), errorl_1(:), errorl_2(:), errorl_3(:)

  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", configFile)
  configFile = TRIM(configFile)//"/gzm_and_laplace_terrain.yaml"
  PRINT *, 'configurefile is', configFile

  ! CALL GET_ENVIRONMENT_VARIABLE("OUTPUT_DIR", outputDir)
  ifile = yaml_get_var(TRIM(configFile), 'IO', 'output_dir', outputDir)
  PRINT *, 'teswt'

  GenContainers = GenContainers_t(TRIM(configFile))
  CALL GenContainers%GenGeometry(geometry1Full)

  vLevel1 = geometry1Full%mg%sg(5)%vLevel - 1

  ALLOCATE (sigma(vLevel1))

  sigma(1:vLevel1) = (geometry1Full%mg%sg(5)%sigma(1:vLevel1) &
                      + geometry1Full%mg%sg(5)%sigma(2:vLevel1 + 1)) / 2.0D0

  CALL GenContainers%GenGeometry(geometry1Half, sigma)

  DEALLOCATE (sigma)

  ASSOCIATE (sg => geometry1Full%mg%sg(5), &
             num_cell => geometry1Full%mg%sg(5)%num_cell, &
             vLevel => geometry1Full%mg%sg(5)%vLevel)

    ALLOCATE (vor(sg%vLevel, sg%num_cell), &
              psi(sg%vLevel, sg%num_cell), &
              div(sg%vLevel, sg%num_cell), &
              chi(sg%vLevel, sg%num_cell), &
              rho(sg%vLevel, sg%num_cell), &
              P(sg%vLevel, sg%num_cell), &
              sigma_3d(sg%vLevel, sg%num_cell), &
              zhght(sg%vLevel, sg%num_cell), &
              z_s(sg%num_cell), &
              Hz(sg%vLevel, sg%num_cell), &
              DM(sg%vLevel, sg%num_cell))

    DO i = 1, sg%num_cell
      sigma_3d(:, i) = sg%sigma
    END DO
    CALL Rossby%initial(sg)
    CALL Rossby%Hght_func(sg%cell_cntr, zhght, z_s, Hz)
    CALL Rossby%RH_psi_func(0.0D0, sg%cell_cntr, sigma_3d, psi, vor)
    CALL Rossby%RH_chi_func(0.0D0, sg%cell_cntr, sigma_3d, chi, div)
    CALL Rossby%DM_func(Hz, DM)
    CALL Rossby%rho_func(0.0D0, sg%cell_cntr, zhght, z_s, rho)
    CALL Rossby%pres_func(0.0D0, sg%cell_cntr, zhght, z_s, P)
    sg%topo = z_s

    CALL GenContainers%updateSg(sg, geometry1Half%mg%sg(5))

    CALL GenContainers%GenStates(sg, X1)

    CALL X1%addVar('psi')
    CALL X1%addVar('chi')
    CALL X1%addVar('rho')
    CALL X1%addVar('div')
    CALL X1%addVar('vor')
    CALL X1%addVar('pres')

    X1%fields(X1%getVarIdx('psi'))%DATA(:, :, 1) = psi
    X1%fields(X1%getVarIdx('chi'))%DATA(:, :, 1) = chi
    X1%fields(X1%getVarIdx('rho'))%DATA(:, :, 1) = rho
    X1%fields(X1%getVarIdx('div'))%DATA(:, :, 1) = div
    X1%fields(X1%getVarIdx('vor'))%DATA(:, :, 1) = vor
    X1%fields(X1%getVarIdx('pres'))%DATA(:, :, 1) = P

    PRINT *, 'ztop in the sg is: ', sg%ztop

    ! CALL Output_NC_State_AV(X1, outputDir, "Test_W")

    CALL Rossby%destroy()

    vLevel1 = (sg%vLevel - 1) * 2 + 1

    ALLOCATE (sigma2(vLevel1))

    sigma2(1:vLevel1:2) = sg%sigma
    sigma2(2:vLevel1:2) = (sg%sigma(1:sg%vLevel - 1) &
                           + sg%sigma(2:sg%vLevel)) / 2.0D0

    DEALLOCATE (zhght, &
                z_s, &
                Hz)

  END ASSOCIATE

  ASSOCIATE (sg => geometry1Half%mg%sg(5), &
             num_cell => geometry1Half%mg%sg(5)%num_cell, &
             vLevel => geometry1Half%mg%sg(5)%vLevel)

    ALLOCATE (tenW1(sg%vLevel, sg%num_cell), &
              tenWt1(sg%vLevel, sg%num_cell), &
              wwnd(sg%vLevel, sg%num_cell), &
              Phalf(sg%vLevel, sg%num_cell), &
              zhght(sg%vLevel, sg%num_cell), &
              z_s(sg%num_cell), &
              Hz(sg%vLevel, sg%num_cell), &
              rhohalf(sg%vLevel, sg%num_cell))

    wwnd = 0.0D0
    CALL Rossby%initial(sg)
    CALL Rossby%Hght_func(sg%cell_cntr, zhght, z_s, Hz)
    CALL Rossby%pres_func(0.0D0, sg%cell_cntr, zhght, z_s, Phalf)
    CALL Rossby%rho_func(0.0D0, sg%cell_cntr, zhght, z_s, rhohalf)
    CALL Rossby%wTen(sg%cell_cntr, rhohalf, Hz, tenWt1)

    sg%topo = z_s

    CALL GenContainers%updateSg(sg, geometry1Full%mg%sg(5))

    TenW = TenW_t(sg)
    CALL TenW%Tendency(wwnd, div, psi, chi, rho, P, tenW1)

    CALL Rossby%destroy()
    CALL TenW%destroy()

    DEALLOCATE (vor, &
                psi, &
                div, &
                chi, &
                rho, rhohalf, &
                P, Phalf, &
                sigma_3d, &
                zhght, &
                z_s, &
                Hz, &
                DM, wwnd)

  END ASSOCIATE

  CALL GenContainers%GenGeometry(geometry2Full, sigma2)
  DEALLOCATE (sigma2)

  vLevel1 = geometry2Full%mg%sg(5)%vLevel - 1
  ALLOCATE (sigma2(vLevel1))
  sigma2(1:vLevel1) = (geometry2Full%mg%sg(5)%sigma(1:vLevel1) &
                       + geometry2Full%mg%sg(5)%sigma(2:vLevel1 + 1)) / 2.0D0

  CALL GenContainers%GenGeometry(geometry2Half, sigma2)
  DEALLOCATE (sigma2)

  ASSOCIATE (sg => geometry2Full%mg%sg(6), &
             num_cell => geometry2Full%mg%sg(6)%num_cell, &
             vLevel => geometry2Full%mg%sg(6)%vLevel)

    ALLOCATE (vor(sg%vLevel, sg%num_cell), &
              psi(sg%vLevel, sg%num_cell), &
              div(sg%vLevel, sg%num_cell), &
              chi(sg%vLevel, sg%num_cell), &
              rho(sg%vLevel, sg%num_cell), &
              P(sg%vLevel, sg%num_cell), &
              sigma_3d(sg%vLevel, sg%num_cell), &
              zhght(sg%vLevel, sg%num_cell), &
              z_s(sg%num_cell), &
              Hz(sg%vLevel, sg%num_cell), &
              DM(sg%vLevel, sg%num_cell))

    DO i = 1, sg%num_cell
      sigma_3d(:, i) = sg%sigma
    END DO
    CALL Rossby%initial(sg)
    CALL Rossby%Hght_func(sg%cell_cntr, zhght, z_s, Hz)
    CALL Rossby%RH_psi_func(0.0D0, sg%cell_cntr, sigma_3d, psi, vor)
    CALL Rossby%RH_chi_func(0.0D0, sg%cell_cntr, sigma_3d, chi, div)
    CALL Rossby%DM_func(Hz, DM)
    CALL Rossby%rho_func(0.0D0, sg%cell_cntr, zhght, z_s, rho)
    CALL Rossby%pres_func(0.0D0, sg%cell_cntr, zhght, z_s, P)
    sg%topo = z_s

    CALL GenContainers%updateSg(sg, geometry2Half%mg%sg(6))

    CALL GenContainers%GenStates(sg, X2)

    CALL X2%addVar('psi')
    CALL X2%addVar('chi')
    CALL X2%addVar('rho')
    CALL X2%addVar('div')
    CALL X2%addVar('vor')
    CALL X2%addVar('pres')

    X2%fields(X2%getVarIdx('psi'))%DATA(:, :, 1) = psi
    X2%fields(X2%getVarIdx('chi'))%DATA(:, :, 1) = chi
    X2%fields(X2%getVarIdx('rho'))%DATA(:, :, 1) = rho
    X2%fields(X2%getVarIdx('div'))%DATA(:, :, 1) = div
    X2%fields(X2%getVarIdx('vor'))%DATA(:, :, 1) = vor
    X2%fields(X2%getVarIdx('pres'))%DATA(:, :, 1) = P

    PRINT *, 'ztop in the sg is: ', sg%ztop

    ! CALL Output_NC_State_AV(X2, outputDir, "Test_W")

    CALL Rossby%destroy()

    vLevel1 = (sg%vLevel - 1) * 2.0D0 + 1

    ALLOCATE (sigma2(vLevel1))

    sigma2(1:vLevel1:2) = sg%sigma
    sigma2(2:vLevel1:2) = (sg%sigma(1:sg%vLevel - 1) &
                           + sg%sigma(2:sg%vLevel)) / 2.0D0

    DEALLOCATE (zhght, &
                z_s, &
                Hz)

  END ASSOCIATE

  ASSOCIATE (sg => geometry2Half%mg%sg(6), &
             num_cell => geometry2Half%mg%sg(6)%num_cell, &
             vLevel => geometry2Half%mg%sg(6)%vLevel)

    ALLOCATE (tenW2(sg%vLevel, sg%num_cell), &
              tenWt2(sg%vLevel, sg%num_cell), &
              wwnd(sg%vLevel, sg%num_cell), &
              Phalf(sg%vLevel, sg%num_cell), &
              zhght(sg%vLevel, sg%num_cell), &
              z_s(sg%num_cell), &
              Hz(sg%vLevel, sg%num_cell), &
              rhohalf(sg%vLevel, sg%num_cell))

    wwnd = 0.0D0
    CALL Rossby%initial(sg)
    CALL Rossby%Hght_func(sg%cell_cntr, zhght, z_s, Hz)
    CALL Rossby%pres_func(0.0D0, sg%cell_cntr, zhght, z_s, Phalf)
    CALL Rossby%rho_func(0.0D0, sg%cell_cntr, zhght, z_s, rhohalf)
    CALL Rossby%wTen(sg%cell_cntr, rhohalf, Hz, tenWt2)

    sg%topo = z_s

    CALL GenContainers%updateSg(sg, geometry2Full%mg%sg(6))

    TenW = TenW_t(sg)
    CALL TenW%Tendency(wwnd, div, psi, chi, rho, P, tenW2)

    CALL Rossby%destroy()
    CALL TenW%destroy()

    DEALLOCATE (vor, &
                psi, &
                div, &
                chi, &
                rho, rhohalf, &
                P, Phalf, &
                sigma_3d, &
                zhght, &
                z_s, &
                Hz, &
                DM, wwnd)

  END ASSOCIATE

  CALL GenContainers%GenGeometry(geometry3Full, sigma2)
  DEALLOCATE (sigma2)

  vLevel1 = geometry3Full%mg%sg(5)%vLevel - 1
  ALLOCATE (sigma2(vLevel1))
  sigma2(1:vLevel1) = (geometry3Full%mg%sg(5)%sigma(1:vLevel1) &
                       + geometry3Full%mg%sg(5)%sigma(2:vLevel1 + 1)) / 2.0D0

  CALL GenContainers%GenGeometry(geometry3Half, sigma2)
  DEALLOCATE (sigma2)

  ASSOCIATE (sg => geometry3Full%mg%sg(7), &
             num_cell => geometry3Full%mg%sg(7)%num_cell, &
             vLevel => geometry3Full%mg%sg(7)%vLevel)

    ALLOCATE (vor(sg%vLevel, sg%num_cell), &
              psi(sg%vLevel, sg%num_cell), &
              div(sg%vLevel, sg%num_cell), &
              chi(sg%vLevel, sg%num_cell), &
              rho(sg%vLevel, sg%num_cell), &
              P(sg%vLevel, sg%num_cell), &
              sigma_3d(sg%vLevel, sg%num_cell), &
              zhght(sg%vLevel, sg%num_cell), &
              z_s(sg%num_cell), &
              Hz(sg%vLevel, sg%num_cell), &
              DM(sg%vLevel, sg%num_cell))

    DO i = 1, sg%num_cell
      sigma_3d(:, i) = sg%sigma
    END DO
    CALL Rossby%initial(sg)
    CALL Rossby%Hght_func(sg%cell_cntr, zhght, z_s, Hz)
    CALL Rossby%RH_psi_func(0.0D0, sg%cell_cntr, sigma_3d, psi, vor)
    CALL Rossby%RH_chi_func(0.0D0, sg%cell_cntr, sigma_3d, chi, div)
    CALL Rossby%DM_func(Hz, DM)
    CALL Rossby%rho_func(0.0D0, sg%cell_cntr, zhght, z_s, rho)
    CALL Rossby%pres_func(0.0D0, sg%cell_cntr, zhght, z_s, P)
    sg%topo = z_s

    CALL GenContainers%updateSg(sg, geometry3Half%mg%sg(7))

    CALL GenContainers%GenStates(sg, X3)

    CALL X3%addVar('psi')
    CALL X3%addVar('chi')
    CALL X3%addVar('rho')
    CALL X3%addVar('div')
    CALL X3%addVar('vor')
    CALL X3%addVar('pres')

    X3%fields(X3%getVarIdx('psi'))%DATA(:, :, 1) = psi
    X3%fields(X3%getVarIdx('chi'))%DATA(:, :, 1) = chi
    X3%fields(X3%getVarIdx('rho'))%DATA(:, :, 1) = rho
    X3%fields(X3%getVarIdx('div'))%DATA(:, :, 1) = div
    X3%fields(X3%getVarIdx('vor'))%DATA(:, :, 1) = vor
    X3%fields(X3%getVarIdx('pres'))%DATA(:, :, 1) = P

    PRINT *, 'ztop in the sg is: ', sg%ztop

    ! CALL Output_NC_State_AV(X3, outputDir, "Test_W")

    CALL Rossby%destroy()

    vLevel1 = (sg%vLevel - 1) * 2.0D0 + 1

    ALLOCATE (sigma2(vLevel1))

    sigma2(1:vLevel1:2) = sg%sigma
    sigma2(2:vLevel1:2) = (sg%sigma(1:sg%vLevel - 1) &
                           + sg%sigma(2:sg%vLevel)) / 2.0D0

    DEALLOCATE (zhght, &
                z_s, &
                Hz)

  END ASSOCIATE

  ASSOCIATE (sg => geometry3Half%mg%sg(7), &
             num_cell => geometry3Half%mg%sg(7)%num_cell, &
             vLevel => geometry3Half%mg%sg(7)%vLevel)

    ALLOCATE (tenW3(sg%vLevel, sg%num_cell), &
              tenWt3(sg%vLevel, sg%num_cell), &
              wwnd(sg%vLevel, sg%num_cell), &
              Phalf(sg%vLevel, sg%num_cell), &
              zhght(sg%vLevel, sg%num_cell), &
              z_s(sg%num_cell), &
              Hz(sg%vLevel, sg%num_cell), &
              rhohalf(sg%vLevel, sg%num_cell))

    wwnd = 0.0D0
    CALL Rossby%initial(sg)
    CALL Rossby%Hght_func(sg%cell_cntr, zhght, z_s, Hz)
    CALL Rossby%pres_func(0.0D0, sg%cell_cntr, zhght, z_s, Phalf)
    CALL Rossby%rho_func(0.0D0, sg%cell_cntr, zhght, z_s, rhohalf)
    CALL Rossby%wTen(sg%cell_cntr, rhohalf, Hz, tenWt3)

    sg%topo = z_s

    CALL GenContainers%updateSg(sg, geometry3Full%mg%sg(7))

    TenW = TenW_t(sg)
    CALL TenW%Tendency(wwnd, div, psi, chi, rho, P, tenW3)

    CALL Rossby%destroy()
    CALL TenW%destroy()

    DEALLOCATE (vor, &
                psi, &
                div, &
                chi, &
                rho, rhohalf, &
                P, Phalf, &
                sigma_3d, &
                zhght, &
                z_s, &
                Hz, &
                DM, wwnd)

  END ASSOCIATE

  ALLOCATE (errorl_1(geometry1half%mg%sg(5)%vLevel), &
            errorl_2(geometry2half%mg%sg(5)%vLevel), &
            errorl_3(geometry3half%mg%sg(5)%vLevel))
  errorw_1 = 0.0D0
  DO k = 1, geometry1half%mg%sg(5)%vLevel
    DO i = 1, geometry1half%mg%sg(5)%num_icell
      IF (geometry1half%mg%sg(5)%bdy_type(i) .LT. 0) THEN
        IF (errorw_1 .LE. DABS(tenW1(k, i) - tenWt1(k, i))) THEN
          errorw_1 = DABS(tenW1(k, i) - tenWt1(k, i))
        END IF
      END IF
    END DO
    errorl_1(k) = errorw_1
  END DO

  errorw_2 = 0.0D0
  DO k = 1, geometry2half%mg%sg(6)%vLevel
    DO i = 1, geometry2half%mg%sg(6)%num_icell
      IF (geometry2half%mg%sg(6)%bdy_type(i) .LT. 0) THEN
        IF (errorw_2 .LE. DABS(tenW2(k, i) - tenWt2(k, i))) THEN
          errorw_2 = DABS(tenW2(k, i) - tenWt2(k, i))
        END IF
      END IF
    END DO
    errorl_2(k) = errorw_2
  END DO

  errorw_3 = 0.0D0
  DO k = 1, geometry3half%mg%sg(7)%vLevel
    DO i = 1, geometry3half%mg%sg(7)%num_icell
      IF (geometry3half%mg%sg(7)%bdy_type(i) .LT. 0) THEN
        IF (errorw_3 .LE. DABS(tenW3(k, i) - tenWt3(k, i))) THEN
          errorw_3 = DABS(tenW3(k, i) - tenWt3(k, i))
        END IF
      END IF
    END DO
    errorl_3(k) = errorw_3
  END DO

  ! errorw_1 = 0.0D0
  ! DO k = 1, geometry1half%mg%sg(5)%vLevel
  !    DO i = 1, geometry1half%mg%sg(5)%num_icell
  !       IF (geometry1half%mg%sg(5)%bdy_type(i) .LT. 0) THEN
  !          IF (errorw_1 .LE. DABS(tenW1(k, i))) THEN
  !             errorw_1 = DABS(tenW1(k, i))
  !          END IF
  !       END IF
  !    END DO
  ! END DO

  ! errorw_2 = 0.0D0
  ! DO k = 1, geometry2half%mg%sg(6)%vLevel
  !    DO i = 1, geometry2half%mg%sg(6)%num_icell
  !       IF (geometry2half%mg%sg(6)%bdy_type(i) .LT. 0) THEN
  !          IF (errorw_2 .LE. DABS(tenW2(k, i))) THEN
  !             errorw_2 = DABS(tenW2(k, i))
  !          END IF
  !       END IF
  !    END DO
  ! END DO

  ! errorw_3 = 0.0D0
  ! DO k = 1, geometry3half%mg%sg(7)%vLevel
  !    DO i = 1, geometry3half%mg%sg(7)%num_icell
  !       IF (geometry3half%mg%sg(7)%bdy_type(i) .LT. 0) THEN
  !          IF (errorw_3 .LE. DABS(tenW3(k, i))) THEN
  !             errorw_3 = DABS(tenW3(k, i))
  !          END IF
  !       END IF
  !    END DO
  ! END DO

  PRINT *, 'error ratio is ', errorw_1, errorw_2, errorw_3, errorw_1 / errorw_2, errorw_2 / errorw_3

  PRINT *, 'error by level '
  DO k = 1, geometry1half%mg%sg(5)%vLevel
    PRINT *, errorl_1(k) / errorl_2(2 * k - 1)
  END DO

  CALL GenContainers%destroy()

END PROGRAM test_tenW
