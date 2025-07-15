PROGRAM test_tenrho
  USE kinds_m, ONLY: i_kind, r_kind
  USE TenRho_m, ONLY: TenRho_t
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
  USE parameters_m, ONLY: Omega
  USE State2NC_m
  USE gzm_m, ONLY: gzm_t

  TYPE(TenRho_t) :: TenRho
  TYPE(State_t) :: X1, X2, X3
  TYPE(PreCal_t) :: PreCal
  TYPE(UnitTestDyn_t) :: Rossby
  TYPE(GenContainers_t) :: GenContainers
  TYPE(mpddGlob_t), TARGET :: mpddGlob
  !  TYPE(SingleGrid_t) :: sg
  TYPE(geometry_t) :: geometry1Full, geometry1Half, geometry2Full, geometry2Half, geometry3Full, geometry3Half
  TYPE(gzm_t) :: gzm

  INTEGER(i_kind) :: i, k, vLevel1, vLevel2, ifile
  REAL(r_kind) :: errorv_1, errorv_2, errorv_3

  CHARACTER(LEN=1024) :: configFile, outputDir

  ! Fields to test:
  REAL(r_kind), ALLOCATABLE :: tenrho1(:, :), tenrho2(:, :), tenrho3(:, :), &
                               tenrhot1(:, :), tenrhot2(:, :), tenrhot3(:, :), &
                               vor(:, :), psi(:, :), &
                               div(:, :), chi(:, :), &
                               wwnd(:, :), rho(:, :), P(:, :), &
                               sigma_3d(:, :), sigma(:), sigma2(:), &
                               zhght(:, :), z_s(:), Hz(:, :), &
                               DM(:, :)

  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", configFile)
  configFile = TRIM(configFile)//"/gzm_and_laplace_terrain.yaml"

  ifile = yaml_get_var(TRIM(configFile), 'IO', 'output_dir', outputDir)
  CALL GET_ENVIRONMENT_VARIABLE("OUTPUT_DIR", outputDir)

  GenContainers = GenContainers_t(TRIM(configFile))
  CALL GenContainers%GenGeometry(geometry1Full)

  vLevel1 = geometry1Full%mg%sg(5)%vLevel - 1

  ALLOCATE (sigma(vLevel1))

  sigma(1:vLevel1) = (geometry1Full%mg%sg(5)%sigma(1:vLevel1) &
                      + geometry1Full%mg%sg(5)%sigma(2:vLevel1 + 1)) / 2.0D0

  CALL GenContainers%GenGeometry(geometry1Half, sigma)

  DEALLOCATE (sigma)

  ASSOCIATE (sg => geometry1Half%mg%sg(5), &
             num_cell => geometry1Half%mg%sg(5)%num_cell, &
             vLevel => geometry1Half%mg%sg(5)%vLevel)

    ALLOCATE (wwnd(sg%vLevel, sg%num_cell))

    wwnd = 0.0D0

  END ASSOCIATE

  ASSOCIATE (sg => geometry1Full%mg%sg(5), &
             num_cell => geometry1Full%mg%sg(5)%num_cell, &
             vLevel => geometry1Full%mg%sg(5)%vLevel)

    ALLOCATE (tenrho1(sg%vLevel, sg%num_cell), &
              !  tenrho2(sg%vLevel, sg%num_cell), &
              !  tenrho3(sg%vLevel, sg%num_cell), &
              tenrhot1(sg%vLevel, sg%num_cell), &
              !  tenrhot2(sg%vLevel, sg%num_cell), &
              !  tenrhot3(sg%vLevel, sg%num_cell), &
              vor(sg%vLevel, sg%num_cell), &
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

    ! P = 0.0D0
    ! rho = 0.5D0
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
    CALL X1%addVar('par_rho_sigma')
    CALL X1%addVar('par_rho_sigmat')
    CALL X1%addVar('F_rho_chit')
    CALL X1%addVar('J_rho_psit')
    CALL X1%addVar('F_rho_chi')
    CALL X1%addVar('J_rho_psi')

    ! CALL X1%addVar('F_eta_chi')
    ! CALL X1%addVar('DM')
    ! CALL X1%addVar('F_DM_psisigma')
    ! CALL X1%addVar('J_DM_chisigma')
    ! CALL X1%addVar('psi_sigma')
    ! CALL X1%addVar('J_eta_psit')
    ! CALL X1%addVar('F_eta_chit')
    ! CALL X1%addVar('DMt')
    ! CALL X1%addVar('F_DM_psisigmat')
    ! CALL X1%addVar('J_DM_chisigmat')
    ! CALL X1%addVar('psi_sigmat')

    CALL Rossby%RhoTen(sg%cell_cntr, rho, div, psi, chi, DM, zhght, Hz, tenrhot1, X1)

    X1%fields(X1%getVarIdx('psi'))%DATA(:, :, 1) = psi
    X1%fields(X1%getVarIdx('chi'))%DATA(:, :, 1) = chi
    X1%fields(X1%getVarIdx('rho'))%DATA(:, :, 1) = rho
    X1%fields(X1%getVarIdx('div'))%DATA(:, :, 1) = div
    X1%fields(X1%getVarIdx('vor'))%DATA(:, :, 1) = vor
    X1%fields(X1%getVarIdx('pres'))%DATA(:, :, 1) = P
    ! X1%fields(X1%getVarIdx('DMt'))%data(:, :, 1) = DM
    ! X1%fields(X1%getVarIdx('psi_sigmat'))%data(:, :, 1) = Rossby%psi1st(1, :, :)
    ! IF (ASSOCIATED(gzm%sg)) NULLIFY (gzm%sg)
    ! gzm = gzm_t(sg)
    ! CALL gzm%Laplacia(psi, vor)

    PreCal = PreCal_t(sg)
    CALL PreCal%PreCals(X1)

    TenRho = TenRho_t(sg)
    CALL TenRho%Tendency(rho, div, psi, chi, wwnd, PreCal, tenrho1, X1)
    PRINT *, 'd_parrho_sigma in X1 is: ', &
      MAXVAL(dabs(X1%fields(X1%getVarIdx('F_rho_chit'))%DATA(:, :, 1) - &
                  X1%fields(X1%getVarIdx('F_rho_chi'))%DATA(:, :, 1))), &
      MAXVAL(dabs(X1%fields(X1%getVarIdx('J_rho_psit'))%DATA(:, :, 1) - &
                  X1%fields(X1%getVarIdx('J_rho_psi'))%DATA(:, :, 1)))

    ! CALL X1%addVar('tenrhot')
    ! CALL X1%addVar('tenrho')
    ! CALL X1%addVar('dtenrho')
    ! CALL X1%addVar('dz')
    ! X1%fields(X1%getVarIdx('tenrho'))%data(:, :, 1) = tenrho1
    ! X1%fields(X1%getVarIdx('tenrhot'))%data(:, :, 1) = tenrhot1
    ! X1%fields(X1%getVarIdx('dtenrho'))%data(:, :, 1) = dabs(tenrhot1 - tenrho1)
    ! X1%fields(X1%getVarIdx('dz'))%data(:, :, 1) = sg%Hz - Hz

    PRINT *, 'ztop in the sg is: ', sg%ztop

    CALL Output_NC_State_AV(X1, outputDir, "Test_rho")

    CALL PreCal%destroy()

    CALL TenRho%destroy()

    CALL Rossby%destroy()

    vLevel1 = (sg%vLevel - 1) * 2 + 1

    ALLOCATE (sigma2(vLevel1))

    sigma2(1:vLevel1:2) = sg%sigma
    sigma2(2:vLevel1:2) = (sg%sigma(1:sg%vLevel - 1) &
                           + sg%sigma(2:sg%vLevel)) / 2.0D0

    DEALLOCATE (vor, &
                psi, &
                div, &
                chi, &
                rho, &
                P, &
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

  ASSOCIATE (sg => geometry2Half%mg%sg(6), &
             num_cell => geometry2Half%mg%sg(6)%num_cell, &
             vLevel => geometry2Half%mg%sg(6)%vLevel)

    ALLOCATE (wwnd(sg%vLevel, sg%num_cell))

    wwnd = 0.0D0

  END ASSOCIATE

  ASSOCIATE (sg => geometry2Full%mg%sg(6), &
             num_cell => geometry2Full%mg%sg(6)%num_cell, &
             vLevel => geometry2Full%mg%sg(6)%vLevel)

    ALLOCATE (tenrho2(sg%vLevel, sg%num_cell), &
              !  tenrho2(sg%vLevel, sg%num_cell), &
              !  tenrho3(sg%vLevel, sg%num_cell), &
              tenrhot2(sg%vLevel, sg%num_cell), &
              !  tenrhot2(sg%vLevel, sg%num_cell), &
              !  tenrhot3(sg%vLevel, sg%num_cell), &
              vor(sg%vLevel, sg%num_cell), &
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

    ! P = 0.0D0
    ! rho = 0.5D0
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
    CALL X2%addVar('par_rho_sigma')
    CALL X2%addVar('par_rho_sigmat')
    CALL X2%addVar('F_rho_chit')
    CALL X2%addVar('J_rho_psit')
    CALL X2%addVar('F_rho_chi')
    CALL X2%addVar('J_rho_psi')
    ! CALL X2%addVar('J_eta_psi')
    ! CALL X2%addVar('F_eta_chi')
    ! CALL X2%addVar('DM')
    ! CALL X2%addVar('F_DM_psisigma')
    ! CALL X2%addVar('J_DM_chisigma')
    ! CALL X2%addVar('psi_sigma')
    ! CALL X2%addVar('J_eta_psit')
    ! CALL X2%addVar('F_eta_chit')
    ! CALL X2%addVar('DMt')
    ! CALL X2%addVar('F_DM_psisigmat')
    ! CALL X2%addVar('J_DM_chisigmat')
    ! CALL X2%addVar('psi_sigmat')

    CALL Rossby%RhoTen(sg%cell_cntr, rho, div, psi, chi, DM, zhght, Hz, tenrhot2, X2)

    X2%fields(X2%getVarIdx('psi'))%DATA(:, :, 1) = psi
    X2%fields(X2%getVarIdx('chi'))%DATA(:, :, 1) = chi
    X2%fields(X2%getVarIdx('rho'))%DATA(:, :, 1) = rho
    X2%fields(X2%getVarIdx('div'))%DATA(:, :, 1) = div
    X2%fields(X2%getVarIdx('vor'))%DATA(:, :, 1) = vor
    X2%fields(X2%getVarIdx('pres'))%DATA(:, :, 1) = P
    ! X2%fields(X2%getVarIdx('DMt'))%data(:, :, 1) = DM
    ! X2%fields(X2%getVarIdx('psi_sigmat'))%data(:, :, 1) = Rossby%psi1st(1, :, :)
    ! IF (ASSOCIATED(gzm%sg)) NULLIFY (gzm%sg)
    ! gzm = gzm_t(sg)
    ! CALL gzm%Laplacia(psi, vor)

    PreCal = PreCal_t(sg)
    TenRho = TenRho_t(sg)
    CALL PreCal%PreCals(X2)

    CALL TenRho%Tendency(rho, div, psi, chi, wwnd, PreCal, tenrho2, X2)

    PRINT *, 'd_parrho_sigma in X2 is: ', &
      MAXVAL(dabs(X2%fields(X2%getVarIdx('F_rho_chit'))%DATA(:, :, 1) - &
                  X2%fields(X2%getVarIdx('F_rho_chi'))%DATA(:, :, 1))), &
      MAXVAL(dabs(X2%fields(X2%getVarIdx('J_rho_psit'))%DATA(:, :, 1) - &
                  X2%fields(X2%getVarIdx('J_rho_psi'))%DATA(:, :, 1)))

    ! CALL X2%addVar('tenrhot')
    ! CALL X2%addVar('tenrho')
    ! CALL X2%addVar('dtenrho')
    ! CALL X2%addVar('dz')
    ! X2%fields(X2%getVarIdx('tenrho'))%data(:, :, 1) = tenrho2
    ! X2%fields(X2%getVarIdx('tenrhot'))%data(:, :, 1) = tenrhot2
    ! X2%fields(X2%getVarIdx('dtenrho'))%data(:, :, 1) = dabs(tenrhot12 - tenrho2)
    ! X2%fields(X2%getVarIdx('dz'))%data(:, :, 1) = sg%Hz - Hz

    PRINT *, 'ztop in the sg is: ', sg%ztop

    CALL Output_NC_State_AV(X2, outputDir, "Test_rho_d2")

    CALL PreCal%destroy()

    CALL TenRho%destroy()

    CALL Rossby%destroy()

    vLevel1 = (sg%vLevel - 1) * 2 + 1

    ALLOCATE (sigma2(vLevel1))

    sigma2(1:vLevel1:2) = sg%sigma
    sigma2(2:vLevel1:2) = (sg%sigma(1:sg%vLevel - 1) &
                           + sg%sigma(2:sg%vLevel)) / 2.0D0

    DEALLOCATE (vor, &
                psi, &
                div, &
                chi, &
                rho, &
                P, &
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

  ASSOCIATE (sg => geometry3Half%mg%sg(7), &
             num_cell => geometry3Half%mg%sg(7)%num_cell, &
             vLevel => geometry3Half%mg%sg(7)%vLevel)

    ALLOCATE (wwnd(sg%vLevel, sg%num_cell))

    wwnd = 0.0D0

  END ASSOCIATE

  ASSOCIATE (sg => geometry3Full%mg%sg(7), &
             num_cell => geometry3Full%mg%sg(7)%num_cell, &
             vLevel => geometry3Full%mg%sg(7)%vLevel)

    ALLOCATE (tenrho3(sg%vLevel, sg%num_cell), &
              !  tenrho2(sg%vLevel, sg%num_cell), &
              !  tenrho3(sg%vLevel, sg%num_cell), &
              tenrhot3(sg%vLevel, sg%num_cell), &
              !  tenrhot2(sg%vLevel, sg%num_cell), &
              !  tenrhot3(sg%vLevel, sg%num_cell), &
              vor(sg%vLevel, sg%num_cell), &
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

    ! P = 0.0D0
    ! rho = 0.5D0
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
    CALL X3%addVar('par_rho_sigma')
    CALL X3%addVar('par_rho_sigmat')
    CALL X3%addVar('F_rho_chit')
    CALL X3%addVar('J_rho_psit')
    CALL X3%addVar('F_rho_chi')
    CALL X3%addVar('J_rho_psi')
    ! CALL X3%addVar('J_eta_psi')
    ! CALL X3%addVar('F_eta_chi')
    ! CALL X3%addVar('DM')
    ! CALL X3%addVar('F_DM_psisigma')
    ! CALL X3%addVar('J_DM_chisigma')
    ! CALL X3%addVar('psi_sigma')
    ! CALL X3%addVar('J_eta_psit')
    ! CALL X3%addVar('F_eta_chit')
    ! CALL X3%addVar('DMt')
    ! CALL X3%addVar('F_DM_psisigmat')
    ! CALL X3%addVar('J_DM_chisigmat')
    ! CALL X3%addVar('psi_sigmat')

    CALL Rossby%RhoTen(sg%cell_cntr, rho, div, psi, chi, DM, zhght, Hz, tenrhot3, X3)

    X3%fields(X3%getVarIdx('psi'))%DATA(:, :, 1) = psi
    X3%fields(X3%getVarIdx('chi'))%DATA(:, :, 1) = chi
    X3%fields(X3%getVarIdx('rho'))%DATA(:, :, 1) = rho
    X3%fields(X3%getVarIdx('div'))%DATA(:, :, 1) = div
    X3%fields(X3%getVarIdx('vor'))%DATA(:, :, 1) = vor
    X3%fields(X3%getVarIdx('pres'))%DATA(:, :, 1) = P
    ! X3%fields(X3%getVarIdx('DMt'))%data(:, :, 1) = DM
    ! X3%fields(X3%getVarIdx('psi_sigmat'))%data(:, :, 1) = Rossby%psi1st(1, :, :)
    ! IF (ASSOCIATED(gzm%sg)) NULLIFY (gzm%sg)
    ! gzm = gzm_t(sg)
    ! CALL gzm%Laplacia(psi, vor)

    PreCal = PreCal_t(sg)
    TenRho = TenRho_t(sg)
    CALL PreCal%PreCals(X3)
    CALL TenRho%Tendency(rho, div, psi, chi, wwnd, PreCal, tenrho3, X3)
    PRINT *, 'd_parrho_sigma in X3 is: ', &
      MAXVAL(dabs(X3%fields(X3%getVarIdx('F_rho_chit'))%DATA(:, :, 1) - &
                  X3%fields(X3%getVarIdx('F_rho_chi'))%DATA(:, :, 1))), &
      MAXVAL(dabs(X3%fields(X3%getVarIdx('J_rho_psit'))%DATA(:, :, 1) - &
                  X3%fields(X3%getVarIdx('J_rho_psi'))%DATA(:, :, 1)))

    ! CALL X3%addVar('tenrhot')
    ! CALL X3%addVar('tenrho')
    ! CALL X3%addVar('dtenrho')
    ! CALL X3%addVar('dz')
    ! X3%fields(X3%getVarIdx('tenrho'))%data(:, :, 1) = tenrho2
    ! X3%fields(X3%getVarIdx('tenrhot'))%data(:, :, 1) = tenrhot2
    ! X3%fields(X3%getVarIdx('dtenrho'))%data(:, :, 1) = dabs(tenrhot12 - tenrho2)
    ! X3%fields(X3%getVarIdx('dz'))%data(:, :, 1) = sg%Hz - Hz

    PRINT *, 'ztop in the sg is: ', sg%ztop

    CALL Output_NC_State_AV(X3, outputDir, "Test_rho_d3")

    CALL PreCal%destroy()

    CALL TenRho%destroy()

    CALL Rossby%destroy()

    DEALLOCATE (vor, &
                psi, &
                div, &
                chi, &
                rho, &
                P, &
                sigma_3d, &
                zhght, &
                z_s, &
                Hz, &
                DM, wwnd)

  END ASSOCIATE

  errorv_1 = 0.0D0
  DO k = 1, geometry1Full%mg%sg(5)%vLevel
    DO i = 1, geometry1Full%mg%sg(5)%num_icell
      IF (geometry1Full%mg%sg(5)%bdy_type(i) .LT. 0) THEN
        IF (errorv_1 .LE. DABS(tenrho1(k, i) - tenrhot1(k, i))) THEN
          errorv_1 = DABS(tenrho1(k, i) - tenrhot1(k, i))
        END IF
      END IF
    END DO
  END DO

  errorv_2 = 0.0D0
  DO k = 1, geometry2Full%mg%sg(6)%vLevel
    DO i = 1, geometry2Full%mg%sg(6)%num_icell
      IF (geometry2Full%mg%sg(6)%bdy_type(i) .LT. 0) THEN
        IF (errorv_2 .LE. DABS(tenrho2(k, i) - tenrhot2(k, i))) THEN
          errorv_2 = DABS(tenrho2(k, i) - tenrhot2(k, i))
        END IF
      END IF
    END DO
  END DO

  errorv_3 = 0.0D0
  DO k = 1, geometry3Full%mg%sg(7)%vLevel
    DO i = 1, geometry3Full%mg%sg(7)%num_icell
      IF (geometry3Full%mg%sg(7)%bdy_type(i) .LT. 0) THEN
        IF (errorv_3 .LE. DABS(tenrho3(k, i) - tenrhot3(k, i))) THEN
          errorv_3 = DABS(tenrho3(k, i) - tenrhot3(k, i))
        END IF
      END IF
    END DO
  END DO

  PRINT *, 'error ratio is ', errorv_1, errorv_2, errorv_3, errorv_1 / errorv_2, errorv_2 / errorv_3

  CALL GenContainers%destroy()

END PROGRAM test_tenrho
