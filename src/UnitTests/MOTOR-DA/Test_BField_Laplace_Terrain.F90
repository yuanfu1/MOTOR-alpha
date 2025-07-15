PROGRAM BField_Laplace_Terrain_test

  USE BFieldLaplace_m, ONLY: BFieldLaplace_t
  USE GenContainers_m, ONLY: GenContainers_t
  USE geometry_m, ONLY: Geometry_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE State_m, ONLY: State_t
  USE kinds_m, ONLY: i_kind, r_kind
  USE IOGrapes_m, ONLY: IOGrapes_t
  USE ModelCoupler_m, ONLY: ModelCoupler_t
  USE parameters_m, ONLY: pi
!    USE Filter_m, ONLY: smoothField
  USE MGOpts_m
  USE YAMLRead_m

  IMPLICIT NONE

  TYPE(BFieldLaplace_t):: BField_Laplace
  TYPE(GenContainers_t) :: GenContainers
  TYPE(mpddGlob_t), TARGET :: mpddGlob
!    TYPE(SingleGrid_t) :: sg
  TYPE(geometry_t) :: geometry, geometry2, geometry3
  TYPE(State_t), ALLOCATABLE :: XbMG(:)
  TYPE(IOGrapes_t) :: ioGrapes
  TYPE(ModelCoupler_t) :: ModelCoupler
  CHARACTER(LEN=1024) :: configFile, varName
  INTEGER(i_kind) :: vLevel, i, k
  REAL(r_kind) :: error_1rst_1, error_1rst_2, error_1rst_3, &
                  error_2end_1, error_2end_2, error_2end_3, &
                  first_A2, first_A1, &
                  error_first_ratio1, error_second_ratio1, &
                  error_first_ratio2, error_second_ratio2
  REAL(r_kind), DIMENSION(:), ALLOCATABLE :: sigma2
  REAL(r_kind), DIMENSION(:, :), ALLOCATABLE:: A, A2, A3, parA_sigma, &
                                               parA_sigma2, parA_sigma3, &
                                               parA_sigmat, &
                                               parA_sigmat2, &
                                               parA_sigmat3, &
                                               parA_sigma_scd, &
                                               parA_sigma_scd2, &
                                               parA_sigma_scd3, &
                                               parA_sigma_scdt, &
                                               parA_sigma_scdt2, &
                                               parA_sigma_scdt3
  REAL(r_kind), ALLOCATABLE :: coef11(:, :, :), coef12(:, :, :), coef13(:, :, :)

  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", configFile)
  configFile = TRIM(configFile)//"/gzm_and_laplace_terrain.yaml"
  varName = 'uwnd'
!    PRINT *, 'configFile in test is', configFile

  GenContainers = GenContainers_t(TRIM(configFile))
  CALL GenContainers%GenGeometry(geometry)

  ALLOCATE (XbMG(geometry%mg%mg_coarsest:geometry%mg%mg_finest))

  ! Initialize the model coupler, to get parameters for topogrophy smoothing
!    CALL ModelCoupler%initialize(configFile, geometry)

  ASSOCIATE (sg => geometry%mg%sg(5))
    ALLOCATE (A(sg%vLevel, sg%num_cell), &
              parA_sigma(sg%vLevel, sg%num_cell), &
              parA_sigma_scd(sg%vLevel, sg%num_cell), &
              parA_sigmat(sg%vLevel, sg%num_cell), &
              parA_sigma_scdt(sg%vLevel, sg%num_cell), &
              coef11(sg%vLevel, sg%num_cell, 3))
    PRINT *, 'sigma is ', sg%sigma
    DO i = 1, sg%vLevel
      A(i, :) = dcos((sg%sigma(i)) / 1.0D04)
      parA_sigmat(i, :) = -1.0D-4 * dsin(sg%sigma(i) / 1.0D04)
      parA_sigma_scdt(i, :) = -1.0D-8 * dcos(sg%sigma(i) / 1.0D04)
    END DO

    !   coef11 = 0.0D0

    sg%topo = 0.0D0
    sg%topo = DSIN(sg%cell_cntr(1, :) / 180.0D0 * pi) + DCOS(sg%cell_cntr(2, :) / 180.0D0 * pi)
    !   Lztopo = -
    !   DO k = 1, sg%vLevel
    !      LzG5(k, :) = (1.0D0 - sg%sigma(k)/sg%ztop)*
    !      Lzt =
    !   END DO

    CALL BField_Laplace%initialize(configFile, sg, varName)
    PRINT *, 'BField initialize finished'

    !   PRINT *, 'size of coef11 and coef_fstdif are', size(coef11), size(sg%coef_fstdif)
    !   coef11 = sg%coef_fstdif

    CALL BField_Laplace%FirstOrderDerivative(A, sg%coef_fstdif, parA_sigma)
    PRINT *, 'Firstorder finished'
    !   PRINT *, sg%coef_scddif
    CALL BField_Laplace%SecondOrderDerivative(A, sg%coef_scddif, parA_sigma_scd)
    PRINT *, 'Secondorder finished'

    vLevel = (sg%vLevel - 1) * 2 + 1

    ALLOCATE (sigma2(vLevel))

    sigma2(1:vLevel:2) = sg%sigma
    sigma2(2:vLevel:2) = (sg%sigma(1:sg%vLevel - 1) &
                          + sg%sigma(2:sg%vLevel)) / 2.0D0

    !   PRINT *, 'dsigma is ', sg%sigma(2:sg%vLevel) - sg%sigma(1:sg%vLevel - 1)
    !   PRINT *, 'dsigma2 is ', sigma2(2:vLevel) - sigma2(1:vLevel - 1)

  END ASSOCIATE

  PRINT *, 'debug here'

!    PRINT *, 'sigma2 is ', sigma2

!    GenContainers = GenContainers_t(TRIM(configFile))
!    PRINT *, 'GenContainers initialize finished'
  CALL GenContainers%GenGeometry(geometry2, sigma2)
  PRINT *, 'GenContainers geometry is finished'
  ASSOCIATE (sg => geometry2%mg%sg(5))

    PRINT *, 'size(sg%sigma) is ', SIZE(sg%sigma), vLevel

    sg%topo = 0.0D0
    sg%topo = DSIN(sg%cell_cntr(1, :)) + DCOS(sg%cell_cntr(2, :))

    ALLOCATE (A2(sg%vLevel, sg%num_cell), &
              parA_sigma2(sg%vLevel, sg%num_cell), &
              parA_sigma_scd2(sg%vLevel, sg%num_cell), &
              parA_sigmat2(sg%vLevel, sg%num_cell), &
              parA_sigma_scdt2(sg%vLevel, sg%num_cell), &
              coef12(sg%vLevel, sg%num_cell, 3))

    !   coef12 = 0.0D0

    DO i = 1, sg%vLevel
      A2(i, :) = dcos((sg%sigma(i)) / 1.0D04)
      parA_sigmat2(i, :) = -1.0D-4 * dsin(sg%sigma(i) / 1.0D04)
      parA_sigma_scdt2(i, :) = -1.0D-8 * dcos(sg%sigma(i) / 1.0D04)
    END DO

    CALL BField_Laplace%initialize(configFile, sg, varName)
    !   coef12 = sg%coef_fstdif

    !   DO i = 1, (sg%vLevel - 1)/2
    !      PRINT *, 'coef_first_1 is ', coef11(i, 100, :)
    !      PRINT *, 'coef_first_2 is ', coef12(i*2 - 1, 100, :)
    !   END DO

    CALL BField_Laplace%FirstOrderDerivative(A2, sg%coef_fstdif, parA_sigma2)
    CALL BField_Laplace%SecondOrderDerivative(A2, sg%coef_scddif, parA_sigma_scd2)

    !   PRINT *, 'parA_sigma and parA_sigma2 are ', size(parA_sigma), size(parA_sigma2)

    vLevel = (sg%vLevel - 1) * 2 + 1

    DEALLOCATE (sigma2)

    ALLOCATE (sigma2(vLevel))

    sigma2(1:vLevel:2) = sg%sigma
    sigma2(2:vLevel:2) = (sg%sigma(1:sg%vLevel - 1) &
                          + sg%sigma(2:sg%vLevel)) / 2.0D0

  END ASSOCIATE

  CALL GenContainers%GenGeometry(geometry3, sigma2)
  PRINT *, 'GenContainers geometry is finished'
  ASSOCIATE (sg => geometry3%mg%sg(5))

    PRINT *, 'size(sg%sigma) is ', SIZE(sg%sigma), vLevel

    sg%topo = 0.0D0
    sg%topo = DSIN(sg%cell_cntr(1, :)) + DCOS(sg%cell_cntr(2, :))

    ALLOCATE (A3(sg%vLevel, sg%num_cell), &
              parA_sigma3(sg%vLevel, sg%num_cell), &
              parA_sigma_scd3(sg%vLevel, sg%num_cell), &
              parA_sigmat3(sg%vLevel, sg%num_cell), &
              parA_sigma_scdt3(sg%vLevel, sg%num_cell), &
              coef13(sg%vLevel, sg%num_cell, 3))

    !   coef12 = 0.0D0

    DO i = 1, sg%vLevel
      A3(i, :) = dcos((sg%sigma(i)) / 1.0D04)
      parA_sigmat3(i, :) = -1.0D-4 * dsin(sg%sigma(i) / 1.0D04)
      parA_sigma_scdt3(i, :) = -1.0D-8 * dcos(sg%sigma(i) / 1.0D04)
    END DO

    CALL BField_Laplace%initialize(configFile, sg, varName)
    !   coef12 = sg%coef_fstdif

    !   DO i = 1, (sg%vLevel - 1)/2
    !      PRINT *, 'coef_first_1 is ', coef11(i, 100, :)
    !      PRINT *, 'coef_first_2 is ', coef12(i*2 - 1, 100, :)
    !   END DO

    CALL BField_Laplace%FirstOrderDerivative(A3, sg%coef_fstdif, parA_sigma3)
    CALL BField_Laplace%SecondOrderDerivative(A3, sg%coef_scddif, parA_sigma_scd3)

    !   PRINT *, 'parA_sigma and parA_sigma2 are ', size(parA_sigma), size(parA_sigma2)

  END ASSOCIATE

  error_1rst_1 = 0.0D0
  error_2end_1 = 0.0D0
  DO k = 1, geometry%mg%sg(5)%vLevel
    DO i = 1, geometry%mg%sg(5)%num_icell
      IF (geometry%mg%sg(5)%cell_type(i) .EQ. 0) THEN
        IF (error_1rst_1 .LE. DABS(parA_sigma(k, i) - parA_sigmat(k, i))) THEN
          error_1rst_1 = DABS(parA_sigma(k, i) - parA_sigmat(k, i))!/size(parA_sigma)
        END IF
        IF (error_2end_1 .LE. DABS(parA_sigma_scd(k, i) - parA_sigma_scdt(k, i))) THEN
          error_2end_1 = DABS(parA_sigma_scd(k, i) - parA_sigma_scdt(k, i))
        END IF
      END IF
    END DO
  END DO
  error_1rst_2 = 0.0D0
  error_2end_2 = 0.0D0

  DO k = 1, geometry2%mg%sg(5)%vLevel
    DO i = 1, geometry2%mg%sg(5)%num_icell
      IF (geometry2%mg%sg(5)%cell_type(i) .EQ. 0) THEN
        IF (error_1rst_2 .LE. DABS(parA_sigma2(k, i) - parA_sigmat2(k, i))) THEN
          error_1rst_2 = DABS(parA_sigma2(k, i) - parA_sigmat2(k, i))!/size(parA_sigma)
        END IF
        IF (error_2end_2 .LE. DABS(parA_sigma_scd2(k, i) - parA_sigma_scdt2(k, i))) THEN
          error_2end_2 = DABS(parA_sigma_scd2(k, i) - parA_sigma_scdt2(k, i))
        END IF
      END IF
    END DO
  END DO
!    error_1rst_2 = error_1rst_2/2.0D0
!    error_2end_2 = error_2end_2/2.0D0

  error_1rst_3 = 0.0D0
  error_2end_3 = 0.0D0
  DO k = 1, geometry3%mg%sg(5)%vLevel
    DO i = 1, geometry3%mg%sg(5)%num_icell
      IF (geometry3%mg%sg(5)%cell_type(i) .EQ. 0) THEN
        IF (error_1rst_3 .LE. DABS(parA_sigma3(k, i) - parA_sigmat3(k, i))) THEN
          error_1rst_3 = DABS(parA_sigma3(k, i) - parA_sigmat3(k, i))!/size(parA_sigma)
        END IF
        IF (error_2end_3 .LE. DABS(parA_sigma_scd3(k, i) - parA_sigma_scdt3(k, i))) THEN
          error_2end_3 = DABS(parA_sigma_scd3(k, i) - parA_sigma_scdt3(k, i))
        END IF
      END IF
    END DO
  END DO
!    error_1rst_3 = error_1rst_3/4.0D0
!    error_2end_3 = error_2end_3/4.0D0

!    ASSOCIATE (num_icell => geometry%mg%sg(5)%num_icell)

!       PRINT *, 'A1-A2 is', SUM(parA_sigmat(50, :) - parA_sigmat2(99, :)), SUM(A(50, :) - A2(99, :))
!       PRINT *, 'errors in the same layer are', &
!          ((parA_sigmat(50, 1011) - parA_sigma(50, 1011))), &
!          ((parA_sigmat2(99, 1011) - parA_sigma2(99, 1011)))
!    END ASSOCIATE

  error_first_ratio1 = error_1rst_1 / error_1rst_2
  error_second_ratio1 = error_2end_1 / error_2end_2

  error_first_ratio2 = error_1rst_1 / error_1rst_3
  error_second_ratio2 = error_2end_1 / error_2end_3

  PRINT *, 'error_1rst_1 and error_1rst_2 are ', error_1rst_1, error_1rst_2, error_1rst_1 / error_1rst_2, error_first_ratio2
  PRINT *, 'error_2end_1 and error_2end_2 are ', error_2end_1, error_2end_2, error_2end_1 / error_2end_2, error_second_ratio2

  IF ((dabs(error_first_ratio1 - 4.0D0) .LE. 0.05D0) .AND. (dabs(error_second_ratio1 - 4.0D0) .LE. 0.05D0) &
      .AND. (dabs(error_first_ratio2 - 16.0D0) .LE. 0.05D0) .AND. (dabs(error_second_ratio2 - 16.0D0) .LE. 0.05D0)) THEN
    PRINT *, 'Vertical difference tested passed'
  ELSE
    PRINT *, 'Vertical difference tested failed'
  END IF

  !  ASSOCIATE (sg => geometry%mg%sg(6))
  !     sg%topo = 0.0D0
  !     sg%topo = DSIN(sg%cell_cntr(1, :)) + DCOS(sg%cell_cntr(2, :))
  !     CALL BField_Laplace%initialize(configFile, sg, varName)
  !  END ASSOCIATE

  !  ASSOCIATE (sg => geometry%mg%sg(7))
  !     sg%topo = 0.0D0
  !     sg%topo = DSIN(sg%cell_cntr(1, :)) + DCOS(sg%cell_cntr(2, :))
  !     CALL BField_Laplace%initialize(configFile, sg, varName)
  !  END ASSOCIATE

END PROGRAM
