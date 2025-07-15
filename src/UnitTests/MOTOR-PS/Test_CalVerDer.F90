PROGRAM SgDerivative_test
  USE GenContainers_m, ONLY: GenContainers_t
  USE geometry_m, ONLY: Geometry_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE GenSgDiffCoef_m, ONLY: GenSgDiffCoef_t
  USE CalVerDer_m, ONLY: CalVerDer_t
  USE kinds_m, ONLY: i_kind, r_kind
  USE YAMLRead_m

  IMPLICIT NONE

  TYPE(GenContainers_t) :: GenContainers
  ! TYPE(SingleGrid_t) :: sg
  TYPE(GenSgDiffCoef_t) :: GenSgDiffCoef
  TYPE(CalVerDer_t) :: CalVerDer
  TYPE(geometry_t) :: geometry, geometry2, geometry3
  CHARACTER(LEN=1024) :: configFile
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

  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", configFile)
  configFile = TRIM(configFile)//"/UnitTest/test_Cum_Grapes.yaml"

  GenContainers = GenContainers_t(TRIM(configFile))
  CALL GenContainers%GenGeometry(geometry)

  ASSOCIATE (sg => geometry%mg%sg(5))
    ALLOCATE (A(sg%vLevel, sg%num_cell), &
              parA_sigma(sg%vLevel, sg%num_cell), &
              parA_sigma_scd(sg%vLevel, sg%num_cell), &
              parA_sigmat(sg%vLevel, sg%num_cell), &
              parA_sigma_scdt(sg%vLevel, sg%num_cell))
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

END PROGRAM
