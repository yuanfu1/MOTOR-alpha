PROGRAM SgDerivative_test
  USE GenContainers_m, ONLY: GenContainers_t
  USE geometry_m, ONLY: Geometry_t
  ! USE SingleGrid_m, ONLY: SingleGrid_t
  USE InterpValue_m, ONLY: InterpValue_t
  USE kinds_m, ONLY: i_kind, r_kind
  USE YAMLRead_m

  IMPLICIT NONE

  TYPE(GenContainers_t) :: GenContainers
  ! TYPE(SingleGrid_t) :: sg
  TYPE(InterpValue_t) :: InterpValue
  TYPE(geometry_t) :: geometry, geometry2, geometry3, geometry4, geometry5, geometry6
  CHARACTER(LEN=1024) :: configFile
  INTEGER(i_kind) :: vLevel, i, k
  REAL(r_kind) :: error_1, error_2, error_3, &
                  error_ratio1, error_ratio2
  REAL(r_kind), DIMENSION(:), ALLOCATABLE :: sigma2, sigma3, sigma4, sigma5, sigma6, errorl_1(:), errorl_2(:), errorl_3(:)
  REAL(r_kind), DIMENSION(:, :), ALLOCATABLE:: A, A2, A3, Apt1, Apt2, Apt3, &
                                               Ap1, Ap2, Ap3

  ! CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", configFile)
  ! configFile = TRIM(configFile)//"/gzm_and_laplace_terrain.yaml"
  configFile = "/Users/jilongchen/research/MOTOR-DA/MOTOR/static/gzm_and_laplace_terrain.yaml"

  GenContainers = GenContainers_t(TRIM(configFile))
  CALL GenContainers%GenGeometry(geometry)

  PRINT *, 'geometry is generated'

  ASSOCIATE (sg => geometry%mg%sg(5))
    ALLOCATE (Ap1(sg%vLevel, sg%num_cell), &
              Apt1(sg%vLevel, sg%num_cell))

    DO i = 1, sg%vLevel
      Apt1(i, :) = dcos((sg%sigma(i)) / 1.0D04)
    END DO

    vLevel = sg%vLevel - 1

    ALLOCATE (sigma2(vLevel))

    sigma2(1:vLevel) = (sg%sigma(1:sg%vLevel - 1) &
                        + sg%sigma(2:sg%vLevel)) / 2.0D0

    vLevel = (sg%vLevel - 1) * 2 + 1

    ALLOCATE (sigma3(vLevel))

    sigma3(1:vLevel:2) = sg%sigma
    sigma3(2:vLevel:2) = (sg%sigma(1:sg%vLevel - 1) &
                          + sg%sigma(2:sg%vLevel)) / 2.0D0

  END ASSOCIATE

  CALL GenContainers%GenGeometry(geometry2, sigma2)

  ASSOCIATE (sg => geometry2%mg%sg(5))
    ALLOCATE (A(sg%vLevel, sg%num_cell))

    ! PRINT *, 'sigma is ', sg%sigma
    DO i = 1, sg%vLevel
      A(i, :) = dcos((sg%sigma(i)) / 1.0D04)
    END DO

    InterpValue = InterpValue_t(geometry%mg%sg(5))
    CALL InterpValue%UpdateSgInterpCoef(sg, geometry%mg%sg(5))

    PRINT *, 'interpvalue%sg%vLevel is', InterpValue%sg%vLevel

    CALL InterpValue%InterpValue(A, Ap1)

    ! vLevel = (sg%vLevel - 1)*2 + 1

    ! ALLOCATE (sigma4(vLevel))

    ! sigma4(1:vLevel:2) = sg%sigma
    ! sigma4(2:vLevel:2) = (sg%sigma(1:sg%vLevel - 1) &
    !                       + sg%sigma(2:sg%vLevel))/2.0D0

  END ASSOCIATE

  CALL GenContainers%GenGeometry(geometry3, sigma3)

  ASSOCIATE (sg => geometry3%mg%sg(5))
    ALLOCATE (Ap2(sg%vLevel, sg%num_cell), &
              Apt2(sg%vLevel, sg%num_cell))

    DO i = 1, sg%vLevel
      Apt2(i, :) = dcos((sg%sigma(i)) / 1.0D04)
    END DO

    vLevel = sg%vLevel - 1

    ALLOCATE (sigma4(vLevel))

    sigma4(1:vLevel) = (sg%sigma(1:sg%vLevel - 1) &
                        + sg%sigma(2:sg%vLevel)) / 2.0D0

    vLevel = (sg%vLevel - 1) * 2 + 1

    ALLOCATE (sigma5(vLevel))

    sigma5(1:vLevel:2) = sg%sigma
    sigma5(2:vLevel:2) = (sg%sigma(1:sg%vLevel - 1) &
                          + sg%sigma(2:sg%vLevel)) / 2.0D0

  END ASSOCIATE

  CALL GenContainers%GenGeometry(geometry4, sigma4)

  ASSOCIATE (sg => geometry4%mg%sg(5))
    ALLOCATE (A2(sg%vLevel, sg%num_cell))

    ! PRINT *, 'sigma is ', sg%sigma
    DO i = 1, sg%vLevel
      A2(i, :) = dcos((sg%sigma(i)) / 1.0D04)
    END DO

    IF (ASSOCIATED(InterpValue%sg)) NULLIFY (InterpValue%sg)

    InterpValue = InterpValue_t(geometry3%mg%sg(5))
    CALL InterpValue%UpdateSgInterpCoef(sg, geometry3%mg%sg(5))

    PRINT *, 'interpvalue%sg%vLevel is', InterpValue%sg%vLevel

    CALL InterpValue%InterpValue(A2, Ap2)

    ! vLevel = (sg%vLevel - 1)*2 + 1

    ! ALLOCATE (sigma6(vLevel))

    ! sigma6(1:vLevel:2) = sg%sigma
    ! sigma6(2:vLevel:2) = (sg%sigma(1:sg%vLevel - 1) &
    !                       + sg%sigma(2:sg%vLevel))/2.0D0

  END ASSOCIATE

  CALL GenContainers%GenGeometry(geometry5, sigma5)

  ASSOCIATE (sg => geometry5%mg%sg(5))
    ALLOCATE (Ap3(sg%vLevel, sg%num_cell), &
              Apt3(sg%vLevel, sg%num_cell))

    DO i = 1, sg%vLevel
      Apt3(i, :) = dcos((sg%sigma(i)) / 1.0D04)
    END DO

    vLevel = sg%vLevel - 1

    ALLOCATE (sigma6(vLevel))

    sigma6(1:vLevel) = (sg%sigma(1:sg%vLevel - 1) &
                        + sg%sigma(2:sg%vLevel)) / 2.0D0

  END ASSOCIATE

  CALL GenContainers%GenGeometry(geometry6, sigma6)

  ASSOCIATE (sg => geometry6%mg%sg(5))
    ALLOCATE (A3(sg%vLevel, sg%num_cell))

    ! PRINT *, 'sigma is ', sg%sigma
    DO i = 1, sg%vLevel
      A3(i, :) = dcos((sg%sigma(i)) / 1.0D04)
    END DO

    IF (ASSOCIATED(InterpValue%sg)) NULLIFY (InterpValue%sg)

    InterpValue = InterpValue_t(geometry5%mg%sg(5))
    CALL InterpValue%UpdateSgInterpCoef(sg, geometry5%mg%sg(5))

    PRINT *, 'interpvalue%sg%vLevel is', InterpValue%sg%vLevel

    CALL InterpValue%InterpValue(A3, Ap3)

  END ASSOCIATE

  ALLOCATE (errorl_1(geometry%mg%sg(5)%vLevel), &
            errorl_2(geometry3%mg%sg(5)%vLevel), &
            errorl_3(geometry5%mg%sg(5)%vLevel))

  error_1 = 0.0D0
  DO k = 1, geometry%mg%sg(5)%vLevel
    DO i = 1, geometry%mg%sg(5)%num_icell
      IF (geometry%mg%sg(5)%cell_type(i) .EQ. 0) THEN
        IF (error_1 .LE. DABS(Ap1(k, i) - Apt1(k, i))) THEN
          error_1 = DABS(Ap1(k, i) - Apt1(k, i))!/size(parA_sigma)
        END IF
      END IF
    END DO
    errorl_1(k) = error_1
  END DO

  error_2 = 0.0D0
  DO k = 1, geometry3%mg%sg(5)%vLevel
    DO i = 1, geometry3%mg%sg(5)%num_icell
      IF (geometry3%mg%sg(5)%cell_type(i) .EQ. 0) THEN
        IF (error_2 .LE. DABS(Ap2(k, i) - Apt2(k, i))) THEN
          error_2 = DABS(Ap2(k, i) - Apt2(k, i))!/size(parA_sigma)
        END IF
      END IF
    END DO
    errorl_2(k) = error_2
  END DO

  error_3 = 0.0D0
  DO k = 1, geometry5%mg%sg(5)%vLevel
    DO i = 1, geometry5%mg%sg(5)%num_icell
      IF (geometry5%mg%sg(5)%cell_type(i) .EQ. 0) THEN
        IF (error_3 .LE. DABS(Ap3(k, i) - Apt3(k, i))) THEN
          error_3 = DABS(Ap3(k, i) - Apt3(k, i))!/size(parA_sigma)
        END IF
      END IF
    END DO
    errorl_3(k) = error_3
  END DO

  PRINT *, 'error is ', error_1, error_2, error_3, error_1 / error_2, error_2 / error_3

  PRINT *, 'error by level '
  DO k = 1, geometry%mg%sg(5)%vLevel
    PRINT *, errorl_1(k) / errorl_2(2 * k - 1)
  END DO
  PRINT *, 'error by level d2'
  DO k = 1, geometry3%mg%sg(5)%vLevel
    PRINT *, errorl_2(k) / errorl_3(2 * k - 1)
  END DO

END PROGRAM
