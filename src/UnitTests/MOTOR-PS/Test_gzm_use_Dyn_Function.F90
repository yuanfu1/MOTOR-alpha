PROGRAM test_gzm
  USE kinds_m, ONLY: i_kind, r_kind
  USE gzm_m, ONLY: gzm_t
  USE UnitTestDyn_m, ONLY: UnitTestDyn_t
  USE GenContainers_m, ONLY: GenContainers_t
  USE geometry_m, ONLY: Geometry_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  !  USE SingleGrid_m, ONLY: SingleGrid_t
  USE MGOpts_m
  USE YAMLRead_m
  USE parameters_m, ONLY: Omega

  TYPE(gzm_t) :: gzm
  TYPE(UnitTestDyn_t) :: Rossby
  TYPE(GenContainers_t) :: GenContainers
  TYPE(mpddGlob_t), TARGET :: mpddGlob
  !  TYPE(SingleGrid_t) :: sg
  TYPE(geometry_t) :: geometry

  INTEGER(i_kind) :: i, k
  REAL(r_kind) :: errorv_1, errorv_2, errorv_3, &
                  errorj_1, errorj_2, errorj_3, &
                  errorf_1, errorf_2, errorf_3

  CHARACTER(LEN=1024) :: configFile

  ! Fields to test:
  REAL(r_kind), ALLOCATABLE :: vor1(:, :), stream1(:, :), &
                               thicks1(:), derivS1(:, :), derivV1(:, :), &
                               jacobi1(:, :), flxdiv1(:, :), &
                               streamtt1(:), jacobitt1(:), &
                               flxdivtt1(:), vorticitytt1(:), &
                               vor2(:, :), stream2(:, :), &
                               thicks2(:), derivS2(:, :), derivV2(:, :), &
                               jacobi2(:, :), flxdiv2(:, :), &
                               streamtt2(:), jacobitt2(:), &
                               flxdivtt2(:), vorticitytt2(:), &
                               vor3(:, :), stream3(:, :), &
                               thicks3(:), derivS3(:, :), derivV3(:, :), &
                               jacobi3(:, :), flxdiv3(:, :), &
                               streamtt3(:), jacobitt3(:), &
                               flxdivtt3(:), vorticitytt3(:), &
                               vort1(:, :), jacobit1(:, :), flxdivt1(:, :), etat1(:, :), &
                               vort2(:, :), jacobit2(:, :), flxdivt2(:, :), etat2(:, :), &
                               vort3(:, :), jacobit3(:, :), flxdivt3(:, :), etat3(:, :)

  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", configFile)
  configFile = TRIM(configFile)//"/gzm_and_laplace_terrain.yaml"

  GenContainers = GenContainers_t(TRIM(configFile))
  CALL GenContainers%GenGeometry(geometry)

  ASSOCIATE (sg => geometry%mg%sg(5), &
             num_cell => geometry%mg%sg(5)%num_cell, &
             vLevel => geometry%mg%sg(5)%vLevel)

    ALLOCATE (vor1(vLevel, num_cell), stream1(vLevel, num_cell), &
              thicks1(num_cell), derivS1(2, num_cell), derivV1(2, num_cell), &
              jacobi1(vLevel, num_cell), flxdiv1(vLevel, num_cell), &
              vort1(vLevel, num_cell), jacobit1(vLevel, num_cell), flxdivt1(vLevel, num_cell), &
              streamtt1(num_cell), jacobitt1(num_cell), flxdivtt1(num_cell), vorticitytt1(num_cell), &
              sg%bdy_type(num_cell), etat1(vLevel, num_cell))

    ALLOCATE (sg%f(sg%vLevel, sg%num_cell))

    DO k = 1, sg%vLevel
      sg%f(k, :) = 2.0D0 * Omega * DSIN(sg%cell_cntr(1, :))
    END DO

    ! CALL Rossby%RH_initial(sg%num_cell)

    DO k = 1, vLevel

      CALL Rossby%RHV_func(0.0D0, sg%cell_cntr, streamtt1, vorticitytt1, &
                           derivS1, derivV1, jacobitt1, flxdivtt1, 'chi')
      stream1(k, :) = streamtt1
      vort1(k, :) = vorticitytt1
      jacobit1(k, :) = jacobitt1
      flxdivt1(k, :) = flxdivtt1
    END DO

    PRINT *, 'max of vort1 is ', MAXVAL(dabs(vort1))
    sg%bdy_type = sg%cell_type
    PRINT *, 'cell type is ', sg%bdy_type(1), sg%cell_type(1)
    vor1 = 0.0D0
    gzm = gzm_t(sg)
    CALL gzm%Laplacia(stream1, vor1)
    PRINT *, 'max of vor1 is ', MAXVAL(dabs(vor1)), MAXLOC(dabs(vor1))
    ! CALL gzm%Divergen(vort1 + sg%f)
    PRINT *, 'gzm laplacia over'
    etat1 = vort1 + sg%f
    CALL gzm%Jacobian(etat1, stream1, jacobi1)
    CALL gzm%Divergen(etat1, stream1, flxdiv1)

    ! CALL Rossby%RH_destroy()

  END ASSOCIATE

  ASSOCIATE (sg => geometry%mg%sg(6), &
             num_cell => geometry%mg%sg(6)%num_cell, &
             vLevel => geometry%mg%sg(6)%vLevel)

    ALLOCATE (vor2(vLevel, num_cell), stream2(vLevel, num_cell), &
              thicks2(num_cell), derivS2(2, num_cell), derivV2(2, num_cell), &
              jacobi2(vLevel, num_cell), flxdiv2(vLevel, num_cell), &
              vort2(vLevel, num_cell), jacobit2(vLevel, num_cell), flxdivt2(vLevel, num_cell), &
              streamtt2(num_cell), jacobitt2(num_cell), flxdivtt2(num_cell), vorticitytt2(num_cell), &
              sg%bdy_type(num_cell))

    ALLOCATE (sg%f(sg%vLevel, sg%num_cell))

    DO k = 1, sg%vLevel
      sg%f(k, :) = 2.0D0 * Omega * DSIN(sg%cell_cntr(1, :))
    END DO

    ! CALL Rossby%RH_initial(sg%num_cell)

    DO k = 1, vLevel

      CALL Rossby%RHV_func(0.0D0, sg%cell_cntr, streamtt2, vorticitytt2, &
                           derivS2, derivV2, jacobitt2, flxdivtt2, 'chi')

      stream2(k, :) = streamtt2
      vort2(k, :) = vorticitytt2
      jacobit2(k, :) = jacobitt2
      flxdivt2(k, :) = flxdivtt2

    END DO

    IF (ASSOCIATED(gzm%sg)) NULLIFY (gzm%sg)
    sg%bdy_type = sg%cell_type
    gzm = gzm_t(sg)
    CALL gzm%Laplacia(stream2, vor2)
    etat2 = vort2 + sg%f
    CALL gzm%Jacobian(etat2, stream2, jacobi2)
    CALL gzm%Divergen(etat2, stream2, flxdiv2)

    ! CALL Rossby%RH_destroy()

  END ASSOCIATE

  ASSOCIATE (sg => geometry%mg%sg(7), &
             num_cell => geometry%mg%sg(7)%num_cell, &
             vLevel => geometry%mg%sg(7)%vLevel)

    ALLOCATE (vor3(vLevel, num_cell), stream3(vLevel, num_cell), &
              thicks3(num_cell), derivS3(2, num_cell), derivV3(2, num_cell), &
              jacobi3(vLevel, num_cell), flxdiv3(vLevel, num_cell), &
              vort3(vLevel, num_cell), jacobit3(vLevel, num_cell), flxdivt3(vLevel, num_cell), &
              streamtt3(num_cell), jacobitt3(num_cell), flxdivtt3(num_cell), vorticitytt3(num_cell), &
              sg%bdy_type(num_cell))

    ALLOCATE (sg%f(sg%vLevel, sg%num_cell))

    DO k = 1, sg%vLevel
      sg%f(k, :) = 2.0D0 * Omega * DSIN(sg%cell_cntr(1, :))
    END DO

    ! CALL Rossby%RH_initial(sg%num_cell)

    DO k = 1, vLevel

      CALL Rossby%RHV_func(0.0D0, sg%cell_cntr, streamtt3, vorticitytt3, &
                           derivS3, derivV3, jacobitt3, flxdivtt3, 'chi')

      stream3(k, :) = streamtt3
      vort3(k, :) = vorticitytt3
      jacobit3(k, :) = jacobitt3
      flxdivt3(k, :) = flxdivtt3

    END DO

    IF (ASSOCIATED(gzm%sg)) NULLIFY (gzm%sg)
    sg%bdy_type = sg%cell_type
    gzm = gzm_t(sg)
    CALL gzm%Laplacia(stream3, vor3)
    etat3 = vort3 + sg%f
    CALL gzm%Jacobian(etat3, stream3, jacobi3)
    CALL gzm%Divergen(etat3, stream3, flxdiv3)

    ! CALL Rossby%RH_destroy()

  END ASSOCIATE

  errorv_1 = 0.0D0
  errorf_1 = 0.0D0
  errorj_1 = 0.0D0
  DO k = 1, geometry%mg%sg(5)%vLevel
    DO i = 1, geometry%mg%sg(5)%num_icell
      IF (geometry%mg%sg(5)%cell_type(i) .EQ. 0) THEN
        IF (errorv_1 .LE. DABS(vor1(k, i) - vort1(k, i))) THEN
          errorv_1 = DABS(vor1(k, i) - vort1(k, i))
        END IF

        IF (errorf_1 .LE. DABS(flxdiv1(k, i) - flxdivt1(k, i))) THEN
          errorf_1 = DABS(flxdiv1(k, i) - flxdivt1(k, i))
        END IF

        IF (errorj_1 .LE. DABS(jacobi1(k, i) - jacobit1(k, i))) THEN
          errorj_1 = DABS(jacobi1(k, i) - jacobit1(k, i))
        END IF

      END IF
    END DO
  END DO

  errorv_2 = 0.0D0
  errorf_2 = 0.0D0
  errorj_2 = 0.0D0
  DO k = 1, geometry%mg%sg(6)%vLevel
    DO i = 1, geometry%mg%sg(6)%num_icell
      IF (geometry%mg%sg(6)%cell_type(i) .EQ. 0) THEN
        IF (errorv_2 .LE. DABS(vor2(k, i) - vort2(k, i))) THEN
          errorv_2 = DABS(vor2(k, i) - vort2(k, i))
        END IF

        IF (errorf_2 .LE. DABS(flxdiv2(k, i) - flxdivt2(k, i))) THEN
          errorf_2 = DABS(flxdiv2(k, i) - flxdivt2(k, i))
        END IF

        IF (errorj_2 .LE. DABS(jacobi2(k, i) - jacobit2(k, i))) THEN
          errorj_2 = DABS(jacobi2(k, i) - jacobit2(k, i))
        END IF

      END IF
    END DO
  END DO

  errorv_3 = 0.0D0
  errorf_3 = 0.0D0
  errorj_3 = 0.0D0
  DO k = 1, geometry%mg%sg(7)%vLevel
    DO i = 1, geometry%mg%sg(7)%num_icell
      IF (geometry%mg%sg(7)%cell_type(i) .EQ. 0) THEN
        IF (errorv_3 .LE. DABS(vor3(k, i) - vort3(k, i))) THEN
          errorv_3 = DABS(vor3(k, i) - vort3(k, i))
        END IF

        IF (errorf_3 .LE. DABS(flxdiv3(k, i) - flxdivt3(k, i))) THEN
          errorf_3 = DABS(flxdiv3(k, i) - flxdivt3(k, i))
        END IF

        IF (errorj_3 .LE. DABS(jacobi3(k, i) - jacobit3(k, i))) THEN
          errorj_3 = DABS(jacobi3(k, i) - jacobit3(k, i))
        END IF

      END IF
    END DO
  END DO

  PRINT *, 'error_v ratio is ', errorv_1, errorv_2, errv_3, errorv_1 / errorv_2, errorv_2 / errorv_3
  PRINT *, 'error_f ratio is ', errorf_1, errorf_2, errf_3, errorf_1 / errorf_2, errorf_2 / errorf_3
  PRINT *, 'error_j ratio is ', errorj_1, errorj_2, errj_3, errorj_1 / errorj_2, errorj_2 / errorj_3

  DEALLOCATE (vor1, stream1, &
              thicks1, derivS1, derivV1, &
              jacobi1, flxdiv1, &
              streamtt1, jacobitt1, &
              flxdivtt1, vorticitytt1, &
              vor2, stream2, &
              thicks2, derivS2, derivV2, &
              jacobi2, flxdiv2, &
              streamtt2, jacobitt2, &
              flxdivtt2, vorticitytt2, &
              vor3, stream3, &
              thicks3, derivS3, derivV3, &
              jacobi3, flxdiv3, &
              streamtt3, jacobitt3, &
              flxdivtt3, vorticitytt3, &
              vort1, jacobit1, flxdivt1, &
              vort2, jacobit2, flxdivt2, &
              vort3, jacobit3, flxdivt3)

END PROGRAM test_gzm
