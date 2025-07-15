!!---------------------------------------------------------------------------------------
! PROJECT           : PRIVATE-PROJECT
! AFFILIATION       : Self-employed
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : V 1
! HISTORY           :
!   Created  by Yuanfu Xie, 2021 Broomfield, CO, USA
!   Modified by
!!---------------------------------------------------------------------------------------

!> @brief
!! Testing single grid structures.
!! @author Yuanfu Xie
!! @copyright (C) 2021 Yuanfu Xie, All rights reserved.
!! @note:
!! @warning
!! @attention

PROGRAM test_gzm
  USE kinds_m, ONLY: i_kind, r_kind
  USE gzm_m, ONLY: gzm_t
  USE RossbyHaurwitz_m, ONLY: RossbyHaurwitzSphere1_t
  ! USE GenContainers_m, ONLY: GenContainers_t
  ! USE geometry_m, ONLY: Geometry_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  !  USE SingleGrid_m, ONLY: SingleGrid_t
  USE State_m, ONLY: State_t
  USE State2NC_m
  USE MGOpts_m
  USE YAMLRead_m

  TYPE(gzm_t) :: gzm
  TYPE(RossbyHaurwitzSphere1_t) :: Rossby
  ! TYPE(GenContainers_t) :: GenContainers
  TYPE(mpddGlob_t), TARGET :: mpddGlob
  !  TYPE(SingleGrid_t) :: sg
  ! TYPE(geometry_t) :: geometry
  TYPE(State_t) :: X1, X2

  INTEGER(i_kind) :: i, k
  REAL(r_kind) :: error_1, error_2

  CHARACTER(LEN=1024) :: configFile, outputDir

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
                               vort1(:, :), jacobit1(:, :), flxdivt1(:, :), &
                               vort2(:, :), jacobit2(:, :), flxdivt2(:, :), &
                               vort3(:, :), jacobit3(:, :), flxdivt3(:, :)

  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", configFile)
  configFile = TRIM(configFile)//"/gzm_and_laplace_terrain.yaml"

  IF (yaml_get_var(configFile, 'IO', 'output_dir', outputDir) /= 0) STOP
  CALL GET_ENVIRONMENT_VARIABLE("OUTPUT_DIR", outputDir)

  CALL gzm%initial(configFile)

  CALL X1%initialize(configFile, gzm%geometry%mg%sg(5))

  CALL X1%addVar('stream1')
  CALL X1%addVar('vor1')

  ASSOCIATE (sg => gzm%geometry%mg%sg(5), &
             num_cell => gzm%geometry%mg%sg(5)%num_cell, &
             vLevel => gzm%geometry%mg%sg(5)%vLevel)

    ALLOCATE (vor1(vLevel, num_cell), stream1(vLevel, num_cell), &
              thicks1(num_cell), derivS1(2, num_cell), derivV1(2, num_cell), &
              jacobi1(vLevel, num_cell), flxdiv1(vLevel, num_cell), &
              vort1(vLevel, num_cell), jacobit1(vLevel, num_cell), flxdivt1(vLevel, num_cell), &
              streamtt1(num_cell), jacobitt1(num_cell), flxdivtt1(num_cell), vorticitytt1(num_cell), &
              sg%bdy_type(num_cell))
    PRINT *, 'size of laton are', SIZE(sg%cell_cntr(1, :)), SIZE(sg%cell_cntr(:, 1))
    PRINT *, 'max of lat and lon are', MINVAL(sg%cell_cntr(1, :)) * 180 / 3.14, MINVAL(sg%cell_cntr(2, :)) * 180 / 3.14
    CALL Rossby%RH_initial(sg%num_cell)

    DO k = 1, vLevel

      PRINT *, 'k is', k

      CALL Rossby%RHV_sphereFunc(0.0D0, sg%cell_cntr, streamtt1, vorticitytt1, &
                                 thicks1, derivS1, derivV1, jacobitt1, flxdivtt1)
      stream1(k, :) = streamtt1
      vort1(k, :) = vorticitytt1
      jacobit1(k, :) = jacobitt1
      flxdivt1(k, :) = flxdivtt1
    END DO

    sg%bdy_type = sg%cell_type
    PRINT *, 'cell type is ', sg%bdy_type(1), sg%cell_type(1)
    ! CALL gzm%Jacobian()

    CALL Rossby%RH_destroy()

  END ASSOCIATE

  PRINT *, 'min of stream1 is ', MAXVAL((stream1)), MINVAL((stream1)), SIZE(stream1)
  PRINT *, 'min of vort1 is ', MAXVAL((vort1)), MINVAL((vort1))

  vor1 = 0.0D0
  !   gzm = gzm_t(sg)
  CALL gzm%Laplacia(5, stream1, vor1)
  PRINT *, 'max of vor1 is ', MAXVAL(dabs(vor1)), MAXLOC(dabs(vor1))
  PRINT *, 'gzm laplacia over'

  X1%Fields(X1%getVarIdx('stream1'))%DATA(:, :, 1) = stream1
  X1%Fields(X1%getVarIdx('vor1'))%DATA(:, :, 1) = vort1
  X1%Fields(X1%getVarIdx('vor1'))%DATA(:, :, 2) = vor1

  ASSOCIATE (sg => gzm%geometry%mg%sg(6), &
             num_cell => gzm%geometry%mg%sg(6)%num_cell, &
             vLevel => gzm%geometry%mg%sg(6)%vLevel)

    ALLOCATE (vor2(vLevel, num_cell), stream2(vLevel, num_cell), &
              thicks2(num_cell), derivS2(2, num_cell), derivV2(2, num_cell), &
              jacobi2(vLevel, num_cell), flxdiv2(vLevel, num_cell), &
              vort2(vLevel, num_cell), jacobit2(vLevel, num_cell), flxdivt2(vLevel, num_cell), &
              streamtt2(num_cell), jacobitt2(num_cell), flxdivtt2(num_cell), vorticitytt2(num_cell), &
              sg%bdy_type(num_cell))

    CALL Rossby%RH_initial(sg%num_cell)

    DO k = 1, vLevel

      CALL Rossby%RHV_func(0.0D0, sg%cell_cntr, streamtt2, vorticitytt2, &
                           thicks2, derivS2, derivV2, jacobitt2, flxdivtt2)

      stream2(k, :) = streamtt2
      vort2(k, :) = vorticitytt2
      jacobit2(k, :) = jacobitt2
      flxdivt2(k, :) = flxdivtt2

    END DO

    !   IF (ASSOCIATED(gzm%sg)) NULLIFY (gzm%sg)
    !   sg%bdy_type = sg%cell_type
    !   gzm = gzm_t(sg)
    ! CALL gzm%Jacobian()

    CALL Rossby%RH_destroy()

  END ASSOCIATE

  CALL gzm%Laplacia(6, stream2, vor2)

  error_1 = 0.0D0
  !  DO k = 1, geometry%mg%sg(5)%vLevel
  DO i = 1, gzm%geometry%mg%sg(5)%num_icell
    IF (gzm%geometry%mg%sg(5)%cell_type(i) .EQ. 0) THEN
      IF (error_1 .LE. DABS(vor1(1, i) - vort1(1, i))) THEN
        error_1 = DABS(vor1(1, i) - vort1(1, i))
      END IF
    END IF
  END DO
  !  END DO

  error_2 = 0.0D0
  !  DO k = 1, geometry%mg%sg(6)%vLevel
  DO i = 1, gzm%geometry%mg%sg(6)%num_icell
    IF (gzm%geometry%mg%sg(6)%cell_type(i) .EQ. 0) THEN
      IF (error_2 .LE. DABS(vor2(1, i) - vort2(1, i))) THEN
        error_2 = DABS(vor2(1, i) - vort2(1, i))
      END IF
    END IF
  END DO
  !  END DO

  PRINT *, 'error ratio is ', error_1, error_2, error_1 / error_2

  CALL Output_NC_State_AV(X1, outputDir, "Test_gzm")

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
