PROGRAM test_laplacian_adjoint
  USE kinds_m, ONLY: i_kind, r_kind
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE geometry_m, ONLY: geometry_t
  USE gzm_Xie_m, ONLY: gzm_Xie_t
  USE State_m, ONLY: State_t
  USE parameters_m, ONLY: pi

  IMPLICIT NONE

  ! Local variables
  CHARACTER(LEN=1024) :: yamlFile
  TYPE(mpddGlob_t) :: mpddGlob
  TYPE(geometry_t) :: geo
  TYPE(gzm_Xie_t) :: gzm
  INTEGER(i_kind) :: i, ic, glevel
  REAL(r_kind), ALLOCATABLE :: u(:,:), v(:,:), w(:,:), Lu(:,:), Lstarv(:,:), Ju(:,:), Jstar1v(:,:), Jstar2u(:,:)
  REAL(r_kind) :: inner_prod1, inner_prod2, rel_diff

  ! Get yaml file from command line
  CALL getarg(1, yamlFile)
  IF (TRIM(yamlFile) .EQ. '') THEN
    PRINT*, 'Usage: <executable> <yaml file>'
    STOP
  END IF

  ! Initialize MPDD and geometry
  CALL mpddGlob%initialize()
  CALL geo%initialize(yamlFile, mpddGlob)

  PRINT*,'Geometry initialized with yaml file: ', TRIM(yamlFile)
  PRINT*, 'Finest grid level: ', geo%mg%mg_finest

  ! Test at finest grid level
  glevel = geo%mg%mg_finest
  
  ASSOCIATE(sg => geo%mg%sg(glevel))
    ! Initialize gzm for this grid level
    PRINT*,'Initializing GZM for grid level: ', glevel,sg%num_cell
    CALL gzm%initial(yamlFile, mpddGlob)
    ! gzm = gzm_t(sg)
    PRINT*, 'GZM initialized for grid level: ', glevel

    ! Allocate arrays
    ALLOCATE(u(sg%vLevel, sg%num_cell))
    ALLOCATE(v(sg%vLevel, sg%num_cell))
    ALLOCATE(w(sg%vLevel, sg%num_cell))
    ALLOCATE(Lu(sg%vLevel, sg%num_cell))
    ALLOCATE(Lstarv(sg%vLevel, sg%num_cell))
    u = 0.0_r_kind
    v = 0.0_r_kind
    w = 0.0_r_kind
    Lu = 0.0_r_kind
    Lstarv = 0.0_r_kind

    ! Initialize test functions (using sinusoidal functions)
    DO ic = 1, sg%num_cell
      IF (sg%cell_type(ic) == 0) THEN  ! Interior points only
        u(:,ic) = SIN(2.0_r_kind*pi*sg%cell_cntr(1,ic)) * &
                  COS(2.0_r_kind*pi*sg%cell_cntr(2,ic))
        v(:,ic) = COS(2.0_r_kind*pi*sg%cell_cntr(1,ic)) * &
                  SIN(2.0_r_kind*pi*sg%cell_cntr(2,ic))
        w(:,ic) = SIN(2.0_r_kind*pi*sg%cell_cntr(1,ic)) * &
                  SIN(2.0_r_kind*pi*sg%cell_cntr(2,ic))
      END IF
    END DO

    ! Apply Laplacian and its adjoint
    CALL gzm%Laplacia(glevel, u, Lu)
    CALL gzm%LaplaciaAdjoint(glevel, v, Lstarv)

    ! Compute inner products <Lu,v> and <u,L*v>
    inner_prod1 = 0.0_r_kind
    inner_prod2 = 0.0_r_kind
    DO ic = 1, sg%num_cell
      IF (sg%cell_type(ic) == 0) THEN
        inner_prod1 = inner_prod1 + DOT_PRODUCT(Lu(:,ic), v(:,ic)) ! * sg%cell_area(ic)
        inner_prod2 = inner_prod2 + DOT_PRODUCT(u(:,ic), Lstarv(:,ic)) ! * sg%cell_area(ic)
      END IF
    END DO

    ! Check adjoint property
    rel_diff = ABS(inner_prod1 - inner_prod2) / &
               (ABS(inner_prod1) + ABS(inner_prod2)) * 2.0_r_kind

    PRINT '(A,E14.6)', 'Relative difference between inner products: ', rel_diff
    IF (rel_diff < 1.0E-10_r_kind) THEN
      PRINT *, 'Adjoint test PASSED'
    ELSE
      PRINT *, 'Adjoint test FAILED'
      PRINT '(A,E14.6)', ' <Lu,v> = ', inner_prod1
      PRINT '(A,E14.6)', '<u,L*v> = ', inner_prod2
    END IF

    ! Jacobian adjoint tests
    ALLOCATE(Ju(sg%vLevel, sg%num_cell))
    ALLOCATE(Jstar1v(sg%vLevel, sg%num_cell))
    ALLOCATE(Jstar2u(sg%vLevel, sg%num_cell))
    Ju = 0.0_r_kind
    Jstar1v = 0.0_r_kind
    Jstar2u = 0.0_r_kind

    ! Test first adjoint: <J(u,v),w> = <u,J*₁(v,w)>
    PRINT '(A)', NEW_LINE('a')//'Testing first Jacobian adjoint (J*₁):'
    CALL gzm%Jacobian(glevel, u, v, Ju)
    CALL gzm%JacobianAdjoint1(glevel, v, w, Jstar1v)

    inner_prod1 = 0.0_r_kind
    inner_prod2 = 0.0_r_kind
    DO ic = 1, sg%num_cell
      IF (sg%cell_type(ic) == 0) THEN
        inner_prod1 = inner_prod1 + DOT_PRODUCT(Ju(:,ic), w(:,ic))
        inner_prod2 = inner_prod2 + DOT_PRODUCT(u(:,ic), Jstar1v(:,ic))
      END IF
    END DO

    rel_diff = ABS(inner_prod1 - inner_prod2) / &
               (ABS(inner_prod1) + ABS(inner_prod2)) * 2.0_r_kind

    PRINT '(A,E14.6)', ' <J(u,v),w> = ', inner_prod1
    PRINT '(A,E14.6)', '<u,J*₁(v,w)> = ', inner_prod2
    PRINT '(A,E14.6)', 'Relative difference: ', rel_diff

    ! Test second adjoint: <J(u,v),w> = <v,J*₂(u,w)>
    PRINT '(A)', NEW_LINE('a')//'Testing second Jacobian adjoint (J*₂):'
    CALL gzm%Jacobian(glevel, u, v, Ju)
    CALL gzm%JacobianAdjoint2(glevel, u, w, Jstar2u)

    inner_prod1 = 0.0_r_kind
    inner_prod2 = 0.0_r_kind
    DO ic = 1, sg%num_cell
      IF (sg%cell_type(ic) == 0) THEN
        inner_prod1 = inner_prod1 + DOT_PRODUCT(Ju(:,ic), w(:,ic))
        inner_prod2 = inner_prod2 + DOT_PRODUCT(v(:,ic), Jstar2u(:,ic))
      END IF
    END DO

    rel_diff = ABS(inner_prod1 - inner_prod2) / &
               (ABS(inner_prod1) + ABS(inner_prod2)) * 2.0_r_kind

    PRINT '(A,E14.6)', ' <J(u,v),w> = ', inner_prod1
    PRINT '(A,E14.6)', '<v,J*₂(u,w)> = ', inner_prod2
    PRINT '(A,E14.6)', 'Relative difference: ', rel_diff

    ! Cleanup
    DEALLOCATE(u, v, w, Lu, Lstarv, Ju, Jstar1v, Jstar2u)
  END ASSOCIATE

  ! Finalize
  CALL mpddGlob%finalize()

END PROGRAM test_laplacian_adjoint
