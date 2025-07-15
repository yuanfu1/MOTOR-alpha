MODULE amg2d_m
  USE, INTRINSIC :: iso_c_binding, ONLY: c_double, c_int
  USE ilu_precond_m, ONLY: ilu_t, create_ilu, apply_ilu  ! Add ILU module
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: amg2d_t, amg_level_t  ! Made amg_level_t public

  !> Parameters for AMG
  INTEGER(c_int), PARAMETER :: MAX_LEVELS = 10
  REAL(c_double), PARAMETER :: COARSE_THRESHOLD = 0.25D0
  REAL(c_double), PARAMETER :: STRONG_THRESHOLD = 0.25D0

  !> AMG grid level type
  TYPE :: amg_level_t
    INTEGER(c_int) :: n                    ! Number of unknowns at this level
    INTEGER(c_int), ALLOCATABLE :: ia(:)   ! Row pointers
    INTEGER(c_int), ALLOCATABLE :: ja(:)   ! Column indices
    REAL(c_double), ALLOCATABLE :: a(:)    ! Matrix values
    REAL(c_double), ALLOCATABLE :: x(:)    ! Solution vector
    REAL(c_double), ALLOCATABLE :: b(:)    ! Right-hand side
    REAL(c_double), ALLOCATABLE :: res(:)  ! Residual vector (renamed from r)
    INTEGER(c_int), ALLOCATABLE :: cf(:)   ! Coarse/Fine splitting
    REAL(c_double), ALLOCATABLE :: P(:)    ! Prolongation operator
    REAL(c_double), ALLOCATABLE :: R(:)    ! Restriction operator
    TYPE(ilu_t) :: ilu_smoother    ! Add ILU smoother to each level
  END TYPE amg_level_t

  !> Main AMG solver type
  TYPE :: amg2d_t
    INTEGER(c_int) :: num_levels           ! Number of multigrid levels
    TYPE(amg_level_t), ALLOCATABLE :: levels(:)
    REAL(c_double) :: tol                  ! Convergence tolerance
    INTEGER(c_int) :: max_iter            ! Maximum iterations
    INTEGER(c_int) :: cycle_type          ! 1=V-cycle, 2=W-cycle
    
    CONTAINS
      PROCEDURE :: initialize => initialize_s
      PROCEDURE :: setup => setup_s
      PROCEDURE :: solve => solve_s
      PROCEDURE :: destroy => destroy_s
      PROCEDURE, PRIVATE :: coarsen => coarsen_s
      PROCEDURE, PRIVATE :: interpolate => interpolate_s
      PROCEDURE, PRIVATE :: smooth => smooth_s
      PROCEDURE, PRIVATE :: cycle => cycle_s
  END TYPE amg2d_t

CONTAINS
  ! Helper subroutines
  SUBROUTINE select_coarse_points(n, is_strong, cf)
    INTEGER(c_int), INTENT(IN) :: n
    LOGICAL, INTENT(IN) :: is_strong(:,:)
    INTEGER(c_int), INTENT(OUT) :: cf(:)
    INTEGER(c_int) :: i, j, max_connections, max_i
    INTEGER(c_int), ALLOCATABLE :: num_connections(:)
    LOGICAL, ALLOCATABLE :: is_processed(:)
    
    ! Initialize arrays
    ALLOCATE(num_connections(n), is_processed(n))
    num_connections = 0
    is_processed = .FALSE.
    cf = 0  ! 0=fine, 1=coarse

    ! Count strong connections for each point
    DO i = 1, n
      DO j = 1, n
        IF (is_strong(i,j)) num_connections(i) = num_connections(i) + 1
      END DO
    END DO

    ! Select coarse points using Ruge-Stuben coarsening
    DO WHILE (.NOT. ALL(is_processed))
      ! Find unprocessed point with maximum connections
      max_connections = -1
      max_i = 0
      DO i = 1, n
        IF (.NOT. is_processed(i) .AND. num_connections(i) > max_connections) THEN
          max_connections = num_connections(i)
          max_i = i
        END IF
      END DO
      
      IF (max_i == 0) EXIT  ! No more points to process

      ! Make this point a coarse point
      cf(max_i) = 1
      is_processed(max_i) = .TRUE.

      ! Make strongly connected neighbors fine points
      DO j = 1, n
        IF (is_strong(max_i,j) .AND. .NOT. is_processed(j)) THEN
          cf(j) = 0
          is_processed(j) = .TRUE.
        END IF
      END DO
    END DO

    DEALLOCATE(num_connections, is_processed)
  END SUBROUTINE select_coarse_points

  SUBROUTINE compute_residual(level)
    TYPE(amg_level_t), INTENT(INOUT) :: level
    INTEGER(c_int) :: i, j
    REAL(c_double) :: temp

    level%res = level%b  ! Start with right-hand side
    
    ! Subtract A*x from b
    DO i = 1, level%n
      temp = 0.0_c_double
      DO j = level%ia(i), level%ia(i+1)-1
        temp = temp + level%a(j) * level%x(level%ja(j))
      END DO
      level%res(i) = level%res(i) - temp
    END DO
  END SUBROUTINE compute_residual

  SUBROUTINE compute_residual_norm(level, norm)
    TYPE(amg_level_t), INTENT(IN) :: level
    REAL(c_double), INTENT(OUT) :: norm
    INTEGER(c_int) :: i
    
    norm = 0.0_c_double
    DO i = 1, level%n
      norm = norm + level%res(i) * level%res(i)
    END DO
    norm = SQRT(norm)
  END SUBROUTINE compute_residual_norm

  SUBROUTINE restrict_residual(fine, coarse)
    TYPE(amg_level_t), INTENT(IN) :: fine
    TYPE(amg_level_t), INTENT(INOUT) :: coarse
    INTEGER(c_int) :: i, j
    
    ! Initialize coarse grid residual to zero
    coarse%b = 0.0_c_double
    
    ! Apply restriction operator R: coarse = R * fine
    DO i = 1, fine%n
      DO j = fine%ia(i), fine%ia(i+1)-1
        IF (fine%cf(fine%ja(j)) == 1) THEN  ! If it's a coarse point
          coarse%b(fine%ja(j)) = coarse%b(fine%ja(j)) + fine%R(j) * fine%res(i)
        END IF
      END DO
    END DO
  END SUBROUTINE restrict_residual

  SUBROUTINE direct_solve(level)
    TYPE(amg_level_t), INTENT(INOUT) :: level
    INTEGER(c_int) :: i, j, k
    REAL(c_double) :: sum
    REAL(c_double), ALLOCATABLE :: L(:,:), U(:,:), y(:)
    
    ! For small systems, use direct LU factorization
    ALLOCATE(L(level%n,level%n), U(level%n,level%n), y(level%n))
    L = 0.0_c_double
    U = 0.0_c_double
    
    ! Construct L and U matrices from CSR format
    DO i = 1, level%n
      DO j = level%ia(i), level%ia(i+1)-1
        IF (level%ja(j) <= i) THEN
          L(i,level%ja(j)) = level%a(j)
        ELSE
          U(i,level%ja(j)) = level%a(j)
        END IF
      END DO
    END DO
    
    ! Forward solve Ly = b
    DO i = 1, level%n
      sum = level%b(i)
      DO j = 1, i-1
        sum = sum - L(i,j) * y(j)
      END DO
      y(i) = sum / L(i,i)
    END DO
    
    ! Backward solve Ux = y
    DO i = level%n, 1, -1
      sum = y(i)
      DO j = i+1, level%n
        sum = sum - U(i,j) * level%x(j)
      END DO
      level%x(i) = sum / U(i,i)
    END DO
    
    DEALLOCATE(L, U, y)
  END SUBROUTINE direct_solve

  SUBROUTINE prolongate_and_correct(fine, coarse)
    TYPE(amg_level_t), INTENT(INOUT) :: fine
    TYPE(amg_level_t), INTENT(IN) :: coarse
    INTEGER(c_int) :: i, j
    REAL(c_double) :: correction
    
    ! For each fine point
    DO i = 1, fine%n
      correction = 0.0_c_double
      
      ! Interpolate from coarse points using prolongation operator P
      DO j = fine%ia(i), fine%ia(i+1)-1
        IF (fine%cf(fine%ja(j)) == 1) THEN  ! If it's a coarse point
          correction = correction + fine%P(j) * coarse%x(fine%ja(j))
        END IF
      END DO
      
      ! Apply correction
      fine%x(i) = fine%x(i) + correction
    END DO
  END SUBROUTINE prolongate_and_correct

  !> Initialize the AMG solver
  SUBROUTINE initialize_s(this, matrix_size, tolerance, max_iterations)
    CLASS(amg2d_t), INTENT(INOUT) :: this
    INTEGER(c_int), INTENT(IN) :: matrix_size
    REAL(c_double), INTENT(IN) :: tolerance
    INTEGER(c_int), INTENT(IN) :: max_iterations
    
    this%tol = tolerance
    this%max_iter = max_iterations
    this%cycle_type = 1  ! Default to V-cycle
    this%num_levels = 1  ! Will be updated in setup

    ! Allocate first level
    ALLOCATE(this%levels(MAX_LEVELS))
    this%levels(1)%n = matrix_size
    ALLOCATE(this%levels(1)%ia(matrix_size + 1))
    ALLOCATE(this%levels(1)%ja(matrix_size * 5))  ! Assuming 5-point stencil
    ALLOCATE(this%levels(1)%a(matrix_size * 5))
    ALLOCATE(this%levels(1)%x(matrix_size))
    ALLOCATE(this%levels(1)%b(matrix_size))
    ALLOCATE(this%levels(1)%res(matrix_size))
    ALLOCATE(this%levels(1)%cf(matrix_size))
  END SUBROUTINE initialize_s

  !> Setup the AMG hierarchy
  SUBROUTINE setup_s(this)
    CLASS(amg2d_t), INTENT(INOUT) :: this
    INTEGER(c_int) :: level
    
    level = 1
    DO WHILE (this%levels(level)%n > 50 .AND. level < MAX_LEVELS)  ! 50 is minimum size
      ! Perform coarsening
      CALL this%coarsen(level)
      
      ! Create interpolation operators
      CALL this%interpolate(level)
      
      ! Setup next level
      level = level + 1
      
      ! Construct coarse grid operator (Galerkin RAP)
      CALL construct_coarse_operator(this%levels(level-1), this%levels(level))
    END DO
    
    this%num_levels = level
  END SUBROUTINE setup_s

  !> Construct coarse grid operator using RAP
  SUBROUTINE construct_coarse_operator(fine, coarse)
    TYPE(amg_level_t), INTENT(IN) :: fine
    TYPE(amg_level_t), INTENT(INOUT) :: coarse
    ! TODO: Implement Galerkin coarse grid operator construction (RAP)
    ! R * A * P where R is restriction and P is prolongation
  END SUBROUTINE construct_coarse_operator

  !> Coarsening routine based on Ruge-Stueben coarsening
  SUBROUTINE coarsen_s(this, level)
    CLASS(amg2d_t), INTENT(INOUT) :: this
    INTEGER(c_int), INTENT(IN) :: level
    INTEGER(c_int) :: i, j, n, nnz
    REAL(c_double) :: diag_i, diag_j
    LOGICAL, ALLOCATABLE :: is_strong(:,:)  ! Strong connection matrix
    
    n = this%levels(level)%n
    nnz = SIZE(this%levels(level)%ja)
    
    ! Allocate temporary arrays
    ALLOCATE(is_strong(n,n))
    is_strong = .FALSE.
    
    ! Step 1: Identify strong connections
    DO i = 1, n
      diag_i = 0.0_c_double
      ! Find diagonal entry
      DO j = this%levels(level)%ia(i), this%levels(level)%ia(i+1)-1
        IF (this%levels(level)%ja(j) == i) THEN
          diag_i = ABS(this%levels(level)%a(j))
          EXIT
        END IF
      END DO
      
      ! Mark strong connections
      DO j = this%levels(level)%ia(i), this%levels(level)%ia(i+1)-1
        IF (ABS(this%levels(level)%a(j)) > STRONG_THRESHOLD * diag_i) THEN
          is_strong(i, this%levels(level)%ja(j)) = .TRUE.
        END IF
      END DO
    END DO
    
    ! Step 2: Select coarse points (C/F splitting)
    CALL select_coarse_points(n, is_strong, this%levels(level)%cf)
    
    DEALLOCATE(is_strong)
  END SUBROUTINE coarsen_s

  !> Interpolation operator construction
  SUBROUTINE interpolate_s(this, level)
    CLASS(amg2d_t), INTENT(INOUT) :: this
    INTEGER(c_int), INTENT(IN) :: level
    INTEGER(c_int) :: i, j, k, nc, nf
    REAL(c_double) :: sum_strong, sum_weak
    
    ! Count coarse points
    nc = COUNT(this%levels(level)%cf == 1)
    nf = this%levels(level)%n - nc
    
    ! Allocate prolongation and restriction operators
    IF (ALLOCATED(this%levels(level)%P)) DEALLOCATE(this%levels(level)%P)
    IF (ALLOCATED(this%levels(level)%R)) DEALLOCATE(this%levels(level)%R)
    ALLOCATE(this%levels(level)%P(this%levels(level)%n * 2))  ! Estimate size
    ALLOCATE(this%levels(level)%R(this%levels(level)%n * 2))  ! Restriction = Transpose(P)
    
    ! Construct interpolation weights
    DO i = 1, this%levels(level)%n
      IF (this%levels(level)%cf(i) == 0) THEN  ! Fine point
        sum_strong = 0.0_c_double
        sum_weak = 0.0_c_double
        
        ! Compute interpolation weights
        DO j = this%levels(level)%ia(i), this%levels(level)%ia(i+1)-1
          k = this%levels(level)%ja(j)
          IF (this%levels(level)%cf(k) == 1) THEN  ! Coarse point
            this%levels(level)%P(j) = -this%levels(level)%a(j)
            sum_strong = sum_strong + ABS(this%levels(level)%a(j))
          ELSE
            sum_weak = sum_weak + this%levels(level)%a(j)
          END IF
        END DO
        
        ! Normalize weights
        IF (sum_strong > 0.0_c_double) THEN
          DO j = this%levels(level)%ia(i), this%levels(level)%ia(i+1)-1
            k = this%levels(level)%ja(j)
            IF (this%levels(level)%cf(k) == 1) THEN
              this%levels(level)%P(j) = this%levels(level)%P(j) / sum_strong
            END IF
          END DO
        END IF
      END IF
    END DO
    
    ! Set restriction as transpose of prolongation
    this%levels(level)%R = this%levels(level)%P
  END SUBROUTINE interpolate_s

  !> Modified smoothing routine to use ILU
  SUBROUTINE smooth_s(this, level, num_sweeps)
    CLASS(amg2d_t), INTENT(INOUT) :: this
    INTEGER(c_int), INTENT(IN) :: level, num_sweeps
    INTEGER(c_int) :: sweep
    REAL(c_double), ALLOCATABLE :: residual(:), correction(:)
    
    ALLOCATE(residual(this%levels(level)%n))
    ALLOCATE(correction(this%levels(level)%n))

    ! Initialize ILU factorization if not already done
    IF (.NOT. ALLOCATED(this%levels(level)%ilu_smoother%d)) THEN
      CALL create_ilu(this%levels(level)%ilu_smoother, &
                     this%levels(level)%ia, &
                     this%levels(level)%ja, &
                     this%levels(level)%a, &
                     this%levels(level)%n)
    END IF

    DO sweep = 1, num_sweeps
      ! Compute residual
      CALL compute_residual(this%levels(level))
      residual = this%levels(level)%res

      ! Apply ILU preconditioner
      CALL apply_ilu(this%levels(level)%ilu_smoother, residual, correction)
      
      ! Update solution
      this%levels(level)%x = this%levels(level)%x + correction
    END DO

    DEALLOCATE(residual, correction)
  END SUBROUTINE smooth_s

  !> Multigrid cycle (V-cycle or W-cycle)
  SUBROUTINE cycle_s(this, level)
    CLASS(amg2d_t), INTENT(INOUT) :: this
    INTEGER(c_int), INTENT(IN) :: level
    INTEGER(c_int) :: i
    
    ! Pre-smoothing
    CALL this%smooth(level, 2)
    
    ! Compute residual
    CALL compute_residual(this%levels(level))
    
    ! Restrict residual
    CALL restrict_residual(this%levels(level), this%levels(level+1))
    
    ! Recursive call or solve directly
    IF (level < this%num_levels-1) THEN
      IF (this%cycle_type == 1) THEN  ! V-cycle
        CALL this%cycle(level+1)
      ELSE  ! W-cycle
        CALL this%cycle(level+1)
        CALL this%cycle(level+1)
      END IF
    ELSE
      ! Direct solve on coarsest level
      CALL direct_solve(this%levels(this%num_levels))
    END IF
    
    ! Prolongate and correct
    CALL prolongate_and_correct(this%levels(level), this%levels(level+1))
    
    ! Post-smoothing
    CALL this%smooth(level, 2)
  END SUBROUTINE cycle_s

  !> Main solve routine
  SUBROUTINE solve_s(this, b, x)
    CLASS(amg2d_t), INTENT(INOUT) :: this
    REAL(c_double), INTENT(IN) :: b(:)
    REAL(c_double), INTENT(OUT) :: x(:)
    INTEGER(c_int) :: iter
    REAL(c_double) :: res_norm
    
    ! Initialize solution
    x = 0.0_c_double
    this%levels(1)%x = x
    this%levels(1)%b = b
    
    ! Main iteration loop
    DO iter = 1, this%max_iter
      CALL this%cycle(1)
      
      ! Check convergence
      CALL compute_residual_norm(this%levels(1), res_norm)
      IF (res_norm < this%tol) EXIT
    END DO
    
    ! Copy solution back
    x = this%levels(1)%x
  END SUBROUTINE solve_s

  !> Deallocate all allocated memory
  SUBROUTINE destroy_s(this)
    CLASS(amg2d_t), INTENT(INOUT) :: this
    INTEGER :: level

    IF (ALLOCATED(this%levels)) THEN
      DO level = 1, this%num_levels
        IF (ALLOCATED(this%levels(level)%ia)) DEALLOCATE(this%levels(level)%ia)
        IF (ALLOCATED(this%levels(level)%ja)) DEALLOCATE(this%levels(level)%ja)
        IF (ALLOCATED(this%levels(level)%a))  DEALLOCATE(this%levels(level)%a)
        IF (ALLOCATED(this%levels(level)%x))  DEALLOCATE(this%levels(level)%x)
        IF (ALLOCATED(this%levels(level)%b))  DEALLOCATE(this%levels(level)%b)
        IF (ALLOCATED(this%levels(level)%res)) DEALLOCATE(this%levels(level)%res)
        IF (ALLOCATED(this%levels(level)%cf)) DEALLOCATE(this%levels(level)%cf)
        IF (ALLOCATED(this%levels(level)%P))  DEALLOCATE(this%levels(level)%P)
        IF (ALLOCATED(this%levels(level)%R))  DEALLOCATE(this%levels(level)%R)
      END DO
      DEALLOCATE(this%levels)
    END IF
  END SUBROUTINE destroy_s

END MODULE amg2d_m