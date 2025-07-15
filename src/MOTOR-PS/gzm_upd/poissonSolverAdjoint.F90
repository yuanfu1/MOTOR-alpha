MODULE poissonSolver_adjoint_m
  USE geometry_m, ONLY: geometry_t
  USE kinds_m, ONLY: i_kind, r_kind
  USE psMatrix_adjoint_m, ONLY: psMatrix_adj_t
  USE psSolutn_adj_m, ONLY: psSolutn_adj_t
  USE namelist_ps_adj_m, ONLY: namelist_ps_adj
  USE parameters_m, ONLY: EarthRadius

  TYPE :: poissonSolver_adj_t
    TYPE(geometry_t), POINTER :: geometry
    TYPE(psMatrix_adj_t), ALLOCATABLE :: psMat_adj(:)
    CHARACTER(LEN=3) :: solver
    INTEGER(i_kind) :: vLevel
    REAL(r_kind), ALLOCATABLE :: omegas(:)
    INTEGER(i_kind) :: &
      nRelax, &
      nCycle, &
      nIterPre, &
      nIterPost, &
      max_band

  CONTAINS
    FINAL :: destructor_adjoint
    PROCEDURE, PUBLIC :: PoissonSol_adjoint
    PROCEDURE, PUBLIC :: PoissonSol_sphere_adjoint
    PROCEDURE, PRIVATE :: fullMcycle_adjoint
    PROCEDURE, PRIVATE :: singleVcyc_adjoint
    PROCEDURE, PRIVATE :: iteration_adjoint

  END TYPE poissonSolver_adj_t

  INTERFACE poissonSolver_adj_t
    PROCEDURE :: constructor_adjoint
  END INTERFACE poissonSolver_adj_t

CONTAINS

  FUNCTION constructor_adjoint(configFile, geometry) RESULT(this)
    TYPE(poissonSolver_adj_t) :: this
    TYPE(geometry_t), TARGET, INTENT(IN) :: geometry
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    INTEGER(i_kind) :: i

    this%geometry => geometry
    ALLOCATE (this%psMat_adj(this%geometry%mg%mg_coarsest:this%geometry%mg%mg_finest))
    ALLOCATE (this%omegas(10))

    CALL namelist_ps_adj(configFile, this%solver, this%nCycle, this%nIterPre, &
                         this%nIterPost, this%nRelax, this%omegas, this%max_band)

    DO i = this%geometry%mg%mg_coarsest, this%geometry%mg%mg_finest
      this%psMat_adj(i) = psMatrix_adj_t(this%max_band, this%geometry%mg%sg(i))
    END DO
  END FUNCTION constructor_adjoint

  SUBROUTINE PoissonSol_adjoint(this, gLevel, vLevel, rights_adj, solutn_adj)
    CLASS(poissonSolver_adj_t) :: this
    INTEGER(i_kind), INTENT(IN) :: vLevel, gLevel
    REAL(r_kind), INTENT(INOUT) :: rights_adj(vlevel, this%geometry%mg%sg(gLevel)%num_cell)
    REAL(r_kind), INTENT(OUT) :: solutn_adj(vlevel, this%geometry%mg%sg(gLevel)%num_cell)

    INTEGER(i_kind) :: icycle, i
    TYPE(psSolutn_adj_t), ALLOCATABLE :: psSolutn_adj(:)

    ALLOCATE (psSolutn_adj(this%geometry%mg%mg_coarsest:gLevel))
    DO i = this%geometry%mg%mg_coarsest, gLevel
      psSolutn_adj(i) = psSolutn_adj_t(vlevel, this%geometry%mg%sg(i)%num_cell)
    END DO

    psSolutn_adj(gLevel)%adj_rights = rights_adj
    psSolutn_adj(gLevel)%adj_solutn = solutn_adj

    DO icycle = 1, this%nCycle
      IF (this%solver .EQ. 'FMV') THEN
        CALL this%fullMcycle_adjoint(gLevel, vlevel, psSolutn_adj)
      ELSE IF (this%solver .EQ. 'VVV') THEN
        CALL this%singleVcyc_adjoint(gLevel, vlevel, psSolutn_adj)
      ELSE
        WRITE (*, *) 'No adjoint solver defined'
      END IF
    END DO

    solutn_adj = psSolutn_adj(gLevel)%adj_solutn
    DEALLOCATE (psSolutn_adj)
  END SUBROUTINE PoissonSol_adjoint

  SUBROUTINE PoissonSol_sphere_adjoint(this, gLevel, vLevel, rights_adj, solutn_adj)
    CLASS(poissonSolver_adj_t) :: this
    INTEGER(i_kind), INTENT(IN) :: vLevel, gLevel
    REAL(r_kind), INTENT(INOUT) :: rights_adj(vlevel, this%geometry%mg%sg(gLevel)%num_cell)
    REAL(r_kind), INTENT(OUT) :: solutn_adj(vlevel, this%geometry%mg%sg(gLevel)%num_cell)

    INTEGER(i_kind) :: k
    REAL(r_kind) :: ratio

    ASSOCIATE (sg => this%geometry%mg%sg(gLevel))
      CALL this%PoissonSol_adjoint(gLevel, vLevel, rights_adj, solutn_adj)
      DO k = 1, sg%vLevel
        ratio = (EarthRadius + sg%sigma(k)) / EarthRadius
        solutn_adj(k, :) = solutn_adj(k, :) * ratio**4.0D0
      END DO
    END ASSOCIATE
  END SUBROUTINE PoissonSol_sphere_adjoint

  RECURSIVE SUBROUTINE fullMcycle_adjoint(this, gLevel, vLevel, psSolutn_adj)
    CLASS(poissonSolver_adj_t) :: this
    INTEGER(i_kind), INTENT(IN) :: vLevel, gLevel
    TYPE(psSolutn_adj_t), INTENT(INOUT) :: psSolutn_adj(:)

    REAL(r_kind), ALLOCATABLE :: tempState(:, :), worksp(:, :)

    ASSOCIATE (sg => this%geometry%mg%sg(gLevel))
      ALLOCATE (tempState(vLevel, sg%num_cell), worksp(vLevel, sg%num_cell))

      IF (gLevel .LE. this%geometry%mg%mg_coarsest) THEN
        CALL this%singleVcyc_adjoint(gLevel, vLevel, psSolutn_adj)
      ELSE
        ! Calculate the adjoint residual
        worksp = 0.0D0
        CALL this%psMat_adj(gLevel)%CalResidus_adjoint(vLevel, psSolutn_adj(gLevel)%adj_rights, &
                                                       psSolutn_adj(gLevel)%adj_solutn, worksp, sg)

        ! Restrict residual to the coarser grid
        psSolutn_adj(gLevel - 1)%adj_rights = 0.0D0
        CALL this%geometry%mg%restrictionAtGLevel(gLevel, vLevel, worksp, psSolutn_adj(gLevel - 1)%adj_rights, .FALSE.)

        ! Recursive call for the coarser grid
        psSolutn_adj(gLevel - 1)%adj_solutn = 0.0D0
        CALL this%fullMcycle_adjoint(gLevel - 1, vLevel, psSolutn_adj)

        ! Prolong the coarse grid correction and apply it to the fine grid
        worksp = 0.0D0
        CALL this%geometry%mg%prolongationAtGLevel(gLevel - 1, vLevel, worksp, psSolutn_adj(gLevel - 1)%adj_solutn, .FALSE.)
        psSolutn_adj(gLevel)%adj_solutn = psSolutn_adj(gLevel)%adj_solutn + worksp

        CALL this%singleVcyc_adjoint(gLevel, vLevel, psSolutn_adj)
      END IF
      DEALLOCATE (tempState, worksp)
    END ASSOCIATE
  END SUBROUTINE fullMcycle_adjoint

  SUBROUTINE iteration_adjoint(this, gLevel, vLevel, psSolutn_adj, nIter)
    CLASS(poissonSolver_adj_t) :: this
    INTEGER(i_kind), INTENT(IN) :: vLevel, gLevel, nIter
    TYPE(psSolutn_adj_t), INTENT(INOUT) :: psSolutn_adj(:)

    REAL(r_kind), ALLOCATABLE :: tempState(:, :)

    ASSOCIATE (sg => this%geometry%mg%sg(gLevel))
      ALLOCATE (tempState(vLevel, sg%num_cell))

      DO iover = 1, nIter
        tempState = psSolutn_adj(gLevel)%adj_solutn
        DO itr = 1, this%nRelax
          IF (sg%mpddInfo_sg%nProc .NE. 1) THEN
            CALL Jacobi(this%psMat_adj(gLevel)%num_eqns, sg%num_cell, vLevel, &
                        this%psMat_adj(gLevel)%max_band, this%psMat_adj(gLevel)%adj_num_elem, &
                        this%psMat_adj(gLevel)%adj_idx_cols, this%psMat_adj(gLevel)%adj_row_elem, &
                        this%omegas(itr), psSolutn_adj(gLevel)%adj_rights, tempState, psSolutn_adj(gLevel)%adj_solutn)
            CALL sg%ExchangeMatOnHalo2D(vLevel, psSolutn_adj(gLevel)%adj_solutn)
          ELSE
            CALL GaussSeidel(this%psMat_adj(gLevel)%num_eqns, vLevel, this%psMat_adj(gLevel)%max_band, &
                             this%psMat_adj(gLevel)%adj_num_elem, this%psMat_adj(gLevel)%adj_idx_cols, &
                             this%psMat_adj(gLevel)%adj_row_elem, this%omegas(itr), psSolutn_adj(gLevel)%adj_rights, &
                             tempState, psSolutn_adj(gLevel)%adj_solutn)
          END IF
          IF (itr .LT. this%nRelax) tempState = psSolutn_adj(gLevel)%adj_solutn
        END DO
      END DO

      DEALLOCATE (tempState)
    END ASSOCIATE
  END SUBROUTINE iteration_adjoint

  RECURSIVE SUBROUTINE singleVcyc_adjoint(this, gLevel, vLevel, psSolutn_adj)
    CLASS(poissonSolver_adj_t) :: this
    INTEGER(i_kind), INTENT(IN) :: vLevel, gLevel
    TYPE(psSolutn_adj_t), INTENT(INOUT) :: psSolutn_adj(:)

    REAL(r_kind), ALLOCATABLE :: worksp(:, :)

    ASSOCIATE (sg => this%geometry%mg%sg(gLevel))
      ALLOCATE (worksp(vLevel, sg%num_cell))

      CALL this%iteration_adjoint(gLevel, vLevel, psSolutn_adj, this%nIterPre)

      IF (gLevel .GT. this%geometry%mg%mg_coarsest) THEN
        ! Calculate the adjoint residual
        worksp = 0.0D0
        CALL this%psMat_adj(gLevel)%CalResidus_adjoint(vLevel, psSolutn_adj(gLevel)%adj_rights, &
                                                       psSolutn_adj(gLevel)%adj_solutn, worksp, sg)

        ! Restrict residual to the coarser grid
        psSolutn_adj(gLevel - 1)%adj_rights = 0.0D0
        CALL this%geometry%mg%restrictionAtGLevel(gLevel, vLevel, worksp, psSolutn_adj(gLevel - 1)%adj_rights, .FALSE.)

        ! Recursive call for the coarser grid
        psSolutn_adj(gLevel - 1)%adj_solutn = 0.0D0
        CALL this%singleVcyc_adjoint(gLevel - 1, vLevel, psSolutn_adj)

        ! Prolong the correction and apply it to the fine grid
        worksp = 0.0D0
        CALL this%geometry%mg%prolongationAtGLevel(gLevel - 1, vLevel, worksp, psSolutn_adj(gLevel - 1)%adj_solutn, .FALSE.)
        psSolutn_adj(gLevel)%adj_solutn = psSolutn_adj(gLevel)%adj_solutn + worksp
      END IF

      CALL this%iteration_adjoint(gLevel, vLevel, psSolutn_adj, this%nIterPost)

      DEALLOCATE (worksp)
    END ASSOCIATE
  END SUBROUTINE singleVcyc_adjoint

  IMPURE ELEMENTAL SUBROUTINE destructor_adjoint(this)
    TYPE(poissonSolver_adj_t), INTENT(INOUT) :: this

    IF (ALLOCATED(this%psMat_adj)) DEALLOCATE (this%psMat_adj)
    IF (ALLOCATED(this%omegas)) DEALLOCATE (this%omegas)
  END SUBROUTINE destructor_adjoint

END MODULE poissonSolver_adjoint_m
