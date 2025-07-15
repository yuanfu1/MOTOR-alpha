MODULE Tenqvapor_m
  USE gzm_m, ONLY: gzm_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE kinds_m, ONLY: i_kind, r_kind
  USE InterpValue_m, ONLY: InterpValue_t
  USE CalVerDer_m, ONLY: CalVerDer_t
  USE PreCal_m, ONLY: PreCal_t
  USE State_m, ONLY: State_t

  IMPLICIT NONE

  TYPE :: Tenqvapor_t
    TYPE(SingleGrid_t), POINTER :: sg
    TYPE(InterpValue_t) :: InterpValue
    TYPE(gzm_t) :: gzm
    TYPE(CalVerDer_t) :: CalVerDer
    ! TYPE(PreCal_t), pointer :: PreCal
  CONTAINS
    PROCEDURE :: Tendency
    PROCEDURE :: destroy
    FINAL :: destructor
  END TYPE Tenqvapor_t

  INTERFACE Tenqvapor_t
    PROCEDURE :: constructor
  END INTERFACE Tenqvapor_t

CONTAINS

  FUNCTION constructor(sg) RESULT(this)
    TYPE(Tenqvapor_t) :: this
    TYPE(SingleGrid_t), TARGET :: sg
    ! TYPE(PreCal_t), target :: PreCal

    this%sg => sg
    ! this%PreCal => PreCal

    this%InterpValue = InterpValue_t(this%sg)
    this%CalVerDer = CalVerDer_t(this%sg)
    this%gzm = gzm_t(this%sg)
  END FUNCTION constructor

  SUBROUTINE Tendency(this, q, div, psi, chi, w, ten_q, X)
    IMPLICIT NONE
    CLASS(Tenqvapor_t) :: this
    TYPE(State_t), INTENT(INOUT), OPTIONAL :: X
    REAL(r_kind), DIMENSION(:, :), INTENT(IN) :: q, psi, chi, w, div
    REAL(r_kind), INTENT(OUT) :: ten_q(:, :)
    REAL(r_kind), DIMENSION(:, :), ALLOCATABLE :: w2, omg, J_q_psi, F_q_chi, &
                                                  parq_parsigma, J_z_psi, F_z_chi

    ASSOCIATE (num_cell => this%sg%num_cell, &
               num_iCell => this%sg%num_icell, &
               vLevel => this%sg%vLevel)

      ALLOCATE (w2(vLevel, num_cell), omg(vLevel, num_cell), &
                J_q_psi(vLevel, num_cell), F_q_chi(vLevel, num_cell), &
                parq_parsigma(vLevel, num_cell), J_z_psi(vLevel, num_cell), &
                F_z_chi(vLevel, num_cell))

      w2 = 0.0D0
      omg = 0.0D0
      J_q_psi = 0.0D0
      F_q_chi = 0.0D0
      parq_parsigma = 0.0D0
      J_z_psi = 0.0D0
      F_z_chi = 0.0D0

      CALL this%InterpValue%InterpValue(w, w2)
      omg(:, 1:num_icell) = this%sg%Hz(:, 1:num_iCell) * w2(:, 1:num_iCell)

      CALL this%CalVerDer%FirstOrder(q, parq_parsigma)

      CALL this%gzm%Jacobian(q, psi, J_q_psi)
      CALL this%gzm%Divergen(q, chi, F_q_chi)
      ! PRINT *, 'max and min of sg%zHght are: ', MAXVAL(this%sg%zHght), MINVAL(this%sg%zHght), MAXVAL(psi), MINVAL(psi)
      CALL this%gzm%Jacobian(this%sg%zHght, psi, J_z_psi)
      CALL this%gzm%Divergen(this%sg%zHght, chi, F_z_chi)

      ten_q(:, 1:num_iCell) = J_q_psi(:, 1:num_iCell) &
                              - F_q_chi(:, 1:num_iCell) &
                              - this%sg%Hz(:, 1:num_iCell) &
                              * parq_parsigma(:, 1:num_iCell) &
                              * (J_z_psi(:, 1:num_iCell) &
                                 - F_z_chi(:, 1:num_iCell)) &
                              + (q(:, 1:num_iCell) &
                                 - this%sg%Hz(:, 1:num_iCell) &
                                 * parq_parsigma(:, 1:num_iCell) &
                                 * this%sg%zHght(:, 1:num_iCell)) * div(:, 1:num_iCell) &
                              - omg(:, 1:num_iCell) * parq_parsigma(:, 1:num_iCell)

      ! IF (PRESENT(X)) THEN
      !    X%fields(X%getVarIdx('J_q_psi'))%data(:, :, 1) = J_q_psi
      !    X%fields(X%getVarIdx('F_q_chi'))%data(:, :, 1) = F_q_chi
      !    X%fields(X%getVarIdx('J_z_psi'))%data(:, :, 1) = J_z_psi
      !    X%fields(X%getVarIdx('F_z_chi'))%data(:, :, 1) = F_z_chi
      !    X%fields(X%getVarIdx('pq_sg'))%data(:, :, 1) = parq_parsigma

      ! END IF

      DEALLOCATE (w2, omg, J_q_psi, F_q_chi, &
                  parq_parsigma, J_z_psi, F_z_chi)
    END ASSOCIATE
  END SUBROUTINE

  SUBROUTINE destroy(this)
    IMPLICIT NONE
    CLASS(Tenqvapor_t), INTENT(INOUT) :: this

    IF (ASSOCIATED(this%sg)) NULLIFY (this%sg)
    ! IF (ASSOCIATED(this%PreCal)) NULLIFY (this%PreCal)

  END SUBROUTINE destroy

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(Tenqvapor_t), INTENT(INOUT) :: this

    IF (ASSOCIATED(this%sg)) NULLIFY (this%sg)
    ! IF (ASSOCIATED(this%PreCal)) NULLIFY (this%PreCal)

  END SUBROUTINE destructor
END MODULE
