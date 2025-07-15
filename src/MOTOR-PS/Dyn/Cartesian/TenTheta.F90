MODULE TenTheta_m
  USE gzm_m, ONLY: gzm_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE kinds_m, ONLY: i_kind, r_kind
  USE InterpValue_m, ONLY: InterpValue_t
  USE CalVerDer_m, ONLY: CalVerDer_t
  USE PreCal_m, ONLY: PreCal_t

  IMPLICIT NONE

  TYPE :: TenTheta_t
    TYPE(SingleGrid_t), POINTER :: sg
    TYPE(InterpValue_t) :: InterpValue
    TYPE(gzm_t) :: gzm
    TYPE(CalVerDer_t) :: CalVerDer
    TYPE(PreCal_t), POINTER :: PreCal
  CONTAINS
    PROCEDURE :: Tendency
    FINAL :: destructor
  END TYPE TenTheta_t

  INTERFACE TenTheta_t
    PROCEDURE :: constructor
  END INTERFACE TenTheta_t

CONTAINS

  FUNCTION constructor(sg, PreCal) RESULT(this)
    TYPE(TenTheta_t) :: this
    TYPE(SingleGrid_t), TARGET :: sg
    TYPE(PreCal_t), TARGET :: PreCal

    this%sg => sg
    this%PreCal => PreCal
    this%InterpValue = InterpValue_t(this%sg)
    this%CalVerDer = CalVerDer_t(this%sg)
    this%gzm = gzm_t(this%sg)

  END FUNCTION constructor

  SUBROUTINE Tendency(this, theta, div, psi, chi, w, ten_theta)
    IMPLICIT NONE

    CLASS(TenTheta_t) :: this
    REAL(r_kind), DIMENSION(:, :), INTENT(IN) :: psi, chi, w, theta, div
    REAL(r_kind), INTENT(OUT) :: ten_theta(:, :)
    REAL(r_kind), DIMENSION(:, :), ALLOCATABLE :: psi2, chi2, omg, J_theta_psi, &
                                                  F_theta_chi, partheta_parsigma, &
                                                  J_z_psi, F_z_chi

    ASSOCIATE (num_cell => this%sg%num_cell, &
               num_iCell => this%sg%num_icell, &
               vLevel => this%sg%vLevel)

      ALLOCATE (psi2(vLevel, num_cell), chi2(vLevel, num_cell), omg(vLevel, num_cell), &
                J_theta_psi(vLevel, num_cell), F_theta_chi(vLevel, num_cell), &
                partheta_parsigma(vLevel, num_cell), J_z_psi(vLevel, num_cell), &
                F_z_chi(vLevel, num_cell))

      omg(:, 1:num_icell) = this%sg%Hz * w(:, 1:num_iCell)

      CALL this%InterpValue%InterpValue(psi, psi2)
      CALL this%InterpValue%InterpValue(chi, chi2)

      CALL this%CalVerDer%FirstOrder(theta, partheta_parsigma)

      CALL this%gzm%Jacobian(theta, psi2, J_theta_psi)
      CALL this%gzm%Divergen(theta, chi2, F_theta_chi)
      CALL this%gzm%Jacobian(this%sg%zHght, psi2, J_z_psi)
      CALL this%gzm%Divergen(this%sg%zHght, chi2, F_z_chi)

      ten_theta(:, 1:num_iCell) = J_theta_psi(:, 1:num_iCell) &
                                  - F_theta_chi(:, 1:num_iCell) &
                                  - this%sg%Hz(:, 1:num_iCell) &
                                  * partheta_parsigma(:, 1:num_iCell) &
                                  * (J_z_psi(:, 1:num_iCell) &
                                     - F_z_chi(:, 1:num_iCell)) &
                                  + (theta(:, 1:num_iCell) &
                                     - this%sg%Hz(:, 1:num_iCell) &
                                     * partheta_parsigma(:, 1:num_iCell) &
                                     * this%sg%zHght(:, 1:num_iCell)) * div(:, 1:num_iCell) &
                                  - omg(:, 1:num_iCell) * partheta_parsigma(:, 1:num_iCell)

      DEALLOCATE (psi2, chi2, omg, J_theta_psi, &
                  F_theta_chi, partheta_parsigma, &
                  J_z_psi, F_z_chi)
    END ASSOCIATE
  END SUBROUTINE

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(TenTheta_t), INTENT(INOUT) :: this

    IF (ASSOCIATED(this%sg)) NULLIFY (this%sg)
    IF (ASSOCIATED(this%PreCal)) NULLIFY (this%PreCal)

  END SUBROUTINE destructor
END MODULE
