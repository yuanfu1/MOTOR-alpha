MODULE TenW_m
  USE gzm_m, ONLY: gzm_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE kinds_m, ONLY: i_kind, r_kind
  USE InterpValue_m, ONLY: InterpValue_t
  USE CalVerDer_m, ONLY: CalVerDer_t
  USE PreCal_m, ONLY: PreCal_t
  USE Parameters_m, ONLY: g

  IMPLICIT NONE

  TYPE :: TenW_t
    TYPE(SingleGrid_t), POINTER :: sg
    TYPE(InterpValue_t) :: InterpValue
    TYPE(gzm_t) :: gzm
    TYPE(CalVerDer_t) :: CalVerDer
  CONTAINS
    PROCEDURE :: Tendency
    FINAL :: destructor
    PROCEDURE :: destroy
  END TYPE TenW_t

  INTERFACE TenW_t
    PROCEDURE :: constructor
  END INTERFACE TenW_t
CONTAINS

  FUNCTION constructor(sg) RESULT(this)
    TYPE(TenW_t) :: this
    TYPE(SingleGrid_t), TARGET :: sg

    this%sg => sg

    this%InterpValue = InterpValue_t(this%sg)
    this%CalVerDer = CalVerDer_t(this%sg)
    this%gzm = gzm_t(this%sg)

  END FUNCTION constructor

  SUBROUTINE Tendency(this, w, div, psi, chi, rho, P, ten_w)

    IMPLICIT NONE

    CLASS(TenW_t) :: this
    REAL(r_kind), DIMENSION(:, :), INTENT(IN) :: div, psi, chi, w, rho, P
    REAL(r_kind), INTENT(INOUT) :: ten_w(:, :)
    REAL(r_kind), DIMENSION(:, :), ALLOCATABLE :: psi2, chi2, div2, J_w_psi, &
                                                  F_w_chi, parw_parsigma, &
                                                  J_z_psi, F_z_chi, &
                                                  parP_parsigma, omg, rho2, P2

    ASSOCIATE (num_cell => this%sg%num_cell, &
               num_iCell => this%sg%num_icell, &
               vLevel => this%sg%vLevel)

      ALLOCATE (psi2(vLevel, num_cell), chi2(vLevel, num_cell), &
                div2(vLevel, num_cell), rho2(vLevel, num_cell), &
                P2(vLevel, num_cell), &
                J_w_psi(vLevel, num_cell), &
                F_w_chi(vLevel, num_cell), &
                parw_parsigma(vLevel, num_cell), &
                F_z_chi(vLevel, num_cell), &
                parP_parsigma(vLevel, num_cell), &
                J_z_psi(vLevel, num_cell), &
                omg(vLevel, num_cell))

      psi2 = 0.0D0
      chi2 = 0.0D0
      div2 = 0.0D0
      rho2 = 0.0D0
      P2 = 0.0D0
      J_w_psi = 0.0D0
      F_w_chi = 0.0D0
      parw_parsigma = 0.0D0
      J_z_psi = 0.0D0
      F_z_chi = 0.0D0
      parP_parsigma = 0.0D0
      omg = 0.0D0

      omg = this%sg%Hz(:, 1:num_iCell) * w(:, 1:num_iCell)

      CALL this%InterpValue%InterpValue(psi, psi2)
      CALL this%InterpValue%InterpValue(chi, chi2)
      CALL this%InterpValue%InterpValue(div, div2)
      CALL this%InterpValue%InterpValue(rho, rho2)
      CALL this%InterpValue%InterpValue(P, P2)

      CALL this%CalVerDer%FirstOrder(w, parw_parsigma)
      CALL this%CalVerDer%FirstOrder(P2, parP_parsigma)

      CALL this%gzm%Jacobian(w, psi2, J_w_psi)
      CALL this%gzm%Divergen(w, chi2, F_w_chi)
      CALL this%gzm%Jacobian(this%sg%zHght, psi2, J_z_psi)
      CALL this%gzm%Divergen(this%sg%zHght, chi2, F_z_chi)

      ten_w(:, 1:num_iCell) = J_w_psi(:, 1:num_iCell) - F_w_chi(:, 1:num_iCell) &
                              - this%sg%Hz(:, 1:num_iCell) * parw_parsigma(:, 1:num_iCell) &
                              * (J_z_psi(:, 1:num_icell) - F_z_chi(:, 1:num_iCell)) &
                              - (w(:, 1:num_iCell) &
                                 - this%sg%Hz(:, 1:num_iCell) * parw_parsigma(:, 1:num_iCell) &
                                 * this%sg%zHght(:, 1:num_icell)) * div2(:, 1:num_iCell) &
                              - omg(:, 1:num_iCell) * parw_parsigma(:, 1:num_iCell) - g &
                              - this%sg%Hz(:, 1:num_iCell) * 1.0D0 / rho2(:, 1:num_iCell) &
                              * parP_parsigma(:, 1:num_iCell)

    END ASSOCIATE

    DEALLOCATE (psi2, chi2, div2, P2, J_w_psi, &
                F_w_chi, parw_parsigma, &
                J_z_psi, F_z_chi, &
                parP_parsigma, omg, rho2)

  END SUBROUTINE

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(TenW_t), INTENT(INOUT) :: this

    IF (ASSOCIATED(this%sg)) NULLIFY (this%sg)

  END SUBROUTINE destructor

  SUBROUTINE destroy(this)
    IMPLICIT NONE
    CLASS(TenW_t), INTENT(INOUT) :: this

    IF (ASSOCIATED(this%sg)) NULLIFY (this%sg)

  END SUBROUTINE destroy

END MODULE

