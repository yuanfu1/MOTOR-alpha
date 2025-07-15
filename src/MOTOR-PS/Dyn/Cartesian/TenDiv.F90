MODULE TenDiv_m
  USE gzm_m, ONLY: gzm_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE kinds_m, ONLY: i_kind, r_kind
  USE InterpValue_m, ONLY: InterpValue_t
  USE CalVerDer_m, ONLY: CalVerDer_t
  USE PreCal_m, ONLY: PreCal_t
  USE State_m, ONLY: State_t

  IMPLICIT NONE

  TYPE :: TenDiv_t
    TYPE(SingleGrid_t), POINTER :: sg
    TYPE(InterpValue_t) :: InterpValue
    TYPE(gzm_t) :: gzm
    TYPE(CalVerDer_t) :: CalVerDer
    TYPE(PreCal_t), POINTER :: PreCal
  CONTAINS
    PROCEDURE :: Tendency
    PROCEDURE :: destroy
    FINAL :: destructor
  END TYPE TenDiv_t

  INTERFACE TenDiv_t
    PROCEDURE :: constructor
  END INTERFACE TenDiv_t

CONTAINS

  FUNCTION constructor(sg, PreCal) RESULT(this)
    TYPE(TenDiv_t) :: this
    TYPE(SingleGrid_t), TARGET :: sg
    TYPE(PreCal_t), TARGET :: PreCal

    this%sg => sg
    this%PreCal => PreCal

    this%gzm = gzm_t(this%sg)
    this%InterpValue = InterpValue_t(this%sg)

  END FUNCTION constructor

  SUBROUTINE Tendency(this, div, vor, psi, chi, w, rho, P, ten_div, X1)
    IMPLICIT NONE

    CLASS(TenDiv_t) :: this
    TYPE(State_t), INTENT(INOUT) :: X1
    REAL(r_kind), DIMENSION(:, :), INTENT(IN) :: div, vor, psi, chi, w, rho, P
    REAL(r_kind), INTENT(INOUT) :: ten_div(:, :)
    REAL(r_kind), DIMENSION(:, :), ALLOCATABLE :: w2, omg, &
                                                  F_vor_psi, J_vor_chi, KE, LKE, &
                                                  F_psi_psi, F_chi_chi, J_psi_chi, &
                                                  J_DM_psisigma, F_DM_chisigma, &
                                                  J_omg_psisigma, F_omg_chisigma, &
                                                  F_invRho_P, F_hzrhoP_z, &
                                                  F_f_psi, J_f_chi

    ASSOCIATE (num_cell => this%sg%num_cell, &
               num_iCell => this%sg%num_icell, &
               vLevel => this%sg%vLevel)

      ALLOCATE (w2(vLevel, num_cell), omg(vLevel, num_cell), &
                F_vor_psi(vLevel, num_cell), J_vor_chi(vLevel, num_cell), &
                KE(vLevel, num_cell), LKE(vLevel, num_cell), &
                F_psi_psi(vLevel, num_cell), F_chi_chi(vLevel, num_cell), &
                J_psi_chi(vLevel, num_cell) &
                J_DM_psisigma(vLevel, num_cell), F_DM_chisigma(vLevel, num_cell), &
                J_omg_psisigma(vLevel, num_cell), F_omg_chisigma(vLevel, num_cell), &
                F_invRho_P(vLevel, num_cell), F_hzrhoP_z(vLevel, num_cell), &
                F_f_psi(vLevel, num_cell), J_f_chi(vLevel, num_cell))

      w2 = 0.0D0
      !  CALL interpvalue(sigma1, sigma2, w(:, 1:num_iCell), w2(:, 1:num_iCell))

      CALL this%InterpValue%InterpValue(w, w2)
      ! w2(:, 1:num_iCell) = (w(1:vLevel, num_iCell) + w(2:vLevel + 1, num_iCell))/2.0D0
      omg(:, 1:num_iCell) = this%sg%Hz(:, 1:num_iCell) * w2(:, 1:num_iCell)

      CALL this%gzm%Divergen(vor, psi, F_vor_psi)
      CALL this%gzm%Jacobian(vor, chi, J_vor_chi)
      CALL this%gzm%Divergen(psi, psi, F_psi_psi)
      CALL this%gzm%Divergen(chi, chi, F_chi_chi)
      CALL this%gzm%Jacobian(psi, chi, J_psi_chi)
      CALL this%gzm%Jacobian(this%PreCal%DM, this%PreCal%par_psi_sigma, J_DM_psisigma)
      CALL this%gzm%Divergen(this%PreCal%DM, this%PreCal%par_chi_sigma, F_DM_chisigma)
      CALL this%gzm%Jacobian(omg, this%PreCal%par_psi_sigma, J_omg_psisigma)
      CALL this%gzm%Divergen(omg, this%PreCal%par_chi_sigma, F_omg_chisigma)
      CALL this%gzm%Divergen(this%PreCal%invRho, P, F_invRho_P)
      CALL this%gzm%Divergen(this%PreCal%HzInvRhoParPsigma, this%sg%zHght, F_hzrhoP_z)
      CALL this%gzm%Divergen(this%sg%f, psi, F_f_psi)
      CALL this%gzm%Jacobian(this%sg%f, chi, J_f_chi)

      KE(:, 1:num_iCell) = (F_psi_psi(:, 1:num_iCell) - psi(:, 1:num_iCell) * vor(:, 1:num_iCell) &
                            + F_chi_chi(:, 1:num_iCell) - chi(:, 1:num_iCell) * div(:, 1:num_iCell)) / 2.0D0 &
                           + J_psi_chi(:, 1:num_iCell)
      CALL this%gzm%Laplacia(KE, LKE)

      ! Calculate Tendency of divergence
      ten_div(:, 1:num_iCell) = F_vor_psi(:, 1:num_iCell) + J_vor_chi(:, 1:num_iCell) &
                                - LKE(:, 1:num_iCell) &
                                - J_DM_psisigma(:, 1:num_iCell) + F_DM_chisigma(:, 1:num_iCell) &
                                + J_omg_psisigma(:, 1:num_iCell) - F_omg_chisigma(:, 1:num_iCell) &
                                - F_invRho_P(:, 1:num_iCell) + F_hzrhoP_z(:, 1:num_iCell) &
                                + F_f_psi(:, 1:num_iCell) + J_f_chi(:, 1:num_iCell)
    END ASSOCIATE

    X1%fields(X1%getVarIdx('J_del_psi'))%DATA(:, :, 1) = J_vor_chi
    X1%fields(X1%getVarIdx('F_del_chi'))%DATA(:, :, 1) = F_vor_psi
    X1%fields(X1%getVarIdx('J_DM_psisigma'))%DATA(:, :, 1) = J_DM_psisigma
    X1%fields(X1%getVarIdx('F_DM_chisigma'))%DATA(:, :, 1) = F_DM_chisigma
    X1%fields(X1%getVarIdx('DM'))%DATA(:, :, 1) = this%PreCal%DM
    X1%fields(X1%getVarIdx('psi_sigma'))%DATA(:, :, 1) = this%PreCal%par_psi_sigma
    X1%fields(X1%getVarIdx('F_invRho_P'))%DATA(:, :, 1) = F_invRho_P
    X1%fields(X1%getVarIdx('F_HzinvRhoPsigma_z'))%DATA(:, :, 1) = F_hzrhoP_z
    X1%fields(X1%getVarIdx('F_f_psi'))%DATA(:, :, 1) = F_f_psi
    X1%fields(X1%getVarIdx('J_f_chi'))%DATA(:, :, 1) = J_f_chi

    DEALLOCATE (w2, omg, &
                F_vor_psi, J_vor_chi, &
                KE, LKE, &
                F_psi_psi, F_chi_chi, &
                J_psi_chi &
                J_DM_psisigma, F_DM_chisigma, &
                J_omg_psisigma, F_omg_chisigma, &
                F_invRho_P, F_hzrhoP_z, &
                F_f_psi, J_f_chi)

  END SUBROUTINE

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(TenDiv_t), INTENT(INOUT) :: this

    IF (ASSOCIATED(this%sg)) NULLIFY (this%sg)
    IF (ASSOCIATED(this%PreCal)) NULLIFY (this%PreCal)

  END SUBROUTINE destructor

  SUBROUTINE destroy(this)
    IMPLICIT NONE
    CLASS(TenDiv_t), INTENT(INOUT) :: this

    IF (ASSOCIATED(this%sg)) NULLIFY (this%sg)
    IF (ASSOCIATED(this%PreCal)) NULLIFY (this%PreCal)

  END SUBROUTINE destroy

END MODULE
