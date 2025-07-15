MODULE TenVor_m
  USE gzm_m, ONLY: gzm_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE kinds_m, ONLY: i_kind, r_kind
  USE InterpValue_m, ONLY: InterpValue_t
  USE CalVerDer_m, ONLY: CalVerDer_t
  USE PreCal_m, ONLY: PreCal_t
  USE State_m, ONLY: State_t

  IMPLICIT NONE

  TYPE :: TenVor_t
    TYPE(SingleGrid_t), POINTER :: sg
    TYPE(InterpValue_t) :: InterpValue
    TYPE(gzm_t) :: gzm
    TYPE(CalVerDer_t) :: CalVerDer
    TYPE(PreCal_t), POINTER :: PreCal
  CONTAINS
    PROCEDURE :: Tendency
    PROCEDURE :: destroy
    FINAL :: destructor
  END TYPE TenVor_t

  INTERFACE TenVor_t
    PROCEDURE :: constructor
  END INTERFACE TenVor_t
CONTAINS

  FUNCTION constructor(sg, PreCal) RESULT(this)
    TYPE(TenVor_t) :: this
    TYPE(SingleGrid_t), TARGET :: sg
    TYPE(PreCal_t), TARGET :: PreCal

    this%sg => sg
    this%PreCal => PreCal

    this%InterpValue = InterpValue_t(this%sg)
    this%gzm = gzm_t(this%sg)

  END FUNCTION constructor

  SUBROUTINE Tendency(this, vor, psi, chi, w, rho, P, ten_vor, X1)

    IMPLICIT NONE

    CLASS(TenVor_t) :: this
    TYPE(State_t), INTENT(INOUT) :: X1
    REAL(r_kind), DIMENSION(:, :), INTENT(IN) :: vor, psi, chi, w, rho, P
    REAL(r_kind), INTENT(INOUT) :: ten_vor(:, :)
    REAL(r_kind), DIMENSION(:, :), ALLOCATABLE :: eta, w2, omg, J_eta_psi, &
                                                  F_eta_chi, &
                                                  F_DM_psisigma, J_DM_chisigma, &
                                                  F_omg_psisigma, J_omg_chisigma, &
                                                  J_invrho_P, J_hzrhoP_z

    ASSOCIATE (num_cell => this%sg%num_cell, &
               num_iCell => this%sg%num_icell, &
               vLevel => this%sg%vLevel)

      ALLOCATE (eta(vLevel, num_cell), w2(vLevel, num_cell), &
                omg(vLevel, num_cell), J_eta_psi(vLevel, num_cell), &
                F_eta_chi(vLevel, num_cell), &
                F_DM_psisigma(vLevel, num_cell), J_DM_chisigma(vLevel, num_cell), &
                F_omg_psisigma(vLevel, num_cell), J_omg_chisigma(vLevel, num_cell), &
                J_invRho_P(vLevel, num_cell), J_hzRhoP_z(vLevel, num_cell))

        !! =================== Calculate Divergence =======================================
      ! Calculate KE(Kinetic Energy)
      w2 = 0.0D0
      PRINT *, 'size of w and w2 are', SIZE(w, 1), SIZE(w2, 1)
      CALL this%InterpValue%InterpValue(w, w2)
      omg(:, 1:num_iCell) = this%sg%Hz(:, 1:num_iCell) * w2(:, 1:num_iCell)
      eta = vor + this%sg%f
      !  CALL interp1d(sigma1, sigma2, w(:, 1:num_iCell), w2(:, 1:num_iCell))
      ! w2(:, 1:num_iCell) = (w(1:vLevel, num_iCell) + w(2:vLevel + 1, num_iCell))/2.0D0

      CALL this%gzm%Jacobian(eta, psi, J_eta_psi)
      CALL this%gzm%Divergen(eta, chi, F_eta_chi)
      CALL this%gzm%Divergen(this%PreCal%DM, this%PreCal%par_psi_sigma, F_DM_psisigma)
      CALL this%gzm%Jacobian(this%PreCal%DM, this%PreCal%par_chi_sigma, J_DM_chisigma)
      CALL this%gzm%Divergen(omg, this%PreCal%par_psi_sigma, F_omg_psisigma)
      CALL this%gzm%Jacobian(omg, this%PreCal%par_chi_sigma, J_omg_chisigma)
      CALL this%gzm%Jacobian(this%PreCal%invRho, P, J_invRho_P)
      CALL this%gzm%Jacobian(this%PreCal%HzInvRhoParPsigma, this%sg%zHght, J_hzrhoP_z)

      ten_vor(:, 1:num_iCell) = J_eta_psi(:, 1:num_iCell) - F_eta_chi(:, 1:num_iCell) &
                                + F_DM_psisigma(:, 1:num_iCell) + J_DM_chisigma(:, 1:num_iCell) &
                                - F_omg_psisigma(:, 1:num_iCell) - J_omg_chisigma(:, 1:num_iCell) &
                                - J_invrho_P(:, 1:num_iCell) + J_hzrhoP_z(:, 1:num_iCell)

      X1%fields(X1%getVarIdx('J_eta_psi'))%DATA(:, :, 1) = J_eta_psi
      X1%fields(X1%getVarIdx('F_eta_chi'))%DATA(:, :, 1) = F_eta_chi
      X1%fields(X1%getVarIdx('F_DM_psisigma'))%DATA(:, :, 1) = F_DM_psisigma
      X1%fields(X1%getVarIdx('J_DM_chisigma'))%DATA(:, :, 1) = J_DM_chisigma
      X1%fields(X1%getVarIdx('DM'))%DATA(:, :, 1) = this%PreCal%DM
      X1%fields(X1%getVarIdx('psi_sigma'))%DATA(:, :, 1) = this%PreCal%par_psi_sigma

      PRINT *, 'max of unused parameters of F and J functions ::', MAXVAL(F_omg_psisigma(:, 1:num_iCell)), &
        MAXVAL(J_omg_chisigma(:, 1:num_iCell)), MAXVAL(J_invRho_P(:, 1:num_iCell)), MAXVAL(J_hzrhoP_z(:, 1:num_iCell))

    END ASSOCIATE

    DEALLOCATE (eta, w2, omg, J_eta_psi, &
                F_eta_chi, &
                F_DM_psisigma, J_DM_chisigma, &
                F_omg_psisigma, J_omg_chisigma, &
                J_invrho_P, J_hzrhoP_z)

  END SUBROUTINE

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(TenVor_t), INTENT(INOUT) :: this

    IF (ASSOCIATED(this%sg)) NULLIFY (this%sg)
    IF (ASSOCIATED(this%PreCal)) NULLIFY (this%PreCal)

  END SUBROUTINE destructor

  SUBROUTINE destroy(this)
    IMPLICIT NONE
    CLASS(TenVor_t), INTENT(INOUT) :: this

    IF (ASSOCIATED(this%sg)) NULLIFY (this%sg)
    IF (ASSOCIATED(this%PreCal)) NULLIFY (this%PreCal)

  END SUBROUTINE destroy

END MODULE
