MODULE TenVor_m
  USE gzm_m, ONLY: gzm_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE kinds_m, ONLY: i_kind, r_kind
  USE InterpValue_m, ONLY: InterpValue_t
  USE CalVerDer_m, ONLY: CalVerDer_t
  USE PreCal_m, ONLY: PreCal_t
  USE State_m, ONLY: State_t
  USE parameters_m, ONLY: EarthRadius

  IMPLICIT NONE

  TYPE :: TenVor_t
    TYPE(SingleGrid_t), POINTER :: sg
    TYPE(InterpValue_t) :: InterpValue
    TYPE(gzm_t) :: gzm
    TYPE(CalVerDer_t) :: CalVerDer
    ! TYPE(PreCal_t), pointer :: PreCal
  CONTAINS
    PROCEDURE :: Tendency
    PROCEDURE :: destroy
    FINAL :: destructor
  END TYPE TenVor_t

  INTERFACE TenVor_t
    PROCEDURE :: constructor
  END INTERFACE TenVor_t
CONTAINS

  FUNCTION constructor(sg) RESULT(this)
    TYPE(TenVor_t) :: this
    TYPE(SingleGrid_t), TARGET :: sg
    ! TYPE(PreCal_t), target :: PreCal

    this%sg => sg
    ! this%PreCal => PreCal

    this%InterpValue = InterpValue_t(this%sg)
    this%gzm = gzm_t(this%sg)

  END FUNCTION constructor

  SUBROUTINE Tendency(this, vor, psi, chi, w, rho, P, PreCal, ten_vor, X1)

    IMPLICIT NONE

    CLASS(TenVor_t) :: this
    TYPE(State_t), INTENT(INOUT), OPTIONAL :: X1
    TYPE(PreCal_t), INTENT(IN) :: PreCal
    REAL(r_kind), DIMENSION(:, :), INTENT(IN) :: vor, psi, chi, w, rho, P
    REAL(r_kind), INTENT(INOUT) :: ten_vor(:, :)
    REAL(r_kind), DIMENSION(:, :), ALLOCATABLE :: eta, w2, omg, J_eta_psi, &
                                                  F_eta_chi, &
                                                  F_DM_psisigma, J_DM_chisigma, &
                                                  F_DM_psi, J_DM_chi, &
                                                  F_omg_psisigma, J_omg_chisigma, &
                                                  F_omg_psi, J_omg_chi, &
                                                  J_invRho_P, J_hzRhoP_z

    REAL(r_kind) :: r_sg
    INTEGER(i_kind) :: k

    ASSOCIATE (num_cell => this%sg%num_cell, &
               num_iCell => this%sg%num_icell, &
               vLevel => this%sg%vLevel, &
               sigma => this%sg%sigma)

      ALLOCATE (eta(vLevel, num_cell), w2(vLevel, num_cell), &
                omg(vLevel, num_cell), J_eta_psi(vLevel, num_cell), &
                F_eta_chi(vLevel, num_cell), &
                F_DM_psisigma(vLevel, num_cell), J_DM_chisigma(vLevel, num_cell), &
                F_omg_psisigma(vLevel, num_cell), J_omg_chisigma(vLevel, num_cell), &
                J_invRho_P(vLevel, num_cell), J_hzRhoP_z(vLevel, num_cell), &
                F_DM_psi(vLevel, num_cell), J_DM_chi(vLevel, num_cell), &
                F_omg_psi(vLevel, num_cell), J_omg_chi(vLevel, num_cell))

      w2 = 0.0D0

      eta = 0.0D0
      w2 = 0.0D0
      omg = 0.0D0
      J_eta_psi = 0.0D0
      F_eta_chi = 0.0D0
      F_DM_psisigma = 0.0D0
      J_DM_chisigma = 0.0D0
      F_DM_psi = 0.0D0
      J_DM_chi = 0.0D0
      F_omg_psisigma = 0.0D0
      J_omg_chisigma = 0.0D0
      F_omg_psi = 0.0D0
      J_omg_chi = 0.0D0
      J_invRho_P = 0.0D0
      J_hzRhoP_z = 0.0D0

      ! print *, 'size of w and w2 are', size(w, 1), size(w2, 1)
      CALL this%InterpValue%InterpValue(w, w2)
      omg(:, 1:num_iCell) = this%sg%Hz(:, 1:num_iCell) * w2(:, 1:num_iCell)
      eta = vor + this%sg%f

      CALL this%gzm%Jacobian(eta, psi, J_eta_psi)
      CALL this%gzm%Divergen(eta, chi, F_eta_chi)
      CALL this%gzm%Divergen(PreCal%DM, PreCal%par_psi_sigma, F_DM_psisigma)
      CALL this%gzm%Jacobian(PreCal%DM, PreCal%par_chi_sigma, J_DM_chisigma)
      CALL this%gzm%Divergen(omg, PreCal%par_psi_sigma, F_omg_psisigma)
      CALL this%gzm%Jacobian(omg, PreCal%par_chi_sigma, J_omg_chisigma)
      CALL this%gzm%Jacobian(PreCal%invRho, P, J_invRho_P)
      CALL this%gzm%Jacobian(PreCal%HzInvRhoParPsigma, this%sg%zHght, J_hzRhoP_z)
      CALL this%gzm%Divergen(omg, psi, F_omg_psi)
      CALL this%gzm%Jacobian(omg, chi, J_omg_chi)
      CALL this%gzm%Divergen(PreCal%DM, psi, F_DM_psi)
      CALL this%gzm%Jacobian(PreCal%DM, chi, J_DM_chi)

      DO k = 1, vLevel
        r_sg = EarthRadius + sigma(k)
        ten_vor(k, 1:num_iCell) = J_eta_psi(k, 1:num_iCell) - F_eta_chi(k, 1:num_iCell) &
                                  + F_DM_psisigma(k, 1:num_iCell) + J_DM_chisigma(k, 1:num_iCell) &
                                  - F_omg_psisigma(k, 1:num_iCell) - J_omg_chisigma(k, 1:num_iCell) &
                                  - 1.0D0 / r_sg * (F_DM_psi(k, 1:num_iCell) + J_DM_chi(k, 1:num_iCell) &
                                                    - PreCal%DM(k, 1:num_iCell) * vor(k, 1:num_iCell) &
                                                    - F_omg_psi(k, 1:num_iCell) - J_omg_chi(k, 1:num_iCell)) &
                                  - J_invRho_P(k, 1:num_iCell) + J_hzRhoP_z(k, 1:num_iCell)
      END DO
      ! IF (PRESENT(X1)) THEN
      !    X1%fields(X1%getVarIdx('J_eta_psi'))%data(:, :, 1) = J_eta_psi
      !    X1%fields(X1%getVarIdx('F_eta_chi'))%data(:, :, 1) = F_eta_chi
      !    X1%fields(X1%getVarIdx('F_DM_psisigma'))%data(:, :, 1) = F_DM_psisigma
      !    X1%fields(X1%getVarIdx('J_DM_chisigma'))%data(:, :, 1) = J_DM_chisigma
      !    X1%fields(X1%getVarIdx('DM'))%data(:, :, 1) = PreCal%DM
      !    X1%fields(X1%getVarIdx('psi_sigma'))%data(:, :, 1) = PreCal%par_psi_sigma
      ! END IF

      ! PRINT *, 'max of unused parameters of F and J functions ::', maxval(F_omg_psisigma(:, 1:num_iCell)), &
      !    maxval(J_omg_chisigma(:, 1:num_iCell)), maxval(J_invRho_P(:, 1:num_iCell)), maxval(J_hzrhoP_z(:, 1:num_iCell))

    END ASSOCIATE

    DEALLOCATE (eta, w2, omg, J_eta_psi, &
                F_eta_chi, &
                F_DM_psisigma, J_DM_chisigma, &
                F_DM_psi, J_DM_chi, &
                F_omg_psisigma, J_omg_chisigma, &
                F_omg_psi, J_omg_chi, &
                J_invRho_P, J_hzRhoP_z)

  END SUBROUTINE

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(TenVor_t), INTENT(INOUT) :: this

    IF (ASSOCIATED(this%sg)) NULLIFY (this%sg)
    ! IF (ASSOCIATED(this%PreCal)) NULLIFY (this%PreCal)

  END SUBROUTINE destructor

  SUBROUTINE destroy(this)
    IMPLICIT NONE
    CLASS(TenVor_t), INTENT(INOUT) :: this

    IF (ASSOCIATED(this%sg)) NULLIFY (this%sg)
    ! IF (ASSOCIATED(this%PreCal)) NULLIFY (this%PreCal)

  END SUBROUTINE destroy

END MODULE
