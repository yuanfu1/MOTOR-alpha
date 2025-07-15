MODULE TenDiv_m
  USE gzm_m, ONLY: gzm_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE kinds_m, ONLY: i_kind, r_kind
  USE InterpValue_m, ONLY: InterpValue_t
  USE CalVerDer_m, ONLY: CalVerDer_t
  USE PreCal_m, ONLY: PreCal_t

  IMPLICIT NONE

  TYPE :: TenDiv_t
    TYPE(SingleGrid_t), POINTER :: sg
    TYPE(InterpValue_t) :: InterpValue
    TYPE(gzm_t) :: gzm
    TYPE(CalVerDer_t) :: CalVerDer
    TYPE(PreCal_t), POINTER :: PreCal
  CONTAINS
    PROCEDURE :: Tendency
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

  SUBROUTINE Tendency(this, div, psi, kai, w, rho, P, ten_div)
    IMPLICIT NONE

    CLASS(TenDiv_t) :: this
    REAL(r_kind), DIMENSION(:, :), INTENT(IN) :: div, psi, kai, w, rho, P
    REAL(r_kind), INTENT(INOUT) :: ten_div(:, :)
    REAL(r_kind), DIMENSION(:, :), ALLOCATABLE :: w2, omg, &
                                                  J_div_psi, F_div_kai, &
                                                  J_DM_psisigma, F_DM_kaisigma, &
                                                  J_ogm_psisigma, F_omg_kaisigma, &
                                                  F_invRho_P, F_hzrhoP_z, &
                                                  F_f_psi, J_f_kai

    ASSOCIATE (num_cell => this%sg%num_cell, &
               num_iCell => this%sg%num_icell, &
               vLevel => this%sg%vLevel)

      ALLOCATE (w2(vLevel, num_cell), omg(vLevel, num_cell), &
                J_div_psi(vLevel, num_cell), F_div_kai(vLevel, num_cell), &
                J_DM_psisigma(vLevel, num_cell), F_DM_kaisigma(vLevel, num_cell), &
                J_ogm_psisigma(vLevel, num_cell), F_omg_kaisigma(vLevel, num_cell), &
                F_invRho_P(vLevel, num_cell), F_hzrhoP_z(vLevel, num_cell), &
                F_f_psi(vLevel, num_cell), J_f_kai(vLevel, num_cell))

      w2 = 0.0D0
      !  CALL interpvalue(sigma1, sigma2, w(:, 1:num_iCell), w2(:, 1:num_iCell))

      CALL this%InterpValue%InterpValue(w, w2)
      ! w2(:, 1:num_iCell) = (w(1:vLevel, num_iCell) + w(2:vLevel + 1, num_iCell))/2.0D0
      omg(:, 1:num_iCell) = this%sg%Hz(:, 1:num_iCell) * w2(:, 1:num_iCell)

      CALL this%gzm%Jacobian(div, psi, J_div_psi)
      CALL this%gzm%Divergen(div, kai, F_div_kai)
      CALL this%gzm%Jacobian(this%PreCal%DM, this%PreCal%par_psi_sigma, J_DM_psisigma)
      CALL this%gzm%Divergen(this%PreCal%DM, this%PreCal%par_kai_sigma, F_DM_kaisigma)
      CALL this%gzm%Jacobian(omg, this%PreCal%par_psi_sigma, J_ogm_psisigma)
      CALL this%gzm%Divergen(omg, this%PreCal%par_kai_sigma, F_omg_kaisigma)
      CALL this%gzm%Divergen(this%PreCal%invRho, P, F_invRho_P)
      CALL this%gzm%Divergen(this%PreCal%HzInvRohParPsigma, this%sg%zHght, F_hzrhoP_z)
      CALL this%gzm%Divergen(this%sg%f, psi, F_f_psi)
      CALL this%gzm%Jacobian(this%sg%f, kai, J_f_kai)

      ! Calculate Tendency of divergence
      ten_div(:, 1:num_iCell) = J_div_psi(:, 1:num_iCell) - F_div_kai(:, 1:num_iCell) &
                                - J_DM_psisigma(:, 1:num_iCell) + F_DM_kaisigma(:, 1:num_iCell) &
                                + J_ogm_psisigma(:, 1:num_iCell) - F_omg_kaisigma(:, 1:num_iCell) &
                                - F_invRho_P(:, 1:num_iCell) + F_hzrhoP_z(:, 1:num_iCell) &
                                + F_f_psi(:, 1:num_iCell) + J_f_kai
    END ASSOCIATE

    DEALLOCATE (w2, omg, &
                J_div_psi, F_div_kai, &
                J_DM_psisigma, F_DM_kaisigma, &
                J_ogm_psisigma, F_omg_kaisigma, &
                F_invRho_P, F_hzrhoP_z, &
                F_f_psi, J_f_kai)

  END SUBROUTINE

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(TenDiv_t), INTENT(INOUT) :: this

    IF (ASSOCIATED(this%sg)) NULLIFY (this%sg)
    IF (ASSOCIATED(this%PreCal)) NULLIFY (this%PreCal)

  END SUBROUTINE destructor

END MODULE
