MODULE TenRho_m
  USE gzm_m, ONLY: gzm_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE kinds_m, ONLY: i_kind, r_kind
  USE InterpValue_m, ONLY: InterpValue_t
  USE CalVerDer_m, ONLY: CalVerDer_t
  USE PreCal_m, ONLY: PreCal_t

  IMPLICIT NONE

  TYPE :: TenRho_t
    TYPE(SingleGrid_t), POINTER :: sg
    TYPE(InterpValue_t) :: InterpValue
    TYPE(gzm_t) :: gzm
    TYPE(CalVerDer_t) :: CalVerDer
    TYPE(PreCal_t), POINTER :: PreCal
  CONTAINS
    PROCEDURE :: Tendency
    FINAL :: destructor
  END TYPE TenRho_t

  INTERFACE TenRho_t
    PROCEDURE :: constructor
  END INTERFACE TenRho_t

CONTAINS

  FUNCTION constructor(sg, PreCal) RESULT(this)
    TYPE(TenRho_t) :: this
    TYPE(SingleGrid_t), TARGET :: sg
    TYPE(PreCal_t), TARGET :: PreCal

    this%sg => sg
    this%PreCal => PreCal
    this%InterpValue = InterpValue_t(this%sg)
    this%CalVerDer = CalVerDer_t(this%sg)
    this%gzm = gzm_t(this%sg)

  END FUNCTION constructor

  SUBROUTINE Tendency(this, rho, div, psi, chi, w, ten_rho)
    IMPLICIT NONE

    CLASS(TenRho_t) :: this
    REAL(r_kind), DIMENSION(:, :), INTENT(IN) :: div, psi, chi, w, rho
    REAL(r_kind), INTENT(INOUT) :: ten_rho(:, :)
    REAL(r_kind), DIMENSION(:, :), ALLOCATABLE :: w2, omg, J_rho_psi, F_rho_chi, &
                                                  pardiv_parsigma, parrho_parsigma, &
                                                  F_z_chisigma, J_z_psisigma, &
                                                  rhomw, parrhomw_parsigma

    ASSOCIATE (num_cell => this%sg%num_cell, &
               num_iCell => this%sg%num_icell, &
               vLevel => this%sg%vLevel)

      ALLOCATE (w2(vLevel, num_cell), omg(vLevel, num_cell), J_rho_psi(vLevel, num_cell), &
                F_rho_chi(vLevel, num_cell), pardiv_parsigma(vLevel, num_cell), &
                parrho_parsigma(vLevel, num_cell), F_z_chisigma(vLevel, num_cell), &
                J_z_psisigma(vLevel, num_cell), rhomw(vLevel, num_cell), &
                parrhomw_parsigma(vLevel, num_cell))

      w2 = 0.0D0
      !  CALL interpvalue(sigma1, sigma2, w(:, 1:num_iCell), w2(:, 1:num_iCell))
      CALL this%InterpValue%InterpValue(w, w2)
      ! w2(:, 1:num_iCell) = (w(1:vLevel, num_iCell) + w(2:vLevel + 1, num_iCell))/2.0D0
      omg(:, 1:num_iCell) = this%sg%Hz(:, 1:num_iCell) * w2(:, 1:num_iCell)
      rhomw(:, 1:num_iCell) = rho(:, 1:num_iCell) * w2(:, 1:num_iCell)

      CALL this%CalVerDer%FirstOrder(div, pardiv_parsigma)
      CALL this%CalVerDer%FirstOrder(rho, parrho_parsigma)
      CALL this%CalVerDer%FirstOrder(rhomw, parrhomw_parsigma)

      CALL this%gzm%Jacobian(rho, psi, J_rho_psi)
      CALL this%gzm%Divergen(rho, chi, F_rho_chi)
      CALL this%gzm%Divergen(this%sg%zHght, this%PreCal%par_chi_sigma, F_z_chisigma)
      CALL this%gzm%Jacobian(this%sg%zHght, this%PreCal%par_psi_sigma, J_z_psisigma)

      ! Calculate Tendency of Rho
      ten_rho(:, 1:num_iCell) = J_rho_psi(:, 1:num_iCell) &
                                - F_rho_chi(:, 1:num_iCell) &
                                + this%sg%Hz(:, 1:num_iCell) * rho(:, 1:num_iCell) &
                                * (F_z_chisigma(:, 1:num_iCell) - J_z_psisigma(:, 1:num_iCell) &
                                   - this%sg%zHght(:, 1:num_iCell) * pardiv_parsigma(:, 1:num_iCell)) &
                                + this%PreCal%DM(:, 1:num_iCell) * parrho_parsigma(:, 1:num_iCell) &
                                - this%sg%Hz(:, 1:num_iCell) * parrhomw_parsigma(:, 1:num_iCell)
    END ASSOCIATE

    DEALLOCATE (w2, J_rho_psi, F_rho_chi, &
                pardiv_parsigma, parrho_parsigma, &
                F_z_chisigma, J_z_psisigma, &
                rhomw, parrhomw_parsigma)

  END SUBROUTINE

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(TenRho_t), INTENT(INOUT) :: this

    IF (ASSOCIATED(this%sg)) NULLIFY (this%sg)
    IF (ASSOCIATED(this%PreCal)) NULLIFY (this%PreCal)

  END SUBROUTINE destructor

END MODULE
