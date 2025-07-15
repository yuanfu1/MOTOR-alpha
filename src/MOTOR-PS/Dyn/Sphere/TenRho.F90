MODULE TenRho_m
  USE gzm_m, ONLY: gzm_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE kinds_m, ONLY: i_kind, r_kind
  USE InterpValue_m, ONLY: InterpValue_t
  USE CalVerDer_m, ONLY: CalVerDer_t
  USE PreCal_m, ONLY: PreCal_t
  USE Parameters_m, ONLY: EarthRadius
  USE State_m, ONLY: State_t

  IMPLICIT NONE

  TYPE :: TenRho_t
    TYPE(SingleGrid_t), POINTER :: sg
    TYPE(InterpValue_t) :: InterpValue
    TYPE(gzm_t) :: gzm
    TYPE(CalVerDer_t) :: CalVerDer
    ! TYPE(PreCal_t), pointer :: PreCal
  CONTAINS
    PROCEDURE :: Tendency
    PROCEDURE :: destroy
    FINAL :: destructor
  END TYPE TenRho_t

  INTERFACE TenRho_t
    PROCEDURE :: constructor
  END INTERFACE TenRho_t

CONTAINS

  FUNCTION constructor(sg) RESULT(this)
    TYPE(TenRho_t) :: this
    TYPE(SingleGrid_t), TARGET :: sg
    ! TYPE(PreCal_t), target :: PreCal

    this%sg => sg
    ! this%PreCal => PreCal
    this%InterpValue = InterpValue_t(this%sg)
    this%CalVerDer = CalVerDer_t(this%sg)
    this%gzm = gzm_t(this%sg)

  END FUNCTION constructor

  SUBROUTINE Tendency(this, rho, div, psi, chi, w, PreCal, ten_rho, X)
    IMPLICIT NONE

    CLASS(TenRho_t) :: this
    TYPE(State_t), INTENT(INOUT), OPTIONAL :: X
    TYPE(PreCal_t), INTENT(IN) :: PreCal
    REAL(r_kind), DIMENSION(:, :), INTENT(IN) :: div, psi, chi, w, rho
    REAL(r_kind), INTENT(INOUT) :: ten_rho(:, :)
    REAL(r_kind), DIMENSION(:, :), ALLOCATABLE :: w2, omg, J_rho_psi, F_rho_chi, &
                                                  L_chisigma, parrho_parsigma, &
                                                  F_z_chisigma, J_z_psisigma, &
                                                  rhomw, parrhomw_parsigma, &
                                                  J_z_psi, F_z_chi, &
                                                  F_z_chisigma_0, J_z_psi_r0, F_z_chi_r0
    REAL(r_kind) :: r_sg
    INTEGER(i_kind) :: k

    ASSOCIATE (num_cell => this%sg%num_cell, &
               num_iCell => this%sg%num_icell, &
               vLevel => this%sg%vLevel, &
               sigma => this%sg%sigma)

      ALLOCATE (w2(vLevel, num_cell), omg(vLevel, num_cell), J_rho_psi(vLevel, num_cell), &
                F_rho_chi(vLevel, num_cell), parrho_parsigma(vLevel, num_cell), &
                L_chisigma(vLevel, num_cell), F_z_chisigma(vLevel, num_cell), &
                J_z_psisigma(vLevel, num_cell), rhomw(vLevel, num_cell), &
                parrhomw_parsigma(vLevel, num_cell), &
                J_z_psi(vLevel, num_cell), F_z_chi(vLevel, num_cell), &
                F_z_chisigma_0(vLevel, num_cell), J_z_psi_r0(vLevel, num_cell), F_z_chi_r0(vLevel, num_cell))

      w2 = 0.0D0
      omg = 0.0D0
      J_rho_psi = 0.0D0
      F_rho_chi = 0.0D0
      L_chisigma = 0.0D0
      parrho_parsigma = 0.0D0
      F_z_chisigma = 0.0D0
      J_z_psisigma = 0.0D0
      rhomw = 0.0D0
      parrhomw_parsigma = 0.0D0
      J_z_psi = 0.0D0
      F_z_chi = 0.0D0

      !  CALL interpvalue(sigma1, sigma2, w(:, 1:num_iCell), w2(:, 1:num_iCell))
      CALL this%InterpValue%InterpValue(w, w2)
      ! w2(:, 1:num_iCell) = (w(1:vLevel, num_iCell) + w(2:vLevel + 1, num_iCell))/2.0D0
      ! omg(:, 1:num_iCell) = this%sg%Hz(:, 1:num_iCell)*w2(:, 1:num_iCell)
      rhomw(:, 1:num_iCell) = rho(:, 1:num_iCell) * w2(:, 1:num_iCell)

      ! CALL this%CalVerDer%FirstOrder(div, pardiv_parsigma)
      CALL this%CalVerDer%FirstOrder(rho, parrho_parsigma)
      CALL this%CalVerDer%FirstOrder(rhomw, parrhomw_parsigma)

      CALL this%gzm%Jacobian(rho, psi, J_rho_psi)
      CALL this%gzm%Divergen(rho, chi, F_rho_chi)
      CALL this%gzm%Divergen(this%sg%zHght, PreCal%par_chi_sigma, F_z_chisigma)
      CALL this%gzm%Jacobian(this%sg%zHght, PreCal%par_psi_sigma, J_z_psisigma)
      CALL this%gzm%Divergen(this%sg%zHght, chi, F_z_chi)
      CALL this%gzm%Jacobian(this%sg%zHght, psi, J_z_psi)
      CALL this%gzm%Laplacia(PreCal%par_chi_sigma, L_chisigma)

      ! Calculate Tendency of Rho
      DO k = 1, vLevel
        r_sg = EarthRadius + sigma(k)
        ten_rho(k, 1:num_iCell) = J_rho_psi(k, 1:num_iCell) &
                                  - F_rho_chi(k, 1:num_iCell) &
                                  + this%sg%Hz(k, 1:num_iCell) * rho(k, 1:num_iCell) &
                                  * (F_z_chisigma(k, 1:num_iCell) - J_z_psisigma(k, 1:num_iCell) &
                                     - this%sg%zHght(k, 1:num_iCell) * L_chisigma(k, 1:num_iCell) &
                                     + 1.0D0 / r_sg * (J_z_psi(k, 1:num_iCell) - F_z_chi(k, 1:num_iCell) &
                                                       + this%sg%zHght(k, 1:num_iCell) * div(k, 1:num_iCell))) &
                                  + PreCal%DM(k, 1:num_iCell) * parrho_parsigma(k, 1:num_iCell) &
                                  - this%sg%Hz(k, 1:num_iCell) * parrhomw_parsigma(k, 1:num_iCell)

        F_z_chisigma_0(k, :) = F_z_chisigma(k, :) - this%sg%zHght(k, :) * L_chisigma(k, :)
        F_z_chi_r0(k, :) = 1.0D0 / r_sg * (F_z_chi(k, :) - this%sg%zHght(k, :) * div(k, :))
        J_z_psi_r0(k, :) = 1.0D0 / r_sg * J_z_psi(k, :)

      END DO
    END ASSOCIATE
    ! IF (PRESENT(X)) THEN
    !    X%fields(X%getVarIdx('prho_sg'))%data(:, :, 1) = parrho_parsigma
    !    X%fields(X%getVarIdx('J_rho_psi'))%data(:, :, 1) = J_rho_psi
    !    X%fields(X%getVarIdx('F_rho_chi'))%data(:, :, 1) = F_rho_chi
    !    X%fields(X%getVarIdx('F_z_chisg'))%data(:, :, 1) = F_z_chisigma_0
    !    X%fields(X%getVarIdx('J_z_psisg'))%data(:, :, 1) = J_z_psisigma
    !    X%fields(X%getVarIdx('F_z_chir'))%data(:, :, 1) = F_z_chi_r0
    !    X%fields(X%getVarIdx('J_z_psir'))%data(:, :, 1) = J_z_psi_r0

    ! END IF

    DEALLOCATE (w2, omg, J_rho_psi, F_rho_chi, &
                L_chisigma, parrho_parsigma, &
                F_z_chisigma, J_z_psisigma, &
                rhomw, parrhomw_parsigma, &
                J_z_psi, F_z_chi)

  END SUBROUTINE

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(TenRho_t), INTENT(INOUT) :: this

    IF (ASSOCIATED(this%sg)) NULLIFY (this%sg)
    ! IF (ASSOCIATED(this%PreCal)) NULLIFY (this%PreCal)

  END SUBROUTINE destructor

  SUBROUTINE destroy(this)
    IMPLICIT NONE
    CLASS(TenRho_t), INTENT(INOUT) :: this

    IF (ASSOCIATED(this%sg)) NULLIFY (this%sg)
    ! IF (ASSOCIATED(this%PreCal)) NULLIFY (this%PreCal)

  END SUBROUTINE destroy

END MODULE
