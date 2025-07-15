MODULE TenDivVor_m
  USE gzmTerr_m, ONLY: gzmTerr_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE kinds_m, ONLY: i_kind, r_kind
  USE InterpValue_m, ONLY: InterpValue_t
  TYPE :: TenDivVor_t
    TYPE(InterpValue_t) :: InterpValue
    TYPE(gzmTerr_t) :: gzmTerr
  CONTAINS
    PROCEDURE :: Tendency
  END TYPE TenDivVor_t

CONTAINS

  SUBROUTINE Tendency(this, sg, vor, phi, kai, w, rho, P, ten_div, ten_vor)
    CLASS(TenDivVor_t) :: this
    TYPE(SingleGrid_t), INTENT(IN) :: sg
    REAL(r_kind), DIMENSION(:, :), INTENT(IN) :: vor, phi, kai, w, rho, P
    REAL(r_kind), INTENT(OUT) :: ten_div(:, :), ten_vor(:, :)
    REAL(r_kind), DIMENSION(:, :), ALLOCATABLE :: eta, KE, LKE, F_phi_phi, F_kai_kai, &
                                                  Lphi, Lkai, J_phi_kai, F_eta_phi, &
                                                  J_eta_kai, J_w_phiz, F_w_kaiz, &
                                                  invrho, F_invrho_P, &
                                                  parphi_parsigma, parphi_parz, &
                                                  parkai_parsigma, parkai_parz, &
                                                  J_eta_phi, F_eta_kai, F_w_phiz, &
                                                  J_w_kaiz, J_rho_P, &
                                                  w2

    ASSOCIATE (num_cell => sg%num_cell, &
               num_iCell => sg%num_icell, &
               vLevel => sg%vLevel, &
               sigma1 => sg%sigma, &
               sigma2 => sg%sigma)

      ALLOCATE (eta(vLevel, num_cell), KE(vLevel, num_cell), &
                F_phi_phi(vLevel, num_cell), &
                F_kai_kai(vLevel, num_cell), &
                Lphi(vLevel, num_cell), Lkai(vLevel, num_cell), &
                J_phi_kai(vLevel, num_cell), LKE(vLevel, num_cell), &
                F_eta_phi(vLevel, num_cell), J_eta_kai(vLevel, num_cell), &
                J_w_phiz(vLevel, num_cell), F_w_kaiz(vLevel, num_cell), &
                invrho(vLevel, num_cell), F_invrho_P(vLevel, num_cell), &
                parphi_parsigma(vLevel, num_cell), parphi_parz(vLevel, num_cell), &
                parkai_parsigma(vLevel, num_cell), parkai_parz(vLevel, num_cell), &
                w2(vLevel, num_cell), J_eta_phi(vLevel, num_cell), &
                F_eta_kai(vLevel, num_cell), F_w_phiz(vLevel, num_cell), &
                J_w_kaiz(vLevel, num_cell), J_rho_P(vLevel, num_cell))

        !! =================== Calculate Divergence =======================================
      ! Calculate KE(Kinetic Energy)
      w2 = 0.0D0
      !  CALL interp1d(sigma1, sigma2, w(:, 1:num_iCell), w2(:, 1:num_iCell))
      w2(:, 1:num_iCell) = (w(1:vLevel, num_iCell) + w(2:vLevel + 1, num_iCell)) / 2.0D0
      CALL this%gzmTerr%initialize(sg)
      CALL this%gzmTerr%Divergen(sg, phi, phi, F_phi_phi)
      CALL this%gzmTerr%Divergen(sg, kai, kai, F_kai_kai)
      CALL this%gzmTerr%Laplacia(sg, phi, Lphi)
      CALL this%gzmTerr%Laplacia(sg, kai, Lkai)
      CALL this%gzmTerr%Jacobian(sg, phi, kai, J_phi_kai)
      KE(:, 1:num_icell) = 0.5D0 * (F_phi_phi(:, 1:num_icell) + F_kai_kai(:, 1:num_icell) &
                                    - phi(:, 1:num_icell) * Lphi(:, 1:num_icell) &
                                    - kai(:, 1:num_icell) * Lkai(:, 1:num_icell)) &
                           + J_phi_kai(:, 1:num_icell)
      CALL this%gzmTerr%Laplacia(sg, KE, LKE)
      ! Calculate eta (vor+f)
      eta(:, 1:num_iCell) = vor(:, 1:num_iCell) + sg%fcor(:, 1:num_iCell)
      CALL this%gzmTerr%Divergen(sg, eta, phi, F_eta_phi)
      CALL this%gzmTerr%Jacobian(sg, eta, kai, J_eta_kai)

      ! Calculate parphi_parz
      CALL this%gzmTerr%FirstOrderDerivative(sg, phi, this%gzmTerr%coef_1rst, parphi_parsigma)
      parphi_parz(:, 1:num_icell) = this%gzmTerr%Hz(:, 1:num_icell) * parphi_parsigma(:, 1:num_icell)
      CALL this%gzmTerr%Jacobian(sg, w2, parphi_parz, J_w_phiz)

      ! Calculate parkai_parz
      CALL this%gzmTerr%FirstOrderDerivative(sg, kai, this%gzmTerr%coef_1rst, parkai_parsigma)
      parkai_parz(:, 1:num_icell) = this%gzmTerr%Hz(:, 1:num_icell) * parkai_parsigma(:, 1:num_icell)
      CALL this%gzmTerr%Divergen(sg, w2, parkai_parz, F_w_kaiz)

      !Calculate invrho
      invrho = 1.0D0 / rho
      CALL this%gzmTerr%Divergen(sg, invrho, P, F_invrho_P)

      ! Calculate Tendency of divergence
      ten_div(:, 1:num_iCell) = F_eta_phi(:, 1:num_iCell) &
                                + J_eta_phi(:, 1:num_iCell) &
                                - LKE(:, 1:num_iCell) &
                                + J_w_phiz(:, 1:num_iCell) &
                                - F_w_kaiz(:, 1:num_iCell) &
                                - F_invrho_P(:, 1:num_iCell)
        !!=========================== End of Divergence calculation ==================================
        !!============================ Calculate Vorticity tendency ==================================
      CALL this%gzmTerr%Divergen(sg, eta, kai, F_eta_kai)
      CALL this%gzmTerr%Jacobian(sg, eta, phi, J_eta_phi)
      CALL this%gzmTerr%Divergen(sg, w2, parphi_parz, F_w_phiz)
      CALL this%gzmTerr%Jacobian(sg, w2, parkai_parz, J_w_kaiz)
      CALL this%gzmTerr%Jacobian(sg, rho, P, J_rho_P)

      ten_vor(:, 1:num_iCell) = J_eta_phi(:, 1:num_iCell) &
                                - F_eta_kai(:, 1:num_iCell) &
                                - F_w_phiz(:, 1:num_iCell) &
                                - J_w_kaiz(:, 1:num_iCell) &
                                + 1.0D0 / rho(:, 1:num_iCell)**2 &
                                * J_rho_P(:, 1:num_iCell)
        !!=========================== End of Vorticity calculation =====================================

      CALL this%gzmTerr%destroy()

    END ASSOCIATE

    DEALLOCATE (eta, KE, F_phi_phi, F_kai_kai, &
                Lphi, Lkai, J_phi_kai, F_eta_phi, &
                J_eta_phi, LKE, J_w_phiz, F_w_kaiz, &
                invrho, F_invrho_P, &
                parphi_parsigma, parphi_parz, &
                parkai_parsigma, parkai_parz, &
                J_eta_phi, F_eta_kai, F_w_phiz, &
                J_w_kaiz, J_rho_P, &
                w2)

  END SUBROUTINE

END MODULE
