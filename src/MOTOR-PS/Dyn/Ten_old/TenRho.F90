MODULE TenRho_m
  USE gzmTerr_m, ONLY: gzmTerr_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE kinds_m, ONLY: i_kind, r_kind
  USE InterpValue_m, ONLY: InterpValue_t

  IMPLICIT NONE

  TYPE :: TenRho_t
    TYPE(InterpValue_t) :: InterpValue
    TYPE(gzmTerr_t) :: gzmTerr
  CONTAINS
    PROCEDURE :: Tendency
  END TYPE TenRho_t

CONTAINS

  SUBROUTINE Tendency(this, sg1, sg2, div, phi, kai, w, rho, ten_rho)
    IMPLICIT NONE
    TYPE(TenRho_t) :: this
    TYPE(SingleGrid_t), INTENT(IN) :: sg1, sg2
    REAL(r_kind), DIMENSION(:, :), INTENT(IN) :: div, phi, kai, w, rho
    REAL(r_kind), INTENT(OUT) :: ten_rho(:, :)
    REAL(r_kind), DIMENSION(:, :), ALLOCATABLE :: J_phi_rho, F_rho_kai, Lkai, &
                                                  parw_parz, parrho_parz, &
                                                  parw_parsigma, parrho_parsigma, &
                                                  w2
    INTEGER(i_kind) :: k

    ASSOCIATE (num_cell => sg1%num_cell, &
               num_iCell => sg1%num_icell, &
               vLevel => sg1%vLevel, &
               sigma1 => sg1%sigma, &
               sigma2 => sg2%sigma, &
               vLevel2 => sg2%vLevel)
      ALLOCATE (J_phi_rho(vLevel, num_cell), F_rho_kai(vLevel, num_cell), &
                Lkai(vLevel, num_cell), &
                parw_parz(vLevel, num_cell), parrho_parz(vLevel, num_cell), &
                parw_parsigma(vLevel, num_cell), parrho_parsigma(vLevel, num_cell), &
                w2(vLevel, num_cell))

      ! calculate w and parw_parz in the half layer
      w2 = 0.0D0
      !  CALL interp1d(sigma1, sigma2, w(:, 1:num_iCell), w2(:, 1:num_iCell))
      w2(:, 1:num_iCell) = (w(1:vLevel, num_iCell) + w(2:vLevel + 1, num_iCell)) / 2.0D0
      DO k = 1, vLevel2 - 1
        parw_parsigma(k, 1:num_iCell) = (w(k + 1, 1:num_iCell) - w(k, 1:num_iCell)) &
                                        / (sigma2(k + 1) - sigma2(k))
      END DO
      parw_parz(:, 1:num_iCell) = this%gzmTerr%Hz(:, 1:num_iCell) * parw_parsigma(k, 1:num_iCell)

      ! calculate parrho_parz
      CALL this%gzmTerr%initialize(sg1)
      CALL this%gzmTerr%FirstOrderDerivative(sg1, rho, this%gzmTerr%coef_1rst, parrho_parsigma)
      parrho_parz(:, 1:num_iCell) = this%gzmTerr%Hz(:, 1:num_iCell) * parrho_parsigma(k, 1:num_iCell)

      ten_rho(:, 1:num_iCell) = -rho(:, 1:num_iCell) * div(:, 1:num_iCell) &
                                - J_phi_rho(:, 1:num_iCell) &
                                - F_rho_kai(:, 1:num_iCell) &
                                + rho(:, 1:num_iCell) * Lkai(:, 1:num_iCell) &
                                - rho(:, 1:num_iCell) * parw_parz(:, 1:num_iCell) &
                                - w2(:, 1:num_iCell) * parrho_parz(:, 1:num_iCell)

      DEALLOCATE (J_phi_rho, F_rho_kai, Lkai, &
                  parw_parz, parrho_parz, &
                  parw_parsigma, parrho_parsigma, &
                  w2)
      CALL this%gzmTerr%destroy
    END ASSOCIATE

  END SUBROUTINE
