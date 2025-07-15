MODULE TenTheta_m
  USE gzmTerr_m, ONLY: gzmTerr_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE kinds_m, ONLY: i_kind, r_kind
  USE InterpValue_m, ONLY: InterpValue_t

  IMPLICIT NONE

  TYPE :: TenTheta_t
    TYPE(InterpValue_t) :: InterpValue
    TYPE(gzmTerr_t) :: gzmTerr
  CONTAINS
    PROCEDURE :: TenTheta
  END TYPE TenTheta_t

CONTAINS

  SUBROUTINE TenTheta(this, sg, phi, kai, w, theta, ten_theta)
    TYPE(TenTheta_t) :: this
    TYPE(SingleGrid_t), INTENT(IN) :: sg
    REAL(r_kind), DIMENSION(:, :), INTENT(IN) :: phi, kai, w, theta
    REAL(r_kind), INTENT(OUT) :: ten_theta(:, :)
    REAL(r_kind), DIMENSION(:, :), ALLOCATABLE :: phi2, kai2, J_phi_theta, &
                                                  F_kai_theta, Ltheta, &
                                                  partheta_parsigma, &
                                                  partheta_parz

    ASSOCIATE (num_cell => sg%num_cell, &
               num_iCell => sg%num_icell, &
               vLevel => sg%vLevel)

      ALLOCATE (phi2(vLevel, num_cell), kai2(vLevel, num_cell), &
                J_phi_theta(vLevel, num_cell), F_kai_theta(vLevel, num_cell), &
                Ltheta(vLevel, num_cell), partheta_parsigma(vLevel, num_cell), &
                partheta_parz(vLevel, num_cell))
      CALL this%InterpValue%InterpValue(sg, phi, phi2)
      CALL this%InterpValue%InterpValue(sg, kai, kai2)
      CALL this%gzmTerr%Divergen(sg, kai2, theta, F_kai_theta)
      CALL this%gzmTerr%Jacobian(sg, phi2, theta, J_phi_theta)
      CALL this%gzmTerr%Laplacia(sg, theta, Ltheta)
      CALL this%gzmTerr%FirstOrderDerivative(sg, theta, this%gzmTerr%coef_1rst, partheta_parsigma)
      partheta_parz(:, 1:num_iCell) = this%gzmTerr%Hz(:, 1:num_iCell) * partheta_parsigma(:, 1:num_iCell)

      ten_theta(:, 1:num_iCell) = -J_phi_theta(:, 1:num_iCell) &
                                  - F_kai_theta(:, 1:num_iCell) &
                                  + kai2(:, 1:num_iCell) * Ltheta(:, 1:num_iCell) &
                                  - w(:, 1:num_iCell) * partheta_parz(:, 1:num_iCell)

      CALL this%gzmTerr%destroy()

      DEALLOCATE (phi2, kai2, J_phi_theta, &
                  F_kai_theta, Ltheta, &
                  partheta_parsigma, &
                  partheta_parz)
    END ASSOCIATE
  END SUBROUTINE
END MODULE
