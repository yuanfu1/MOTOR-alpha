!!--------------------------------------------------------------------------------------------------
! PROJECT           : gzmTerr - generalized Z-grid model in the Terrain following system
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Jilong Chen
! VERSION           : V 0.0
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module provides the set of operators needed for a Z-grid model

MODULE gzmTerr_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE geometry_m, ONLY: geometry_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE Multigrid_m, ONLY: MultiGrid_t
  USE state_m, ONLY: state_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE gzm_m, ONLY: gzm_t
  USE Field_m, ONLY: Field_t
  ! USE vertexSpace_m, ONLY : vertexSpace_t

  IMPLICIT NONE

  TYPE :: gzmTerr_t
    TYPE(gzm_t) :: gzm
    REAL(r_kind), POINTER :: z(:, :), sigma(:), &
                             coef_1rst(:, :, :), coef_2end(:, :, :), &
                             Hz(:, :), invHz(:, :), &
                             parHz_sigma(:, :), &
                             F_z_z(:, :), Lz(:, :), &
                             F_Hz_z(:, :), F_invHz_z(:, :)
    ! REAL(r_kind), ALLOCATABLE :: parA_parsigma_1rst(:, :), &
    !                              !   parA_parsigma_2end(:, :), &
    !                              parB_parsigma_1rst(:, :)
    !   parB_parsigma_2end(:, :)

  CONTAINS
    PROCEDURE :: initialize
    PROCEDURE :: destroy
    PROCEDURE :: Divergen
    PROCEDURE :: Jacobian
    PROCEDURE :: Laplacia
    PROCEDURE :: FirstOrderDerivative
    PROCEDURE :: SecondOrderDerivative
  END TYPE gzmTerr_t

CONTAINS

  SUBROUTINE initialize(this, sg)
    IMPLICIT NONE
    CLASS(gzmTerr_t), INTENT(INOUT) :: this
    TYPE(SingleGrid_t), INTENT(IN), TARGET :: sg
    ! REAL(r_kind), INTENT(IN) :: A(:, :)
    ! REAL(r_kind), INTENT(IN), optional :: B(:, :)
    ! CHARACTER(LEN=1024), INTENT(IN) :: configFile

    ! mpdd initialization:
    ! CALL this%gzm%initialize(configFile)
    ! load coefs

    this%z => sg%zHght
    this%sigma => sg%sigma
    this%coef_1rst => sg%coef_fstdif
    this%coef_2end => sg%coef_scddif
    this%Hz => sg%Hz
    this%invHz => sg%parz_parsigma
    this%parHz_sigma => sg%parHz_parsigma
    this%F_z_z => sg%F_z_z
    this%Lz => sg%Lz
    this%F_Hz_z => sg%F_Hz_z
    this%F_invHz_z => sg%F_invHz_z

    ! calculate first and second order

    ! IF (PRESENT(B)) THEN

    !    ASSOCIATE (num_cell => sg%num_cell, &
    !               vLevel => sg%vLevel)
    !       ALLOCATE (this%parA_parsigma_1rst(vLevel, num_cell), &
    !                 this%parB_parsigma_1rst(vLevel, num_cell))
    !       this%parA_parsigma_1rst = 0.0D0
    !       ! this%parA_parsigma_2end = 0.0D0
    !       this%parB_parsigma_1rst = 0.0D0
    !       ! this%parB_parsigma_2end = 0.0D0
    !       CALL this%FirstOrderDerivative(sg, A, this%coef_1rst, this%parA_parsigma_1rst)
    !       ! CALL this%SecondOrderDerivative(sg, A, this%coef_2end, this%parA_parsigma_2end)
    !       CALL this%FirstOrderDerivative(sg, B, this%coef_1rst, this%parB_parsigma_1rst)
    !       ! CALL this%SecondOrderDerivative(sg, B, this%coef_2end, this%parB_parsigma_2end)
    !    END ASSOCIATE
    ! END IF

  END SUBROUTINE initialize

  SUBROUTINE destroy(this)

    IMPLICIT NONE
    CLASS(gzmTerr_t) :: this

    PRINT *, 'gzmTerr destructor works'

    IF (ASSOCIATED(this%sigma)) THEN
      NULLIFY (this%sigma)
      NULLIFY (this%z)
      NULLIFY (this%coef_1rst)
      NULLIFY (this%coef_2end)
      NULLIFY (this%invHz)
      NULLIFY (this%Hz)
      NULLIFY (this%parHz_sigma)
      NULLIFY (this%Lz)
      NULLIFY (this%F_z_z)
      NULLIFY (this%F_Hz_z)
      NULLIFY (this%F_invHz_z)
    END IF
    ! IF (ALLOCATED(this%parA_parsigma_1rst)) THEN
    !    DEALLOCATE (this%parA_parsigma_1rst, this%parB_parsigma_1rst)
    !    ! this%parA_parsigma_2end, , this%parB_parsigma_2end)
    ! END IF

  END SUBROUTINE destroy

  SUBROUTINE Divergen(this, sg, A, B, value_out)

    IMPLICIT NONE

    CLASS(gzmTerr_t), INTENT(IN) :: this
    TYPE(SingleGrid_t), INTENT(IN) :: sg
    REAL(r_kind), INTENT(IN) :: A(:, :), B(:, :)
    REAL(r_kind), INTENT(OUT) :: value_out(:, :)
    INTEGER(i_kind) :: i
    REAL(r_kind), ALLOCATABLE :: F_A_B(:, :), LprimeB(:, :), &
                                 F_A_z(:, :), F_B_z(:, :), LB(:, :), &
                                 parA_parsigma_1rst(:, :), &
                                 parB_parsigma_1rst(:, :)

    ASSOCIATE (num_icell => sg%num_icell, &
               num_cell => sg%num_cell, &
               cell_type => sg%cell_type, &
               vLevel => sg%vLevel)
      ALLOCATE (F_A_B(vLevel, num_cell), LprimeB(vLevel, num_cell), &
                F_A_z(vLevel, num_cell), F_B_z(vLevel, num_cell), &
                LB(vLevel, num_cell), &
                parA_parsigma_1rst(vLevel, num_cell), &
                parB_parsigma_1rst(vLevel, num_cell))
      F_A_B = 0.0D0
      LprimeB = 0.0D0
      F_A_z = 0.0D0
      F_B_z = 0.0D0
      LB = 0.0D0
      parA_parsigma_1rst = 0.0D0
      parB_parsigma_1rst = 0.0D0
      CALL this%FirstOrderDerivative(sg, A, this%coef_1rst, parA_parsigma_1rst)
      CALL this%FirstOrderDerivative(sg, B, this%coef_1rst, parB_parsigma_1rst)
      CALL this%gzm%Divergen(sg, A, B, F_A_B)
      CALL this%gzm%Laplacia(sg, B, LprimeB)
      CALL this%gzm%Divergen(sg, B, this%z, F_B_z)
      CALL this%gzm%Divergen(sg, A, this%z, F_A_z)
      CALL this%Laplacia(sg, B, LB)

      DO i = 1, num_icell
        IF (cell_type(i) .EQ. 0) THEN
          value_out(:, i) = F_A_B(:, i) - A(:, i) * LprimeB(:, i) &
                            - this%Hz(:, i) * parB_parsigma_1rst(:, i) * F_A_z(:, i) &
                            - this%Hz(:, i) * parA_parsigma_1rst(:, i) * F_B_z(:, i) &
                            + this%Hz(:, i)**2 * parB_parsigma_1rst(:, i) &
                            * parA_parsigma_1rst(:, i) * this%F_z_z(:, i) &
                            + this%Hz(:, i) * (A(:, i) * parB_parsigma_1rst(:, i) &
                                               + B(:, i) * parA_parsigma_1rst(:, i) &
                                               - this%z(:, i) * this%Hz(:, i) &
                                               * parA_parsigma_1rst(:, i) &
                                               * parB_parsigma_1rst(:, i)) * this%Lz(:, i) &
                            + A(:, i) * LB(:, i)
        END IF
      END DO
    END ASSOCIATE

    DEALLOCATE (F_A_B, LprimeB, F_A_z, F_B_z, LB)

  END SUBROUTINE Divergen

  SUBROUTINE Jacobian(this, sg, A, B, value_out)

    IMPLICIT NONE

    CLASS(gzmTerr_t), INTENT(IN) :: this
    TYPE(SingleGrid_t), INTENT(IN) :: sg
    REAL(r_kind), INTENT(IN) :: A(:, :), B(:, :)
    REAL(r_kind), INTENT(OUT) :: value_out(:, :)
    INTEGER :: i
    REAL(r_kind), ALLOCATABLE :: J_A_B(:, :), J_A_z(:, :), &
                                 J_z_B(:, :), &
                                 parA_parsigma_1rst(:, :), &
                                 parB_parsigma_1rst(:, :)

    ASSOCIATE (num_icell => sg%num_icell, &
               num_cell => sg%num_cell, &
               cell_type => sg%cell_type, &
               vLevel => sg%vLevel)

      ALLOCATE (J_A_B(vLevel, num_cell), J_A_z(vLevel, num_cell), &
                J_z_B(vLevel, num_cell), &
                parA_parsigma_1rst(vLevel, num_cell), &
                parB_parsigma_1rst(vLevel, num_cell))
      J_A_B = 0.0D0
      J_A_z = 0.0D0
      J_z_B = 0.0D0
      parA_parsigma_1rst = 0.0D0
      parB_parsigma_1rst = 0.0D0
      CALL this%FirstOrderDerivative(sg, A, this%coef_1rst, parA_parsigma_1rst)
      CALL this%FirstOrderDerivative(sg, B, this%coef_1rst, parB_parsigma_1rst)
      CALL this%gzm%Jacobian(sg, A, B, J_A_B)
      CALL this%gzm%Jacobian(sg, A, this%z, J_A_z)
      CALL this%gzm%Jacobian(sg, this%z, B, J_z_B)
      DO i = 1, num_icell
        value_out(:, i) = J_A_B(:, i) &
                          - this%Hz(:, i) * parB_parsigma_1rst(:, i) * J_A_z(:, i) &
                          - this%Hz(:, i) * parA_parsigma_1rst(:, i) * J_z_B(:, i)
      END DO
      DEALLOCATE (J_A_B, J_A_z, J_z_B)
      ! CALL sg%ExchangeMatOnHalo2D(vLevel, value_out)
    END ASSOCIATE

  END SUBROUTINE Jacobian

  SUBROUTINE Laplacia(this, sg, opr, value_out)

    IMPLICIT NONE

    CLASS(gzmTerr_t) :: this
    TYPE(SingleGrid_t), INTENT(IN) :: sg
    REAL(r_kind), INTENT(IN) :: opr(:, :)
    ! These operands are in cell space:
    REAL(r_kind), INTENT(OUT) :: value_out(:, :)
    REAL(r_kind), ALLOCATABLE :: Lopr(:, :), &
                                 F_oprsigma_z(:, :), &
                                 paropr_sigma_1rst(:, :), &
                                 paropr_sigma_2end(:, :)
    INTEGER(i_kind) :: i

    ASSOCIATE (num_icell => sg%num_icell, &
               num_cell => sg%num_cell, &
               cell_type => sg%cell_type, &
               vLevel => sg%vLevel)

      ALLOCATE (Lopr(vLevel, num_cell), &
                F_oprsigma_z(vLevel, num_cell), &
                paropr_sigma_1rst(vLevel, num_cell), &
                paropr_sigma_2end(vLevel, num_cell))

      Lopr = 0.0D0; 
      F_oprsigma_z = 0.0D0
      paropr_sigma_1rst = 0.0D0
      paropr_sigma_2end = 0.0D0

      CALL this%FirstOrderDerivative(sg, opr, this%coef_1rst, paropr_sigma_1rst)
      CALL this%SecondOrderDerivative(sg, opr, this%coef_2end, paropr_sigma_2end)

      CALL this%gzm%Divergen(sg, paropr_sigma_1rst, this%z, F_oprsigma_z)

      DO i = 1, num_icell
        IF (cell_type(i) .EQ. 0) THEN
          value_out(:, i) = Lopr(:, i) &
                            - this%Hz(:, i) * (this%z(:, i) * this%Hz(:, i) &
                                               * paropr_sigma_2end(:, i) &
                                               - paropr_sigma_1rst(:, i) &
                                               + this%z(:, i) &
                                               * paropr_sigma_1rst(:, i) &
                                               * this%parHz_sigma(:, i)) * this%Lz(:, i) &
                            + this%Hz(:, i) * (this%Hz(:, i) * paropr_sigma_2end(:, i) &
                                               + paropr_sigma_1rst(:, i) &
                                               * this%parHz_sigma(:, i)) * this%F_z_z(:, i) &
                            - 2 * this%Hz(:, i) * F_oprsigma_z(:, i) &
                            + this%Hz(:, i)**2 * paropr_sigma_1rst(:, i) * this%F_invHz_z(:, i) &
                            - paropr_sigma_1rst(:, i) * this%F_Hz_z(:, i)
        END IF
      END DO
      ! CALL sg%ExchangeMatOnHalo2D(vLevel, value_out)
    END ASSOCIATE
    DEALLOCATE (Lopr, &
                F_oprsigma_z, &
                paropr_sigma_1rst, &
                paropr_sigma_2end)
  END SUBROUTINE Laplacia

  SUBROUTINE FirstOrderDerivative(this, sg, A, coef, parA_sigma)
    CLASS(gzmTerr_t), INTENT(IN) :: this
    TYPE(SingleGrid_t), INTENT(IN) :: sg
    REAL(r_kind), INTENT(IN) :: coef(:, :, :), A(:, :)
    REAL(r_kind), INTENT(OUT) :: parA_sigma(:, :)
    INTEGER(i_kind) :: i, k

    ASSOCIATE (num_icell => sg%num_icell, &
               cell_type => sg%cell_type, &
               vLevel => sg%vLevel)
      parA_sigma = 0.0D0
      DO i = 1, num_icell
        IF (cell_type(i) .EQ. 0) THEN
          DO k = 1, vLevel
            IF (k .EQ. 1) THEN
              parA_sigma(k, i) = SUM(coef(k, i, :) * A(k:k + 2, i))
            ELSEIF (k .EQ. vLevel) THEN
              parA_sigma(k, i) = SUM(coef(k, i, :) * A(k - 2:k, i))
            ELSE
              parA_sigma(k, i) = SUM(coef(k, i, :) * A(k - 1:k + 1, i))
            END IF
          END DO
        END IF
      END DO
      ! CALL sg%ExchangeMatOnHalo2D(vLevel, parA_sigma)
    END ASSOCIATE
  END SUBROUTINE FirstOrderDerivative

  SUBROUTINE SecondOrderDerivative(this, sg, A, coef, parA_sigma_2nd)
    CLASS(gzmTerr_t), INTENT(IN) :: this
    TYPE(SingleGrid_t), INTENT(IN) :: sg
    REAL(r_kind), INTENT(IN) :: coef(:, :, :), A(:, :)
    REAL(r_kind), INTENT(OUT) :: parA_sigma_2nd(:, :)
    INTEGER(i_kind) :: i, k

    ASSOCIATE (num_icell => sg%num_icell, &
               cell_type => sg%cell_type, &
               vLevel => sg%vLevel)

      parA_sigma_2nd = 0.0D0
      DO i = 1, num_icell
        IF (cell_type(i) .EQ. 0) THEN
          DO k = 1, vLevel
            IF (k .EQ. 1) THEN
              parA_sigma_2nd(k, i) = SUM(coef(k, i, :) * A(k:k + 4, i))
            ELSEIF (k .EQ. 2) THEN
              parA_sigma_2nd(k, i) = SUM(coef(k, i, :) * A(k - 1:k + 3, i))
            ELSEIF (k .EQ. vLevel) THEN
              parA_sigma_2nd(k, i) = SUM(coef(k, i, :) * A(k - 4:k, i))
            ELSEIF (k .EQ. vLevel - 1) THEN
              parA_sigma_2nd(k, i) = SUM(coef(k, i, :) * A(k - 3:k + 1, i))
            ELSE
              parA_sigma_2nd(k, i) = SUM(coef(k, i, :) * A(k - 2:k + 2, i))
            END IF
          END DO
        END IF
      END DO
    END ASSOCIATE
    ! CALL sg%ExchangeMatOnHalo2D(vLevel, parA_sigma_2nd)
  END SUBROUTINE

END MODULE gzmTerr_m
