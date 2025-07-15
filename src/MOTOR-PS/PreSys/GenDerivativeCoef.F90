MODULE GenSgDiffCoef_m
  USE Derivative_m, ONLY: FirstOrder_t, FirstOrderHalfLayerShortestDist_t, SecOrderMinNorm_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE kinds_m, ONLY: i_kind, r_kind
  USE gzm_m, ONLY: gzm_t

  IMPLICIT NONE

  TYPE :: GenSgDiffCoef_t

    TYPE(FirstOrder_t) :: FirstDiff
    TYPE(FirstOrderHalfLayerShortestDist_t) :: FirstDiffHalf
    TYPE(SecOrderMinNorm_t) :: SecDiff

  CONTAINS
    PROCEDURE :: GenSgVerDiffCoefs
    PROCEDURE:: GenFirstOrderFull
    PROCEDURE :: GenFirstOrderHalf
    PROCEDURE :: GenSecondOrder
    PROCEDURE :: GenHzParams
    PROCEDURE :: GenLapaceParams

  END TYPE GenSgDiffCoef_t

CONTAINS

  SUBROUTINE GenFirstOrderFull(this, sg)

    IMPLICIT NONE
    CLASS(GenSgDiffCoef_t) :: this
    TYPE(SingleGrid_t), INTENT(INOUT) :: sg

    INTEGER(i_kind) :: i, k, ks

    IF (.NOT. ALLOCATED(sg%coef_fstdif)) ALLOCATE (sg%coef_fstdif(sg%vLevel, sg%num_cell, 3))

    DO i = 1, sg%num_cell
      DO k = 1, sg%vLevel
        IF (k .EQ. 1) THEN
          ks = 1
          CALL this%FirstDiff%FirstOrder(sg%coef_fstdif(k, i, :), &
                                         sg%sigma(k:k + 2), ks)
        ELSEIF (k .EQ. sg%vLevel) THEN
          ks = 3
          CALL this%FirstDiff%FirstOrder(sg%coef_fstdif(k, i, :), &
                                         sg%sigma(k - 2:k), ks)
        ELSE
          ks = 2
          CALL this%FirstDiff%FirstOrder(sg%coef_fstdif(k, i, :), &
                                         sg%sigma(k - 1:k + 1), ks)
        END IF
      END DO
    END DO

  END SUBROUTINE

  SUBROUTINE GenFirstOrderHalf(this, sg_in, sg_out)
    IMPLICIT NONE

    CLASS(GenSgDiffCoef_t) :: this
    TYPE(SingleGrid_t), INTENT(IN) :: sg_in
    TYPE(SingleGrid_t), INTENT(INOUT) :: sg_out
    INTEGER(i_kind) :: i, k

    ALLOCATE (sg_out%coef_fstdif_half(sg_out%vLevel, sg_out%num_cell, 4))

    IF (sg_out%vLevel .LT. sg_in%vLevel) THEN
      DO i = 1, sg_out%num_cell
        DO k = 1, sg_out%vLevel
          IF (k .EQ. 1) THEN
            CALL this%FirstDiffHalf%FirstOrder(sg_in%sigma(k:k + 3), sg_out%sigma(k), sg_out%coef_fstdif_half(k, i, :))
          ELSEIF (k .EQ. sg_out%vLevel) THEN
            CALL this%FirstDiffHalf%FirstOrder(sg_in%sigma(k - 2:k + 1), sg_out%sigma(k), sg_out%coef_fstdif_half(k, i, :))
          ELSE
            CALL this%FirstDiffHalf%FirstOrder(sg_in%sigma(k - 1:k + 2), sg_out%sigma(k), sg_out%coef_fstdif_half(k, i, :))
          END IF
        END DO
      END DO
    ELSE
      DO i = 1, sg_out%num_cell
        DO k = 1, sg_out%vLevel
          IF (k .EQ. 1) THEN
            CALL this%FirstDiffHalf%FirstOrder(sg_in%sigma(k:k + 3), sg_out%sigma(k), sg_out%coef_fstdif_half(k, i, :))
          ELSEIF (k .EQ. 2) THEN
            CALL this%FirstDiffHalf%FirstOrder(sg_in%sigma(k - 1:k + 2), sg_out%sigma(k), sg_out%coef_fstdif_half(k, i, :))
          ELSEIF (k .EQ. sg_out%vLevel) THEN
            CALL this%FirstDiffHalf%FirstOrder(sg_in%sigma(k - 4:k - 1), sg_out%sigma(k), sg_out%coef_fstdif_half(k, i, :))
          ELSEIF (k .EQ. sg_out%vLevel - 1) THEN
            CALL this%FirstDiffHalf%FirstOrder(sg_in%sigma(k - 3:k), sg_out%sigma(k), sg_out%coef_fstdif_half(k, i, :))
          ELSE
            CALL this%FirstDiffHalf%FirstOrder(sg_in%sigma(k - 2:k + 1), sg_out%sigma(k), sg_out%coef_fstdif_half(k, i, :))
          END IF
        END DO
      END DO
    END IF
  END SUBROUTINE

  SUBROUTINE GenSecondOrder(this, sg)

    IMPLICIT NONE
    CLASS(GenSgDiffCoef_t) :: this
    TYPE(SingleGrid_t), INTENT(INOUT) :: sg

    INTEGER(i_kind) :: i, k, ks

    IF (.NOT. ALLOCATED(sg%coef_scddif))  ALLOCATE (sg%coef_scddif(sg%vLevel, sg%num_cell, 5))

    DO i = 1, sg%num_cell
      DO k = 1, sg%vLevel
        IF (k .EQ. 1) THEN
          ks = 1
          CALL this%SecDiff%SecondOrder(sg%coef_scddif(k, i, :), &
                                        sg%sigma(k:k + 4), ks)
        ELSEIF (k .EQ. 2) THEN
          ks = 2
          CALL this%SecDiff%SecondOrder(sg%coef_scddif(k, i, :), &
                                        sg%sigma(k - 1:k + 3), ks)
        ELSEIF (k .EQ. sg%vLevel) THEN
          ks = 5
          CALL this%SecDiff%SecondOrder(sg%coef_scddif(k, i, :), &
                                        sg%sigma(k - 4:k), ks)
        ELSEIF (k .EQ. sg%vLevel - 1) THEN
          ks = 4
          CALL this%SecDiff%SecondOrder(sg%coef_scddif(k, i, :), &
                                        sg%sigma(k - 3:k + 1), ks)
        ELSE
          ks = 3
          CALL this%SecDiff%SecondOrder(sg%coef_scddif(k, i, :), &
                                        sg%sigma(k - 2:k + 2), ks)

        END IF
      END DO
    END DO
  END SUBROUTINE

  SUBROUTINE GenHzParams(this, sg)

    IMPLICIT NONE
    CLASS(GenSgDiffCoef_t) :: this
    TYPE(SingleGrid_t), INTENT(INOUT) :: sg
    INTEGER(i_kind) :: i, k

    IF (ALLOCATED(sg%zHght)) THEN
      DEALLOCATE (sg%zHght)
      ALLOCATE (sg%zHght(sg%vLevel, sg%num_cell))
    END IF
    IF (.NOT. ALLOCATED(sg%parz_parsigma)) ALLOCATE (sg%parz_parsigma(sg%vLevel, sg%num_cell))
    IF (.NOT. ALLOCATED(sg%Hz)) ALLOCATE (sg%Hz(sg%vLevel, sg%num_cell))
    IF (.NOT. ALLOCATED(sg%parHz_parsigma)) ALLOCATE (sg%parHz_parsigma(sg%vLevel, sg%num_cell))

    IF (sg%vLevel > 1) THEN
      ! sg%ztop = sg%sigma(sg%vLevel)
      DO i = 1, sg%vLevel
        sg%zHght(i, :) = sg%sigma(i) * (sg%ztop - sg%topo) / sg%ztop + sg%topo
        sg%Hz(i, :) = sg%ztop / (sg%ztop - sg%topo)
      END DO

      sg%parz_parsigma = 1.0D0 / sg%Hz
      sg%parHz_parsigma = 0.0D0

    ELSE
      sg%zHght = 0.0D0
      sg%sigma = 0.0D0
      sg%parHz_parsigma = 0.0D0
      sg%parz_parsigma = 0.0D0
      sg%Hz = 0.0D0
    END IF

    ! ALLOCATE (sg%parz_parsigma(sg%vLevel, sg%num_cell), &
    !           sg%Hz(sg%vLevel, sg%num_cell))
    ! ! sg%parHz_parsigma(sg%vLevel, sg%num_cell), &
    ! ! parz_parsigma_2nd(sg%vLevel, sg%num_cell))

    ! !   REAL(r_kind), ALLOCATABLE :: parz_parsigma_2nd(:, :)
    ! DO i = 1, sg%num_cell
    !    DO k = 1, sg%vLevel
    !       IF (k .EQ. 1) THEN
    !          sg%parz_parsigma(k, i) = SUM(sg%coef_fstdif(k, i, :)*sg%zHght(k:k + 2, i))
    !       ELSEIF (k .EQ. sg%vLevel) THEN
    !          sg%parz_parsigma(k, i) = SUM(sg%coef_fstdif(k, i, :)*sg%zHght(k - 2:k, i))
    !       ELSE
    !          sg%parz_parsigma(k, i) = SUM(sg%coef_fstdif(k, i, :)*sg%zHght(k - 1:k + 1, i))
    !       END IF
    !    END DO
    ! END DO

    ! !   DO i = 1, sg%num_cell
    ! !      DO k = 1, sg%vLevel
    ! !         IF (k .EQ. 1) THEN
    ! !            parz_parsigma_2nd(k, i) = SUM(sg%coef_scddif(k, i, :)*sg%zHght(k:k + 4, i))
    ! !         ELSEIF (k .EQ. 2) THEN
    ! !            parz_parsigma_2nd(k, i) = SUM(sg%coef_scddif(k, i, :)*sg%zHght(k - 1:k + 3, i))
    ! !         ELSEIF (k .EQ. sg%vLevel) THEN
    ! !            parz_parsigma_2nd(k, i) = SUM(sg%coef_scddif(k, i, :)*sg%zHght(k - 4:k, i))
    ! !         ELSEIF (k .EQ. sg%vLevel - 1) THEN
    ! !            parz_parsigma_2nd(k, i) = SUM(sg%coef_scddif(k, i, :)*sg%zHght(k - 3:k + 1, i))
    ! !         ELSE
    ! !            parz_parsigma_2nd(k, i) = SUM(sg%coef_scddif(k, i, :)*sg%zHght(k - 2:k + 2, i))
    ! !         END IF
    ! !      END DO
    ! !   END DO

    ! sg%Hz = 1.0D0/sg%parz_parsigma

    ! !   sg%parHz_parsigma = -1.0D0*sg%Hz**2*parz_parsigma_2nd

    ! ! DO i = 1, sg%num_cell
    ! !    DO k = 1, sg%vLevel
    ! !       IF (k .EQ. 1) THEN
    ! !          sg%parHz_parsigma(k, i) = SUM(sg%coef_fstdif(k, i, :)*sg%Hz(k:k + 2, i))
    ! !       ELSEIF (k .EQ. sg%vLevel) THEN
    ! !          sg%parHz_parsigma(k, i) = SUM(sg%coef_fstdif(k, i, :)*sg%Hz(k - 2:k, i))
    ! !       ELSE
    ! !          sg%parHz_parsigma(k, i) = SUM(sg%coef_fstdif(k, i, :)*sg%Hz(k - 1:k + 1, i))
    ! !       END IF
    ! !    END DO
    ! ! END DO

    ! !   DEALLOCATE (parz_parsigma_2nd)
  END SUBROUTINE

  SUBROUTINE GenLapaceParams(this, sg, DA_flag)
    IMPLICIT NONE

    CLASS(GenSgDiffCoef_t) :: this
    TYPE(SingleGrid_t), INTENT(INOUT) :: sg
    TYPE(gzm_t) :: gzm
    CHARACTER(*), INTENT(IN), OPTIONAL :: DA_flag

    ! Yuanfu Xie added the deallocations to prevent redundant allocation 2025-01-10:
    IF (ALLOCATED(sg%F_z_z)) DEALLOCATE(sg%F_z_z)
    IF (ALLOCATED(sg%Lz)) DEALLOCATE(sg%Lz)
    IF (ALLOCATED(sg%F_Hz_z)) DEALLOCATE(sg%F_Hz_z)
    IF (ALLOCATED(sg%F_invHz_z)) DEALLOCATE(sg%F_invHz_z)
    ALLOCATE (sg%F_z_z(sg%vLevel, sg%num_cell), &
              sg%Lz(sg%vLevel, sg%num_cell), &
              sg%F_Hz_z(sg%vLevel, sg%num_cell), &
              sg%F_invHz_z(sg%vLevel, sg%num_cell))

    sg%F_Hz_z = 0.0D0
    sg%F_z_z = 0.0D0
    sg%Lz = 0.0D0
    sg%F_invHz_z = 0.0D0

    ! IF (PRESENT(DA_flag)) THEN
    ! IF (.NOT. ALLOCATED(sg%bdy_type)) ALLOCATE (sg%bdy_type(sg%num_cell))
    ! sg%bdy_type = sg%cell_type
    ! END IF

    gzm = gzm_t(sg)

    CALL gzm%Divergen(sg%Hz, sg%zHght, sg%F_Hz_z)
    CALL gzm%Divergen(sg%zHght, sg%zHght, sg%F_z_z)
    CALL gzm%Divergen(sg%parz_parsigma, sg%zHght, sg%F_invHz_z)
    CALL gzm%Laplacia(sg%zHght, sg%Lz)

    CALL sg%ExchangeMatOnHalo2D(sg%vLevel, sg%F_Hz_z)
    CALL sg%ExchangeMatOnHalo2D(sg%vLevel, sg%F_z_z)
    CALL sg%ExchangeMatOnHalo2D(sg%vLevel, sg%F_invHz_z)
    CALL sg%ExchangeMatOnHalo2D(sg%vLevel, sg%Lz)

  END SUBROUTINE

  SUBROUTINE GenSgVerDiffCoefs(this, sg_other, sg_this)

    IMPLICIT NONE

    CLASS(GenSgDiffCoef_t) :: this
    TYPE(SingleGrid_t), INTENT(IN) :: sg_other
    TYPE(SingleGrid_t), INTENT(INOUT) :: sg_this

    CALL this%GenFirstOrderFull(sg_this)
    CALL this%GenFirstOrderHalf(sg_other, sg_this)
    CALL this%GenSecondOrder(sg_this)
    CALL this%GenHzParams(sg_this)
    CALL this%GenLapaceParams(sg_this)

  END SUBROUTINE

END MODULE GenSgDiffCoef_m
