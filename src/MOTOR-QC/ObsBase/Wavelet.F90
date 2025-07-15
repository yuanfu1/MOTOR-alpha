!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-QC.Wavelet
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Ting Shu
! FUNCTION          : 2D symlet-2 wavelet decomposition and reconstruction
! VERSION           : V 0.0
! HISTORY           :
!   Created by Ting Shu (shuting@gbamwf.com), 2023/08/08, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module implements 2D symlet-2 wavelet decomposition and reconstruction for satellite data.
MODULE wavelet_m
    USE kinds_m, ONLY: r_kind, i_kind
    IMPLICIT NONE
    TYPE coeff_cd
        REAL(r_kind), DIMENSION(:, :), ALLOCATABLE :: cd1, cd2, cd3
    END TYPE coeff_cd
    ! Cofficient type
    ! Shapes increase, size of [cd1, cd2, cd3] = level
    ! [ca, [cd1, cd2, cd3], [cd1, cd2, cd3], ...]
    TYPE wavelet_t
        REAL(r_kind), DIMENSION(:, :), ALLOCATABLE :: ca
        TYPE(coeff_cd), DIMENSION(:), ALLOCATABLE :: cds

        CONTAINS
        PROCEDURE, PUBLIC :: decom => decomWave2d_m
        PROCEDURE, PUBLIC :: recon => reconWave2d_m
        PROCEDURE, PUBLIC :: coeff2matrix => coeff_2d
        PROCEDURE, PUBLIC :: matrix2coeff => d2_coeff
        PROCEDURE, PUBLIC :: release => deallocate_coeff
    END TYPE wavelet_t
    PRIVATE :: decomWave2d_m, reconWave2d_m, coeff_2d, d2_coeff, deallocate_coeff, &
               decomWave, decomWave2d, reconWave, reconWave2d, reconWave2d_two
CONTAINS
! Multiple-level two-dimensional wavelet decomposition
    ! Based on decomWave2d
    SUBROUTINE decomWave2d_m(this, x, level)
        IMPLICIT NONE
        CLASS(wavelet_t), INTENT(INOUT) :: this
        REAL(r_kind), DIMENSION(:, :), INTENT(IN) :: x
        INTEGER(i_kind), INTENT(IN) :: level

        ! Local variables
        INTEGER(i_kind) :: i, status, a, b, a2, b2
        REAL(r_kind), DIMENSION(:, :), ALLOCATABLE :: temp_ca, temp_x
        
        a = UBOUND(x, 1)
        b = UBOUND(x, 2)
        ! Allocate coefficients
        ALLOCATE(this%cds(level), temp_x(a, b), STAT=status)
        temp_x = x

        DO i = 1, level
            a2 = FLOOR((a-1)/2.0)+2
            b2 = FLOOR((b-1)/2.0)+2
            ! Allocate coefficients
            ALLOCATE(temp_ca(a2, b2), STAT=status)
            ALLOCATE(this%cds(level-i+1)%cd1(a2, b2), &
                    this%cds(level-i+1)%cd2(a2, b2), &
                    this%cds(level-i+1)%cd3(a2, b2), STAT=status)
            CALL decomWave2d(temp_x, temp_ca, this%cds(level-i+1)%cd1, &
                            this%cds(level-i+1)%cd2, this%cds(level-i+1)%cd3)
            IF (i .EQ. level) THEN
                ALLOCATE(this%ca(a2, b2), STAT=status)
                this%ca = temp_ca
            END IF
            IF(ALLOCATED(temp_x)) DEALLOCATE(temp_x, STAT=status)
            a = a2
            b = b2
            ALLOCATE(temp_x(a, b), STAT=status)
            temp_x = temp_ca
            IF(ALLOCATED(temp_ca)) DEALLOCATE(temp_ca, STAT=status)
        END DO
    END SUBROUTINE decomWave2d_m

    ! Multiple-level two-dimension
    ! Based on reconWave2d, corresponding to decomWave2d_m
    SUBROUTINE reconWave2d_m(this, re_level, re_x)
        IMPLICIT NONE
        CLASS(wavelet_t), INTENT(INOUT) :: this
        INTEGER(i_kind), INTENT(IN) :: re_level
        REAL(r_kind), DIMENSION(:, :), ALLOCATABLE, INTENT(OUT) :: re_x
        ! Local variables
        INTEGER(i_kind) :: c_level, i, start_i, c_a, c_b, status
        REAL(r_kind), DIMENSION(:, :), ALLOCATABLE :: temp_x, temp_ca

        c_level = SIZE(this%cds)
        ! If re_level > level, then exist
        IF (re_level .GT. c_level) THEN
            WRITE(*, *) "Please enter a reconstruction level not larger than coefficient level: ", c_level
            RETURN
        END IF

        c_a = UBOUND(this%ca, 1)
        c_b = UBOUND(this%ca, 2)
        ALLOCATE(temp_ca(c_a, c_b), STAT=status)
        temp_ca = this%ca
        IF (re_level .EQ. 0) THEN
            CALL reconWave2d_two(this, 0, temp_ca, .TRUE., temp_x)
            start_i = 2
            IF (c_level .LT. start_i) THEN
                ALLOCATE(re_x(SIZE(temp_x, 1), SIZE(temp_x, 2)), STAT=status)
                re_x = temp_x
                IF(ALLOCATED(temp_ca)) DEALLOCATE(temp_ca, STAT=status)
                IF(ALLOCATED(temp_x)) DEALLOCATE(temp_x, STAT=status)
                RETURN
            END IF
        ELSE
            ALLOCATE(temp_x(c_a, c_b), STAT=status)
            temp_x = temp_ca
            start_i = 1
        END IF
        DO i = start_i, c_level
            c_a = UBOUND(temp_x, 1)
            c_b = UBOUND(temp_x, 2)
            IF(ALLOCATED(temp_ca)) DEALLOCATE(temp_ca, STAT=status)
            ALLOCATE(temp_ca(c_a, c_b), STAT=status)
            temp_ca = temp_x
            IF(ALLOCATED(temp_x)) DEALLOCATE(temp_x, STAT=status)
            ! For i = 0 - re_level: zero_flag=False, use original cds
            ! For i = re_level+1 - c_level: zero_flag=True, use zero cds
            CALL reconWave2d_two(this, i, temp_ca, i.GT.re_level, temp_x)
        END DO
        ALLOCATE(re_x(SIZE(temp_x, 1), SIZE(temp_x, 2)), STAT=status)
        re_x = temp_x
        IF(ALLOCATED(temp_ca)) DEALLOCATE(temp_ca, STAT=status)
        IF(ALLOCATED(temp_x)) DEALLOCATE(temp_x, STAT=status)
    END SUBROUTINE reconWave2d_m

    ! Convert wavelet_coeff type to 2d data
    SUBROUTINE coeff_2d(this, coeffMatrix, coeffIdx)
        IMPLICIT NONE
        CLASS(wavelet_t), INTENT(INOUT) :: this
        REAL(r_kind), DIMENSION(:, :), ALLOCATABLE, INTENT(OUT) :: coeffMatrix
        INTEGER(i_kind), DIMENSION(:, :), ALLOCATABLE, INTENT(OUT) :: coeffIdx
        ! Local variables
        INTEGER(i_kind) :: i, level, status, iS, iE, jS, jE

        level = SIZE(this%cds)
        ! Allocate coefficient matrix index vector
        ALLOCATE(coeffIdx(level+1, 2), STAT=status)
        coeffIdx(1, 1) = UBOUND(this%ca, 1)
        coeffIdx(1, 2) = UBOUND(this%ca, 2)
        DO i = 2, level+1
            coeffIdx(i, 1) = UBOUND(this%cds(i-1)%cd2, 1)
            coeffIdx(i, 2) = UBOUND(this%cds(i-1)%cd1, 2)
        END DO

        ! Allocate and assign coefficient matrix
        ALLOCATE(coeffMatrix(SUM(coeffIdx(:, 1)), SUM(coeffIdx(:, 2))), STAT=status)
        coeffMatrix = 0.0D0
        coeffMatrix(1:coeffIdx(1,1), 1:coeffIdx(1,2)) = this%ca
        DO i = 2, level+1
            iS = SUM(coeffIdx(1:i-1, 1)) + 1
            iE = SUM(coeffIdx(1:i, 1))
            jS = SUM(coeffIdx(1:i-1, 2)) + 1
            jE = SUM(coeffIdx(1:i, 2))
            ! cd1 at right top
            coeffMatrix(1:coeffIdx(i, 1), jS:jE) = this%cds(i-1)%cd1
            ! cd2 at left bottom
            coeffMatrix(iS:iE, 1:coeffIdx(i, 2)) = this%cds(i-1)%cd2
            ! cd3 at right bottom
            coeffMatrix(iS:iE, jS:jE) = this%cds(i-1)%cd3
        END DO
    END SUBROUTINE coeff_2d

    ! Convert 2d coefficient to type
    SUBROUTINE d2_coeff(this, coeffMatrix, coeffIdx)
        IMPLICIT NONE
        CLASS(wavelet_t), INTENT(INOUT) :: this
        REAL(r_kind), DIMENSION(:, :), INTENT(IN) :: coeffMatrix
        INTEGER(i_kind), DIMENSION(:, :), INTENT(IN) :: coeffIdx
        ! Local variables
        INTEGER(i_kind) :: i, level, a, b, status, iS, iE, jS, jE

        ! Get level
        level = UBOUND(coeffIdx, 1) - 1
        ! Allocate and assign ca
        a = coeffIdx(1, 1)
        b = coeffIdx(1, 2)
        ALLOCATE(this%ca(a, b), this%cds(level), STAT=status)
        this%ca = coeffMatrix(1:a, 1:b)
        DO i = 2, level+1
            iS = SUM(coeffIdx(1:i-1, 1)) + 1
            iE = SUM(coeffIdx(1:i, 1))
            jS = SUM(coeffIdx(1:i-1, 2)) + 1
            jE = SUM(coeffIdx(1:i, 2))
            a = coeffIdx(i, 1)
            b = coeffIdx(i, 2)
            ALLOCATE(this%cds(i-1)%cd1(a, b), this%cds(i-1)%cd2(a, b), this%cds(i-1)%cd3(a, b), STAT=status)
            ! cd1 at right top
            this%cds(i-1)%cd1 = coeffMatrix(1:a, jS:jE)
            ! cd2 at left bottom
            this%cds(i-1)%cd2 = coeffMatrix(iS:iE, 1:b)
            ! cd3 at right bottom
            this%cds(i-1)%cd3 = coeffMatrix(iS:iE, jS:jE)
        END DO
    END SUBROUTINE d2_coeff

    ! Deallocate wavelet_coefficient data
    SUBROUTINE deallocate_coeff(this)
        IMPLICIT NONE
        CLASS(wavelet_t), INTENT(INOUT) :: this
        ! Local variables
        INTEGER(i_kind) :: i, status, level

        ! Get coefficient level
        level = SIZE(this%cds)
        ! Deallocate ca
        IF(ALLOCATED(this%ca)) DEALLOCATE(this%ca, STAT=status)
        DO i = 1, level
            IF(ALLOCATED(this%cds(i)%cd1)) DEALLOCATE(this%cds(i)%cd1, STAT=status)
            IF(ALLOCATED(this%cds(i)%cd2)) DEALLOCATE(this%cds(i)%cd2, STAT=status)
            IF(ALLOCATED(this%cds(i)%cd3)) DEALLOCATE(this%cds(i)%cd3, STAT=status)
        END DO
        ! Deallocate cds
        IF(ALLOCATED(this%cds)) DEALLOCATE(this%cds, STAT=status)
    END SUBROUTINE deallocate_coeff

    ! One-dimensional wavelet decomposition
    ! Symlet2, mode='zero'
    SUBROUTINE decomWave(x, ca, cd)
        IMPLICIT NONE
        REAL(r_kind), DIMENSION(:), INTENT(IN) :: x
        REAL(r_kind), DIMENSION(:), INTENT(OUT) :: ca, cd
        ! Local variables, low and high symlet2 coefficients
        REAL(r_kind), PARAMETER, DIMENSION(4) :: lf = (/-0.12940952255092145, 0.22414386804185735, &
                                                        0.836516303737469, 0.48296291314469025/), &
                                                 hf = (/-0.48296291314469025, 0.836516303737469, &
                                                        -0.22414386804185735, -0.12940952255092145/)
        INTEGER(i_kind) :: len_x, i, status, idx
        REAL(r_kind), DIMENSION(:), ALLOCATABLE :: conv_x

        ! Get length of original signal x
        len_x = SIZE(x)
        ! Allocate variables
        ALLOCATE(conv_x(len_x+6), STAT=status)
        conv_x = 0.0D0
        conv_x(4:len_x+3) = x

        idx = 0
        DO i = 1, len_x+3
            IF (MOD(i, 2) .EQ. 0) THEN
                idx = idx + 1
                ca(idx) = conv_x(i)*lf(4) + conv_x(i+1)*lf(3) + conv_x(i+2)*lf(2) + conv_x(i+3)*lf(1)
                cd(idx) = conv_x(i)*hf(4) + conv_x(i+1)*hf(3) + conv_x(i+2)*hf(2) + conv_x(i+3)*hf(1)
            END IF
        END DO

        IF(ALLOCATED(conv_x)) DEALLOCATE(conv_x, STAT=status)
    END SUBROUTINE decomWave

    ! Two-dimensional wavelet decomposition
    ! Based on decomWave
    SUBROUTINE decomWave2d(x, ca, cd1, cd2, cd3)
        IMPLICIT NONE
        REAL(r_kind), DIMENSION(:, :), INTENT(IN) :: x
        REAL(r_kind), DIMENSION(:, :), INTENT(OUT) :: ca, cd1, cd2, cd3
        ! Local variables
        INTEGER(i_kind) :: n_rows, n_columns, n_new_rows, n_new_columns, i, status
        REAL(r_kind), DIMENSION(:, :), ALLOCATABLE :: r_ca, r_cd

        n_rows = UBOUND(x, 1)
        n_columns = UBOUND(x, 2)
        n_new_rows = FLOOR((n_rows-1)/2.0)+2
        n_new_columns = FLOOR((n_columns-1)/2.0)+2
        ! Allocate r_ca and r_cd
        ALLOCATE(r_ca(n_rows, n_new_columns), r_cd(n_rows, n_new_columns), STAT=status)
        ! Get row convolution arrays
        DO i = 1, n_rows
            CALL decomWave(x(i, :), r_ca(i, :), r_cd(i, :))
        END DO
        DO i = 1, n_new_columns
            CALL decomWave(r_ca(:, i), ca(:, i), cd1(:, i))
            CALL decomWave(r_cd(:, i), cd2(:, i), cd3(:, i))
        END DO

        IF(ALLOCATED(r_ca)) DEALLOCATE(r_ca, STAT=status)
        IF(ALLOCATED(r_cd)) DEALLOCATE(r_cd)
    END SUBROUTINE decomWave2d

    ! Wavelet reconstruction
    ! One-dimension, corresponding to decomWave
    SUBROUTINE reconWave(ca, cd, re_x)
        IMPLICIT NONE
        ! Declare variables
        REAL(r_kind), DIMENSION(:), INTENT(IN) :: ca, cd
        REAL(r_kind), DIMENSION(:), INTENT(OUT) :: re_x
        ! Coefficient of wavelet reconstruction
        REAL(r_kind), PARAMETER, DIMENSION(4) :: re_ll = (/0.48296291314469025, 0.836516303737469, &
                                                   0.22414386804185735, -0.12940952255092145/)
        REAL(r_kind), PARAMETER, DIMENSION(4) :: re_lh = (/-0.12940952255092145, -0.22414386804185735, &
                                                   0.836516303737469, -0.48296291314469025/)
        INTEGER(i_kind) :: i, idx

        ! The sizes of ca and cd should be the same
        IF (SIZE(ca) .NE. SIZE(cd)) THEN
            WRITE(*, *) "The sizes of ca and cd shoud be the same!"
            RETURN
        END IF

        idx = 0
        DO i = 1, SIZE(ca)-1
            idx = idx + 1
            re_x(idx) = re_ll(3)*ca(i) + re_ll(1)*ca(i+1) + re_lh(3)*cd(i) + re_lh(1) * cd(i+1)
            idx = idx + 1
            re_x(idx) = re_ll(4)*ca(i) + re_ll(2)*ca(i+1) + re_lh(4)*cd(i) + re_lh(2) * cd(i+1)
        END DO
    END SUBROUTINE reconWave

    ! Two-dimension
    ! Based on reconWave, corresponding to decomWave2d
    SUBROUTINE reconWave2d(ca, cd1, cd2, cd3, re_x)
        IMPLICIT NONE
        ! declare variables
        REAL(r_kind), DIMENSION(:, :), INTENT(IN) :: ca, cd1, cd2, cd3
        REAL(r_kind), DIMENSION(:, :), INTENT(OUT) :: re_x
        INTEGER(i_kind) :: i
        REAL(r_kind), DIMENSION(UBOUND(re_x, 2)) :: one_row
        REAL(r_kind), DIMENSION(UBOUND(re_x, 1)) :: one_column
        REAL(r_kind), DIMENSION(UBOUND(re_x, 1), UBOUND(ca, 2)) :: low_arr, high_arr

        DO i = 1, UBOUND(ca, 2)
            CALL reconWave(ca(:, i), cd1(:, i), one_column)
            low_arr(:, i) = one_column
            CALL reconWave(cd2(:, i), cd3(:, i), one_column)
            high_arr(:, i) = one_column
        END DO

        DO i = 1, UBOUND(re_x, 1)
            CALL reconWave(low_arr(i, :), high_arr(i, :), one_row)
            re_x(i, :) = one_row
        END DO
    END SUBROUTINE reconWave2d

    SUBROUTINE reconWave2d_two(coeff, cur_level, ca, zero_flag, re_x)
        IMPLICIT NONE
        TYPE(wavelet_t), INTENT(IN) :: coeff
        INTEGER(i_kind), INTENT(IN) :: cur_level
        REAL(r_kind), DIMENSION(:, :), INTENT(IN) :: ca
        LOGICAL, INTENT(IN) :: zero_flag
        REAL(r_kind), DIMENSION(:, :), ALLOCATABLE, INTENT(OUT) :: re_x
        ! Local variables
        INTEGER(i_kind) :: c_a, c_b, a, b, status
        REAL(r_kind), DIMENSION(:, :), ALLOCATABLE :: cd1, cd2, cd3

        IF (cur_level .EQ. 0) THEN
            c_a = UBOUND(coeff%ca, 1)
            c_b = UBOUND(coeff%ca, 2)
        ELSE
            c_a = UBOUND(coeff%cds(cur_level)%cd1, 1)
            c_b = UBOUND(coeff%cds(cur_level)%cd1, 2)
        END IF

        ALLOCATE(cd1(c_a, c_b), cd2(c_a, c_b), cd3(c_a, c_b), STAT=status)
        IF (zero_flag .OR. cur_level .EQ. 0) THEN
            cd1 = 0.0D0
            cd2 = 0.0D0
            cd3 = 0.0D0
        ELSE
            cd1 = coeff%cds(cur_level)%cd1
            cd2 = coeff%cds(cur_level)%cd2
            cd3 = coeff%cds(cur_level)%cd3
        END IF

        a = (c_a - 1)*2
        b = (c_b - 1)*2
        ALLOCATE(re_x(a, b), STAT=status)

        CALL reconWave2d(ca(1:c_a, 1:c_b), cd1, cd2, cd3, re_x)

        IF(ALLOCATED(cd1)) DEALLOCATE(cd1, STAT=status)
        IF(ALLOCATED(cd2)) DEALLOCATE(cd2, STAT=status)
        IF(ALLOCATED(cd3)) DEALLOCATE(cd3, STAT=status)
    END SUBROUTINE reconWave2d_two
END MODULE wavelet_m