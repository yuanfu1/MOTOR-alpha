!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA.possionSolver_Kp.psMatrix
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Yuanfu Xie, Zilong Qin
! VERSION           : V 0.0
! HISTORY           :
!   Created by Yuanfu Xie, in grids/gzGrid.F90
!     Created by Yuanfu Xie, 2019/04, @GBA-MWF, Shenzhen
!     Modified by Yuanfu Xie 2020-4 for adding boundary conditions
!   Reforged by Zilong Qin (zilong.qin@gmail.com), 2020/12/17, @GBA-MWF, Shenzhen, add parallelization
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module contains the data type for model_states
!! method.
!! @copyright (C) 2020 GBA-MWF, All rights reserved.
! @note
! @warning
! @attention
MODULE psMatrix_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE singleGrid_m, ONLY: singleGrid_t

  TYPE psMatrix_t
    ! New data structure of the coefficient matrix:
    INTEGER(i_kind)              :: num_eqns = 0 &
                                    , max_band = 0
    INTEGER(i_kind), ALLOCATABLE :: num_elem(:), idx_cols(:, :)
    REAL(r_kind), ALLOCATABLE    :: row_elem(:, :)

  CONTAINS
    FINAL :: destructor

    PROCEDURE, PUBLIC :: loadMatrix
    PROCEDURE, PUBLIC :: CalResidus

  END TYPE psMatrix_t
  INTERFACE psMatrix_t
    PROCEDURE :: constructor
  END INTERFACE psMatrix_t
CONTAINS

  SUBROUTINE loadMatrix(this, sg)
    CLASS(psMatrix_t) :: this
    TYPE(singleGrid_t), INTENT(IN) :: sg

    ! Local variables:
    INTEGER(i_kind) :: i, j, istcl, ir

    ! Debug variables:
    INTEGER(i_kind) :: idebug

    ! Initialize: New data structure of the coefficient matrix
    this%num_eqns = sg%num_icell

    ! For all cells:
    DO i = 1, this%num_eqns

      ! New data structure of the coefficient matrix:
      this%idx_cols(:, i) = 0
      this%row_elem(:, i) = 0.0D0

      ! Assume the diagonal element always exist and saved in the first slot:
      this%num_elem(i) = 1
      this%idx_cols(1, i) = i

      IF (sg%cell_type(i) .NE. 0) THEN ! The boundary cell
        this%row_elem(1, i) = 1.0D0
      ELSE              ! inner cell
        ! For all edge surrounding a cell:
        DO j = 1, sg%num_edge
          ! For each edge, go through all stencil:
          DO istcl = 1, sg%num_stcl
            ! New data structure of the coefficient matrix ordered by rows:
            DO ir = 1, this%num_elem(i)
              IF (sg%edge_stcl(istcl, j, i) .EQ. this%idx_cols(ir, i)) THEN
                this%row_elem(ir, i) = this%row_elem(ir, i) + &
                                       sg%coef_norm(istcl, 1, j, i) * sg%edge_leng(j, i) / &
                                       sg%cell_area(i)
                EXIT
              END IF
            END DO
            ! a new element:
            IF (ir .GT. this%num_elem(i)) THEN
              IF (sg%coef_norm(istcl, 1, j, i) .NE. 0.0D0) THEN
                this%num_elem(i) = this%num_elem(i) + 1
                this%row_elem(ir, i) = this%row_elem(ir, i) + &
                                       sg%coef_norm(istcl, 1, j, i) * sg%edge_leng(j, i) / &
                                       sg%cell_area(i)
                this%idx_cols(ir, i) = sg%edge_stcl(istcl, j, i)
                IF (this%idx_cols(ir, i) .LE. 0) THEN
                  PRINT *, 'LoadMatrix: Negative or zero index: ', sg%edge_stcl(istcl, j, i)
                  STOP
                END IF
              END IF
            END IF
          END DO
        END DO
      END IF

    END DO
    CALL sg%ExchangeMatOnHalo2D(this%max_band, this%row_elem)

    ! ! Make the transpose for the PSMatrix structure:
    ! ALLOCATE (this%num_elem_adj(this%num_eqns))
    ! ALLOCATE (this%idx_cols_adj(max_band, this%num_eqns))
    ! ALLOCATE (this%row_elem_adj(max_band, this%num_eqns))
    ! this%num_elem_adj = 1
    ! this%idx_cols_adj = 0
    ! this%row_elem_adj = 0.0D0

    ! DO i = 1, this%num_eqns
    !   this%idx_cols_adj(1, i) = i
    !   this%row_elem_adj(1, i) = this%row_elem(1, i)
    !   DO j = 2, this%num_elem(i)
    !     this%num_elem_adj(this%idx_cols(j, i)) = this%num_elem_adj(this%idx_cols(j, i)) + 1
    !     this%idx_cols_adj(this%num_elem_adj(this%idx_cols(j, i)), this%idx_cols(j, i)) = i
    !     this%row_elem_adj(this%num_elem_adj(this%idx_cols(j, i)), this%idx_cols(j, i)) = this%row_elem(j, i)
    !   END DO
    ! END DO

    ! ! PRINT*, 'Matrix:'
    ! ! BLOCK
    ! !   REAL(r_kind) :: rowValue(this%num_eqns)

    ! !   DO i = 1, this%num_eqns
    ! !     rowValue = 0.0
    ! !     DO j = 1, this%num_elem(i)
    ! !       rowValue((this%idx_cols(j, i))) = this%row_elem(j, i)
    ! !     ENDDO

    ! !     write(*,"(16E15.7)") rowValue

    ! !   END DO
    ! ! END BLOCK

    ! ! PRINT*, 'Transpose matrix:'
    ! ! BLOCK
    ! !   REAL(r_kind) :: rowValue(this%num_eqns)

    ! !   DO i = 1, this%num_eqns
    ! !     rowValue = 0.0
    ! !     DO j = 1, this%num_elem_adj(i)
    ! !       rowValue((this%idx_cols_adj(j, i))) = this%row_elem_adj(j, i)
    ! !     ENDDO

    ! !     write(*,"(16E15.7)") rowValue

    ! !   END DO
    ! ! END BLOCK

    ! ! STOP

    ! this%num_elem = 0
    ! this%idx_cols = 0
    ! this%row_elem = 0.0D0

    ! this%num_elem = this%num_elem_adj
    ! this%idx_cols = this%idx_cols_adj
    ! this%row_elem = this%row_elem_adj
    IF (ALLOCATED(this%num_elem)) PRINT *, 'this%num_elem max, gLevel', MAXVAL(this%num_elem), sg%gLevel, this%num_eqns
  END SUBROUTINE loadMatrix

!> @brief
!! This routine loads the coefficient matrix of Poisson equation
!! History
!!   Created by Yuanfu Xie May 2019
!! Modified by Zilong Qin (zilong.qin@gmail.com), 2020/12/18, @GBA-MWF, Shenzhen! @see
!!@Author Yuanfu Xie
! @note
! @warning
! @attention
  SUBROUTINE CalResidus(this, vlevel, rhs, sol, rsd, sg, oprType)
    CLASS(psMatrix_t) :: this
    TYPE(singleGrid_t), INTENT(IN) :: sg

    INTEGER(i_kind), INTENT(IN) :: vlevel ! Number of vertical levels
    REAL(r_kind), INTENT(IN) :: rhs(vlevel, sg%num_cell), &
                                sol(vlevel, sg%num_cell)
    REAL(r_kind), INTENT(OUT) :: rsd(vlevel, sg%num_cell)
    INTEGER(i_kind) :: i, j
    INTEGER(i_kind) :: oprType

    ! ! debugging:
    ! INTEGER(i_kind) :: imx
    ! REAL(r_kind) :: emx

    ! All elements:

    rsd = rhs

    IF(oprType == 1) THEN
      DO i = 1, sg%num_icell
        DO j = 1, this%num_elem(i) ! Off diagonal: starts from 2
          rsd(:, i) = rsd(:, i) - this%row_elem(j, i) * sol(:, this%idx_cols(j, i))
        END DO
      END DO
    END IF

    IF(oprType == 2) THEN
      rsd(:, sg%num_icell+1:sg%num_cell) = 0.0D0
      DO i = 1, sg%num_icell
        DO j = 1, this%num_elem(i) ! Off diagonal: starts from 2
          rsd(:, this%idx_cols(j, i)) = rsd(:, this%idx_cols(j, i)) - this%row_elem(j, i) * sol(:, i)
        END DO
      END DO
      CALL sg%ExchangeMatOnHaloReverseSum2D(vLevel, rsd)
    END IF

    ! ! Debugging:
    ! imx = 0
    ! emx = 0.0D0
    ! DO i = 1, sg%num_icell
    !   if (emx .LT. ABS(rsd(1, i))) then
    !     imx = i
    !     emx = ABS(rsd(1, i))
    !   end if
    ! END DO

    CALL sg%ExchangeMatOnHalo2D(vLevel, rsd)

    ! comm rsd
  END SUBROUTINE CalResidus

  FUNCTION constructor(max_band, sg) RESULT(this)
    IMPLICIT NONE
    TYPE(psMatrix_t) :: this
    TYPE(singleGrid_t), INTENT(IN) :: sg
    INTEGER(i_kind), INTENT(IN) :: max_band

    this%num_eqns = sg%num_icell
    this%max_band = max_band

    ALLOCATE (this%num_elem(this%num_eqns))
    ALLOCATE (this%idx_cols(max_band, this%num_eqns))
    ALLOCATE (this%row_elem(max_band, sg%num_cell))

    ! ALLOCATE (this%rHandSide(vLevel, sg%num_cell), &
    !           this%solutions(vLevel, sg%num_cell))

    CALL this%loadMatrix(sg)

  END FUNCTION constructor

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    TYPE(psMatrix_t), INTENT(INOUT) :: this

    IF (ALLOCATED(this%num_elem)) DEALLOCATE (this%num_elem)
    IF (ALLOCATED(this%idx_cols)) DEALLOCATE (this%idx_cols)
    IF (ALLOCATED(this%row_elem)) DEALLOCATE (this%row_elem)
    ! IF (ALLOCATED(this%rHandSide)) DEALLOCATE (this%rHandSide)
    ! IF (ALLOCATED(this%solutions)) DEALLOCATE (this%solutions)
  END SUBROUTINE destructor

END MODULE psMatrix_m
