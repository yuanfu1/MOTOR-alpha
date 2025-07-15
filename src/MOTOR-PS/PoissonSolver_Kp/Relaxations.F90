!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA.utility.relaxation
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Yuanfu Xie, 2019
! Modified by Zilong Qin (zilong.qin@gmail.com), 2020/12/23, @GBA-MWF, Shenzhen
! Modified by Zilong Qin (zilong.qin@gmail.com), 2024/12/03, @GBA-MWF, Shenzhen, add adjoint for each iteration
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module contains the data type for model_states
!! method.
!! @author Yuanfu Xie, Zilong Qin
!! @copyright (C) 2020 GBA-MWF, All rights reserved.
! @note
! @warning
! @attention

SUBROUTINE Jacobi(numOfRows, numOfCells, numVLevel, maxRowWid, numRowEle, indexColm, &
                  rowMatrix, relaxPara, rightHand, tempState, solutions, iterType)

  !>@brief
  !!=================================================================
  !!  This is a low level routine for a Jacobi relaxation iterative
  !!  method for solving a line system derived from a Poisson eqn on
  !!  a given grid, either icosahedron or latlon or rectangle. The
  !!  Laplace operator is implemented following finite volume scheme
  !!
  !!  int_A nabla^2 f = int_partial A (df)/(dn) dl
  !!
  !!  This implementation follows the Relaxed Jacobi (RJ) iteration
  !!  method with scheduled relaxation Jacobi factor proposed by
  !!  Yang and Mittal (2014)
  !!    iterType:   1: Forward, 2: Adjoint
  !!  References:
  !!    Yang, X., & Mittal, R. (2017). Efficient relaxed-Jacobi
  !!      smoothers for multigrid on parallel computers. Journal of
  !!      Computational Physics, 332, 135–142.
  !!      doi:10.1016/j.jcp.2016.12.010
  !!
  !!  uthor Yuanfu Xie
  !!=================================================================
  !
  USE kinds_m, ONLY: i_kind, r_kind, r_double
  USE parameters_m, ONLY: machineEps

  IMPLICIT NONE
  !!  INPUT:
  INTEGER(i_kind), INTENT(IN) :: numOfRows, &                          !< number of rows
                                 numOfCells, &                         !< number of cells indexed by each row (include halos in MPI)
                                 numVLevel, &                          !< number of vertical levels
                                 maxRowWid, &                          !< maximum width of the banded matrix
                                 numRowEle(numOfRows), &               !< numbers of row elements
                                 indexColm(maxRowWid, numOfRows)       !< colume index at a given row

  REAL(r_kind), INTENT(IN) :: rowMatrix(maxRowWid, numOfCells), &       !< matrix stored by row, diagonal element at first
                              relaxPara, &                             !< the relaxation parameters of each relaxation
                              righthand(numVLevel, numOfCells), &      !< values of right hand side of Poission
                              tempstate(numVLevel, numOfCells)         !< initial guesses

  !!  OUTPUT:
  REAL(r_kind), INTENT(OUT) :: solutions(numVLevel, numOfCells)        !< output solution after Jacobi iterations

  ! Local variables:
  INTEGER(i_kind) :: irow, iele, &
                     iterType                 !< 1: Forward, 2: Adjoint

  ! Check diagonal elements: quit if one is too small
  IF (MINVAL(ABS(rowMatrix(1, 1:numOfRows))) .LE. machineEps) THEN
    WRITE (*, 1) MINVAL(ABS(rowMatrix(1, 1:numOfRows)))
1   FORMAT('Relaxations->Jacobi: some diagonal elements are too small--', F16.8)
  END IF

  ! Starting iteration:
  ! IF (ioOptions .GE. 1) PRINT*,'Jacobi iteration starts...'

  ! Off diagonal elements: Jacobi iteration
  solutions = 0

  IF (iterType == 1) THEN ! Forward
    DO irow = 1, numOfRows
      solutions(:, irow) = righthand(:, irow) / rowMatrix(1, irow)
      DO iele = 2, numRowEle(irow)
        solutions(:, irow) = solutions(:, irow) - &
                             rowMatrix(iele, irow) * tempState(:, indexColm(iele, irow)) / rowMatrix(1, irow)
      END DO
    END DO
  ELSE IF (iterType == 2) THEN ! Adjoint
    DO irow = 1, numOfRows
      solutions(:, irow) = righthand(:, irow) / rowMatrix(1, irow)
    END DO

    DO irow = 1, numOfRows
      DO iele = 2, numRowEle(irow)
        solutions(:, indexColm(iele, irow)) = solutions(:, indexColm(iele, irow)) - &
                                              rowMatrix(iele, irow) * tempState(:, irow) / rowMatrix(1, indexColm(iele, irow))
      END DO
    END DO
  END IF

END SUBROUTINE Jacobi

SUBROUTINE Relaxation(solutions, tempState, relaxPara, numOfCells, numOfRows, numVLevel)

  !>@brief
    !!=================================================================
    !!  This is a low level routine for a Jacobi relaxation iterative
    !!  method for solving a line system derived from a Poisson eqn on
    !!  a given grid, either icosahedron or latlon or rectangle. The
    !!  Laplace operator is implemented following finite volume scheme
    !!
    !!  int_A nabla^2 f = int_partial A (df)/(dn) dl
    !!
    !!  This implementation follows the Relaxed Jacobi (RJ) iteration
    !!  method with scheduled relaxation Jacobi factor proposed by
    !!  Yang and Mittal (2014)
    !!
    !!  References:
    !!    Yang, X., & Mittal, R. (2017). Efficient relaxed-Jacobi
    !!      smoothers for multigrid on parallel computers. Journal of
    !!      Computational Physics, 332, 135–142.
    !!      doi:10.1016/j.jcp.2016.12.010
    !!
    !!  uthor Yuanfu Xie
    !!=================================================================
  !
  USE kinds_m, ONLY: i_kind, r_kind, r_double
  USE parameters_m, ONLY: machineEps

  IMPLICIT NONE
    !!  INPUT:
  REAL(r_kind), INTENT(IN) :: relaxPara, &
                              tempState(numVLevel, numOfCells)
  INTEGER(i_kind), INTENT(IN) :: numOfCells, numOfRows, numVLevel

    !!  OUTPUT:
  REAL(r_kind), INTENT(INOUT) :: solutions(numVLevel, numOfCells)        !< output solution after Jacobi iterations

  ! Local variables:
  INTEGER(i_kind) :: irow, iele

  ! Starting iteration:
  ! IF (ioOptions .GE. 1) PRINT*,'Jacobi iteration starts...'
  !Finish Jacobi iteration - divide the diagonal element; + Relaxation:

  DO irow = 1, numOfRows
    solutions(:, irow) = (1.0D0 - relaxPara) * tempState(:, irow) + relaxPara * solutions(:, irow)
  END DO
END SUBROUTINE

SUBROUTINE GaussSeidel(numOfRows, numVLevel, maxRowWid, numRowEle, indexColm, &
                       rowMatrix, relaxPara, rightHand, tempState, solutions, iterType)

  !>
  !!=================================================================
  !!  This is a low level routine for a Gauss-Seidel relaxation
  !!  method for solving a line system derived from a Poisson eqn on
  !!  a given grid, either icosahedron or latlon or rectangle. The
  !!  Laplace operator is implemented following finite volume scheme
  !!
  !!  int_A nabla^2 f = int_partial A (df)/(dn) dl
  !!
  !!  This implementation follows the Relaxed Jacobi (RJ) iteration
  !!  method with scheduled relaxation Jacobi factor proposed by
  !!  Yang and Mittal (2014)
  !!
  !!  INPUT:
  !!    numOfRows:  number of rows
  !!    numVLevel:  number of vertical levels
  !!    maxRowWid:  maximum width of the banded matrix
  !!    numRowEle:  numbers of row elements
  !!    indexColm:  colume index at a given row
  !!    ioOptions:  output options: 0 no text write out
  !!    rowMatrix:  matrix stored by row, diagonal element at first
  !!    numRelaxs:  number of relaxations
  !!    relaxPara:  the relaxation parameters of each relaxation
  !!    rightHand:  values of right hand side of Poission
  !!    tempState:  initial guesses
  !!    iterType:   1: Forward, 2: Adjoint
  !!
  !!  OUTPUT:
  !!    solutions:  output solution after Gauss-Seidel iterations
  !!
  !!  References:
  !!    Yang, X., & Mittal, R. (2017). Efficient relaxed-Jacobi
  !!      smoothers for multigrid on parallel computers. Journal of
  !!      Computational Physics, 332, 135–142.
  !!      doi:10.1016/j.jcp.2016.12.010
  !!
  !!  uthor Yuanfu Xie
  !!=================================================================
  !
  USE kinds_m, ONLY: i_kind, r_kind, r_double
  USE parameters_m, ONLY: machineEps

  IMPLICIT NONE

  INTEGER(i_kind), INTENT(IN) :: numOfRows, numVLevel, maxRowWid, &
                                 numRowEle(numOfRows), &
                                 indexColm(maxRowWid, numOfRows), iterType
  REAL(r_kind), INTENT(IN) :: rowMatrix(maxRowWid, numOfRows), &
                              relaxPara, &
                              righthand(numVLevel, numOfRows), &
                              tempState(numVLevel, numOfRows)
  REAL(r_kind), INTENT(OUT) :: solutions(numVLevel, numOfRows)

  ! Local variables:
  INTEGER(i_kind) :: ir, ie
  REAL(r_kind), ALLOCATABLE :: row(:)

  ! Debugging:
  INTEGER :: mxidx
  REAL(r_kind) :: amx

  ! Allocate temporary variables:
  ALLOCATE (row(numVLevel))

  ! Check diagonal elements: quit if one is too small
  IF (MINVAL(ABS(rowMatrix(1, :))) .LE. machineEps) THEN
    WRITE (*, 1) MINVAL(ABS(rowMatrix(1, :)))
1   FORMAT('Relaxations->Gauss-Seidel: some diagonal elements are too small--', F16.8)
  END IF

  ! Starting iterations:
  ! PRINT *, 'Gauss-Seidel iteration starts...'
  solutions = 0.0D0
  IF (iterType == 1) THEN ! Forward
    DO ir = 1, numOfRows
      solutions(:, ir) = solutions(:, ir) + righthand(:, ir)

      DO ie = 2, numRowEle(ir)
        solutions(:, ir) = solutions(:, ir) - rowMatrix(ie, ir) * &
                           tempState(:, indexColm(ie, ir))
      END DO

    END DO
  ELSEIF (iterType == 2) THEN ! Adjoint
    DO ir = 1, numOfRows
      solutions(:, ir) = solutions(:, ir) + righthand(:, ir)

      DO ie = 2, numRowEle(ir)
        solutions(:, indexColm(ie, ir)) = solutions(:, indexColm(ie, ir)) - rowMatrix(ie, ir) * &
                                          tempState(:, ir)
      END DO
    END DO
  END IF

  DO ir = 1, numOfRows
    solutions(:, ir) = solutions(:, ir) / rowMatrix(1, ir)
  END DO

END SUBROUTINE GaussSeidel

! SUBROUTINE Mixed(numOfRows, numOfCells, numVLevel, maxRowWid, numRowEle, indexColm, &
!                  rowMatrix, relaxPara, rightHand, tempState, solutions)

!   !>@brief
!   !!=================================================================
!   !!  This is a Mixed solver that use the GS on the inner points, and Jacobi on the boundaries.
!   !!=================================================================
!   !
!   USE kinds_m, ONLY: i_kind, r_kind, r_double
!   USE parameters_m, ONLY: machineEps

!   IMPLICIT NONE
!   !!  INPUT:
!   INTEGER(i_kind), INTENT(IN) :: numOfRows, &                          !< number of rows
!                                  numOfCells, &                         !< number of cells indexed by each row (include halos in MPI)
!                                  numVLevel, &                          !< number of vertical levels
!                                  maxRowWid, &                          !< maximum width of the banded matrix
!                                  numRowEle(numOfRows), &               !< numbers of row elements
!                                  indexColm(maxRowWid, numOfRows)       !< colume index at a given row

!   REAL(r_kind), INTENT(IN) :: rowMatrix(maxRowWid, numOfRows), &       !< matrix stored by row, diagonal element at first
!                               relaxPara, &                             !< the relaxation parameters of each relaxation
!                               righthand(numVLevel, numOfCells), &      !< values of right hand side of Poission
!                               tempstate(numVLevel, numOfCells)         !< initial guesses

!   !!  OUTPUT:
!   REAL(r_kind), INTENT(OUT) :: solutions(numVLevel, numOfCells)        !< output solution after Jacobi iterations

!   ! Local variables:
!   INTEGER(i_kind) :: irow, iele, it, ir, ie
!   REAL(r_kind), ALLOCATABLE :: row(:)

!   ! Check diagonal elements: quit if one is too small
!   ! IF (MINVAL(ABS(rowMatrix(1,:))) .LE. machineEps) THEN
!   !    WRITE(*,1) MINVAL(ABS(rowMatrix(1,:)))
!   !1   FORMAT('Relaxations->Jacobi: some diagonal elements are too small--',F16.8)
!   ! END IF

!   ! Starting iteration:
!   ! IF (ioOptions .GE. 1) PRINT*,'Jacobi iteration starts...'
!   ! Allocate temporary variables:
!   ALLOCATE (row(numVLevel))

!   solutions = tempState
!   DO ir = 1, numOfRows
!     row = righthand(:, ir)
!     DO ie = 2, numRowEle(ir)
!       row = row - rowMatrix(ie, ir) * &
!             solutions(:, indexColm(ie, ir))
!     END DO
!     solutions(:, ir) = &
!       (1.0D0 - 1) * solutions(:, ir) + &
!       1 * row / rowMatrix(1, ir)
!   END DO

! END SUBROUTINE Mixed
