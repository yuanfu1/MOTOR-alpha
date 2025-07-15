!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2021/1/26, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!! @copyright (C) 2020 GBA-MWF, All rights reserved.
MODULE SolverLBFGS_m
  USE State_m, ONLY: State_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE JFunc_m, ONLY: JFunc_t
  USE kinds_m, ONLY: i_kind, r_kind
  USE MiniSolver_m, ONLY: MiniSolver_t

  TYPE, EXTENDS(MiniSolver_t) :: SolverLBFGS_t

  CONTAINS
    PROCEDURE, PUBLIC :: initialize
    FINAL :: destructor
    PROCEDURE, PUBLIC, PASS :: run
    PROCEDURE, PUBLIC, NOPASS :: aggr_states_to_base_proc
    PROCEDURE, PUBLIC, NOPASS :: scatter_states_to_each_proc
  END TYPE SolverLBFGS_t

CONTAINS

  SUBROUTINE initialize(this, configFile)
    IMPLICIT NONE
    CLASS(SolverLBFGS_t) :: this
    CHARACTER(LEN=1024), INTENT(IN) :: configFile

    CALL this%MiniSolver_t%initialize(configFile)
  END SUBROUTINE

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(SolverLBFGS_t), INTENT(INOUT) :: this

  END SUBROUTINE destructor

  SUBROUTINE aggr_states_to_base_proc(X, x1D, sg)
    TYPE(State_t), INTENT(INOUT) :: X
    REAL(r_kind), ALLOCATABLE, INTENT(INOUT)  :: x1D(:)
    TYPE(SingleGrid_t), INTENT(IN) :: sg

    REAL(r_kind), ALLOCATABLE :: dataSwap(:, :, :)
    INTEGER(i_kind) :: idxStart, idxEnd, nGlobSizePerVar

    idxStart = 1

    IF (sg%mpddInfo_sg%isBaseProc()) THEN
      ALLOCATE (dataSwap(sg%vLevel, X%sg%num_icell_global, sg%tSlots))
    END IF

    DO i = LBOUND(X%fields, 1), UBOUND(X%fields, 1)
      IF (sg%mpddInfo_sg%isBaseProc()) THEN
        nGlobSizePerVar = X%sg%tSlots * X%sg%vLevel * X%sg%num_icell_global
      END IF

      CALL X%sg%aggrGridRealForFieldGrid(X%fields(i)%DATA, dataSwap, [X%sg%vLevel, X%sg%num_icell_global, X%sg%tSlots])

      IF (sg%mpddInfo_sg%isBaseProc()) THEN
        idxEnd = idxStart + nGlobSizePerVar - 1
        x1D(idxStart:idxEnd) = RESHAPE(dataSwap, (/nGlobSizePerVar/))
        idxStart = idxEnd + 1
      END IF
    END DO

    IF (sg%mpddInfo_sg%isBaseProc()) DEALLOCATE (dataSwap)
  END SUBROUTINE aggr_states_to_base_proc

  SUBROUTINE scatter_states_to_each_proc(X, x1D, sg)
    TYPE(State_t), INTENT(INOUT) :: X
    REAL(r_kind), ALLOCATABLE, INTENT(INOUT)  :: x1D(:)
    TYPE(SingleGrid_t), INTENT(IN) :: sg

    INTEGER(i_kind) :: idxStart, idxEnd, nGlobSizePerVar
    REAL(r_kind), ALLOCATABLE :: dataSwap(:, :, :)

    IF (sg%mpddInfo_sg%isBaseProc()) idxStart = 1

    IF (sg%mpddInfo_sg%isBaseProc()) THEN
      ALLOCATE (dataSwap(X%sg%vLevel, X%sg%num_icell_global, X%sg%tSlots))
    END IF

    DO i = LBOUND(X%fields, 1), UBOUND(X%fields, 1)
      IF (sg%mpddInfo_sg%isBaseProc()) THEN
        nGlobSizePerVar = X%sg%tSlots * X%sg%vLevel * X%sg%dimCell_global(1) * X%sg%dimCell_global(2)
        idxEnd = idxStart + nGlobSizePerVar - 1
      END IF

      IF (sg%mpddInfo_sg%isBaseProc()) THEN
        dataSwap = RESHAPE(x1D(idxStart:idxEnd), (/X%sg%vLevel, X%sg%num_icell_global, X%sg%tSlots/))
      END IF
      CALL X%sg%DistGridRealWithHaloExForFieldGrid(dataSwap, X%fields(i)%DATA, [X%sg%vLevel, X%sg%num_icell_global, X%sg%tSlots])

      IF (sg%mpddInfo_sg%isBaseProc()) idxStart = idxEnd + 1
    END DO
    IF (sg%mpddInfo_sg%isBaseProc()) DEALLOCATE (dataSwap)

  END SUBROUTINE scatter_states_to_each_proc

  SUBROUTINE run(this, X, JFunc, sg, iters)
    IMPLICIT NONE
    CLASS(SolverLBFGS_t) :: this
    TYPE(State_t), INTENT(INOUT) :: X
    TYPE(JFunc_t), INTENT(INOUT) :: JFunc
    TYPE(SingleGrid_t), INTENT(INOUT) :: sg
    ! TYPE(State_t) :: Xb
    INTEGER(i_kind), OPTIONAL :: iters

    INTEGER    :: i, j, k, loopIdx
    INTEGER    :: n, m = 5, iprint = 1

    REAL(r_kind), PARAMETER    :: factr = 1.0D+7, pgtol = 1.0D-5
    CHARACTER(len=60)          :: task, csave
    LOGICAL                    :: lsave(4)
    INTEGER                    :: isave(44)
    REAL(r_kind)                 :: f
    REAL(r_kind)                 :: dsave(29)
    INTEGER, ALLOCATABLE         :: nbd(:), iwa(:)
    REAL(r_kind), ALLOCATABLE    :: x1D(:), l(:), u(:), g(:), wa(:)

    !     Declare a few additional variables for this sample problem
    IF (.NOT. sg%mpddInfo_sg%isActiveProc()) RETURN

    ! Xb = this%X
    PRINT *, 'Start LBFGS minimization...'

    n = 0
    DO i = LBOUND(X%fields, 1), UBOUND(X%fields, 1)
      n = n + sg%num_icell_global * sg%vLevel * sg%tSlots
    END DO

    PRINT *, 'n is: ', n, 'myrank is: ', sg%mpddInfo_sg%myrank, sg%num_icell_global, sg%vLevel, sg%tSlots

    !   Allocate dynamic arrays
    IF (sg%mpddInfo_sg%isBaseProc()) THEN
      ALLOCATE (nbd(2 * n), x1D(n), g(n), l(n), u(n))
      ALLOCATE (iwa(3 * n))
      ALLOCATE (wa(2 * m * n + 5 * n + 11 * m * m + 8 * m))
    END IF

    IF (sg%mpddInfo_sg%isBaseProc()) THEN
      nbd = 0
    END IF

    !   We now define the starting point.
    CALL this%aggr_states_to_base_proc(X, x1D, sg)

    IF (sg%isBaseProc()) PRINT *, 'Solving the simple variational problem by LBFGS:'

    !   We start the iteration by initializing task.
    task = 'START'
    !   The beginning of the loop
    DO WHILE (task(1:2) .EQ. 'FG' .OR. task .EQ. 'NEW_X' .OR. &
              task .EQ. 'START')

      !     This is the call to the L-BFGS-B code.
      IF (sg%mpddInfo_sg%isBaseProc()) THEN
        CALL setulb(n, m, x1D, l, u, nbd, f, g, factr, pgtol, &
                    wa, iwa, task, iprint, &
                    csave, lsave, isave, dsave)
      END IF

      ! Boardcast the calculation state to all procs.
      CALL sg%mpddInfo_sg%bCastString(task)
      CALL sg%mpddInfo_sg%bCastInt1D(isave)

      IF (task(1:2) .EQ. 'FG') THEN

        ! Gather to the base proc
        BLOCK
          TYPE(State_t) :: GG

          ! Scatter x1d to X
          CALL this%scatter_states_to_each_proc(X, x1D, sg)

          ! Compute gradient g
          CALL JFunc%get_J_and_grad_J_vec(X, f, GG)

          ! Aggregate GG to g
          CALL this%aggr_states_to_base_proc(GG, g, sg)
        END BLOCK
      END IF

      IF (iSave(30) .EQ. this%maxOptStep) EXIT
      ! IF (f < 2030.0D0) EXIT
    END DO

    CALL this%scatter_states_to_each_proc(X, x1D, sg)

    ! Save the f value at the last.
    this%J = f

    !   End of loop do while
    IF (sg%mpddInfo_sg%isBaseProc()) THEN
      DEALLOCATE (nbd, x1D, g)
      DEALLOCATE (iwa)
      DEALLOCATE (wa)
    END IF
  END SUBROUTINE run

END MODULE SolverLBFGS_m
