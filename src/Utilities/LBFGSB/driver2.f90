!
!  L-BFGS-B is released under the “New BSD License” (aka “Modified BSD License”
!  or “3-clause license”)
!  Please read attached file License.txt
!
!
!                  DRIVER 2  in Fortran 90
!    --------------------------------------------------------------
!             CUSTOMIZED DRIVER FOR L-BFGS-B
!    --------------------------------------------------------------
!
!       L-BFGS-B is a code for solving large nonlinear optimization
!            problems with simple bounds on the variables.
!
!       The code can also be used for unconstrained problems and is
!       as efficient for these problems as the earlier limited memory
!                         code L-BFGS.
!
!       This driver illustrates how to control the termination of the
!       run and how to design customized output.
!
!    References:
!
!       [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
!       memory algorithm for bound constrained optimization'',
!       SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.
!
!       [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: FORTRAN
!       Subroutines for Large Scale Bound Constrained Optimization''
!       Tech. Report, NAM-11, EECS Department, Northwestern University,
!       1994.
!
!
!         (Postscript files of these papers are available via anonymous
!          ftp to eecs.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.)
!
!                             *  *  *
!
!         February 2011   (latest revision)
!         Optimization Center at Northwestern University
!         Instituto Tecnologico Autonomo de Mexico
!
!         Jorge Nocedal and Jose Luis Morales
!
!    **************
PROGRAM driver

!     This driver shows how to replace the default stopping test
!     by other termination criteria. It also illustrates how to
!     print the values of several parameters during the course of
!     the iteration. The sample problem used here is the same as in
!     DRIVER1 (the extended Rosenbrock function with bounds on the
!     variables).

  IMPLICIT NONE

!     Declare variables and parameters needed by the code.
!
!     Note that we suppress the default output (iprint = -1)
!     We suppress both code-supplied stopping tests because the
!     user is providing his/her own stopping criteria.

  INTEGER, PARAMETER    :: n = 25, m = 5, iprint = -1
  INTEGER, PARAMETER    :: dp = KIND(1.0D0)
  REAL(dp), PARAMETER    :: factr = 0.0D0, pgtol = 0.0D0

  CHARACTER(len=60)      :: task, csave
  LOGICAL                :: lsave(4)
  INTEGER                :: isave(44)
  REAL(dp)               :: f
  REAL(dp)               :: dsave(29)
  INTEGER, ALLOCATABLE  :: nbd(:), iwa(:)
  REAL(dp), ALLOCATABLE  :: x(:), l(:), u(:), g(:), wa(:)
!
  REAL(dp)               :: t1, t2
  INTEGER                :: i

  ALLOCATE (nbd(n), x(n), l(n), u(n), g(n))
  ALLOCATE (iwa(3 * n))
  ALLOCATE (wa(2 * m * n + 5 * n + 11 * m * m + 8 * m))
!
!    This driver shows how to replace the default stopping test
!    by other termination criteria. It also illustrates how to
!    print the values of several parameters during the course of
!    the iteration. The sample problem used here is the same as in
!    DRIVER1 (the extended Rosenbrock function with bounds on the
!    variables).
!    We now specify nbd which defines the bounds on the variables:
!                   l   specifies the lower bounds,
!                   u   specifies the upper bounds.

!    First set bounds on the odd numbered variables.

  DO 10 i = 1, n, 2
    nbd(i) = 2
    l(i) = 1.0D0
    u(i) = 1.0D2
10  CONTINUE

!    Next set bounds on the even numbered variables.

    DO 12 i = 2, n, 2
      nbd(i) = 2
      l(i) = -1.0D2
      u(i) = 1.0D2
12    CONTINUE

!     We now define the starting point.

      DO 14 i = 1, n
        x(i) = 3.0D0
14      CONTINUE

!     We now write the heading of the output.

        WRITE (6, 16)
16      FORMAT(/, 5X, 'Solving sample problem.', &
                /, 5X, ' (f = 0.0 at the optimal solution.)',/)

!     We start the iteration by initializing task.
!
        task = 'START'

!        ------- the beginning of the loop ----------

        DO WHILE (task(1:2) .EQ. 'FG' .OR. task .EQ. 'NEW_X' .OR. &
                  task .EQ. 'START')

!     This is the call to the L-BFGS-B code.

          CALL setulb(n, m, x, l, u, nbd, f, g, factr, pgtol, wa, iwa, task, iprint, &
                      csave, lsave, isave, dsave)

          IF (task(1:2) .EQ. 'FG') THEN

!        the minimization routine has returned to request the
!        function f and gradient g values at the current x.

!        Compute function value f for the sample problem.

            f = .25D0 * (x(1) - 1.D0)**2
            DO 20 i = 2, n
              f = f + (x(i) - x(i - 1)**2)**2
20            CONTINUE
              f = 4.D0 * f

!        Compute gradient g for the sample problem.

              t1 = x(2) - x(1)**2
              g(1) = 2.D0 * (x(1) - 1.D0) - 1.6D1 * x(1) * t1
              DO 22 i = 2, n - 1
                t2 = t1
                t1 = x(i + 1) - x(i)**2
                g(i) = 8.D0 * t2 - 1.6D1 * x(i) * t1
22              CONTINUE
                g(n) = 8.D0 * t1
!
                ELSE
!
                IF (task(1:5) .EQ. 'NEW_X') THEN
!
!       the minimization routine has returned with a new iterate.
!       At this point have the opportunity of stopping the iteration
!       or observing the values of certain parameters
!
!       First are two examples of stopping tests.

!       Note: task(1:4) must be assigned the value 'STOP' to terminate
!         the iteration and ensure that the final results are
!         printed in the default format. The rest of the character
!         string TASK may be used to store other information.

!       1) Terminate if the total number of f and g evaluations
!            exceeds 99.

                  IF (isave(34) .GE. 99) &
                    task = 'STOP: TOTAL NO. of f AND g EVALUATIONS EXCEEDS LIMIT'

!       2) Terminate if  |proj g|/(1+|f|) < 1.0d-10, where
!          "proj g" denoted the projected gradient

                  IF (dsave(13) .LE. 1.D-10 * (1.0D0 + ABS(f))) &
                    task = 'STOP: THE PROJECTED GRADIENT IS SUFFICIENTLY SMALL'

!       We now wish to print the following information at each
!       iteration:
!
!         1) the current iteration number, isave(30),
!         2) the total number of f and g evaluations, isave(34),
!         3) the value of the objective function f,
!         4) the norm of the projected gradient,  dsve(13)
!
!       See the comments at the end of driver1 for a description
!       of the variables isave and dsave.

                  WRITE (6, '(2(a,i5,4x),a,1p,d12.5,4x,a,1p,d12.5)') 'Iterate' &
                    , isave(30), 'nfg =', isave(34), 'f =', f, '|proj g| =', dsave(13)

!        If the run is to be terminated, we print also the information
!        contained in task as well as the final value of x.

                  IF (task(1:4) .EQ. 'STOP') THEN
                    WRITE (6, *) task
                    WRITE (6, *) 'Final X='
                    WRITE (6, '((1x,1p, 6(1x,d11.4)))') (x(i), i=1, n)
                  END IF

                END IF
                END IF

              END DO
!           ---------- the end of the loop -------------

!     If task is neither FG nor NEW_X we terminate execution.

              END PROGRAM driver

!======================= The end of driver2 ============================

