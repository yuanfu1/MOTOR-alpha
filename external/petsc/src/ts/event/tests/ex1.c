static char help[] = "Solves the trivial ODE du/dt = 1, u(0) = 0. \n\n";

#include <petscts.h>
#include <petscpc.h>

static PetscErrorCode RHSFunction(TS, PetscReal, Vec, Vec, void *);
static PetscErrorCode RHSJacobian(TS, PetscReal, Vec, Mat, Mat, void *);

static PetscErrorCode PreStep(TS);
static PetscErrorCode PostStep(TS);
static PetscErrorCode Monitor(TS, PetscInt, PetscReal, Vec, void *);
static PetscErrorCode Event(TS, PetscReal, Vec, PetscReal *, void *);
static PetscErrorCode PostEvent(TS, PetscInt, PetscInt[], PetscReal, Vec, PetscBool, void *);
static PetscErrorCode TransferSetUp(TS, PetscInt, PetscReal, Vec, PetscBool *, void *);
static PetscErrorCode Transfer(TS, PetscInt, Vec[], Vec[], void *);

int main(int argc, char **argv)
{
  TS              ts;
  PetscInt        n, ntransfer[] = {2, 2};
  const PetscInt  n_end = 11;
  PetscReal       t;
  const PetscReal t_end = 11;
  Vec             x;
  Vec             f;
  Mat             A;

  PetscFunctionBeginUser;
  PetscCall(PetscInitialize(&argc, &argv, NULL, help));

  PetscCall(TSCreate(PETSC_COMM_WORLD, &ts));

  PetscCall(VecCreate(PETSC_COMM_WORLD, &f));
  PetscCall(VecSetSizes(f, 1, PETSC_DECIDE));
  PetscCall(VecSetFromOptions(f));
  PetscCall(VecSetUp(f));
  PetscCall(TSSetRHSFunction(ts, f, RHSFunction, NULL));
  PetscCall(VecDestroy(&f));

  PetscCall(MatCreate(PETSC_COMM_WORLD, &A));
  PetscCall(MatSetSizes(A, 1, 1, PETSC_DECIDE, PETSC_DECIDE));
  PetscCall(MatSetFromOptions(A));
  PetscCall(MatSetUp(A));
  PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
  /* ensure that the Jacobian matrix has diagonal entries since that is required by TS */
  PetscCall(MatShift(A, (PetscReal)1));
  PetscCall(MatShift(A, (PetscReal)-1));
  PetscCall(TSSetRHSJacobian(ts, A, A, RHSJacobian, NULL));
  PetscCall(MatDestroy(&A));

  PetscCall(VecCreate(PETSC_COMM_WORLD, &x));
  PetscCall(VecSetSizes(x, 1, PETSC_DECIDE));
  PetscCall(VecSetFromOptions(x));
  PetscCall(VecSetUp(x));
  PetscCall(TSSetSolution(ts, x));
  PetscCall(VecDestroy(&x));

  PetscCall(TSMonitorSet(ts, Monitor, NULL, NULL));
  PetscCall(TSSetPreStep(ts, PreStep));
  PetscCall(TSSetPostStep(ts, PostStep));

  {
    TSAdapt adapt;
    PetscCall(TSGetAdapt(ts, &adapt));
    PetscCall(TSAdaptSetType(adapt, TSADAPTNONE));
  }
  {
    PetscInt  direction[3];
    PetscBool terminate[3];
    direction[0] = +1;
    terminate[0] = PETSC_FALSE;
    direction[1] = -1;
    terminate[1] = PETSC_FALSE;
    direction[2] = 0;
    terminate[2] = PETSC_FALSE;
    PetscCall(TSSetTimeStep(ts, 1));
    PetscCall(TSSetEventHandler(ts, 3, direction, terminate, Event, PostEvent, NULL));
  }
  PetscCall(TSSetExactFinalTime(ts, TS_EXACTFINALTIME_STEPOVER));
  PetscCall(TSSetResize(ts, PETSC_TRUE, TransferSetUp, Transfer, ntransfer));
  PetscCall(TSSetFromOptions(ts));

  /* --- First Solve --- */

  PetscCall(TSSetStepNumber(ts, 0));
  PetscCall(TSSetTimeStep(ts, 1));
  PetscCall(TSSetTime(ts, 0));
  PetscCall(TSSetMaxTime(ts, PETSC_MAX_REAL));
  PetscCall(TSSetMaxSteps(ts, 3));

  PetscCall(TSGetTime(ts, &t));
  PetscCall(TSGetSolution(ts, &x));
  PetscCall(VecSet(x, t));
  while (t < t_end) {
    PetscCall(PetscPrintf(PetscObjectComm((PetscObject)ts), "TSSolve: Begin\n"));
    PetscCall(TSSolve(ts, NULL));
    PetscCall(PetscPrintf(PetscObjectComm((PetscObject)ts), "TSSolve: End\n\n"));
    PetscCall(TSGetTime(ts, &t));
    PetscCall(TSGetStepNumber(ts, &n));
    PetscCall(TSSetMaxSteps(ts, PetscMin(n + 3, n_end)));
  }
  PetscCall(PetscPrintf(PetscObjectComm((PetscObject)ts), "TSSolve: Begin\n"));
  PetscCall(TSSolve(ts, NULL));
  PetscCall(PetscPrintf(PetscObjectComm((PetscObject)ts), "TSSolve: End\n\n"));

  /* --- Second Solve --- */

  PetscCall(TSSetStepNumber(ts, 0));
  PetscCall(TSSetTimeStep(ts, 1));
  PetscCall(TSSetTime(ts, 0));
  PetscCall(TSSetMaxTime(ts, 3));
  PetscCall(TSSetMaxSteps(ts, PETSC_INT_MAX));

  PetscCall(TSGetTime(ts, &t));
  PetscCall(TSGetSolution(ts, &x));
  PetscCall(VecSet(x, t));
  while (t < t_end) {
    PetscCall(PetscPrintf(PetscObjectComm((PetscObject)ts), "TSSolve: Begin\n"));
    PetscCall(TSSolve(ts, NULL));
    PetscCall(PetscPrintf(PetscObjectComm((PetscObject)ts), "TSSolve: End\n\n"));
    PetscCall(TSGetTime(ts, &t));
    PetscCall(TSSetMaxTime(ts, PetscMin(t + 3, t_end)));
  }
  PetscCall(PetscPrintf(PetscObjectComm((PetscObject)ts), "TSSolve: Begin\n"));
  PetscCall(TSSolve(ts, NULL));
  PetscCall(PetscPrintf(PetscObjectComm((PetscObject)ts), "TSSolve: End\n\n"));

  /* --- */

  PetscCall(TSDestroy(&ts));

  PetscCall(PetscFinalize());
  return 0;
}

/* -------------------------------------------------------------------*/

PetscErrorCode RHSFunction(TS ts, PetscReal t, Vec x, Vec f, void *ctx)
{
  PetscFunctionBeginUser;
  PetscCall(VecSet(f, (PetscReal)1));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode RHSJacobian(TS ts, PetscReal t, Vec x, Mat A, Mat B, void *ctx)
{
  PetscFunctionBeginUser;
  PetscCall(MatZeroEntries(B));
  if (B != A) PetscCall(MatZeroEntries(A));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode PreStep(TS ts)
{
  PetscInt           n;
  PetscReal          t;
  Vec                x;
  const PetscScalar *a;
  PetscBool          flg;

  PetscFunctionBeginUser;
  PetscCall(TSGetStepNumber(ts, &n));
  PetscCall(TSGetTime(ts, &t));
  PetscCall(TSGetSolution(ts, &x));
  PetscCall(VecGetArrayRead(x, &a));
  PetscCall(TSGetStepResize(ts, &flg));
  PetscCall(PetscPrintf(PetscObjectComm((PetscObject)ts), "%-10s-> step %" PetscInt_FMT " time %g value %g%s\n", PETSC_FUNCTION_NAME, n, (double)t, (double)PetscRealPart(a[0]), flg ? " resized" : ""));
  PetscCall(VecRestoreArrayRead(x, &a));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode PostStep(TS ts)
{
  PetscInt           n;
  PetscReal          t;
  Vec                x;
  const PetscScalar *a;

  PetscFunctionBeginUser;
  PetscCall(TSGetStepNumber(ts, &n));
  PetscCall(TSGetTime(ts, &t));
  PetscCall(TSGetSolution(ts, &x));
  PetscCall(VecGetArrayRead(x, &a));
  PetscCall(PetscPrintf(PetscObjectComm((PetscObject)ts), "%-10s-> step %" PetscInt_FMT " time %g value %g\n", PETSC_FUNCTION_NAME, n, (double)t, (double)PetscRealPart(a[0])));
  PetscCall(VecRestoreArrayRead(x, &a));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode Monitor(TS ts, PetscInt n, PetscReal t, Vec x, void *ctx)
{
  const PetscScalar *a;

  PetscFunctionBeginUser;
  PetscCall(VecGetArrayRead(x, &a));
  PetscCall(PetscPrintf(PetscObjectComm((PetscObject)ts), "%-10s-> step %" PetscInt_FMT " time %g value %g\n", PETSC_FUNCTION_NAME, n, (double)t, (double)PetscRealPart(a[0])));
  PetscCall(VecRestoreArrayRead(x, &a));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode Event(TS ts, PetscReal t, Vec x, PetscReal *fvalue, void *ctx)
{
  PetscFunctionBeginUser;
  fvalue[0] = t - 5;
  fvalue[1] = 7 - t;
  fvalue[2] = t - 9;
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode PostEvent(TS ts, PetscInt nevents, PetscInt event_list[], PetscReal t, Vec x, PetscBool forwardsolve, void *ctx)
{
  PetscInt           i;
  const PetscScalar *a;

  PetscFunctionBeginUser;
  PetscCall(TSGetStepNumber(ts, &i));
  PetscCall(VecGetArrayRead(x, &a));
  PetscCall(PetscPrintf(PetscObjectComm((PetscObject)ts), "%-10s-> step %" PetscInt_FMT " time %g value %g\n", PETSC_FUNCTION_NAME, i, (double)t, (double)PetscRealPart(a[0])));
  PetscCall(VecRestoreArrayRead(x, &a));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode TransferSetUp(TS ts, PetscInt n, PetscReal t, Vec x, PetscBool *flg, void *ctx)
{
  PetscInt *nt = (PetscInt *)ctx;

  PetscFunctionBeginUser;
  if (n == 1) {
    nt[0] = 2;
    nt[1] = 2;
  }
  *flg = (PetscBool)(nt[0] && n && n % (nt[0]) == 0);
  if (*flg) nt[0] += nt[1];
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode Transfer(TS ts, PetscInt nv, Vec vin[], Vec vout[], void *ctx)
{
  PetscInt i;

  PetscFunctionBeginUser;
  PetscCall(TSGetStepNumber(ts, &i));
  PetscCall(PetscPrintf(PetscObjectComm((PetscObject)ts), "%-10s-> step %" PetscInt_FMT " nv %" PetscInt_FMT "\n", PETSC_FUNCTION_NAME, i, nv));
  for (i = 0; i < nv; i++) {
    PetscCall(VecDuplicate(vin[i], &vout[i]));
    PetscCall(VecCopy(vin[i], vout[i]));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*TEST

    test:
      suffix: euler
      diff_args: -j
      args: -ts_type euler
      output_file: output/ex1.out

    test:
      suffix: ssp
      diff_args: -j
      args: -ts_type ssp
      output_file: output/ex1.out

    test:
      suffix: rk
      diff_args: -j
      args: -ts_type rk
      output_file: output/ex1.out

    test:
      suffix: beuler
      diff_args: -j
      args: -ts_type beuler
      output_file: output/ex1_theta.out

    test:
      suffix: cn
      diff_args: -j
      args: -ts_type cn
      output_file: output/ex1_theta.out

    test:
      suffix: theta
      args: -ts_type theta
      diff_args: -j
      output_file: output/ex1_theta.out

    test:
      suffix: bdf
      diff_args: -j
      args: -ts_type bdf
      output_file: output/ex1_bdf.out

    test:
      suffix: alpha
      diff_args: -j
      args: -ts_type alpha
      output_file: output/ex1_alpha.out

    test:
      suffix: rosw
      diff_args: -j
      args: -ts_type rosw
      output_file: output/ex1.out

    test:
      suffix: arkimex
      diff_args: -j
      args: -ts_type arkimex
      output_file: output/ex1_arkimex.out

TEST*/
