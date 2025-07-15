#include <petsc/private/ftnimpl.h>
#include <petsc/private/kspimpl.h>

#if defined(PETSC_HAVE_FORTRAN_CAPS)
  #define dmkspsetcomputerhs_          DMKSPSETCOMPUTERHS
  #define dmkspsetcomputeinitialguess_ DMKSPSETCOMPUTEINITIALGUESS
  #define dmkspsetcomputeoperators_    DMKSPSETCOMPUTEOPERATORS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE)
  #define dmkspsetcomputerhs_          dmkspsetcomputerhs
  #define dmkspsetcomputeinitialguess_ dmkspsetcomputeinitialguess
  #define dmkspsetcomputeoperators_    dmkspsetcomputeoperators
#endif

static PetscErrorCode ourkspcomputerhs(KSP ksp, Vec b, void *ctx)
{
  DM    dm;
  DMKSP kdm;
  PetscCall(KSPGetDM(ksp, &dm));
  PetscCall(DMGetDMKSP(dm, &kdm));
  PetscCallFortranVoidFunction((*(void (*)(KSP *, Vec *, void *, PetscErrorCode *))kdm->fortran_func_pointers[0])(&ksp, &b, ctx, &ierr));
  return PETSC_SUCCESS;
}

static PetscErrorCode ourkspcomputeinitialguess(KSP ksp, Vec b, void *ctx)
{
  DM    dm;
  DMKSP kdm;
  PetscCall(KSPGetDM(ksp, &dm));
  PetscCall(DMGetDMKSP(dm, &kdm));
  PetscCallFortranVoidFunction((*(void (*)(KSP *, Vec *, void *, PetscErrorCode *))kdm->fortran_func_pointers[2])(&ksp, &b, ctx, &ierr));
  return PETSC_SUCCESS;
}

static PetscErrorCode ourkspcomputeoperators(KSP ksp, Mat A, Mat B, void *ctx)
{
  DM    dm;
  DMKSP kdm;
  PetscCall(KSPGetDM(ksp, &dm));
  PetscCall(DMGetDMKSP(dm, &kdm));
  PetscCallFortranVoidFunction((*(void (*)(KSP *, Mat *, Mat *, void *, PetscErrorCode *))kdm->fortran_func_pointers[1])(&ksp, &A, &B, ctx, &ierr));
  return PETSC_SUCCESS;
}

/* The counting for fortran_func_pointers is insanely brittle. We're putting these inside the base DM, but we have no
 * way to be sure there is room other than to grep the sources from src/dm (and any other possible client). Fortran
 * function pointers need an overhaul.
 */

PETSC_EXTERN void dmkspsetcomputerhs_(DM *dm, void (*func)(KSP *, Vec *, void *, PetscErrorCode *), void *ctx, PetscErrorCode *ierr)
{
  DMKSP kdm;
  *ierr = DMGetDMKSP(*dm, &kdm);
  if (!*ierr) {
    kdm->fortran_func_pointers[0] = (PetscVoidFn *)func;
    *ierr                         = DMKSPSetComputeRHS(*dm, ourkspcomputerhs, ctx);
  }
}

PETSC_EXTERN void dmkspsetcomputeinitialguess_(DM *dm, void (*func)(KSP *, Vec *, void *, PetscErrorCode *), void *ctx, PetscErrorCode *ierr)
{
  DMKSP kdm;
  *ierr = DMGetDMKSP(*dm, &kdm);
  if (!*ierr) {
    kdm->fortran_func_pointers[2] = (PetscVoidFn *)func;

    *ierr = DMKSPSetComputeInitialGuess(*dm, ourkspcomputeinitialguess, ctx);
  }
}

PETSC_EXTERN void dmkspsetcomputeoperators_(DM *dm, void (*func)(KSP *, Vec *, void *, PetscErrorCode *), void *ctx, PetscErrorCode *ierr)
{
  DMKSP kdm;
  *ierr = DMGetDMKSP(*dm, &kdm);
  if (!*ierr) {
    kdm->fortran_func_pointers[1] = (PetscVoidFn *)func;
    *ierr                         = DMKSPSetComputeOperators(*dm, ourkspcomputeoperators, ctx);
  }
}
