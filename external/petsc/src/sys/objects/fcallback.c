#include <petsc/private/petscimpl.h> /*I   "petscsys.h"    I*/

typedef struct _FortranCallbackLink *FortranCallbackLink;
struct _FortranCallbackLink {
  char                  *type_name;
  PetscFortranCallbackId max;
  FortranCallbackLink    next;
};

typedef struct {
  PetscFortranCallbackId basecount;
  PetscFortranCallbackId maxsubtypecount;
  FortranCallbackLink    subtypes;
} FortranCallbackBase;

static FortranCallbackBase *_classbase;
static PetscClassId         _maxclassid = PETSC_SMALLEST_CLASSID;

static PetscErrorCode PetscFortranCallbackFinalize(void)
{
  PetscFunctionBegin;
  for (PetscInt i = PETSC_SMALLEST_CLASSID; i < _maxclassid; i++) {
    FortranCallbackBase *base = &_classbase[i - PETSC_SMALLEST_CLASSID];
    FortranCallbackLink  next, link = base->subtypes;
    for (; link; link = next) {
      next = link->next;
      PetscCall(PetscFree(link->type_name));
      PetscCall(PetscFree(link));
    }
  }
  PetscCall(PetscFree(_classbase));
  _maxclassid = PETSC_SMALLEST_CLASSID;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
  PetscFortranCallbackRegister - register a type+subtype callback. This is used by the PETSc Fortran stubs to allow the use of user Fortran functions
  as arguments to PETSc functions that take function pointers

  Not Collective, No Fortran Support

  Input Parameters:
+ classid - ID of class on which to register callback
- subtype - subtype string, or `NULL` for class ids

  Output Parameter:
. id - callback id

  Level: developer

.seealso: `PetscFortranCallbackGetSizes()`, `PetscObjectCopyFortranFunctionPointers()`, `PetscObjectSetFortranCallback()`, `PetscObjectGetFortranCallback()`
@*/
PetscErrorCode PetscFortranCallbackRegister(PetscClassId classid, const char *subtype, PetscFortranCallbackId *id)
{
  FortranCallbackBase *base;
  FortranCallbackLink  link;

  PetscFunctionBegin;
  if (subtype) PetscAssertPointer(subtype, 2);
  PetscAssertPointer(id, 3);
  PetscCheck(classid >= PETSC_SMALLEST_CLASSID && classid <= PETSC_LARGEST_CLASSID, PETSC_COMM_SELF, PETSC_ERR_ARG_CORRUPT, "ClassId %d corrupt", classid);
  *id = 0;
  if (classid >= _maxclassid) {
    PetscClassId         newmax = PETSC_SMALLEST_CLASSID + 2 * (PETSC_LARGEST_CLASSID - PETSC_SMALLEST_CLASSID);
    FortranCallbackBase *newbase;
    if (!_classbase) PetscCall(PetscRegisterFinalize(PetscFortranCallbackFinalize));
    PetscCall(PetscCalloc1(newmax - PETSC_SMALLEST_CLASSID, &newbase));
    PetscCall(PetscArraycpy(newbase, _classbase, _maxclassid - PETSC_SMALLEST_CLASSID));
    PetscCall(PetscFree(_classbase));

    _classbase  = newbase;
    _maxclassid = newmax;
  }
  base = &_classbase[classid - PETSC_SMALLEST_CLASSID];
  if (!subtype) *id = PETSC_SMALLEST_FORTRAN_CALLBACK + base->basecount++;
  else {
    for (link = base->subtypes; link; link = link->next) { /* look for either both NULL or matching values (implies both non-NULL) */
      PetscBool match;
      PetscCall(PetscStrcmp(subtype, link->type_name, &match));
      if (match) { /* base type or matching subtype */
        goto found;
      }
    }
    /* Not found. Create node and prepend to class' subtype list */
    PetscCall(PetscNew(&link));
    PetscCall(PetscStrallocpy(subtype, &link->type_name));

    link->max      = PETSC_SMALLEST_FORTRAN_CALLBACK;
    link->next     = base->subtypes;
    base->subtypes = link;

  found:
    *id = link->max++;

    base->maxsubtypecount = PetscMax(base->maxsubtypecount, link->max - PETSC_SMALLEST_FORTRAN_CALLBACK);
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*@C
  PetscFortranCallbackGetSizes - get sizes of class and subtype pointer arrays

  Collective, No Fortran Support

  Input Parameter:
. classid - class Id

  Output Parameters:
+ numbase    - number of registered class callbacks
- numsubtype - max number of registered subtype callbacks

  Level: developer

.seealso: `PetscFortranCallbackRegister()`, `PetscObjectCopyFortranFunctionPointers()`, `PetscObjectSetFortranCallback()`, `PetscObjectGetFortranCallback()`
@*/
PetscErrorCode PetscFortranCallbackGetSizes(PetscClassId classid, PetscFortranCallbackId *numbase, PetscFortranCallbackId *numsubtype)
{
  PetscFunctionBegin;
  PetscAssertPointer(numbase, 2);
  PetscAssertPointer(numsubtype, 3);
  if (classid < _maxclassid) {
    FortranCallbackBase *base = &_classbase[classid - PETSC_SMALLEST_CLASSID];
    *numbase                  = base->basecount;
    *numsubtype               = base->maxsubtypecount;
  } else { /* nothing registered */
    *numbase    = 0;
    *numsubtype = 0;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}
