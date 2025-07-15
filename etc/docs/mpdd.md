# mpdd (Multiprocess for domain decomposition)

`mpdd` is a class saving the global MPI communicators and some basic MPI functions. In the MOTOR-DA, mpdd is always required, and has to be initialized firstly:

```fortran
! link the lib `mpdd` is required
  USE mpddGlob_m, only: mpddGlob_t 
  TYPE(mpddGlob_t), target :: mpddGlob   ! Defined as the target for sharing with other classes.

  CALL mpdd%initialize        ! Initialize
  ...                         ! Other codes
  CALL mpdd%destroy           ! Always destory in the end
```
### members in mpdd

```fortran
  TYPE, PUBLIC :: mpddGlob_t
    INTEGER(i_kind) :: nProc                  !< Number of total processes running.
    INTEGER(i_kind) :: rank                   !< Rank number of current process.
    INTEGER(i_kind) :: group                  !< World group number of current process
    INTEGER(i_kind) :: comm                   !< MPI_COMM_WORLD
    INTEGER(i_kind) :: myrank                        !< myrank = rank + 1
    INTEGER(i_kind) :: ierr                          !< Error hanlder of MPI
    INTEGER(i_kind) :: STATUS(MPI_STATUS_SIZE)       !< Status info of MPI

  CONTAINS
    PROCEDURE, PUBLIC, PASS :: initialize            !< Initialization
    PROCEDURE, PUBLIC, PASS :: destroy               !< Deconstructor
    PROCEDURE, PUBLIC, PASS :: isBaseProc            !< Check if the proc is base proc (rank = 0)
    PROCEDURE, PUBLIC, PASS :: isNProcMatch          !< Reserved

    PROCEDURE, PUBLIC, PASS :: barrier               !< MPI barrier in global

    GENERIC, PUBLIC :: bCast&                        !< Boardcast the varibales in single int, single real and 1D int
    ...
```
