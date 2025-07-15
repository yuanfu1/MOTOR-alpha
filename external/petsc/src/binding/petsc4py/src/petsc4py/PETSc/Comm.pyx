# --------------------------------------------------------------------

cdef class Comm:
    """Communicator object.

    Predefined instances:

    `COMM_NULL`
        The *null* (or invalid) communicator.
    `COMM_SELF`
        The *self* communicator.
    `COMM_WORLD`
        The *world* communicator.

    See Also
    --------
    Sys.setDefaultComm, Sys.getDefaultComm

    """

    #

    def __cinit__(self, comm=None):
        self.comm = def_Comm(comm, MPI_COMM_NULL)
        self.isdup = 0
        if self.comm != MPI_COMM_NULL:
            self.base = comm
        else:
            self.base = None

    def __dealloc__(self):
        if self.isdup:
            CHKERR(PetscCommDEALLOC(&self.comm))
        self.comm = MPI_COMM_NULL
        self.isdup = 0
        self.base = None

    def __richcmp__(self, other, int op):
        if not isinstance(self,  Comm): return NotImplemented
        if not isinstance(other, Comm): return NotImplemented
        if op!=2 and op!=3: raise TypeError("only '==' and '!='")
        cdef Comm s = self
        cdef Comm o = other
        cdef int eq = (op == 2)
        cdef MPI_Comm comm1 = s.comm
        cdef MPI_Comm comm2 = o.comm
        cdef int flag = 0
        if comm1 != MPI_COMM_NULL and comm2 != MPI_COMM_NULL:
            CHKERR(<PetscErrorCode>MPI_Comm_compare(comm1, comm2, &flag))
            if eq: return (flag==<int>MPI_IDENT or  flag==<int>MPI_CONGRUENT)
            else:  return (flag!=<int>MPI_IDENT and flag!=<int>MPI_CONGRUENT)
        else:
            if eq: return (comm1 == comm2)
            else:  return (comm1 != comm2)

    def __bool__(self) -> bool:
        return self.comm != MPI_COMM_NULL

    #

    def destroy(self) -> None:
        """Destroy the communicator.

        Collective.

        See Also
        --------
        petsc.PetscCommDestroy

        """
        if self.comm == MPI_COMM_NULL: return
        if not self.isdup:
            raise ValueError("communicator not owned")
        CHKERR(PetscCommDestroy(&self.comm))
        self.comm = MPI_COMM_NULL
        self.isdup = 0
        self.base = None

    def duplicate(self) -> Self:
        """Duplicate the communicator.

        Collective.

        See Also
        --------
        petsc.PetscCommDuplicate

        """
        if self.comm == MPI_COMM_NULL:
            raise ValueError("null communicator")
        cdef MPI_Comm newcomm = MPI_COMM_NULL
        CHKERR(PetscCommDuplicate(self.comm, &newcomm, NULL))
        cdef Comm comm = type(self)()
        comm.comm  = newcomm
        comm.isdup = 1
        comm.base = self.base
        return comm

    def getSize(self) -> int:
        """Return the number of processes in the communicator.

        Not collective.

        """
        if self.comm == MPI_COMM_NULL:
            raise ValueError("null communicator")
        cdef int size=0
        CHKERRMPI(MPI_Comm_size(self.comm, &size))
        return size

    def getRank(self) -> int:
        """Return the rank of the calling processes in the communicator.

        Not collective.

        """
        if self.comm == MPI_COMM_NULL:
            raise ValueError("null communicator")
        cdef int rank=0
        CHKERRMPI(MPI_Comm_rank(self.comm, &rank))
        return rank

    def barrier(self) -> None:
        """Barrier synchronization.

        Collective.

        """
        if self.comm == MPI_COMM_NULL:
            raise ValueError("null communicator")
        CHKERRMPI(MPI_Barrier(self.comm))

    # --- properties ---

    property size:
        """Communicator size."""
        def __get__(self) -> int:
            return self.getSize()

    property rank:
        """Communicator rank."""
        def __get__(self) -> int:
            return self.getRank()

    # --- Fortran support ---

    property fortran:
        """Fortran handle."""
        def __get__(self) -> int:
            cdef MPI_Comm comm = self.comm
            return MPI_Comm_c2f(comm)

    # --- mpi4py support ---

    def tompi4py(self) -> Intracomm:
        """Convert communicator to `mpi4py`.

        Not collective.

        See Also
        --------
        mpi4py.MPI.Comm, mpi4py.MPI.Intracomm

        """
        cdef MPI_Comm comm = self.comm
        return mpi4py_Comm_New(comm)

    # --- mpi4py compatibility API ---

    Free     = destroy
    Clone    = duplicate
    Dup      = duplicate
    Get_size = getSize
    Get_rank = getRank
    Barrier  = barrier

# --------------------------------------------------------------------

cdef Comm __COMM_NULL__  = Comm()
cdef Comm __COMM_SELF__  = Comm()
cdef Comm __COMM_WORLD__ = Comm()

COMM_NULL  = __COMM_NULL__
COMM_SELF  = __COMM_SELF__
COMM_WORLD = __COMM_WORLD__

# --------------------------------------------------------------------

cdef MPI_Comm PETSC_COMM_DEFAULT = MPI_COMM_NULL

cdef MPI_Comm GetComm(
    object comm, MPI_Comm defv,
) except? MPI_COMM_NULL:
    return def_Comm(comm, defv)

cdef MPI_Comm GetCommDefault():
    return PETSC_COMM_DEFAULT

# --------------------------------------------------------------------
