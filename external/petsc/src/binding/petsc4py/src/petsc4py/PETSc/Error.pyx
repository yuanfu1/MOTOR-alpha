# --------------------------------------------------------------------

class Error(RuntimeError):
    """PETSc Error.

    Attributes
    ----------
    ierr : int
        PETSc error code.

    """

    def __init__(self, int ierr=0):
        self.ierr = ierr
        self._traceback = []
        RuntimeError.__init__(self, self.ierr)

    def __bool__(self):
        cdef int ierr = self.ierr
        return ierr != 0

    def __repr__(self):
        return 'PETSc.Error(%d)' % self.ierr

    def __str__(self):
        cdef int csize=1, crank=0
        if not (<int>PetscFinalizeCalled):
            MPI_Comm_size(PETSC_COMM_WORLD, &csize)
            MPI_Comm_rank(PETSC_COMM_WORLD, &crank)
        width, rank = len(str(csize-1)), crank
        tblist = ['error code %d' % self.ierr]
        for entry in self._traceback:
            tbline = '[%*d] %s' % (width, rank, entry)
            tblist.append(tbline)
        return '\n'.join(tblist)

# --------------------------------------------------------------------
