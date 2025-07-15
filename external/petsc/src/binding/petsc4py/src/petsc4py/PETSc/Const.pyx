# ------------------------------------------------------------------------------

DECIDE    = PETSC_DECIDE
DEFAULT   = PETSC_DEFAULT
DETERMINE = PETSC_DETERMINE
CURRENT   = PETSC_CURRENT
UNLIMITED = PETSC_UNLIMITED

__doc__ += """
Basic constants:

`DECIDE`
    Use a default value for an `int` or `float` parameter.
`DEFAULT`
    Use a default value chosen by PETSc.
`DETERMINE`
    Compute a default value for an `int` or `float` parameter.
    For tolerances this uses the default value from when
    the object's type was set.
`CURRENT`
    Do not change the current value that is set.
`UNLIMITED`
    For a parameter that is a bound, such as the maximum
    number of iterations, do not bound the value.
"""

# ------------------------------------------------------------------------------

INFINITY  = toReal(PETSC_INFINITY)
NINFINITY = toReal(PETSC_NINFINITY)
PINFINITY = toReal(PETSC_INFINITY)

__doc__ += """
More constants:

`INFINITY`
    Very large real value.
`NINFINITY`
    Very large negative real value.
`PINFINITY`
    Very large positive real value, same as `INFINITY`.
"""

# ------------------------------------------------------------------------------


class InsertMode(object):
    """Insertion mode.

    Most commonly used insertion modes are:

    `INSERT`
        Insert provided value/s discarding previous value/s.
    `ADD`
        Add provided value/s to current value/s.
    `MAX`
        Insert the maximum of provided value/s and current value/s.

    See Also
    --------
    petsc.InsertMode

    """
    # native
    NOT_SET_VALUES    = PETSC_NOT_SET_VALUES
    INSERT_VALUES     = PETSC_INSERT_VALUES
    ADD_VALUES        = PETSC_ADD_VALUES
    MAX_VALUES        = PETSC_MAX_VALUES
    INSERT_ALL_VALUES = PETSC_INSERT_ALL_VALUES
    ADD_ALL_VALUES    = PETSC_ADD_ALL_VALUES
    INSERT_BC_VALUES  = PETSC_INSERT_BC_VALUES
    ADD_BC_VALUES     = PETSC_ADD_BC_VALUES
    # aliases
    INSERT     = INSERT_VALUES
    ADD        = ADD_VALUES
    MAX        = MAX_VALUES
    INSERT_ALL = INSERT_ALL_VALUES
    ADD_ALL    = ADD_ALL_VALUES
    INSERT_BC  = INSERT_BC_VALUES
    ADD_BC     = ADD_BC_VALUES

# ------------------------------------------------------------------------------


class ScatterMode(object):
    """Scatter mode.

    Most commonly used scatter modes are:

    `FORWARD`
        Scatter values in the forward direction.
    `REVERSE`
        Scatter values in the reverse direction.

    See Also
    --------
    Scatter.create, Scatter.begin, Scatter.end
    petsc.ScatterMode

    """
    # native
    SCATTER_FORWARD       = PETSC_SCATTER_FORWARD
    SCATTER_REVERSE       = PETSC_SCATTER_REVERSE
    SCATTER_FORWARD_LOCAL = PETSC_SCATTER_FORWARD_LOCAL
    SCATTER_REVERSE_LOCAL = PETSC_SCATTER_REVERSE_LOCAL
    # aliases
    FORWARD       = SCATTER_FORWARD
    REVERSE       = SCATTER_REVERSE
    FORWARD_LOCAL = SCATTER_FORWARD_LOCAL
    REVERSE_LOCAL = SCATTER_REVERSE_LOCAL

# ------------------------------------------------------------------------------


class NormType(object):
    """Norm type.

    Commonly used norm types:

    `N1`
        The one norm.
    `N2`
        The two norm.
    `FROBENIUS`
        The Frobenius norm.
    `INFINITY`
        The infinity norm.

    See Also
    --------
    petsc.NormType

    """
    # native
    NORM_1         = PETSC_NORM_1
    NORM_2         = PETSC_NORM_2
    NORM_1_AND_2   = PETSC_NORM_1_AND_2
    NORM_FROBENIUS = PETSC_NORM_FROBENIUS
    NORM_INFINITY  = PETSC_NORM_INFINITY
    NORM_MAX       = PETSC_NORM_MAX
    # aliases
    N1        = NORM_1
    N2        = NORM_2
    N12       = NORM_1_AND_2
    MAX       = NORM_MAX
    FROBENIUS = NORM_FROBENIUS
    INFINITY  = NORM_INFINITY
    # extra aliases
    FRB = FROBENIUS
    INF = INFINITY

# ------------------------------------------------------------------------------
