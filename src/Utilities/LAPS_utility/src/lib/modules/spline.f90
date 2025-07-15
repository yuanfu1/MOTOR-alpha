
! Downloaded by Steve Albers from an Iowa State Math Dept. Website produced
! by John Burkhardt:

! http://orion.math.iastate.edu/burkardt/f_src/spline/spline.html

SUBROUTINE basis_function_b_val(tdata, tval, yval)
!
!*******************************************************************************
!
!! BASIS_FUNCTION_B_VAL evaluates the B spline basis function.
!
!
!  Discussion:
!
!    The B spline basis function is a piecewise cubic which
!    has the properties that:
!
!    * it equals 2/3 at TDATA(3), 1/6 at TDATA(2) and TDATA(4);
!    * it is 0 for TVAL <= TDATA(1) or TDATA(5) <= TVAL;
!    * it is strictly increasing from TDATA(1) to TDATA(3),
!      and strictly decreasing from TDATA(3) to TDATA(5);
!    * the function and its first two derivatives are continuous
!      at each node TDATA(I).
!
!  Reference:
!
!    Alan Davies and Philip Samuels,
!    An Introduction to Computational Geometry for Curves and Surfaces,
!    Clarendon Press, 1996.
!
!  Modified:
!
!    06 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real TDATA(5), the nodes associated with the basis function.
!    The entries of TDATA are assumed to be distinct and increasing.
!
!    Input, real TVAL, a point at which the B spline basis function is
!    to be evaluated.
!
!    Output, real YVAL, the value of the function at TVAL.
!
  IMPLICIT NONE
!
  INTEGER, PARAMETER :: ndata = 5
!
  INTEGER left
  INTEGER right
  REAL tdata(ndata)
  REAL tval
  REAL u
  REAL yval
!
  IF (tval <= tdata(1) .OR. tval >= tdata(ndata)) THEN
    yval = 0.0E+00
    RETURN
  END IF
!
!  Find the interval [ TDATA(LEFT), TDATA(RIGHT) ] containing TVAL.
!
  CALL rvec_bracket(ndata, tdata, tval, left, right)
!
!  U is the normalized coordinate of TVAL in this interval.
!
  u = (tval - tdata(left)) / (tdata(right) - tdata(left))
!
!  Now evaluate the function.
!
  IF (tval < tdata(2)) THEN
    yval = u**3 / 6.0E+00
  ELSE IF (tval < tdata(3)) THEN
    yval = (1.0E+00 + 3.0E+00 * u + 3.0E+00 * u**2 - 3.0E+00 * u**3) / 6.0E+00
  ELSE IF (tval < tdata(4)) THEN
    yval = (4.0E+00 - 6.0E+00 * u**2 + 3.0E+00 * u**3) / 6.0E+00
  ELSE IF (tval < tdata(5)) THEN
    yval = (1.0E+00 - u)**3 / 6.0E+00
  END IF

  RETURN
END
SUBROUTINE basis_function_beta_val(beta1, beta2, tdata, tval, yval)
!
!*******************************************************************************
!
!! BASIS_FUNCTION_BETA_VAL evaluates the beta spline basis function.
!
!
!  Discussion:
!
!    With BETA1 = 1 and BETA2 = 0, the beta spline basis function
!    equals the B spline basis function.
!
!    With BETA1 large, and BETA2 = 0, the beta spline basis function
!    skews to the right, that is, its maximum increases, and occurs
!    to the right of the center.
!
!    With BETA1 = 1 and BETA2 large, the beta spline becomes more like
!    a linear basis function; that is, its value in the outer two intervals
!    goes to zero, and its behavior in the inner two intervals approaches
!    a piecewise linear function that is 0 at the second node, 1 at the
!    third, and 0 at the fourth.
!
!  Reference:
!
!    Alan Davies and Philip Samuels,
!    An Introduction to Computational Geometry for Curves and Surfaces,
!    Clarendon Press, 1996, page 129.
!
!  Modified:
!
!    09 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real BETA1, the skew or bias parameter.
!    BETA1 = 1 for no skew or bias.
!
!    Input, real BETA2, the tension parameter.
!    BETA2 = 0 for no tension.
!
!    Input, real TDATA(5), the nodes associated with the basis function.
!    The entries of TDATA are assumed to be distinct and increasing.
!
!    Input, real TVAL, a point at which the B spline basis function is
!    to be evaluated.
!
!    Output, real YVAL, the value of the function at TVAL.
!
  IMPLICIT NONE
!
  INTEGER, PARAMETER :: ndata = 5
!
  REAL a
  REAL b
  REAL beta1
  REAL beta2
  REAL c
  REAL d
  INTEGER left
  INTEGER right
  REAL tdata(ndata)
  REAL tval
  REAL u
  REAL yval
!
  IF (tval <= tdata(1) .OR. tval >= tdata(ndata)) THEN
    yval = 0.0E+00
    RETURN
  END IF
!
!  Find the interval [ TDATA(LEFT), TDATA(RIGHT) ] containing TVAL.
!
  CALL rvec_bracket(ndata, tdata, tval, left, right)
!
!  U is the normalized coordinate of TVAL in this interval.
!
  u = (tval - tdata(left)) / (tdata(right) - tdata(left))
!
!  Now evaluate the function.
!
  IF (tval < tdata(2)) THEN

    yval = 2.0E+00 * u**3

  ELSE IF (tval < tdata(3)) THEN

    a = beta2 + 4.0E+00 * beta1 + 4.0E+00 * beta1**2 &
        + 6.0E+00 * (1.0E+00 - beta1**2) &
        - 3.0E+00 * (2.0E+00 + beta2 + 2.0E+00 * beta1) &
        + 2.0E+00 * (1.0E+00 + beta2 + beta1 + beta1**2)

    b = -6.0E+00 * (1.0E+00 - beta1**2) &
        + 6.0E+00 * (2.0E+00 + beta2 + 2.0E+00 * beta1) &
        - 6.0E+00 * (1.0E+00 + beta2 + beta1 + beta1**2)

    c = -3.0E+00 * (2.0E+00 + beta2 + 2.0E+00 * beta1) &
        + 6.0E+00 * (1.0E+00 + beta2 + beta1 + beta1**2)

    d = -2.0E+00 * (1.0E+00 + beta2 + beta1 + beta1**2)

    yval = a + b * u + c * u**2 + d * u**3

  ELSE IF (tval < tdata(4)) THEN

    a = beta2 + 4.0E+00 * beta1 + 4.0E+00 * beta1**2

    b = -6.0E+00 * (beta1 - beta1**3)

    c = -3.0E+00 * (beta2 + 2.0E+00 * beta1**2 + 2.0E+00 * beta1**3)

    d = 2.0E+00 * (beta2 + beta1 + beta1**2 + beta1**3)

    yval = a + b * u + c * u**2 + d * u**3

  ELSE IF (tval < tdata(5)) THEN

    yval = 2.0E+00 * beta1**3 * (1.0E+00 - u)**3

  END IF

  yval = yval / (2.0E+00 + beta2 + 4.0E+00 * beta1 + 4.0E+00 * beta1**2 &
                 + 2.0E+00 * beta1**3)

  RETURN
END
SUBROUTINE basis_matrix_b_uni(mbasis)
!
!*******************************************************************************
!
!! BASIS_MATRIX_B_UNI sets up the uniform B spline basis matrix.
!
!
!  Reference:
!
!    Foley, van Dam, Feiner, Hughes,
!    Computer Graphics: Principles and Practice,
!    page 493.
!
!  Modified:
!
!    05 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real MBASIS(4,4), the basis matrix.
!
  IMPLICIT NONE
!
  REAL mbasis(4, 4)
!
  mbasis(1, 1) = -1.0E+00 / 6.0E+00
  mbasis(1, 2) = 3.0E+00 / 6.0E+00
  mbasis(1, 3) = -3.0E+00 / 6.0E+00
  mbasis(1, 4) = 1.0E+00 / 6.0E+00

  mbasis(2, 1) = 3.0E+00 / 6.0E+00
  mbasis(2, 2) = -6.0E+00 / 6.0E+00
  mbasis(2, 3) = 3.0E+00 / 6.0E+00
  mbasis(2, 4) = 0.0E+00

  mbasis(3, 1) = -3.0E+00 / 6.0E+00
  mbasis(3, 2) = 0.0E+00
  mbasis(3, 3) = 3.0E+00 / 6.0E+00
  mbasis(3, 4) = 0.0E+00

  mbasis(4, 1) = 1.0E+00 / 6.0E+00
  mbasis(4, 2) = 4.0E+00 / 6.0E+00
  mbasis(4, 3) = 1.0E+00 / 6.0E+00
  mbasis(4, 4) = 0.0E+00

  RETURN
END
SUBROUTINE basis_matrix_beta_uni(beta1, beta2, mbasis)
!
!*******************************************************************************
!
!! BASIS_MATRIX_BETA_UNI sets up the uniform beta spline basis matrix.
!
!
!  Discussion:
!
!    If BETA1 = 1 and BETA2 = 0, then the beta spline reduces to
!    the B spline.
!
!  Reference:
!
!    Foley, van Dam, Feiner, Hughes,
!    Computer Graphics: Principles and Practice,
!    page 505.
!
!  Modified:
!
!    03 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real BETA1, the skew or bias parameter.
!    BETA1 = 1 for no skew or bias.
!
!    Input, real BETA2, the tension parameter.
!    BETA2 = 0 for no tension.
!
!    Output, real MBASIS(4,4), the basis matrix.
!
  IMPLICIT NONE
!
  REAL beta1
  REAL beta2
  REAL delta
  INTEGER i
  INTEGER j
  REAL mbasis(4, 4)
!
  mbasis(1, 1) = -2.0E+00 * beta1**3
  mbasis(1, 2) = 2.0E+00 * (beta2 + beta1**3 + beta1**2 + beta1)
  mbasis(1, 3) = -2.0E+00 * (beta2 + beta1**2 + beta1 + 1.0E+00)
  mbasis(1, 4) = 2.0E+00

  mbasis(2, 1) = 6.0E+00 * beta1**3
  mbasis(2, 2) = -3.0E+00 * (beta2 + 2.0E+00 * beta1**3 + 2.0E+00 * beta1**2)
  mbasis(2, 3) = 3.0E+00 * (beta2 + 2.0E+00 * beta1**2)
  mbasis(2, 4) = 0.0E+00

  mbasis(3, 1) = -6.0E+00 * beta1**3
  mbasis(3, 2) = 6.0E+00 * beta1 * (beta1 - 1.0E+00) * (beta1 + 1.0E+00)
  mbasis(3, 3) = 6.0E+00 * beta1
  mbasis(3, 4) = 0.0E+00

  mbasis(4, 1) = 2.0E+00 * beta1**3
  mbasis(4, 2) = 4.0E+00 * beta1 * (beta1 + 1.0E+00) + beta2
  mbasis(4, 3) = 2.0E+00
  mbasis(4, 4) = 0.0E+00

  delta = beta2 + 2.0E+00 * beta1**3 + 4.0E+00 * beta1**2 &
          + 4.0E+00 * beta1 + 2.0E+00

  mbasis(1:4, 1:4) = mbasis(1:4, 1:4) / delta

  RETURN
END
SUBROUTINE basis_matrix_bezier(mbasis)
!
!*******************************************************************************
!
!! BASIS_MATRIX_BEZIER_UNI sets up the cubic Bezier spline basis matrix.
!
!
!  Discussion:
!
!    This basis matrix assumes that the data points are stored as
!    ( P1, P2, P3, P4 ).  P1 is the function value at T = 0, while
!    P2 is used to approximate the derivative at T = 0 by
!    dP/dt = 3 * ( P2 - P1 ).  Similarly, P4 is the function value
!    at T = 1, and P3 is used to approximate the derivative at T = 1
!    by dP/dT = 3 * ( P4 - P3 ).
!
!  Reference:
!
!    Foley, van Dam, Feiner, Hughes,
!    Computer Graphics: Principles and Practice,
!    page 489.
!
!  Modified:
!
!    05 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real MBASIS(4,4), the basis matrix.
!
  IMPLICIT NONE
!
  REAL mbasis(4, 4)
!
  mbasis(1, 1) = -1.0E+00
  mbasis(1, 2) = 3.0E+00
  mbasis(1, 3) = -3.0E+00
  mbasis(1, 4) = 1.0E+00

  mbasis(2, 1) = 3.0E+00
  mbasis(2, 2) = -6.0E+00
  mbasis(2, 3) = 3.0E+00
  mbasis(2, 4) = 0.0E+00

  mbasis(3, 1) = -3.0E+00
  mbasis(3, 2) = 3.0E+00
  mbasis(3, 3) = 0.0E+00
  mbasis(3, 4) = 0.0E+00

  mbasis(4, 1) = 1.0E+00
  mbasis(4, 2) = 0.0E+00
  mbasis(4, 3) = 0.0E+00
  mbasis(4, 4) = 0.0E+00

  RETURN
END
SUBROUTINE basis_matrix_hermite(mbasis)
!
!*******************************************************************************
!
!! BASIS_MATRIX_HERMITE sets up the Hermite spline basis matrix.
!
!
!  Discussion:
!
!    This basis matrix assumes that the data points are stored as
!    ( P1, P2, P1', P2' ), with P1 and P1' being the data value and
!    the derivative dP/dT at T = 0, while P2 and P2' apply at T = 1.
!
!  Reference:
!
!    Foley, van Dam, Feiner, Hughes,
!    Computer Graphics: Principles and Practice,
!    page 484.
!
!  Modified:
!
!    06 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real MBASIS(4,4), the basis matrix.
!
  IMPLICIT NONE
!
  REAL mbasis(4, 4)
!
  mbasis(1, 1) = 2.0E+00
  mbasis(1, 2) = -2.0E+00
  mbasis(1, 3) = 1.0E+00
  mbasis(1, 4) = 1.0E+00

  mbasis(2, 1) = -3.0E+00
  mbasis(2, 2) = 3.0E+00
  mbasis(2, 3) = -2.0E+00
  mbasis(2, 4) = -1.0E+00

  mbasis(3, 1) = 0.0E+00
  mbasis(3, 2) = 0.0E+00
  mbasis(3, 3) = 1.0E+00
  mbasis(3, 4) = 0.0E+00

  mbasis(4, 1) = 1.0E+00
  mbasis(4, 2) = 0.0E+00
  mbasis(4, 3) = 0.0E+00
  mbasis(4, 4) = 0.0E+00

  RETURN
END
SUBROUTINE basis_matrix_overhauser_nonuni(alpha, beta, mbasis)
!
!*******************************************************************************
!
!! BASIS_MATRIX_OVERHAUSER_NONUNI sets up the nonuniform Overhauser spline basis matrix.
!
!
!  Discussion:
!
!    This basis matrix assumes that the data points P1, P2, P3 and
!    P4 are not uniformly spaced in T, and that P2 corresponds to T = 0,
!    and P3 to T = 1.
!
!  Modified:
!
!    05 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ALPHA, BETA.
!    ALPHA = | P2 - P1 | / ( | P3 - P2 | + | P2 - P1 | )
!    BETA  = | P3 - P2 | / ( | P4 - P3 | + | P3 - P2 | ).
!
!    Output, real MBASIS(4,4), the basis matrix.
!
  IMPLICIT NONE
!
  REAL alpha
  REAL beta
  REAL mbasis(4, 4)
!
  mbasis(1, 1) = -(1.0E+00 - alpha)**2 / alpha
  mbasis(1, 2) = beta + (1.0E+00 - alpha) / alpha
  mbasis(1, 3) = alpha - 1.0E+00 / (1.0E+00 - beta)
  mbasis(1, 4) = beta**2 / (1.0E+00 - beta)

  mbasis(2, 1) = 2.0E+00 * (1.0E+00 - alpha)**2 / alpha
  mbasis(2, 2) = (-2.0E+00 * (1.0E+00 - alpha) - alpha * beta) / alpha
  mbasis(2, 3) = (2.0E+00 * (1.0E+00 - alpha) &
                  - beta * (1.0E+00 - 2.0E+00 * alpha)) / (1.0E+00 - beta)
  mbasis(2, 4) = -beta**2 / (1.0E+00 - beta)

  mbasis(3, 1) = -(1.0E+00 - alpha)**2 / alpha
  mbasis(3, 2) = (1.0E+00 - 2.0E+00 * alpha) / alpha
  mbasis(3, 3) = alpha
  mbasis(3, 4) = 0.0E+00

  mbasis(4, 1) = 0.0E+00
  mbasis(4, 2) = 1.0E+00
  mbasis(4, 3) = 0.0E+00
  mbasis(4, 4) = 0.0E+00

  RETURN
END
SUBROUTINE basis_matrix_overhauser_nul(alpha, mbasis)
!
!*******************************************************************************
!
!! BASIS_MATRIX_OVERHAUSER_NUL sets up the left nonuniform Overhauser spline basis matrix.
!
!
!  Discussion:
!
!    This basis matrix assumes that the data points P1, P2, and
!    P3 are not uniformly spaced in T, and that P1 corresponds to T = 0,
!    and P2 to T = 1. (???)
!
!  Modified:
!
!    27 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ALPHA.
!    ALPHA = | P2 - P1 | / ( | P3 - P2 | + | P2 - P1 | )
!
!    Output, real MBASIS(3,3), the basis matrix.
!
  IMPLICIT NONE
!
  REAL alpha
  REAL mbasis(3, 3)
!
  mbasis(1, 1) = 1.0E+00 / alpha
  mbasis(1, 2) = -1.0E+00 / (alpha * (1.0E+00 - alpha))
  mbasis(1, 3) = 1.0E+00 / (1.0E+00 - alpha)

  mbasis(2, 1) = -(1.0E+00 + alpha) / alpha
  mbasis(2, 2) = 1.0E+00 / (alpha * (1.0E+00 - alpha))
  mbasis(2, 3) = -alpha / (1.0E+00 - alpha)

  mbasis(3, 1) = 1.0E+00
  mbasis(3, 2) = 0.0E+00
  mbasis(3, 3) = 0.0E+00

  RETURN
END
SUBROUTINE basis_matrix_overhauser_nur(beta, mbasis)
!
!*******************************************************************************
!
!! BASIS_MATRIX_OVERHAUSER_NUR sets up the right nonuniform Overhauser spline basis matrix.
!
!
!  Discussion:
!
!    This basis matrix assumes that the data points PN-2, PN-1, and
!    PN are not uniformly spaced in T, and that PN-1 corresponds to T = 0,
!    and PN to T = 1. (???)
!
!  Modified:
!
!    27 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real BETA.
!    BETA = | P(N) - P(N-1) | / ( | P(N) - P(N-1) | + | P(N-1) - P(N-2) | )
!
!    Output, real MBASIS(3,3), the basis matrix.
!
  IMPLICIT NONE
!
  REAL beta
  REAL mbasis(3, 3)
!
  mbasis(1, 1) = 1.0E+00 / beta
  mbasis(1, 2) = -1.0E+00 / (beta * (1.0E+00 - beta))
  mbasis(1, 3) = 1.0E+00 / (1.0E+00 - beta)

  mbasis(2, 1) = -(1.0E+00 + beta) / beta
  mbasis(2, 2) = 1.0E+00 / (beta * (1.0E+00 - beta))
  mbasis(2, 3) = -beta / (1.0E+00 - beta)

  mbasis(3, 1) = 1.0E+00
  mbasis(3, 2) = 0.0E+00
  mbasis(3, 3) = 0.0E+00

  RETURN
END
SUBROUTINE basis_matrix_overhauser_uni(mbasis)
!
!*******************************************************************************
!
!! BASIS_MATRIX_OVERHAUSER_UNI sets up the uniform Overhauser spline basis matrix.
!
!
!  Discussion:
!
!    This basis matrix assumes that the data points P1, P2, P3 and
!    P4 are uniformly spaced in T, and that P2 corresponds to T = 0,
!    and P3 to T = 1.
!
!  Reference:
!
!    Foley, van Dam, Feiner, Hughes,
!    Computer Graphics: Principles and Practice,
!    page 505.
!
!  Modified:
!
!    05 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real MBASIS(4,4), the basis matrix.
!
  IMPLICIT NONE
!
  REAL mbasis(4, 4)
!
  mbasis(1, 1) = -1.0E+00 / 2.0E+00
  mbasis(1, 2) = 3.0E+00 / 2.0E+00
  mbasis(1, 3) = -3.0E+00 / 2.0E+00
  mbasis(1, 4) = 1.0E+00 / 2.0E+00

  mbasis(2, 1) = 2.0E+00 / 2.0E+00
  mbasis(2, 2) = -5.0E+00 / 2.0E+00
  mbasis(2, 3) = 4.0E+00 / 2.0E+00
  mbasis(2, 4) = -1.0E+00 / 2.0E+00

  mbasis(3, 1) = -1.0E+00 / 2.0E+00
  mbasis(3, 2) = 0.0E+00
  mbasis(3, 3) = 1.0E+00 / 2.0E+00
  mbasis(3, 4) = 0.0E+00

  mbasis(4, 1) = 0.0E+00
  mbasis(4, 2) = 2.0E+00 / 2.0E+00
  mbasis(4, 3) = 0.0E+00
  mbasis(4, 4) = 0.0E+00

  RETURN
END
SUBROUTINE basis_matrix_overhauser_uni_l(mbasis)
!
!*******************************************************************************
!
!! BASIS_MATRIX_OVERHAUSER_UNI_L sets up the left uniform Overhauser spline basis matrix.
!
!
!  Discussion:
!
!    This basis matrix assumes that the data points P1, P2, and P3
!    are not uniformly spaced in T, and that P1 corresponds to T = 0,
!    and P2 to T = 1.
!
!  Modified:
!
!    05 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real MBASIS(3,3), the basis matrix.
!
  IMPLICIT NONE
!
  REAL mbasis(3, 3)
!
  mbasis(1, 1) = 2.0E+00
  mbasis(1, 2) = -4.0E+00
  mbasis(1, 3) = 2.0E+00

  mbasis(2, 1) = -3.0E+00
  mbasis(2, 2) = 4.0E+00
  mbasis(2, 3) = -1.0E+00

  mbasis(3, 1) = 1.0E+00
  mbasis(3, 2) = 0.0E+00
  mbasis(3, 3) = 0.0E+00

  RETURN
END
SUBROUTINE basis_matrix_overhauser_uni_r(mbasis)
!
!*******************************************************************************
!
!! BASIS_MATRIX_OVERHAUSER_UNI_R sets up the right uniform Overhauser spline basis matrix.
!
!
!  Discussion:
!
!    This basis matrix assumes that the data points P(N-2), P(N-1),
!    and P(N) are uniformly spaced in T, and that P(N-1) corresponds to
!    T = 0, and P(N) to T = 1.
!
!  Modified:
!
!    05 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real MBASIS(3,3), the basis matrix.
!
  IMPLICIT NONE
!
  REAL mbasis(3, 3)
!
  mbasis(1, 1) = 2.0E+00
  mbasis(1, 2) = -4.0E+00
  mbasis(1, 3) = 2.0E+00

  mbasis(2, 1) = -3.0E+00
  mbasis(2, 2) = 4.0E+00
  mbasis(2, 3) = -1.0E+00

  mbasis(3, 1) = 1.0E+00
  mbasis(3, 2) = 0.0E+00
  mbasis(3, 3) = 0.0E+00

  RETURN
END
SUBROUTINE basis_matrix_tmp(left, n, mbasis, ndata, tdata, ydata, tval, yval)
!
!*******************************************************************************
!
!! BASIS_MATRIX_TMP computes Q = T * MBASIS * P
!
!
!  Discussion:
!
!    YDATA is a vector of data values, most frequently the values of some
!    function sampled at uniformly spaced points.  MBASIS is the basis
!    matrix for a particular kind of spline.  T is a vector of the
!    powers of the normalized difference between TVAL and the left
!    endpoint of the interval.
!
!  Modified:
!
!    06 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LEFT, indicats that TVAL is in the interval
!    [ TDATA(LEFT), TDATA(LEFT+1) ], or that this is the "nearest"
!    interval to TVAL.
!    For TVAL < TDATA(1), use LEFT = 1.
!    For TVAL > TDATA(NDATA), use LEFT = NDATA - 1.
!
!    Input, integer N, the order of the basis matrix.
!
!    Input, real MBASIS(N,N), the basis matrix.
!
!    Input, integer NDATA, the dimension of the vectors TDATA and YDATA.
!
!    Input, real TDATA(NDATA), the abscissa values.  This routine
!    assumes that the TDATA values are uniformly spaced, with an
!    increment of 1.0.
!
!    Input, real YDATA(NDATA), the data values to be interpolated or
!    approximated.
!
!    Input, real TVAL, the value of T at which the spline is to be
!    evaluated.
!
!    Output, real YVAL, the value of the spline at TVAL.
!
  IMPLICIT NONE
!
  INTEGER, PARAMETER :: maxn = 4
  INTEGER n
  INTEGER ndata
!
  REAL arg
  INTEGER first
  INTEGER i
  INTEGER j
  INTEGER left
  REAL mbasis(n, n)
  REAL tdata(ndata)
  REAL temp
  REAL tval
  REAL tvec(maxn)
  REAL ydata(ndata)
  REAL yval
!
  IF (left == 1) THEN
    arg = 0.5 * (tval - tdata(left))
    first = left
  ELSE IF (left < ndata - 1) THEN
    arg = tval - tdata(left)
    first = left - 1
  ELSE IF (left == ndata - 1) THEN
    arg = 0.5E+00 * (1.0E+00 + tval - tdata(left))
    first = left - 1
  END IF

  DO i = 1, n - 1
    tvec(i) = arg**(n - i)
  END DO
  tvec(n) = 1.0E+00

  yval = 0.0E+00
  DO j = 1, n
    yval = yval + DOT_PRODUCT(tvec(1:n), mbasis(1:n, j)) &
           * ydata(first - 1 + j)
  END DO

  RETURN
END
SUBROUTINE bc_val(n, t, xcon, ycon, xval, yval)
!
!*******************************************************************************
!
!! BC_VAL evaluates a parameterized Bezier curve.
!
!
!  Discussion:
!
!    BC_VAL(T) is the value of a vector function of the form
!
!      BC_VAL(T) = ( X(T), Y(T) )
!
!    where
!
!      X(T) = Sum (i = 0 to N) XCON(I) * BERN(I,N)(T)
!      Y(T) = Sum (i = 0 to N) YCON(I) * BERN(I,N)(T)
!
!    where BERN(I,N)(T) is the I-th Bernstein polynomial of order N
!    defined on the interval [0,1], and where XCON(*) and YCON(*)
!    are the coordinates of N+1 "control points".
!
!  Modified:
!
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the Bezier curve, which
!    must be at least 0.
!
!    Input, real T, the point at which the Bezier curve should
!    be evaluated.  The best results are obtained within the interval
!    [0,1] but T may be anywhere.
!
!    Input, real XCON(0:N), YCON(0:N), the X and Y coordinates
!    of the control points.  The Bezier curve will pass through
!    the points ( XCON(0), YCON(0) ) and ( XCON(N), YCON(N) ), but
!    generally NOT through the other control points.
!
!    Output, real XVAL, YVAL, the X and Y coordinates of the point
!    on the Bezier curve corresponding to the given T value.
!
  IMPLICIT NONE
!
  INTEGER n
!
  REAL bval(0:n)
  INTEGER i
  REAL t
  REAL xcon(0:n)
  REAL xval
  REAL ycon(0:n)
  REAL yval
!
  CALL bp01(n, bval, t)

  xval = DOT_PRODUCT(xcon(0:n), bval(0:n))
  yval = DOT_PRODUCT(ycon(0:n), bval(0:n))

  RETURN
END
FUNCTION bez_val(n, x, a, b, y)
!
!*******************************************************************************
!
!! BEZ_VAL evaluates a Bezier function at a point.
!
!
!  Discussion:
!
!    The Bezier function has the form:
!
!      BEZ(X) = Sum (i = 0 to N) Y(I) * BERN(N,I)( (X-A)/(B-A) )
!
!    where BERN(N,I)(X) is the I-th Bernstein polynomial of order N
!    defined on the interval [0,1], and Y(*) is a set of
!    coefficients.
!
!    If we define the N+1 equally spaced
!
!      X(I) = ( (N-I)*A + I*B)/ N,
!
!    for I = 0 to N, then the pairs ( X(I), Y(I) ) can be regarded as
!    "control points".
!
!  Modified:
!
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the Bezier function, which
!    must be at least 0.
!
!    Input, real X, the point at which the Bezier function should
!    be evaluated.  The best results are obtained within the interval
!    [A,B] but X may be anywhere.
!
!    Input, real A, B, the interval over which the Bezier function
!    has been defined.  This is the interval in which the control
!    points have been set up.  Note BEZ(A) = Y(0) and BEZ(B) = Y(N),
!    although BEZ will not, in general pass through the other
!    control points.  A and B must not be equal.
!
!    Input, real Y(0:N), a set of data defining the Y coordinates
!    of the control points.
!
!    Output, real BEZ_VAL, the value of the Bezier function at X.
!
  IMPLICIT NONE
!
  INTEGER n
!
  REAL a
  REAL b
  REAL bez_val
  REAL bval(0:n)
  INTEGER i
  INTEGER ncopy
  REAL x
  REAL x01
  REAL y(0:n)
!
  IF (b - a == 0.0E+00) THEN
    WRITE (*, '(a)') ' '
    WRITE (*, '(a)') 'BEZ_VAL - Fatal error!'
    WRITE (*, '(a,g14.6)') '  Null interval, A = B = ', a
    STOP
  END IF

  x01 = (x - a) / (b - a)

  ncopy = n
  CALL bp01(ncopy, bval, x01)

  bez_val = DOT_PRODUCT(y(0:n), bval(0:n))

  RETURN
END
SUBROUTINE bp01(n, bern, x)
!
!*******************************************************************************
!
!! BP01 evaluates the Bernstein basis polynomials for [01,1] at a point X.
!
!
!  Formula:
!
!    BERN(N,I,X) = [N!/(I!*(N-I)!)] * (1-X)**(N-I) * X**I
!
!  First values:
!
!    B(0,0,X) = 1
!
!    B(1,0,X) =      1-X
!    B(1,1,X) =                X
!
!    B(2,0,X) =     (1-X)**2
!    B(2,1,X) = 2 * (1-X)    * X
!    B(2,2,X) =                X**2
!
!    B(3,0,X) =     (1-X)**3
!    B(3,1,X) = 3 * (1-X)**2 * X
!    B(3,2,X) = 3 * (1-X)    * X**2
!    B(3,3,X) =                X**3
!
!    B(4,0,X) =     (1-X)**4
!    B(4,1,X) = 4 * (1-X)**3 * X
!    B(4,2,X) = 6 * (1-X)**2 * X**2
!    B(4,3,X) = 4 * (1-X)    * X**3
!    B(4,4,X) =                X**4
!
!  Special values:
!
!    B(N,I,1/2) = C(N,K) / 2**N
!
!  Modified:
!
!    12 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the degree of the Bernstein basis polynomials.
!    For any N greater than or equal to 0, there is a set of N+1 Bernstein
!    basis polynomials, each of degree N, which form a basis for
!    all polynomials on [0,1].
!
!    Output, real BERN(0:N), the values of the N+1 Bernstein basis
!    polynomials at X.
!
!    Input, real X, the point at which the polynomials are to be
!    evaluated.
!
  IMPLICIT NONE
!
  INTEGER n
!
  REAL bern(0:n)
  INTEGER i
  INTEGER j
  REAL x
!
  IF (n == 0) THEN

    bern(0) = 1.0E+00

  ELSE IF (n > 0) THEN

    bern(0) = 1.0E+00 - x
    bern(1) = x

    DO i = 2, n
      bern(i) = x * bern(i - 1)
      DO j = i - 1, 1, -1
        bern(j) = x * bern(j - 1) + (1.0E+00 - x) * bern(j)
      END DO
      bern(0) = (1.0E+00 - x) * bern(0)
    END DO

  END IF

  RETURN
END
SUBROUTINE bp_approx(n, a, b, ydata, xval, yval)
!
!*******************************************************************************
!
!! BP_APPROX evaluates the Bernstein polynomial for F(X) on [A,B].
!
!
!  Formula:
!
!    BERN(F)(X) = SUM ( 0 <= I <= N ) F(X(I)) * B_BASE(I,X)
!
!    where
!
!      X(I) = ( ( N - I ) * A + I * B ) / N
!      B_BASE(I,X) is the value of the I-th Bernstein basis polynomial at X.
!
!  Discussion:
!
!    The Bernstein polynomial BERN(F) for F(X) is an approximant, not an
!    interpolant; in other words, its value is not guaranteed to equal
!    that of F at any particular point.  However, for a fixed interval
!    [A,B], if we let N increase, the Bernstein polynomial converges
!    uniformly to F everywhere in [A,B], provided only that F is continuous.
!    Even if F is not continuous, but is bounded, the polynomial converges
!    pointwise to F(X) at all points of continuity.  On the other hand,
!    the convergence is quite slow compared to other interpolation
!    and approximation schemes.
!
!  Modified:
!
!    10 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the degree of the Bernstein polynomial to be used.
!
!    Input, real A, B, the endpoints of the interval on which the
!    approximant is based.  A and B should not be equal.
!
!    Input, real YDATA(0:N), the data values at N+1 equally spaced points
!    in [A,B].  If N = 0, then the evaluation point should be 0.5 * ( A + B).
!    Otherwise, evaluation point I should be ( (N-I)*A + I*B ) / N ).
!
!    Input, real XVAL, the point at which the Bernstein polynomial
!    approximant is to be evaluated.  XVAL does not have to lie in the
!    interval [A,B].
!
!    Output, real YVAL, the value of the Bernstein polynomial approximant
!    for F, based in [A,B], evaluated at XVAL.
!
  IMPLICIT NONE
!
  INTEGER n
!
  REAL a
  REAL b
  REAL bvec(0:n)
  INTEGER i
  REAL xval
  REAL ydata(0:n)
  REAL yval
!
!  Evaluate the Bernstein basis polynomials at XVAL.
!
  CALL bpab(n, bvec, xval, a, b)
!
!  Now compute the sum of YDATA(I) * BVEC(I).
!
  yval = DOT_PRODUCT(ydata(0:n), bvec(0:n))

  RETURN
END
SUBROUTINE bpab(n, bern, x, a, b)
!
!*******************************************************************************
!
!! BPAB evaluates the Bernstein basis polynomials for [A,B] at a point X.
!
!
!  Formula:
!
!    BERN(N,I,X) = [N!/(I!*(N-I)!)] * (B-X)**(N-I) * (X-A)**I / (B-A)**N
!
!  First values:
!
!    B(0,0,X) =   1
!
!    B(1,0,X) = (      B-X                ) / (B-A)
!    B(1,1,X) = (                 X-A     ) / (B-A)
!
!    B(2,0,X) = (     (B-X)**2            ) / (B-A)**2
!    B(2,1,X) = ( 2 * (B-X)    * (X-A)    ) / (B-A)**2
!    B(2,2,X) = (                (X-A)**2 ) / (B-A)**2
!
!    B(3,0,X) = (     (B-X)**3            ) / (B-A)**3
!    B(3,1,X) = ( 3 * (B-X)**2 * (X-A)    ) / (B-A)**3
!    B(3,2,X) = ( 3 * (B-X)    * (X-A)**2 ) / (B-A)**3
!    B(3,3,X) = (                (X-A)**3 ) / (B-A)**3
!
!    B(4,0,X) = (     (B-X)**4            ) / (B-A)**4
!    B(4,1,X) = ( 4 * (B-X)**3 * (X-A)    ) / (B-A)**4
!    B(4,2,X) = ( 6 * (B-X)**2 * (X-A)**2 ) / (B-A)**4
!    B(4,3,X) = ( 4 * (B-X)    * (X-A)**3 ) / (B-A)**4
!    B(4,4,X) = (                (X-A)**4 ) / (B-A)**4
!
!  Modified:
!
!    12 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the degree of the Bernstein basis polynomials.
!    For any N greater than or equal to 0, there is a set of N+1
!    Bernstein basis polynomials, each of degree N, which form a basis
!    for polynomials on [A,B].
!
!    Output, real BERN(0:N), the values of the N+1 Bernstein basis
!    polynomials at X.
!
!    Input, real X, the point at which the polynomials are to be
!    evaluated.  X need not lie in the interval [A,B].
!
!    Input, real A, B, the endpoints of the interval on which the
!    polynomials are to be based.  A and B should not be equal.
!
  IMPLICIT NONE
!
  INTEGER n
!
  REAL a
  REAL b
  REAL bern(0:n)
  INTEGER i
  INTEGER j
  REAL x
!
  IF (b == a) THEN
    WRITE (*, '(a)') ' '
    WRITE (*, '(a)') 'BPAB - Fatal error!'
    WRITE (*, '(a,g14.6)') '  A = B = ', a
    STOP
  END IF

  IF (n == 0) THEN

    bern(0) = 1.0E+00

  ELSE IF (n > 0) THEN

    bern(0) = (b - x) / (b - a)
    bern(1) = (x - a) / (b - a)

    DO i = 2, n
      bern(i) = (x - a) * bern(i - 1) / (b - a)
      DO j = i - 1, 1, -1
        bern(j) = ((b - x) * bern(j) + (x - a) * bern(j - 1)) / (b - a)
      END DO
      bern(0) = (b - x) * bern(0) / (b - a)
    END DO

  END IF

  RETURN
END
SUBROUTINE data_to_dif(diftab, ntab, xtab, ytab)
!
!*******************************************************************************
!
!! DATA_TO_DIF sets up a divided difference table from raw data.
!
!
!  Discussion:
!
!    Space can be saved by using a single array for both the DIFTAB and
!    YTAB dummy parameters.  In that case, the divided difference table will
!    overwrite the Y data without interfering with the computation.
!
!  Modified:
!
!    11 April 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real DIFTAB(NTAB), the divided difference coefficients
!    corresponding to the input (XTAB,YTAB).
!
!    Input, integer NTAB, the number of pairs of points
!    (XTAB(I),YTAB(I)) which are to be used as data.  The
!    number of entries to be used in DIFTAB, XTAB and YTAB.
!
!    Input, real XTAB(NTAB), the X values at which data was taken.
!    These values must be distinct.
!
!    Input, real YTAB(NTAB), the corresponding Y values.
!
  IMPLICIT NONE
!
  INTEGER ntab
!
  REAL diftab(ntab)
  INTEGER i
  INTEGER j
  LOGICAL rvec_distinct
  REAL xtab(ntab)
  REAL ytab(ntab)
!
  IF (.NOT. rvec_distinct(ntab, xtab)) THEN
    WRITE (*, '(a)') ' '
    WRITE (*, '(a)') 'DATA_TO_DIF - Fatal error!'
    WRITE (*, '(a)') '  Two entries of XTAB are equal!'
    RETURN
  END IF
!
!  Copy the data values into DIFTAB.
!
  diftab(1:ntab) = ytab(1:ntab)
!
!  Compute the divided differences.
!
  DO i = 2, ntab
    DO j = ntab, i, -1

      diftab(j) = (diftab(j) - diftab(j - 1)) / (xtab(j) - xtab(j + 1 - i))

    END DO
  END DO

  RETURN
END
SUBROUTINE dif_val(diftab, ntab, xtab, xval, yval)
!
!*******************************************************************************
!
!! DIF_VAL evaluates a divided difference polynomial at a point.
!
!
!  Modified:
!
!    11 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real DIFTAB(NTAB), the divided difference polynomial coefficients.
!
!    Input, integer NTAB, the number of divided difference
!    coefficients in DIFTAB, and the number of points XTAB.
!
!    Input, real XTAB(NTAB), the X values upon which the
!    divided difference polynomial is based.
!
!    Input, real XVAL, a value of X at which the polynomial
!    is to be evaluated.
!
!    Output, real YVAL, the value of the polynomial at XVAL.
!
  IMPLICIT NONE
!
  INTEGER ntab
!
  REAL diftab(ntab)
  INTEGER i
  REAL xtab(ntab)
  REAL xval
  REAL yval
!
  yval = diftab(ntab)
  DO i = 1, ntab - 1
    yval = diftab(ntab - i) + (xval - xtab(ntab - i)) * yval
  END DO

  RETURN
END
SUBROUTINE least_set(ntab, xtab, ytab, ndeg, ptab, array, eps, ierror)
!
!*******************************************************************************
!
!! LEAST_SET constructs the least squares polynomial approximation to data.
!
!
!  Discussion:
!
!    The routine LEAST_EVAL must be used to evaluate the approximation at a
!    point.
!
!  Modified:
!
!    20 November 2000
!
!  Parameters:
!
!    Input, integer NTAB, the number of data points.
!
!    Input, real XTAB(NTAB), the X data.  The values in XTAB
!    should be distinct, and in increasing order.
!
!    Input, real YTAB(NTAB), the Y data values corresponding
!    to the X data in XTAB.
!
!    Input, integer NDEG, the degree of the polynomial which the
!    program is to use.  NDEG must be at least 1, and less than or
!    equal to NTAB-1.
!
!    Output, real PTAB(NTAB).  PTAB(I) is the value of the
!    least squares polynomial at the point XTAB(I).
!
!    Output, real ARRAY(2*NTAB+3*NDEG), an array containing data about
!    the polynomial.
!
!    Output, real EPS, the root-mean-square discrepancy of the
!    polynomial fit.
!
!    Output, integer IERROR, error flag.
!    zero, no error occurred;
!    nonzero, an error occurred, and the polynomial could not be computed.
!
  IMPLICIT NONE
!
  INTEGER ndeg
  INTEGER ntab
!
  REAL array(2 * ntab + 3 * ndeg)
  REAL eps
  REAL error
  INTEGER i
  INTEGER i0l1
  INTEGER i1l1
  INTEGER ierror
  INTEGER it
  INTEGER k
  INTEGER mdeg
  REAL ptab(ntab)
  REAL rn0
  REAL rn1
  REAL s
  REAL sum2
  REAL xtab(ntab)
  REAL y_sum
  REAL ytab(ntab)
!
  ierror = 0
!
!  Check NDEG.
!
  IF (ndeg < 1) THEN
    ierror = 1
    WRITE (*, '(a)') ' '
    WRITE (*, '(a)') 'LEAST_SET - Fatal error!'
    WRITE (*, '(a)') '  NDEG < 1.'
    RETURN
  ELSE IF (ndeg >= ntab) THEN
    ierror = 1
    WRITE (*, '(a)') ' '
    WRITE (*, '(a)') 'LEAST_SET - Fatal error!'
    WRITE (*, '(a)') '  NDEG >= NTAB.'
    RETURN
  END IF
!
!  Check that the abscissas are strictly increasing.
!
  DO i = 1, ntab - 1
    IF (xtab(i) >= xtab(i + 1)) THEN
      ierror = 1
      WRITE (*, '(a)') ' '
      WRITE (*, '(a)') 'LEAST_SET - Fatal error!'
      WRITE (*, '(a)') '  XTAB must be strictly increasing, but'
      WRITE (*, '(a,i6,a,g14.6)') '  XTAB(', i, ') = ', xtab(i)
      WRITE (*, '(a,i6,a,g14.6)') '  XTAB(', i + 1, ') = ', xtab(i + 1)
      RETURN
    END IF
  END DO

  i0l1 = 3 * ndeg
  i1l1 = 3 * ndeg + ntab

  y_sum = SUM(ytab)
  rn0 = ntab
  array(2 * ndeg) = y_sum / REAL(ntab)

  ptab(1:ntab) = y_sum / REAL(ntab)

  error = 0.0E+00
  DO i = 1, ntab
    error = error + (y_sum / REAL(ntab) - ytab(i))**2
  END DO

  IF (ndeg == 0) THEN
    eps = SQRT(error / REAL(ntab))
    RETURN
  END IF

  array(1) = SUM(xtab) / REAL(ntab)

  s = 0.0E+00
  sum2 = 0.0E+00
  DO i = 1, ntab
    array(i1l1 + i) = xtab(i) - array(1)
    s = s + array(i1l1 + i)**2
    sum2 = sum2 + array(i1l1 + i) * (ytab(i) - ptab(i))
  END DO

  rn1 = s
  array(2 * ndeg + 1) = sum2 / s

  DO i = 1, ntab
    ptab(i) = ptab(i) + sum2 * array(i1l1 + i) / s
  END DO

  error = SUM((ptab(1:ntab) - ytab(1:ntab))**2)

  IF (ndeg == 1) THEN
    eps = SQRT(error / REAL(ntab))
    RETURN
  END IF

  DO i = 1, ntab
    array(3 * ndeg + i) = 1.0E+00
  END DO

  mdeg = 2
  k = 2

  DO

    array(ndeg - 1 + k) = rn1 / rn0

    sum2 = 0.0E+00
    DO i = 1, ntab
      sum2 = sum2 + xtab(i) * array(i1l1 + i)**2
    END DO

    array(k) = sum2 / rn1

    s = 0.0E+00
    sum2 = 0.0E+00
    DO i = 1, ntab
      array(i0l1 + i) = (xtab(i) - array(k)) * array(i1l1 + i) &
                        - array(ndeg - 1 + k) * array(i0l1 + i)
      s = s + array(i0l1 + i)**2
      sum2 = sum2 + array(i0l1 + i) * (ytab(i) - ptab(i))
    END DO

    rn0 = rn1
    rn1 = s
    it = i0l1
    i0l1 = i1l1
    i1l1 = it
    array(2 * ndeg + k) = sum2 / rn1

    DO i = 1, ntab
      ptab(i) = ptab(i) + sum2 * array(i1l1 + i) / rn1
    END DO

    error = SUM((ptab(1:ntab) - ytab(1:ntab))**2)

    IF (mdeg >= ndeg) THEN
      EXIT
    END IF

    mdeg = mdeg + 1
    k = k + 1

  END DO

  eps = SQRT(error / REAL(ntab))

  RETURN
END
SUBROUTINE least_val(x, ndeg, array, VALUE)
!
!*******************************************************************************
!
!! LEAST_VAL evaluates a least squares polynomial defined by LEAST_SET.
!
!
!  Modified:
!
!    01 March 1999
!
!  Parameters:
!
!    Input, real X, the point at which the polynomial is to be evaluated.
!
!    Input, integer NDEG, the degree of the polynomial fit used.
!    This is the value of NDEG as returned from LEAST_SET.
!
!    Input, real ARRAY(*), an array of a certain dimension.
!    See LEAST_SET for details on the size of ARRAY.
!    ARRAY contains information about the polynomial, as set up by LEAST_SET.
!
!    Output, real VALUE, the value of the polynomial at X.
!
  IMPLICIT NONE
!
  REAL array(*)
  REAL dk
  REAL dkp1
  REAL dkp2
  INTEGER k
  INTEGER l
  INTEGER ndeg
  REAL VALUE
  REAL x
!
  IF (ndeg <= 0) THEN

    VALUE = array(2 * ndeg)

  ELSE IF (ndeg == 1) THEN

    VALUE = array(2 * ndeg) + array(2 * ndeg + 1) * (x - array(1))

  ELSE

    dkp2 = array(3 * ndeg)
    dkp1 = array(3 * ndeg - 1) + (x - array(ndeg)) * dkp2

    DO l = 1, ndeg - 2

      k = ndeg - 1 - l

      dk = array(2 * ndeg + k) + (x - array(k + 1)) * dkp1 - array(ndeg + 1 + k) * dkp2

      dkp2 = dkp1

      dkp1 = dk

    END DO

    VALUE = array(2 * ndeg) + (x - array(1)) * dkp1 - array(ndeg + 1) * dkp2

  END IF

  RETURN
END
SUBROUTINE parabola_val2(ndim, ndata, tdata, ydata, left, tval, yval)
!
!*******************************************************************************
!
!! PARABOLA_VAL2 evaluates a parabolic interpolant through tabular data.
!
!
!  Discussion:
!
!    This routine is a utility routine used by OVERHAUSER_SPLINE_VAL.
!    It constructs the parabolic interpolant through the data in
!    3 consecutive entries of a table and evaluates this interpolant
!    at a given abscissa value.
!
!  Modified:
!
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NDIM, the dimension of a single data point.
!    NDIM must be at least 1.
!
!    Input, integer NDATA, the number of data points.
!    NDATA must be at least 3.
!
!    Input, real TDATA(NDATA), the abscissas of the data points.  The
!    values in TDATA must be in strictly ascending order.
!
!    Input, real YDATA(NDIM,NDATA), the data points corresponding to
!    the abscissas.
!
!    Input, integer LEFT, the location of the first of the three
!    consecutive data points through which the parabolic interpolant
!    must pass.  1 <= LEFT <= NDATA - 2.
!
!    Input, real TVAL, the value of T at which the parabolic interpolant
!    is to be evaluated.  Normally, TDATA(1) <= TVAL <= T(NDATA), and
!    the data will be interpolated.  For TVAL outside this range,
!    extrapolation will be used.
!
!    Output, real YVAL(NDIM), the value of the parabolic interpolant at TVAL.
!
  IMPLICIT NONE
!
  INTEGER ndata
  INTEGER ndim
!
  REAL dif1
  REAL dif2
  INTEGER i
  INTEGER left
  REAL t1
  REAL t2
  REAL t3
  REAL tval
  REAL tdata(ndata)
  REAL ydata(ndim, ndata)
  REAL y1
  REAL y2
  REAL y3
  REAL yval(ndim)
!
!  Check.
!
  IF (left < 1 .OR. left > ndata - 2) THEN
    WRITE (*, '(a)') ' '
    WRITE (*, '(a)') 'PARABOLA_VAL2 - Fatal error!'
    WRITE (*, '(a)') '  LEFT < 1 or LEFT > NDATA-2.'
    STOP
  END IF

  IF (ndim < 1) THEN
    WRITE (*, '(a)') ' '
    WRITE (*, '(a)') 'PARABOLA_VAL2 - Fatal error!'
    WRITE (*, '(a)') '  NDIM < 1.'
    STOP
  END IF
!
!  Copy out the three abscissas.
!
  t1 = tdata(left)
  t2 = tdata(left + 1)
  t3 = tdata(left + 2)

  IF (t1 >= t2 .OR. t2 >= t3) THEN
    WRITE (*, '(a)') ' '
    WRITE (*, '(a)') 'PARABOLA_VAL2 - Fatal error!'
    WRITE (*, '(a)') '  T1 >= T2 or T2 >= T3.'
    STOP
  END IF
!
!  Construct and evaluate a parabolic interpolant for the data
!  in each dimension.
!
  DO i = 1, ndim

    y1 = ydata(i, left)
    y2 = ydata(i, left + 1)
    y3 = ydata(i, left + 2)

    dif1 = (y2 - y1) / (t2 - t1)
    dif2 = ((y3 - y1) / (t3 - t1) &
            - (y2 - y1) / (t2 - t1)) / (t3 - t2)

    yval(i) = y1 + (tval - t1) * (dif1 + (tval - t2) * dif2)

  END DO

  RETURN
END
SUBROUTINE r_random(rlo, rhi, r)
!
!*******************************************************************************
!
!! R_RANDOM returns a random real in a given range.
!
!
!  Modified:
!
!    06 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real RLO, RHI, the minimum and maximum values.
!
!    Output, real R, the randomly chosen value.
!
  IMPLICIT NONE
!
  REAL r
  REAL rhi
  REAL rlo
  REAL t
!
!  Pick T, a random number in (0,1).
!
  CALL RANDOM_NUMBER(harvest=t)
!
!  Set R in ( RLO, RHI ).
!
  r = (1.0E+00 - t) * rlo + t * rhi

  RETURN
END
SUBROUTINE r_swap(x, y)
!
!*******************************************************************************
!
!! R_SWAP swaps two real values.
!
!
!  Modified:
!
!    01 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real X, Y.  On output, the values of X and
!    Y have been interchanged.
!
  IMPLICIT NONE
!
  REAL x
  REAL y
  REAL z
!
  z = x
  x = y
  y = z

  RETURN
END
SUBROUTINE rvec_bracket(n, x, xval, left, right)
!
!*******************************************************************************
!
!! RVEC_BRACKET searches a sorted array for successive brackets of a value.
!
!
!  Discussion:
!
!    If the values in the vector are thought of as defining intervals
!    on the real line, then this routine searches for the interval
!    nearest to or containing the given value.
!
!  Modified:
!
!    06 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, length of input array.
!
!    Input, real X(N), an array sorted into ascending order.
!
!    Input, real XVAL, a value to be bracketed.
!
!    Output, integer LEFT, RIGHT, the results of the search.
!    Either:
!      XVAL < X(1), when LEFT = 1, RIGHT = 2;
!      XVAL > X(N), when LEFT = N-1, RIGHT = N;
!    or
!      X(LEFT) <= XVAL <= X(RIGHT).
!
  IMPLICIT NONE
!
  INTEGER n
!
  INTEGER i
  INTEGER left
  INTEGER right
  REAL x(n)
  REAL xval
!
  DO i = 2, n - 1

    IF (xval < x(i)) THEN
      left = i - 1
      right = i
      RETURN
    END IF

  END DO

  left = n - 1
  right = n

  RETURN
END
SUBROUTINE rvec_bracket3(n, t, tval, left)
!
!*******************************************************************************
!
!! RVEC_BRACKET3 finds the interval containing or nearest a given value.
!
!
!  Discussion:
!
!    The routine always returns the index LEFT of the sorted array
!    T with the property that either
!    *  T is contained in the interval [ T(LEFT), T(LEFT+1) ], or
!    *  T < T(LEFT) = T(1), or
!    *  T > T(LEFT+1) = T(N).
!
!    The routine is useful for interpolation problems, where
!    the abscissa must be located within an interval of data
!    abscissas for interpolation, or the "nearest" interval
!    to the (extreme) abscissa must be found so that extrapolation
!    can be carried out.
!
!  Modified:
!
!    05 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, length of the input array.
!
!    Input, real T(N), an array sorted into ascending order.
!
!    Input, real TVAL, a value to be bracketed by entries of T.
!
!    Input/output, integer LEFT.
!
!    On input, if 1 <= LEFT <= N-1, LEFT is taken as a suggestion for the
!    interval [ T(LEFT), T(LEFT+1) ] in which TVAL lies.  This interval
!    is searched first, followed by the appropriate interval to the left
!    or right.  After that, a binary search is used.
!
!    On output, LEFT is set so that the interval [ T(LEFT), T(LEFT+1) ]
!    is the closest to TVAL; it either contains TVAL, or else TVAL
!    lies outside the interval [ T(1), T(N) ].
!
  IMPLICIT NONE
!
  INTEGER n
!
  INTEGER high
  INTEGER left
  INTEGER low
  INTEGER mid
  REAL t(n)
  REAL tval
!
!  Check the input data.
!
  IF (n < 2) THEN
    WRITE (*, '(a)') ' '
    WRITE (*, '(a)') 'RVEC_BRACKET3 - Fatal error!'
    WRITE (*, '(a)') '  N must be at least 2.'
    STOP
  END IF
!
!  If LEFT is not between 1 and N-1, set it to the middle value.
!
  IF (left < 1 .OR. left > n - 1) THEN
    left = (n + 1) / 2
  END IF
!
!  CASE 1: TVAL < T(LEFT):
!  Search for TVAL in [T(I), T(I+1)] for intervals I = 1 to LEFT-1.
!
  IF (tval < t(left)) THEN

    IF (left == 1) THEN
      RETURN
    ELSE IF (left == 2) THEN
      left = 1
      RETURN
    ELSE IF (tval >= t(left - 1)) THEN
      left = left - 1
      RETURN
    ELSE IF (tval <= t(2)) THEN
      left = 1
      RETURN
    END IF
!
!  ...Binary search for TVAL in [T(I), T(I+1)] for intervals I = 2 to LEFT-2.
!
    low = 2
    high = left - 2

    DO

      IF (low == high) THEN
        left = low
        RETURN
      END IF

      mid = (low + high + 1) / 2

      IF (tval >= t(mid)) THEN
        low = mid
      ELSE
        high = mid - 1
      END IF

    END DO
!
!  CASE2: T(LEFT+1) < TVAL:
!  Search for TVAL in {T(I),T(I+1)] for intervals I = LEFT+1 to N-1.
!
  ELSE IF (tval > t(left + 1)) THEN

    IF (left == n - 1) THEN
      RETURN
    ELSE IF (left == n - 2) THEN
      left = left + 1
      RETURN
    ELSE IF (tval <= t(left + 2)) THEN
      left = left + 1
      RETURN
    ELSE IF (tval >= t(n - 1)) THEN
      left = n - 1
      RETURN
    END IF
!
!  ...Binary search for TVAL in [T(I), T(I+1)] for intervals I = LEFT+2 to N-2.
!
    low = left + 2
    high = n - 2

    DO

      IF (low == high) THEN
        left = low
        RETURN
      END IF

      mid = (low + high + 1) / 2

      IF (tval >= t(mid)) THEN
        low = mid
      ELSE
        high = mid - 1
      END IF

    END DO
!
!  CASE3: T(LEFT) <= TVAL <= T(LEFT+1):
!  T is in [T(LEFT), T(LEFT+1)], as the user said it might be.
!
  ELSE

  END IF

  RETURN
END
FUNCTION rvec_distinct(n, x)
!
!*******************************************************************************
!
!! RVEC_DISTINCT is true if the entries in a real vector are distinct.
!
!
!  Modified:
!
!    16 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, real X(N), the vector to be checked.
!
!    Output, logical RVEC_DISTINCT is .TRUE. if all N elements of X
!    are distinct.
!
  IMPLICIT NONE
!
  INTEGER n
!
  INTEGER i
  INTEGER j
  LOGICAL rvec_distinct
  REAL x(n)
!
  rvec_distinct = .FALSE.

  DO i = 2, n
    DO j = 1, i - 1
      IF (x(i) == x(j)) THEN
        RETURN
      END IF
    END DO
  END DO

  rvec_distinct = .TRUE.

  RETURN
END
SUBROUTINE rvec_even(alo, ahi, n, a)
!
!*******************************************************************************
!
!! RVEC_EVEN returns N real values, evenly spaced between ALO and AHI.
!
!
!  Modified:
!
!    31 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ALO, AHI, the low and high values.
!
!    Input, integer N, the number of values.
!
!    Output, real A(N), N evenly spaced values.
!    Normally, A(1) = ALO and A(N) = AHI.
!    However, if N = 1, then A(1) = 0.5*(ALO+AHI).
!
  IMPLICIT NONE
!
  INTEGER n
!
  REAL a(n)
  REAL ahi
  REAL alo
  INTEGER i
!
  IF (n == 1) THEN

    a(1) = 0.5E+00 * (alo + ahi)

  ELSE

    DO i = 1, n
      a(i) = (REAL(n - i) * alo + REAL(i - 1) * ahi) / REAL(n - 1)
    END DO

  END IF

  RETURN
END
SUBROUTINE rvec_order_type(n, a, order)
!
!*******************************************************************************
!
!! RVEC_ORDER_TYPE determines if a real array is (non)strictly ascending/descending.
!
!
!  Modified:
!
!    20 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries of the array.
!
!    Input, real A(N), the array to be checked.
!
!    Output, integer ORDER, order indicator:
!    -1, no discernable order;
!    0, all entries are equal;
!    1, ascending order;
!    2, strictly ascending order;
!    3, descending order;
!    4, strictly descending order.
!
  IMPLICIT NONE
!
  INTEGER n
!
  REAL a(n)
  INTEGER i
  INTEGER order
!
!  Search for the first value not equal to A(1).
!
  i = 1

  DO

    i = i + 1

    IF (i > n) THEN
      order = 0
      RETURN
    END IF

    IF (a(i) > a(1)) THEN

      IF (i == 2) THEN
        order = 2
      ELSE
        order = 1
      END IF

      EXIT

    ELSE IF (a(i) < a(1)) THEN

      IF (i == 2) THEN
        order = 4
      ELSE
        order = 3
      END IF

      EXIT

    END IF

  END DO
!
!  Now we have a "direction".  Examine subsequent entries.
!
  DO

    i = i + 1
    IF (i > n) THEN
      EXIT
    END IF

    IF (order == 1) THEN

      IF (a(i) < a(i - 1)) THEN
        order = -1
        EXIT
      END IF

    ELSE IF (order == 2) THEN

      IF (a(i) < a(i - 1)) THEN
        order = -1
        EXIT
      ELSE IF (a(i) == a(i - 1)) THEN
        order = 1
      END IF

    ELSE IF (order == 3) THEN

      IF (a(i) > a(i - 1)) THEN
        order = -1
        EXIT
      END IF

    ELSE IF (order == 4) THEN

      IF (a(i) > a(i - 1)) THEN
        order = -1
        EXIT
      ELSE IF (a(i) == a(i - 1)) THEN
        order = 3
      END IF

    END IF

  END DO

  RETURN
END
SUBROUTINE rvec_print(n, a, title)
!
!*******************************************************************************
!
!! RVEC_PRINT prints a real vector.
!
!
!  Modified:
!
!    16 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components of the vector.
!
!    Input, real A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  IMPLICIT NONE
!
  INTEGER n
!
  REAL a(n)
  INTEGER i
  CHARACTER(len=*) title
!
  IF (title /= ' ') THEN
    WRITE (*, '(a)') ' '
    WRITE (*, '(a)') TRIM(title)
  END IF

  WRITE (*, '(a)') ' '
  DO i = 1, n
    WRITE (*, '(i6,g14.6)') i, a(i)
  END DO

  RETURN
END
SUBROUTINE rvec_random(alo, ahi, n, a)
!
!*******************************************************************************
!
!! RVEC_RANDOM returns a random real vector in a given range.
!
!
!  Modified:
!
!    04 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ALO, AHI, the range allowed for the entries.
!
!    Input, integer N, the number of entries in the vector.
!
!    Output, real A(N), the vector of randomly chosen values.
!
  IMPLICIT NONE
!
  INTEGER n
!
  REAL a(n)
  REAL ahi
  REAL alo
  INTEGER i
!
  DO i = 1, n
    CALL r_random(alo, ahi, a(i))
  END DO

  RETURN
END
SUBROUTINE rvec_sort_bubble_a(n, a)
!
!*******************************************************************************
!
!! RVEC_SORT_BUBBLE_A ascending sorts a real array using bubble sort.
!
!
!  Discussion:
!
!    Bubble sort is simple to program, but inefficient.  It should not
!    be used for large arrays.
!
!  Modified:
!
!    01 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the array.
!
!    Input/output, real A(N).
!    On input, an unsorted array.
!    On output, the array has been sorted.
!
  IMPLICIT NONE
!
  INTEGER n
!
  REAL a(n)
  INTEGER i
  INTEGER j
!
  DO i = 1, n - 1
    DO j = i + 1, n
      IF (a(i) > a(j)) THEN
        CALL r_swap(a(i), a(j))
      END IF
    END DO
  END DO

  RETURN
END
SUBROUTINE s3_fs(a1, a2, a3, n, b, x)
!
!*******************************************************************************
!
!! S3_FS factors and solves a tridiagonal linear system.
!
!
!  Note:
!
!    This algorithm requires that each diagonal entry be nonzero.
!
!  Modified:
!
!    05 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real A1(2:N), A2(1:N), A3(1:N-1).
!    On input, the nonzero diagonals of the linear system.
!    On output, the data in these vectors has been overwritten
!    by factorization information.
!
!    Input, integer N, the order of the linear system.
!
!    Input/output, real B(N).
!    On input, B contains the right hand side of the linear system.
!    On output, B has been overwritten by factorization information.
!
!    Output, real X(N), the solution of the linear system.
!
  IMPLICIT NONE
!
  INTEGER n
!
  REAL a1(2:n)
  REAL a2(1:n)
  REAL a3(1:n - 1)
  REAL b(n)
  INTEGER i
  REAL x(n)
  REAL xmult
!
!  The diagonal entries can't be zero.
!
  DO i = 1, n
    IF (a2(i) == 0.0E+00) THEN
      WRITE (*, '(a)') ' '
      WRITE (*, '(a)') 'S3_FS - Fatal error!'
      WRITE (*, '(a,i6,a)') '  A2(', i, ') = 0.'
      RETURN
    END IF
  END DO

  DO i = 2, n - 1

    xmult = a1(i) / a2(i - 1)
    a2(i) = a2(i) - xmult * a3(i - 1)

    b(i) = b(i) - xmult * b(i - 1)

  END DO

  xmult = a1(n) / a2(n - 1)
  a2(n) = a2(n) - xmult * a3(n - 1)

  x(n) = (b(n) - xmult * b(n - 1)) / a2(n)
  DO i = n - 1, 1, -1
    x(i) = (b(i) - a3(i) * x(i + 1)) / a2(i)
  END DO

  RETURN
END
SUBROUTINE sgtsl(n, c, d, e, b, info)
!
!*******************************************************************************
!
!! SGTSL solves a general tridiagonal linear system.
!
!
!  Reference:
!
!    Dongarra, Moler, Bunch and Stewart,
!    LINPACK User's Guide,
!    SIAM, (Society for Industrial and Applied Mathematics),
!    3600 University City Science Center,
!    Philadelphia, PA, 19104-2688.
!    ISBN 0-89871-172-X
!
!  Modified:
!
!    31 October 2001
!
!  Parameters:
!
!    Input, integer N, the order of the tridiagonal matrix.
!
!    Input/output, real C(N), contains the subdiagonal of the tridiagonal
!    matrix in entries C(2:N).  On output, C is destroyed.
!
!    Input/output, real D(N).  On input, the diagonal of the matrix.
!    On output, D is destroyed.
!
!    Input/output, real E(N), contains the superdiagonal of the tridiagonal
!    matrix in entries E(1:N-1).  On output E is destroyed.
!
!    Input/output, real B(N).  On input, the right hand side.  On output,
!    the solution.
!
!    Output, integer INFO, error flag.
!    0, normal value.
!    K, the K-th element of the diagonal becomes exactly zero.  The
!       subroutine returns if this error condition is detected.
!
  IMPLICIT NONE
!
  INTEGER n
!
  REAL b(n)
  REAL c(n)
  REAL d(n)
  REAL e(n)
  INTEGER info
  INTEGER k
  REAL t
!
  info = 0
  c(1) = d(1)

  IF (n >= 2) THEN

    d(1) = e(1)
    e(1) = 0.0E+00
    e(n) = 0.0E+00

    DO k = 1, n - 1
!
!  Find the larger of the two rows, and interchange if necessary.
!
      IF (ABS(c(k + 1)) >= ABS(c(k))) THEN
        CALL r_swap(c(k), c(k + 1))
        CALL r_swap(d(k), d(k + 1))
        CALL r_swap(e(k), e(k + 1))
        CALL r_swap(b(k), b(k + 1))
      END IF
!
!  Fail if no nonzero pivot could be found.
!
      IF (c(k) == 0.0E+00) THEN
        info = k
        RETURN
      END IF
!
!  Zero elements.
!
      t = -c(k + 1) / c(k)
      c(k + 1) = d(k + 1) + t * d(k)
      d(k + 1) = e(k + 1) + t * e(k)
      e(k + 1) = 0.0E+00
      b(k + 1) = b(k + 1) + t * b(k)

    END DO

  END IF

  IF (c(n) == 0.0E+00) THEN
    info = n
    RETURN
  END IF
!
!  Back solve.
!
  b(n) = b(n) / c(n)

  IF (n > 1) THEN

    b(n - 1) = (b(n - 1) - d(n - 1) * b(n)) / c(n - 1)

    DO k = n - 2, 1, -1
      b(k) = (b(k) - d(k) * b(k + 1) - e(k) * b(k + 2)) / c(k)
    END DO

  END IF

  RETURN
END
SUBROUTINE spline_b_val(ndata, tdata, ydata, tval, yval)
!
!*******************************************************************************
!
!! SPLINE_B_VAL evaluates a cubic B spline approximant.
!
!
!  Discussion:
!
!    The cubic B spline will approximate the data, but is not
!    designed to interpolate it.
!
!    In effect, two "phantom" data values are appended to the data,
!    so that the spline will interpolate the first and last data values.
!
!  Modified:
!
!    07 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NDATA, the number of data values.
!
!    Input, real TDATA(NDATA), the abscissas of the data.
!
!    Input, real YDATA(NDATA), the data values.
!
!    Input, real TVAL, a point at which the spline is to be evaluated.
!
!    Output, real YVAL, the value of the function at TVAL.
!
  IMPLICIT NONE
!
  INTEGER ndata
!
  REAL bval
  INTEGER left
  INTEGER right
  REAL tdata(ndata)
  REAL tval
  REAL u
  REAL ydata(ndata)
  REAL yval
!
!  Find the nearest interval [ TDATA(LEFT), TDATA(RIGHT) ] to TVAL.
!
  CALL rvec_bracket(ndata, tdata, tval, left, right)
!
!  Evaluate the 5 nonzero B spline basis functions in the interval,
!  weighted by their corresponding data values.
!
  u = (tval - tdata(left)) / (tdata(right) - tdata(left))
  yval = 0.0E+00
!
!  B function associated with node LEFT - 1, (or "phantom node"),
!  evaluated in its 4th interval.
!
  bval = (1.0E+00 - 3.0E+00 * u + 3.0E+00 * u**2 - u**3) / 6.0E+00
  IF (left - 1 > 0) THEN
    yval = yval + ydata(left - 1) * bval
  ELSE
    yval = yval + (2.0E+00 * ydata(1) - ydata(2)) * bval
  END IF
!
!  B function associated with node LEFT,
!  evaluated in its third interval.
!
  bval = (4.0E+00 - 6.0E+00 * u**2 + 3.0E+00 * u**3) / 6.0E+00
  yval = yval + ydata(left) * bval
!
!  B function associated with node RIGHT,
!  evaluated in its second interval.
!
  bval = (1.0E+00 + 3.0E+00 * u + 3.0E+00 * u**2 - 3.0E+00 * u**3) / 6.0E+00
  yval = yval + ydata(right) * bval
!
!  B function associated with node RIGHT+1, (or "phantom node"),
!  evaluated in its first interval.
!
  bval = u**3 / 6.0E+00
  IF (right + 1 <= ndata) THEN
    yval = yval + ydata(right + 1) * bval
  ELSE
    yval = yval + (2.0E+00 * ydata(ndata) - ydata(ndata - 1)) * bval
  END IF

  RETURN
END
SUBROUTINE spline_beta_val(beta1, beta2, ndata, tdata, ydata, tval, yval)
!
!*******************************************************************************
!
!! SPLINE_BETA_VAL evaluates a cubic beta spline approximant.
!
!
!  Discussion:
!
!    The cubic beta spline will approximate the data, but is not
!    designed to interpolate it.
!
!    If BETA1 = 1 and BETA2 = 0, the cubic beta spline will be the
!    same as the cubic B spline approximant.
!
!    With BETA1 = 1 and BETA2 large, the beta spline becomes more like
!    a linear spline.
!
!    In effect, two "phantom" data values are appended to the data,
!    so that the spline will interpolate the first and last data values.
!
!  Modified:
!
!    12 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real BETA1, the skew or bias parameter.
!    BETA1 = 1 for no skew or bias.
!
!    Input, real BETA2, the tension parameter.
!    BETA2 = 0 for no tension.
!
!    Input, integer NDATA, the number of data values.
!
!    Input, real TDATA(NDATA), the abscissas of the data.
!
!    Input, real YDATA(NDATA), the data values.
!
!    Input, real TVAL, a point at which the spline is to be evaluated.
!
!    Output, real YVAL, the value of the function at TVAL.
!
  IMPLICIT NONE
!
  INTEGER ndata
!
  REAL a
  REAL b
  REAL beta1
  REAL beta2
  REAL bval
  REAL c
  REAL d
  REAL delta
  INTEGER left
  INTEGER right
  REAL tdata(ndata)
  REAL tval
  REAL u
  REAL ydata(ndata)
  REAL yval
!
!  Find the nearest interval [ TDATA(LEFT), TDATA(RIGHT) ] to TVAL.
!
  CALL rvec_bracket(ndata, tdata, tval, left, right)
!
!  Evaluate the 5 nonzero beta spline basis functions in the interval,
!  weighted by their corresponding data values.
!
  u = (tval - tdata(left)) / (tdata(right) - tdata(left))

  delta = 2.0E+00 + beta2 + 4.0E+00 * beta1 + 4.0E+00 * beta1**2 &
          + 2.0E+00 * beta1**3

  yval = 0.0E+00
!
!  Beta function associated with node LEFT - 1, (or "phantom node"),
!  evaluated in its 4th interval.
!
  bval = (2.0E+00 * beta1**3 * (1.0E+00 - u)**3) / delta

  IF (left - 1 > 0) THEN
    yval = yval + ydata(left - 1) * bval
  ELSE
    yval = yval + (2.0E+00 * ydata(1) - ydata(2)) * bval
  END IF
!
!  Beta function associated with node LEFT,
!  evaluated in its third interval.
!
  a = beta2 + 4.0E+00 * beta1 + 4.0E+00 * beta1**2

  b = -6.0E+00 * beta1 * (1.0E+00 - beta1) * (1.0E+00 + beta1)

  c = -3.0E+00 * (beta2 + 2.0E+00 * beta1**2 + 2.0E+00 * beta1**3)

  d = 2.0E+00 * (beta2 + beta1 + beta1**2 + beta1**3)

  bval = (a + u * (b + u * (c + u * d))) / delta

  yval = yval + ydata(left) * bval
!
!  Beta function associated with node RIGHT,
!  evaluated in its second interval.
!
  a = 2.0E+00

  b = 6.0E+00 * beta1

  c = 3.0E+00 * beta2 + 6.0E+00 * beta1**2

  d = -2.0E+00 * (1.0E+00 + beta2 + beta1 + beta1**2)

  bval = (a + u * (b + u * (c + u * d))) / delta

  yval = yval + ydata(right) * bval
!
!  Beta function associated with node RIGHT+1, (or "phantom node"),
!  evaluated in its first interval.
!
  bval = 2.0E+00 * u**3 / delta

  IF (right + 1 <= ndata) THEN
    yval = yval + ydata(right + 1) * bval
  ELSE
    yval = yval + (2.0E+00 * ydata(ndata) - ydata(ndata - 1)) * bval
  END IF

  RETURN
END
SUBROUTINE spline_constant_val(ndata, tdata, ydata, tval, yval)
!
!*******************************************************************************
!
!! SPLINE_CONSTANT_VAL evaluates a piecewise constant spline at a point.
!
!
!  Discussion:
!
!    NDATA-1 points TDATA define NDATA intervals, with the first
!    and last being semi-infinite.
!
!    The value of the spline anywhere in interval I is YDATA(I).
!
!  Modified:
!
!    16 November 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NDATA, the number of data points defining the spline.
!
!    Input, real TDATA(NDATA-1), the breakpoints.  The values of TDATA should
!    be distinct and increasing.
!
!    Input, YDATA(NDATA), the values of the spline in the intervals
!    defined by the breakpoints.
!
!    Input, real TVAL, the point at which the spline is to be evaluated.
!
!    Output, real YVAL, the value of the spline at TVAL.
!
  IMPLICIT NONE
!
  INTEGER ndata
!
  INTEGER i
  REAL tdata(ndata - 1)
  REAL tval
  REAL ydata(ndata)
  REAL yval
!
  DO i = 1, ndata - 1
    IF (tval <= tdata(i)) THEN
      yval = ydata(i)
      RETURN
    END IF
  END DO

  yval = ydata(ndata)

  RETURN
END
SUBROUTINE spline_cubic_set(n, t, y, ibcbeg, ybcbeg, ibcend, ybcend, ypp)
!
!*******************************************************************************
!
!! SPLINE_CUBIC_SET computes the second derivatives of a cubic spline.
!
!
!  Discussion:
!
!    For data interpolation, the user must call SPLINE_CUBIC_SET to
!    determine the second derivative data, passing in the data to be
!    interpolated, and the desired boundary conditions.
!
!    The data to be interpolated, plus the SPLINE_CUBIC_SET output,
!    defines the spline.  The user may then call SPLINE_CUBIC_VAL to
!    evaluate the spline at any point.
!
!    The cubic spline is a piecewise cubic polynomial.  The intervals
!    are determined by the "knots" or abscissas of the data to be
!    interpolated.  The cubic spline has continous first and second
!    derivatives over the entire interval of interpolation.
!
!    For any point T in the interval T(IVAL), T(IVAL+1), the form of
!    the spline is
!
!      SPL(T) = A(IVAL)
!             + B(IVAL) * ( T - T(IVAL) )
!             + C(IVAL) * ( T - T(IVAL) )**2
!             + D(IVAL) * ( T - T(IVAL) )**3
!
!    If we assume that we know the values Y(*) and YPP(*), which represent
!    the values and second derivatives of the spline at each knot, then
!    the coefficients can be computed as:
!
!      A(IVAL) = Y(IVAL)
!      B(IVAL) = ( Y(IVAL+1) - Y(IVAL) ) / ( T(IVAL+1) - T(IVAL) )
!        - ( YPP(IVAL+1) + 2 * YPP(IVAL) ) * ( T(IVAL+1) - T(IVAL) ) / 6
!      C(IVAL) = YPP(IVAL) / 2
!      D(IVAL) = ( YPP(IVAL+1) - YPP(IVAL) ) / ( 6 * ( T(IVAL+1) - T(IVAL) ) )
!
!    Since the first derivative of the spline is
!
!      SPL'(T) =     B(IVAL)
!              + 2 * C(IVAL) * ( T - T(IVAL) )
!              + 3 * D(IVAL) * ( T - T(IVAL) )**2,
!
!    the requirement that the first derivative be continuous at interior
!    knot I results in a total of N-2 equations, of the form:
!
!      B(IVAL-1) + 2 C(IVAL-1) * (T(IVAL)-T(IVAL-1))
!      + 3 * D(IVAL-1) * (T(IVAL) - T(IVAL-1))**2 = B(IVAL)
!
!    or, setting H(IVAL) = T(IVAL+1) - T(IVAL)
!
!      ( Y(IVAL) - Y(IVAL-1) ) / H(IVAL-1)
!      - ( YPP(IVAL) + 2 * YPP(IVAL-1) ) * H(IVAL-1) / 6
!      + YPP(IVAL-1) * H(IVAL-1)
!      + ( YPP(IVAL) - YPP(IVAL-1) ) * H(IVAL-1) / 2
!      =
!      ( Y(IVAL+1) - Y(IVAL) ) / H(IVAL)
!      - ( YPP(IVAL+1) + 2 * YPP(IVAL) ) * H(IVAL) / 6
!
!    or
!
!      YPP(IVAL-1) * H(IVAL-1) + 2 * YPP(IVAL) * ( H(IVAL-1) + H(IVAL) )
!      + YPP(IVAL) * H(IVAL)
!      =
!      6 * ( Y(IVAL+1) - Y(IVAL) ) / H(IVAL)
!      - 6 * ( Y(IVAL) - Y(IVAL-1) ) / H(IVAL-1)
!
!    Boundary conditions must be applied at the first and last knots.
!    The resulting tridiagonal system can be solved for the YPP values.
!
!  Modified:
!
!    20 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of data points; N must be at least 2.
!
!    Input, real T(N), the points where data is specified.
!    The values should be distinct, and increasing.
!
!    Input, real Y(N), the data values to be interpolated.
!
!    Input, integer IBCBEG, the left boundary condition flag:
!
!      0: the spline should be a quadratic over the first interval;
!      1: the first derivative at the left endpoint should be YBCBEG;
!      2: the second derivative at the left endpoint should be YBCBEG.
!
!    Input, real YBCBEG, the left boundary value, if needed.
!
!    Input, integer IBCEND, the right boundary condition flag:
!
!      0: the spline should be a quadratic over the last interval;
!      1: the first derivative at the right endpoint should be YBCEND;
!      2: the second derivative at the right endpoint should be YBCEND.
!
!    Input, real YBCEND, the right boundary value, if needed.
!
!    Output, real YPP(N), the second derivatives of the cubic spline.
!
  IMPLICIT NONE
!
  INTEGER n
!
  REAL diag(n)
  INTEGER i
  INTEGER ibcbeg
  INTEGER ibcend
  REAL sub(2:n)
  REAL sup(1:n - 1)
  REAL t(n)
  REAL y(n)
  REAL ybcbeg
  REAL ybcend
  REAL ypp(n)
!
!  Check.
!
  IF (n <= 1) THEN
    WRITE (*, '(a)') ' '
    WRITE (*, '(a)') 'SPLINE_CUBIC_SET - Fatal error!'
    WRITE (*, '(a)') '  The number of knots must be at least 2.'
    WRITE (*, '(a,i6)') '  The input value of N = ', n
    STOP
  END IF

  DO i = 1, n - 1
    IF (t(i) >= t(i + 1)) THEN
      WRITE (*, '(a)') ' '
      WRITE (*, '(a)') 'SPLINE_CUBIC_SET - Fatal error!'
      WRITE (*, '(a)') '  The knots must be strictly increasing, but'
      WRITE (*, '(a,i6,a,g14.6)') '  T(', i, ') = ', t(i)
      WRITE (*, '(a,i6,a,g14.6)') '  T(', i + 1, ') = ', t(i + 1)
      STOP
    END IF
  END DO
!
!  Set the first equation.
!
  IF (ibcbeg == 0) THEN
    ypp(1) = 0.0E+00
    diag(1) = 1.0E+00
    sup(1) = -1.0E+00
  ELSE IF (ibcbeg == 1) THEN
    ypp(1) = (y(2) - y(1)) / (t(2) - t(1)) - ybcbeg
    diag(1) = (t(2) - t(1)) / 3.0E+00
    sup(1) = (t(2) - t(1)) / 6.0E+00
  ELSE IF (ibcbeg == 2) THEN
    ypp(1) = ybcbeg
    diag(1) = 1.0E+00
    sup(1) = 0.0E+00
  ELSE
    WRITE (*, '(a)') ' '
    WRITE (*, '(a)') 'SPLINE_CUBIC_SET - Fatal error!'
    WRITE (*, '(a)') '  The boundary flag IBCBEG must be 0, 1 or 2.'
    WRITE (*, '(a,i6)') '  The input value is IBCBEG = ', ibcbeg
    STOP
  END IF
!
!  Set the intermediate equations.
!
  DO i = 2, n - 1
    ypp(i) = (y(i + 1) - y(i)) / (t(i + 1) - t(i)) &
             - (y(i) - y(i - 1)) / (t(i) - t(i - 1))
    sub(i) = (t(i) - t(i - 1)) / 6.0E+00
    diag(i) = (t(i + 1) - t(i - 1)) / 3.0E+00
    sup(i) = (t(i + 1) - t(i)) / 6.0E+00
  END DO
!
!  Set the last equation.
!
  IF (ibcend == 0) THEN
    ypp(n) = 0.0E+00
    sub(n) = -1.0E+00
    diag(n) = 1.0E+00
  ELSE IF (ibcend == 1) THEN
    ypp(n) = ybcend - (y(n) - y(n - 1)) / (t(n) - t(n - 1))
    sub(n) = (t(n) - t(n - 1)) / 6.0E+00
    diag(n) = (t(n) - t(n - 1)) / 3.0E+00
  ELSE IF (ibcend == 2) THEN
    ypp(n) = ybcend
    sub(n) = 0.0E+00
    diag(n) = 1.0E+00
  ELSE
    WRITE (*, '(a)') ' '
    WRITE (*, '(a)') 'SPLINE_CUBIC_SET - Fatal error!'
    WRITE (*, '(a)') '  The boundary flag IBCEND must be 0, 1 or 2.'
    WRITE (*, '(a,i6)') '  The input value is IBCEND = ', ibcend
    STOP
  END IF
!
!  Special case:
!    N = 2, IBCBEG = IBCEND = 0.
!
  IF (n == 2 .AND. ibcbeg == 0 .AND. ibcend == 0) THEN

    ypp(1) = 0.0E+00
    ypp(2) = 0.0E+00
!
!  Solve the linear system.
!
  ELSE

    CALL s3_fs(sub, diag, sup, n, ypp, ypp)

  END IF

  RETURN
END
SUBROUTINE spline_cubic_val(n, t, y, ypp, tval, yval, ypval, yppval)
!
!*******************************************************************************
!
!! SPLINE_CUBIC_VAL evaluates a cubic spline at a specific point.
!
!
!  Discussion:
!
!    SPLINE_CUBIC_SET must have already been called to define the
!    values of YPP.
!
!    For any point T in the interval T(IVAL), T(IVAL+1), the form of
!    the spline is
!
!      SPL(T) = A
!             + B * ( T - T(IVAL) )
!             + C * ( T - T(IVAL) )**2
!             + D * ( T - T(IVAL) )**3
!
!    Here:
!      A = Y(IVAL)
!      B = ( Y(IVAL+1) - Y(IVAL) ) / ( T(IVAL+1) - T(IVAL) )
!        - ( YPP(IVAL+1) + 2 * YPP(IVAL) ) * ( T(IVAL+1) - T(IVAL) ) / 6
!      C = YPP(IVAL) / 2
!      D = ( YPP(IVAL+1) - YPP(IVAL) ) / ( 6 * ( T(IVAL+1) - T(IVAL) ) )
!
!  Modified:
!
!    20 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of data values.
!
!    Input, real T(N), the knot values.
!
!    Input, real Y(N), the data values at the knots.
!
!    Input, real YPP(N), the second derivatives of the spline at the knots.
!
!    Input, real TVAL, a point, typically between T(1) and T(N), at
!    which the spline is to be evalulated.  If TVAL lies outside
!    this range, extrapolation is used.
!
!    Output, real YVAL, YPVAL, YPPVAL, the value of the spline, and
!    its first two derivatives at TVAL.
!
  IMPLICIT NONE
!
  INTEGER n
!
  REAL dt
  REAL h
  INTEGER left
  INTEGER right
  REAL t(n)
  REAL tval
  REAL y(n)
  REAL ypp(n)
  REAL yppval
  REAL ypval
  REAL yval
!
!  Determine the interval [T(LEFT), T(RIGHT)] that contains TVAL.
!  Values below T(1) or above T(N) use extrapolation.
!
  CALL rvec_bracket(n, t, tval, left, right)
!
!  Evaluate the polynomial.
!
  dt = tval - t(left)
  h = t(right) - t(left)

  yval = y(left) &
         + dt * ((y(right) - y(left)) / h &
                 - (ypp(right) / 6.0E+00 + ypp(left) / 3.0E+00) * h &
                 + dt * (0.5E+00 * ypp(left) &
                         + dt * ((ypp(right) - ypp(left)) / (6.0E+00 * h))))

  ypval = (y(right) - y(left)) / h &
          - (ypp(right) / 6.0E+00 + ypp(left) / 3.0E+00) * h &
          + dt * (ypp(left) &
                  + dt * (0.5E+00 * (ypp(right) - ypp(left)) / h))

  yppval = ypp(left) + dt * (ypp(right) - ypp(left)) / h

  RETURN
END
SUBROUTINE spline_cubic_val2(n, t, y, ypp, left, tval, yval, ypval, yppval)
!
!*******************************************************************************
!
!! SPLINE_CUBIC_VAL2 evaluates a cubic spline at a specific point.
!
!
!  Discussion:
!
!    This routine is a modification of SPLINE_CUBIC_VAL; it allows the
!    user to speed up the code by suggesting the appropriate T interval
!    to search first.
!
!    SPLINE_CUBIC_SET must have already been called to define the
!    values of YPP.
!
!    In the LEFT interval, let RIGHT = LEFT+1.  The form of the spline is
!
!      SPL(T) =
!          A
!        + B * ( T - T(LEFT) )
!        + C * ( T - T(LEFT) )**2
!        + D * ( T - T(LEFT) )**3
!
!    Here:
!      A = Y(LEFT)
!      B = ( Y(RIGHT) - Y(LEFT) ) / ( T(RIGHT) - T(LEFT) )
!        - ( YPP(RIGHT) + 2 * YPP(LEFT) ) * ( T(RIGHT) - T(LEFT) ) / 6
!      C = YPP(LEFT) / 2
!      D = ( YPP(RIGHT) - YPP(LEFT) ) / ( 6 * ( T(RIGHT) - T(LEFT) ) )
!
!  Modified:
!
!    20 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of knots.
!
!    Input, real T(N), the knot values.
!
!    Input, real Y(N), the data values at the knots.
!
!    Input, real YPP(N), the second derivatives of the spline at
!    the knots.
!
!    Input/output, integer LEFT, the suggested T interval to search.
!    LEFT should be between 1 and N-1.  If LEFT is not in this range,
!    then its value will be ignored.  On output, LEFT is set to the
!    actual interval in which TVAL lies.
!
!    Input, real TVAL, a point, typically between T(1) and T(N), at
!    which the spline is to be evalulated.  If TVAL lies outside
!    this range, extrapolation is used.
!
!    Output, real YVAL, YPVAL, YPPVAL, the value of the spline, and
!    its first two derivatives at TVAL.
!
  IMPLICIT NONE
!
  INTEGER n
!
  REAL dt
  REAL h
  INTEGER left
  INTEGER right
  REAL t(n)
  REAL tval
  REAL y(n)
  REAL ypp(n)
  REAL yppval
  REAL ypval
  REAL yval
!
!  Determine the interval [T(LEFT), T(RIGHT)] that contains TVAL.
!
!  What you want from RVEC_BRACKET3 is that TVAL is to be computed
!  by the data in interval {T(LEFT), T(RIGHT)].
!
  left = 0
  CALL rvec_bracket3(n, t, tval, left)
  right = left + 1
!
!  In the interval LEFT, the polynomial is in terms of a normalized
!  coordinate  ( DT / H ) between 0 and 1.
!
  dt = tval - t(left)
  h = t(right) - t(left)

  yval = y(left) + dt * ((y(right) - y(left)) / h &
                         - (ypp(right) / 6.0E+00 + ypp(left) / 3.0E+00) * h &
                         + dt * (0.5E+00 * ypp(left) &
                                 + dt * ((ypp(right) - ypp(left)) / (6.0E+00 * h))))

  ypval = (y(right) - y(left)) / h &
          - (ypp(right) / 6.0E+00 + ypp(left) / 3.0E+00) * h &
          + dt * (ypp(left) &
                  + dt * (0.5E+00 * (ypp(right) - ypp(left)) / h))

  yppval = ypp(left) + dt * (ypp(right) - ypp(left)) / h

  RETURN
END
SUBROUTINE spline_hermite_set(ndata, tdata, c)
!
!*************************************************************************
!
!! SPLINE_HERMITE_SET sets up a piecewise cubic Hermite interpolant.
!
!
!  Reference:
!
!    Conte and de Boor,
!    Algorithm CALCCF,
!    Elementary Numerical Analysis,
!    1973, page 235.
!
!  Modified:
!
!    06 April 1999
!
!  Parameters:
!
!    Input, integer NDATA, the number of data points.
!    NDATA must be at least 2.
!
!    Input, real TDATA(NDATA), the abscissas of the data points.
!    The entries of TDATA are assumed to be strictly increasing.
!
!    Input/output, real C(4,NDATA).
!
!    On input, C(1,I) and C(2,I) should contain the value of the
!    function and its derivative at TDATA(I), for I = 1 to NDATA.
!    These values will not be changed by this routine.
!
!    On output, C(3,I) and C(4,I) contain the quadratic
!    and cubic coefficients of the Hermite polynomial
!    in the interval (TDATA(I), TDATA(I+1)), for I=1 to NDATA-1.
!    C(3,NDATA) and C(4,NDATA) are set to 0.
!
!    In the interval (TDATA(I), TDATA(I+1)), the interpolating Hermite
!    polynomial is given by
!
!    SVAL(TVAL) =                 C(1,I)
!       + ( TVAL - TDATA(I) ) * ( C(2,I)
!       + ( TVAL - TDATA(I) ) * ( C(3,I)
!       + ( TVAL - TDATA(I) ) *   C(4,I) ) )
!
  IMPLICIT NONE
!
  INTEGER ndata
!
  REAL c(4, ndata)
  REAL divdif1
  REAL divdif3
  REAL dt
  INTEGER i
  REAL tdata(ndata)
!
  DO i = 1, ndata - 1
    dt = tdata(i + 1) - tdata(i)
    divdif1 = (c(1, i + 1) - c(1, i)) / dt
    divdif3 = c(2, i) + c(2, i + 1) - 2.0E+00 * divdif1
    c(3, i) = (divdif1 - c(2, i) - divdif3) / dt
    c(4, i) = divdif3 / dt**2
  END DO

  c(3, ndata) = 0.0E+00
  c(4, ndata) = 0.0E+00

  RETURN
END
SUBROUTINE spline_hermite_val(ndata, tdata, c, tval, sval)
!
!*************************************************************************
!
!! SPLINE_HERMITE_VAL evaluates a piecewise cubic Hermite interpolant.
!
!
!  Discussion:
!
!    SPLINE_HERMITE_SET must be called first, to set up the
!    spline data from the raw function and derivative data.
!
!  Reference:
!
!    Conte and de Boor,
!    Algorithm PCUBIC,
!    Elementary Numerical Analysis,
!    1973, page 234.
!
!  Modified:
!
!    06 April 1999
!
!  Parameters:
!
!    Input, integer NDATA, the number of data points.
!    NDATA is assumed to be at least 2.
!
!    Input, real TDATA(NDATA), the abscissas of the data points.
!    The entries of TDATA are assumed to be strictly increasing.
!
!    Input, real C(4,NDATA), contains the data computed by
!    SPLINE_HERMITE_SET.
!
!    Input, real TVAL, the point where the interpolant is to
!    be evaluated.
!
!    Output, real SVAL, the value of the interpolant at TVAL.
!
  IMPLICIT NONE
!
  INTEGER ndata
!
  REAL c(4, ndata)
  REAL dt
  INTEGER left
  INTEGER right
  REAL sval
  REAL tdata(ndata)
  REAL tval
!
!  Find the interval [ TDATA(LEFT), TDATA(RIGHT) ] that contains
!  or is nearest to TVAL.
!
  CALL rvec_bracket(ndata, tdata, tval, left, right)
!
!  Evaluate the cubic polynomial.
!
  dt = tval - tdata(left)

  sval = c(1, left) + dt * (c(2, left) + dt * (c(3, left) + dt * c(4, left)))

  RETURN
END
SUBROUTINE spline_linear_int(ndata, tdata, ydata, a, b, int_val)
!
!*******************************************************************************
!
!! SPLINE_LINEAR_INT evaluates the integral of a linear spline.
!
!
!  Modified:
!
!    01 November 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NDATA, the number of data points defining the spline.
!
!    Input, real TDATA(NDATA), YDATA(NDATA), the values of the independent
!    and dependent variables at the data points.  The values of TDATA should
!    be distinct and increasing.
!
!    Input, real A, B, the interval over which the integral is desired.
!
!    Output, real INT_VAL, the value of the integral.
!
  IMPLICIT NONE
!
  INTEGER ndata
!
  REAL a
  REAL a_copy
  INTEGER a_left
  INTEGER a_right
  REAL b
  REAL b_copy
  INTEGER b_left
  INTEGER b_right
  INTEGER i_left
  REAL int_val
  REAL tdata(ndata)
  REAL tval
  REAL ydata(ndata)
  REAL yp
  REAL yval
!
  int_val = 0.0E+00

  IF (a == b) THEN
    RETURN
  END IF

  a_copy = MIN(a, b)
  b_copy = MAX(a, b)
!
!  Find the interval [ TDATA(A_LEFT), TDATA(A_RIGHT) ] that contains, or is
!  nearest to, A.
!
  CALL rvec_bracket(ndata, tdata, a_copy, a_left, a_right)
!
!  Find the interval [ TDATA(B_LEFT), TDATA(B_RIGHT) ] that contains, or is
!  nearest to, B.
!
  CALL rvec_bracket(ndata, tdata, b_copy, b_left, b_right)
!
!  If A and B are in the same interval...
!
  IF (a_left == b_left) THEN

    tval = (a_copy + b_copy) / 2.0E+00

    yp = (ydata(a_right) - ydata(a_left)) / &
         (tdata(a_right) - tdata(a_left))

    yval = ydata(a_left) + (tval - tdata(a_left)) * yp

    int_val = yval * (b_copy - a_copy)

    RETURN
  END IF
!
!  Otherwise, integrate from:
!
!  A               to TDATA(A_RIGHT),
!  TDATA(A_RIGHT)  to TDATA(A_RIGHT+1),...
!  TDATA(B_LEFT-1) to TDATA(B_LEFT),
!  TDATA(B_LEFT)   to B.
!
!  Use the fact that the integral of a linear function is the
!  value of the function at the midpoint times the width of the interval.
!
  tval = (a_copy + tdata(a_right)) / 2.0E+00

  yp = (ydata(a_right) - ydata(a_left)) / &
       (tdata(a_right) - tdata(a_left))

  yval = ydata(a_left) + (tval - tdata(a_left)) * yp

  int_val = int_val + yval * (tdata(a_right) - a_copy)

  DO i_left = a_right, b_left - 1

    tval = (tdata(i_left + 1) + tdata(i_left)) / 2.0E+00

    yp = (ydata(i_left + 1) - ydata(i_left)) / &
         (tdata(i_left + 1) - tdata(i_left))

    yval = ydata(i_left) + (tval - tdata(i_left)) * yp

    int_val = int_val + yval * (tdata(i_left + 1) - tdata(i_left))

  END DO

  tval = (tdata(b_left) + b_copy) / 2.0E+00

  yp = (ydata(b_right) - ydata(b_left)) / &
       (tdata(b_right) - tdata(b_left))

  yval = ydata(b_left) + (tval - tdata(b_left)) * yp

  int_val = int_val + yval * (b_copy - tdata(b_left))

  IF (b < a) THEN
    int_val = -int_val
  END IF

  RETURN
END
SUBROUTINE spline_linear_intset(int_n, int_x, int_v, data_n, data_x, data_y)
!
!*******************************************************************************
!
!! SPLINE_LINEAR_INTSET sets a linear spline with given integral properties.
!
!
!  Discussion:
!
!    The user has in mind an interval, divided by INT_N+1 points into
!    INT_N intervals.  A linear spline is to be constructed,
!    with breakpoints at the centers of each interval, and extending
!    continuously to the left of the first and right of the last
!    breakpoints.  The constraint on the linear spline is that it is
!    required that it have a given integral value over each interval.
!
!    A tridiagonal linear system of equations is solved for the
!    values of the spline at the breakpoints.
!
!  Modified:
!
!    02 November 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer INT_N, the number of intervals.
!
!    Input, real INT_X(INT_N+1), the points that define the intervals.
!    Interval I lies between INT_X(I) and INT_X(I+1).
!
!    Input, real INT_V(INT_N), the desired value of the integral of the
!    linear spline over each interval.
!
!    Output, integer DATA_N, the number of data points defining the spline.
!    (This is the same as INT_N).
!
!    Output, real DATA_X(DATA_N), DATA_Y(DATA_N), the values of the independent
!    and dependent variables at the data points.  The values of DATA_X are
!    the interval midpoints.  The values of DATA_Y are determined in such
!    a way that the exact integral of the linear spline over interval I
!    is equal to INT_V(I).
!
  IMPLICIT NONE
!
  INTEGER int_n
!
  REAL c(int_n)
  REAL d(int_n)
  INTEGER data_n
  REAL data_x(int_n)
  REAL data_y(int_n)
  REAL e(int_n)
  INTEGER info
  REAL int_v(int_n)
  REAL int_x(int_n + 1)
!
!  Set up the easy stuff.
!
  data_n = int_n
  data_x(1:data_n) = 0.5E+00 * (int_x(1:data_n) + int_x(2:data_n + 1))
!
!  Set up C, D, E, the coefficients of the linear system.
!
  c(1) = 0.0E+00
  c(2:data_n - 1) = 1.0 &
                    - (0.5 * (data_x(2:data_n - 1) + int_x(2:data_n - 1)) &
                       - data_x(1:data_n - 2)) &
                    / (data_x(2:data_n - 1) - data_x(1:data_n - 2))
  c(data_n) = 0.0E+00

  d(1) = int_x(2) - int_x(1)

  d(2:data_n - 1) = 1.0 &
                    + (0.5 * (data_x(2:data_n - 1) + int_x(2:data_n - 1)) &
                       - data_x(1:data_n - 2)) &
                    / (data_x(2:data_n - 1) - data_x(1:data_n - 2)) &
                    - (0.5 * (data_x(2:data_n - 1) + int_x(3:data_n)) - data_x(2:data_n - 1)) &
                    / (data_x(3:data_n) - data_x(2:data_n - 1))

  d(data_n) = int_x(data_n + 1) - int_x(data_n)

  e(1) = 0.0E+00

  e(2:data_n - 1) = (0.5 * (data_x(2:data_n - 1) + int_x(3:data_n)) &
                     - data_x(2:data_n - 1)) / (data_x(3:data_n) - data_x(2:data_n - 1))

  e(data_n) = 0.0E+00
!
!  Set up DATA_Y, which begins as the right hand side of the linear system.
!
  data_y(1) = int_v(1)
  data_y(2:data_n - 1) = 2.0E+00 * int_v(2:data_n - 1) &
                         / (int_x(3:int_n) - int_x(2:int_n - 1))
  data_y(data_n) = int_v(data_n)
!
!  Solve the linear system.
!
  CALL sgtsl(data_n, c, d, e, data_y, info)

  IF (info /= 0) THEN
    WRITE (*, '(a)') ' '
    WRITE (*, '(a)') 'SPLINE_LINEAR_INTSET - Fatal error!'
    WRITE (*, '(a)') '  The linear system is singular.'
    STOP
  END IF

  RETURN
END
SUBROUTINE spline_linear_val(ndata, tdata, ydata, tval, yval, ypval)
!
!*******************************************************************************
!
!! SPLINE_LINEAR_VAL evaluates a linear spline at a specific point.
!
!
!  Discussion:
!
!    Because of the extremely simple form of the linear spline,
!    the raw data points ( TDATA(I), YDATA(I)) can be used directly to
!    evaluate the spline at any point.  No processing of the data
!    is required.
!
!  Modified:
!
!    06 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NDATA, the number of data points defining the spline.
!
!    Input, real TDATA(NDATA), YDATA(NDATA), the values of the independent
!    and dependent variables at the data points.  The values of TDATA should
!    be distinct and increasing.
!
!    Input, real TVAL, the point at which the spline is to be evaluated.
!
!    Output, real YVAL, YPVAL, the value of the spline and its first
!    derivative dYdT at TVAL.  YPVAL is not reliable if TVAL is exactly
!    equal to TDATA(I) for some I.
!
  IMPLICIT NONE
!
  INTEGER ndata
!
  INTEGER left
  INTEGER right
  REAL tdata(ndata)
  REAL tval
  REAL ydata(ndata)
  REAL ypval
  REAL yval
!
!  Find the interval [ TDATA(LEFT), TDATA(RIGHT) ] that contains, or is
!  nearest to, TVAL.
!
  CALL rvec_bracket(ndata, tdata, tval, left, right)
!
!  Now evaluate the piecewise linear function.
!
  ypval = (ydata(right) - ydata(left)) / (tdata(right) - tdata(left))

  yval = ydata(left) + (tval - tdata(left)) * ypval

  RETURN
END
SUBROUTINE spline_overhauser_nonuni_val(ndata, tdata, ydata, tval, yval)
!
!*******************************************************************************
!
!! SPLINE_OVERHAUSER_NONUNI_VAL evaluates the nonuniform Overhauser spline.
!
!
!  Discussion:
!
!    The nonuniformity refers to the fact that the abscissas values
!    need not be uniformly spaced.
!
!  Diagnostics:
!
!    The values of ALPHA and BETA have to be properly assigned.
!    The basis matrices for the first and last interval have to
!    be computed.
!
!  Modified:
!
!    06 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NDATA, the number of data points.
!
!    Input, real TDATA(NDATA), the abscissas of the data points.
!    The values of TDATA are assumed to be distinct and increasing.
!
!    Input, real YDATA(NDATA), the data values.
!
!    Input, real TVAL, the value where the spline is to
!    be evaluated.
!
!    Output, real YVAL, the value of the spline at TVAL.
!
  IMPLICIT NONE
!
  INTEGER ndata
!
  REAL alpha
  REAL beta
  INTEGER left
  REAL mbasis(4, 4)
  REAL mbasis_l(3, 3)
  REAL mbasis_r(3, 3)
  INTEGER right
  REAL tdata(ndata)
  REAL tval
  REAL ydata(ndata)
  REAL yval
!
!  Find the nearest interval [ TDATA(LEFT), TDATA(RIGHT) ] to TVAL.
!
  CALL rvec_bracket(ndata, tdata, tval, left, right)
!
!  Evaluate the spline in the given interval.
!
  IF (left == 1) THEN

    alpha = 1.0E+00
    CALL basis_matrix_overhauser_nul(alpha, mbasis_l)

    CALL basis_matrix_tmp(1, 3, mbasis_l, ndata, tdata, ydata, tval, yval)

  ELSE IF (left < ndata - 1) THEN

    alpha = 1.0E+00
    beta = 1.0E+00
    CALL basis_matrix_overhauser_nonuni(alpha, beta, mbasis)

    CALL basis_matrix_tmp(left, 4, mbasis, ndata, tdata, ydata, tval, yval)

  ELSE IF (left == ndata - 1) THEN

    beta = 1.0E+00
    CALL basis_matrix_overhauser_nur(beta, mbasis_r)

    CALL basis_matrix_tmp(left, 3, mbasis_r, ndata, tdata, ydata, tval, yval)

  END IF

  RETURN
END
SUBROUTINE spline_overhauser_uni_val(ndata, tdata, ydata, tval, yval)
!
!*******************************************************************************
!
!! SPLINE_OVERHAUSER_UNI_VAL evaluates the uniform Overhauser spline.
!
!
!  Modified:
!
!    06 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NDATA, the number of data points.
!
!    Input, real TDATA(NDATA), the abscissas of the data points.
!    The values of TDATA are assumed to be distinct and increasing.
!    This routine also assumes that the values of TDATA are uniformly
!    spaced; for instance, TDATA(1) = 10, TDATA(2) = 11, TDATA(3) = 12...
!
!    Input, real YDATA(NDATA), the data values.
!
!    Input, real TVAL, the value where the spline is to
!    be evaluated.
!
!    Output, real YVAL, the value of the spline at TVAL.
!
  IMPLICIT NONE
!
  INTEGER ndata
!
  INTEGER left
  REAL mbasis(4, 4)
  REAL mbasis_l(3, 3)
  REAL mbasis_r(3, 3)
  INTEGER right
  REAL tdata(ndata)
  REAL tval
  REAL ydata(ndata)
  REAL yval
!
!  Find the nearest interval [ TDATA(LEFT), TDATA(RIGHT) ] to TVAL.
!
  CALL rvec_bracket(ndata, tdata, tval, left, right)
!
!  Evaluate the spline in the given interval.
!
  IF (left == 1) THEN

    CALL basis_matrix_overhauser_uni_l(mbasis_l)

    CALL basis_matrix_tmp(1, 3, mbasis_l, ndata, tdata, ydata, tval, yval)

  ELSE IF (left < ndata - 1) THEN

    CALL basis_matrix_overhauser_uni(mbasis)

    CALL basis_matrix_tmp(left, 4, mbasis, ndata, tdata, ydata, tval, yval)

  ELSE IF (left == ndata - 1) THEN

    CALL basis_matrix_overhauser_uni_r(mbasis_r)

    CALL basis_matrix_tmp(left, 3, mbasis_r, ndata, tdata, ydata, tval, yval)

  END IF

  RETURN
END
SUBROUTINE spline_overhauser_val(ndim, ndata, tdata, ydata, tval, yval)
!
!*******************************************************************************
!
!! SPLINE_OVERHAUSER_VAL evaluates an Overhauser spline.
!
!
!  Discussion:
!
!    Over the first and last intervals, the Overhauser spline is a
!    quadratic.  In the intermediate intervals, it is a piecewise cubic.
!    The Overhauser spline is also known as the Catmull-Rom spline.
!
!  Reference:
!
!    H Brewer and D Anderson,
!    Visual Interaction with Overhauser Curves and Surfaces,
!    SIGGRAPH 77, pages 132-137.
!
!    E Catmull and R Rom,
!    A Class of Local Interpolating Splines,
!    in Computer Aided Geometric Design,
!    edited by R Barnhill and R Reisenfeld,
!    Academic Press, 1974, pages 317-326.
!
!    David Rogers and Alan Adams,
!    Mathematical Elements of Computer Graphics,
!    McGraw Hill, 1990, Second Edition, pages 278-289.
!
!  Modified:
!
!   08 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NDIM, the dimension of a single data point.
!    NDIM must be at least 1.  There is an internal limit on NDIM,
!    called MAXDIM, which is presently set to 5.
!
!    Input, integer NDATA, the number of data points.
!    NDATA must be at least 3.
!
!    Input, real TDATA(NDATA), the abscissas of the data points.  The
!    values in TDATA must be in strictly ascending order.
!
!    Input, real YDATA(NDIM,NDATA), the data points corresponding to
!    the abscissas.
!
!    Input, real TVAL, the abscissa value at which the spline
!    is to be evaluated.  Normally, TDATA(1) <= TVAL <= T(NDATA), and
!    the data will be interpolated.  For TVAL outside this range,
!    extrapolation will be used.
!
!    Output, real YVAL(NDIM), the value of the spline at TVAL.
!
  IMPLICIT NONE
!
  INTEGER, PARAMETER :: MAXDIM = 5
  INTEGER ndata
  INTEGER ndim
!
  INTEGER i
  INTEGER left
  INTEGER order
  INTEGER right
  REAL tdata(ndata)
  REAL tval
  REAL ydata(ndim, ndata)
  REAL yl(MAXDIM)
  REAL yr(MAXDIM)
  REAL yval(ndim)
!
!  Check.
!
  CALL rvec_order_type(ndata, tdata, order)

  IF (order /= 2) THEN
    WRITE (*, '(a)') ' '
    WRITE (*, '(a)') 'SPLINE_OVERHAUSER_VAL - Fatal error!'
    WRITE (*, '(a)') '  The data abscissas are not strictly ascending.'
    STOP
  END IF

  IF (ndata < 3) THEN
    WRITE (*, '(a)') ' '
    WRITE (*, '(a)') 'SPLINE_OVERHAUSER_VAL - Fatal error!'
    WRITE (*, '(a)') '  NDATA < 3.'
    STOP
  END IF

  IF (ndim < 1) THEN
    WRITE (*, '(a)') ' '
    WRITE (*, '(a)') 'SPLINE_OVERHAUSER_VAL - Fatal error!'
    WRITE (*, '(a)') '  NDIM < 1.'
    STOP
  END IF

  IF (ndim > maxdim) THEN
    WRITE (*, '(a)') ' '
    WRITE (*, '(a)') 'SPLINE_OVERHAUSER_VAL - Fatal error!'
    WRITE (*, '(a)') '  NDIM > MAXDIM.'
    STOP
  END IF
!
!  Locate the abscissa interval T(LEFT), T(LEFT+1) nearest to or
!  containing TVAL.
!
  CALL rvec_bracket(ndata, tdata, tval, left, right)
!
!  Evaluate the "left hand" quadratic defined at T(LEFT-1), T(LEFT), T(RIGHT).
!
  IF (left - 1 > 0) THEN
    CALL parabola_val2(ndim, ndata, tdata, ydata, left - 1, tval, yl)
  END IF
!
!  Evaluate the "right hand" quadratic defined at T(LEFT), T(RIGHT), T(RIGHT+1).
!
  IF (right + 1 <= ndata) THEN
    CALL parabola_val2(ndim, ndata, tdata, ydata, left, tval, yr)
  END IF
!
!  Average the quadratics.
!
  IF (left == 1) THEN

    yval(1:ndim) = yr(1:ndim)

  ELSE IF (right < ndata) THEN

    yval(1:ndim) = ((tdata(right) - tval) * yl(1:ndim) &
                    + (tval - tdata(left)) * yr(1:ndim)) / (tdata(right) - tdata(left))

  ELSE

    yval(1:ndim) = yl(1:ndim)

  END IF

  RETURN
END
SUBROUTINE spline_quadratic_val(ndata, tdata, ydata, tval, yval, ypval)
!
!*******************************************************************************
!
!! SPLINE_QUADRATIC_VAL evaluates a quadratic spline at a specific point.
!
!
!  Discussion:
!
!    Because of the simple form of a piecewise quadratic spline,
!    the raw data points ( TDATA(I), YDATA(I)) can be used directly to
!    evaluate the spline at any point.  No processing of the data
!    is required.
!
!  Modified:
!
!    24 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NDATA, the number of data points defining the spline.
!    NDATA should be odd.
!
!    Input, real TDATA(NDATA), YDATA(NDATA), the values of the independent
!    and dependent variables at the data points.  The values of TDATA should
!    be distinct and increasing.
!
!    Input, real TVAL, the point at which the spline is to be evaluated.
!
!    Output, real YVAL, YPVAL, the value of the spline and its first
!    derivative dYdT at TVAL.  YPVAL is not reliable if TVAL is exactly
!    equal to TDATA(I) for some I.
!
  IMPLICIT NONE
!
  INTEGER ndata
!
  REAL dif1
  REAL dif2
  INTEGER left
  INTEGER right
  REAL t1
  REAL t2
  REAL t3
  REAL tdata(ndata)
  REAL tval
  REAL y1
  REAL y2
  REAL y3
  REAL ydata(ndata)
  REAL ypval
  REAL yval
!
  IF (MOD(ndata, 3) == 0) THEN
    WRITE (*, '(a)') ' '
    WRITE (*, '(a)') 'SPLINE_QUADRATIC_VAL - Fatal error!'
    WRITE (*, '(a)') '  NDATA must be odd.'
    STOP
  END IF
!
!  Find the interval [ TDATA(LEFT), TDATA(RIGHT) ] that contains, or is
!  nearest to, TVAL.
!
  CALL rvec_bracket(ndata, tdata, tval, left, right)
!
!  Force LEFT to be odd.
!
  IF (MOD(left, 2) == 0) THEN
    left = left - 1
  END IF
!
!  Copy out the three abscissas.
!
  t1 = tdata(left)
  t2 = tdata(left + 1)
  t3 = tdata(left + 2)

  IF (t1 >= t2 .OR. t2 >= t3) THEN
    WRITE (*, '(a)') ' '
    WRITE (*, '(a)') 'SPLINE_QUADRATIC_VAL - Fatal error!'
    WRITE (*, '(a)') '  T1 >= T2 or T2 >= T3.'
    STOP
  END IF
!
!  Construct and evaluate a parabolic interpolant for the data
!  in each dimension.
!
  y1 = ydata(left)
  y2 = ydata(left + 1)
  y3 = ydata(left + 2)

  dif1 = (y2 - y1) / (t2 - t1)

  dif2 = ((y3 - y1) / (t3 - t1) &
          - (y2 - y1) / (t2 - t1)) / (t3 - t2)

  yval = y1 + (tval - t1) * (dif1 + (tval - t2) * dif2)
  ypval = dif1 + dif2 * (2.0E+00 * tval - t1 - t2)

  RETURN
END
SUBROUTINE timestamp()
!
!*******************************************************************************
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!
!  Example:
!
!    May 31 2001   9:45:54.872 AM
!
!  Modified:
!
!    31 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  IMPLICIT NONE
!
  CHARACTER(len=8) ampm
  INTEGER d
  CHARACTER(len=8) date
  INTEGER h
  INTEGER m
  INTEGER mm
  CHARACTER(len=9), PARAMETER, DIMENSION(12) :: month = (/ &
                                                'January  ', 'February ', 'March    ', 'April    ', &
                                                'May      ', 'June     ', 'July     ', 'August   ', &
                                                'September', 'October  ', 'November ', 'December '/)
  INTEGER n
  INTEGER s
  CHARACTER(len=10) time
  INTEGER values(8)
  INTEGER y
  CHARACTER(len=5) zone
!
  CALL DATE_AND_TIME(date, time, zone, values)

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  IF (h < 12) THEN
    ampm = 'AM'
  ELSE IF (h == 12) THEN
    IF (n == 0 .AND. s == 0) THEN
      ampm = 'Noon'
    ELSE
      ampm = 'PM'
    END IF
  ELSE
    h = h - 12
    IF (h < 12) THEN
      ampm = 'PM'
    ELSE IF (h == 12) THEN
      IF (n == 0 .AND. s == 0) THEN
        ampm = 'Midnight'
      ELSE
        ampm = 'AM'
      END IF
    END IF
  END IF

  WRITE (*, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)') &
    TRIM(month(m)), d, y, h, ':', n, ':', s, '.', mm, TRIM(ampm)

  RETURN
END
