MODULE generalHarmonicFunctions_m
  USE, INTRINSIC :: iso_c_binding, ONLY: c_double
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: harmonicFunction, sphericalHarmonic, planeHarmonic
  PUBLIC :: dSphericalHarmonicDTheta, dSphericalHarmonicDPhi
  PUBLIC :: dPlaneHarmonicDX, dPlaneHarmonicDY

  TYPE, ABSTRACT :: harmonics_t
    INTEGER :: l, m
    REAL(c_double) :: value

    CONTAINS
      PROCEDURE(forward), DEFERRED :: harmonicFunction
  END TYPE harmonics_t

  ! Abstract interface for general harmonic functions
  ABSTRACT INTERFACE
    FUNCTION harmonicFunction(l, m, x, y) RESULT(value)
      IMPORT :: c_double
      INTEGER, INTENT(IN) :: l, m
      REAL(c_double), INTENT(IN) :: x, y
      REAL(c_double) :: value
    END FUNCTION harmonicFunction
  END INTERFACE

CONTAINS

  !***********************************************************************
  ! Function: sphericalHarmonic
  ! Description: Computes the real-valued spherical harmonic function Y_lm(theta, phi).
  ! Input: l - integer, the degree of the harmonic
  !        m - integer, the order of the harmonic
  !        theta - real, the polar angle in radians
  !        phi - real, the azimuthal angle in radians
  ! Output: result - real, the value of the spherical harmonic function
  !***********************************************************************
  REAL(c_double) FUNCTION sphericalHarmonic(l, m, theta, phi)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: l, m
    REAL(c_double), INTENT(IN) :: theta, phi
    REAL(c_double) :: P_lm, normalization

    ! Compute the associated Legendre polynomial P_lm
    P_lm = legendreP(l, m, cos(theta))

    ! Compute the normalization factor
    normalization = sqrt((2.0 * l + 1.0) / (4.0 * 3.14159265358979323846) * &
                         factorial(l - m) / factorial(l + m))

    ! Compute the real-valued spherical harmonic
    IF (m > 0) THEN
      sphericalHarmonic = sqrt(2.0) * normalization * P_lm * cos(m * phi)
    ELSE IF (m < 0) THEN
      sphericalHarmonic = sqrt(2.0) * normalization * P_lm * sin(abs(m) * phi)
    ELSE
      sphericalHarmonic = normalization * P_lm
    END IF
  END FUNCTION sphericalHarmonic

  !***********************************************************************
  ! Function: dSphericalHarmonicDTheta
  ! Description: Computes the derivative of spherical harmonic with respect to theta
  !***********************************************************************
  REAL(c_double) FUNCTION dSphericalHarmonicDTheta(l, m, theta, phi)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: l, m
    REAL(c_double), INTENT(IN) :: theta, phi
    REAL(c_double) :: P_lm, dP_lm, normalization

    ! Get the associated Legendre polynomial and its derivative
    P_lm = legendreP(l, m, cos(theta))
    dP_lm = dLegendrePDTheta(l, m, theta)
    
    ! Compute the normalization factor
    normalization = sqrt((2.0 * l + 1.0) / (4.0 * 3.14159265358979323846) * &
                         factorial(l - m) / factorial(l + m))

    IF (m > 0) THEN
      dSphericalHarmonicDTheta = sqrt(2.0) * normalization * dP_lm * cos(m * phi)
    ELSE IF (m < 0) THEN
      dSphericalHarmonicDTheta = sqrt(2.0) * normalization * dP_lm * sin(abs(m) * phi)
    ELSE
      dSphericalHarmonicDTheta = normalization * dP_lm
    END IF
  END FUNCTION dSphericalHarmonicDTheta

  !***********************************************************************
  ! Function: dSphericalHarmonicDPhi
  ! Description: Computes the derivative of spherical harmonic with respect to phi
  !***********************************************************************
  REAL(c_double) FUNCTION dSphericalHarmonicDPhi(l, m, theta, phi)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: l, m
    REAL(c_double), INTENT(IN) :: theta, phi
    REAL(c_double) :: P_lm, normalization

    P_lm = legendreP(l, m, cos(theta))
    normalization = sqrt((2.0 * l + 1.0) / (4.0 * 3.14159265358979323846) * &
                         factorial(l - m) / factorial(l + m))

    IF (m > 0) THEN
      dSphericalHarmonicDPhi = -m * sqrt(2.0) * normalization * P_lm * sin(m * phi)
    ELSE IF (m < 0) THEN
      dSphericalHarmonicDPhi = abs(m) * sqrt(2.0) * normalization * P_lm * cos(abs(m) * phi)
    ELSE
      dSphericalHarmonicDPhi = 0.0_c_double
    END IF
  END FUNCTION dSphericalHarmonicDPhi

  !***********************************************************************
  ! Function: planeHarmonic
  ! Description: Computes the real-valued plane harmonic function.
  ! Input: l - integer, the degree of the harmonic
  !        m - integer, the order of the harmonic (not used in this implementation)
  !        x - real, the x-coordinate
  !        y - real, the y-coordinate
  ! Output: result - real, the value of the plane harmonic function
  !***********************************************************************
  REAL(c_double) FUNCTION planeHarmonic(l, m, x, y)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: l, m
    REAL(c_double), INTENT(IN) :: x, y
    REAL(c_double) :: r, theta

    ! Compute polar coordinates
    r = sqrt(x*x + y*y)
    theta = atan2(y, x)

    ! Compute the real-valued plane harmonic
    ! Using r^l * cos(l*theta) as our harmonic function
    IF (r > 0.0_c_double) THEN
      planeHarmonic = r**l * cos(l*theta)
    ELSE
      planeHarmonic = 0.0_c_double
    END IF
  END FUNCTION planeHarmonic

  !***********************************************************************
  ! Function: dPlaneHarmonicDX
  ! Description: Computes the derivative of plane harmonic with respect to x
  !***********************************************************************
  REAL(c_double) FUNCTION dPlaneHarmonicDX(l, m, x, y)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: l, m
    REAL(c_double), INTENT(IN) :: x, y
    REAL(c_double) :: r, theta, r_pow

    r = sqrt(x*x + y*y)
    theta = atan2(y, x)

    IF (r > 0.0_c_double) THEN
      r_pow = r**(l-1)
      dPlaneHarmonicDX = l*x/r * r_pow * cos(l*theta) - &
                         l*y/(r*r) * r_pow * sin(l*theta)
    ELSE
      dPlaneHarmonicDX = 0.0_c_double
    END IF
  END FUNCTION dPlaneHarmonicDX

  !***********************************************************************
  ! Function: dPlaneHarmonicDY
  ! Description: Computes the derivative of plane harmonic with respect to y
  !***********************************************************************
  REAL(c_double) FUNCTION dPlaneHarmonicDY(l, m, x, y)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: l, m
    REAL(c_double), INTENT(IN) :: x, y
    REAL(c_double) :: r, theta, r_pow

    r = sqrt(x*x + y*y)
    theta = atan2(y, x)

    IF (r > 0.0_c_double) THEN
      r_pow = r**(l-1)
      dPlaneHarmonicDY = l*y/r * r_pow * cos(l*theta) + &
                         l*x/(r*r) * r_pow * sin(l*theta)
    ELSE
      dPlaneHarmonicDY = 0.0_c_double
    END IF
  END FUNCTION dPlaneHarmonicDY

  !***********************************************************************
  ! Function: legendreP
  ! Description: Computes the associated Legendre polynomial P_lm(x).
  ! Input: l - integer, the degree of the polynomial
  !        m - integer, the order of the polynomial
  !        x - real, the argument of the polynomial
  ! Output: result - real, the value of the associated Legendre polynomial
  !***********************************************************************
  REAL(c_double) FUNCTION legendreP(l, m, x)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: l, m
    REAL(c_double), INTENT(IN) :: x
    INTEGER :: i
    REAL(c_double) :: pmm, pmmp1, pll

    IF (m == 0 .AND. l == 0) THEN
      legendreP = 1.0
      RETURN
    END IF

    pmm = 1.0
    IF (m > 0) THEN
      DO i = 1, m
        pmm = pmm * (-sqrt(1.0 - x * x)) * (2.0 * i - 1.0) / i
      END DO
    END IF

    IF (l == m) THEN
      legendreP = pmm
      RETURN
    END IF

    pmmp1 = x * (2.0 * m + 1.0) * pmm
    IF (l == m + 1) THEN
      legendreP = pmmp1
      RETURN
    END IF

    DO i = m + 2, l
      pll = ((2.0 * i - 1.0) * x * pmmp1 - (i + m - 1.0) * pmm) / (i - m)
      pmm = pmmp1
      pmmp1 = pll
    END DO

    legendreP = pll
  END FUNCTION legendreP

  !***********************************************************************
  ! Function: dLegendrePDTheta
  ! Description: Computes the derivative of associated Legendre polynomial
  !              with respect to theta
  !***********************************************************************
  REAL(c_double) FUNCTION dLegendrePDTheta(l, m, theta)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: l, m
    REAL(c_double), INTENT(IN) :: theta
    REAL(c_double) :: sin_theta, cos_theta

    sin_theta = sin(theta)
    cos_theta = cos(theta)
    
    ! The derivative with respect to theta is related to the derivative
    ! with respect to x = cos(theta) by the chain rule:
    ! dP/dtheta = -(sin theta) * dP/dx
    dLegendrePDTheta = -sin_theta * dLegendrePDX(l, m, cos_theta)
  END FUNCTION dLegendrePDTheta

  !***********************************************************************
  ! Function: dLegendrePDX
  ! Description: Computes the derivative of associated Legendre polynomial
  !              with respect to x
  !***********************************************************************
  REAL(c_double) FUNCTION dLegendrePDX(l, m, x)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: l, m
    REAL(c_double), INTENT(IN) :: x
    
    IF (l == m) THEN
      dLegendrePDX = l * x * legendreP(l, m, x) / (x * x - 1.0_c_double)
    ELSE
      dLegendrePDX = l * x * legendreP(l, m, x) / (x * x - 1.0_c_double) - &
                     (l + m) * legendreP(l-1, m, x) / (x * x - 1.0_c_double)
    END IF
  END FUNCTION dLegendrePDX

  !***********************************************************************
  ! Function: factorial
  ! Description: Computes the factorial of a non-negative integer n.
  ! Input: n - integer, the number to compute the factorial of
  ! Output: result - real, the factorial of n
  !***********************************************************************
  REAL(c_double) FUNCTION factorial(n)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    INTEGER :: i

    factorial = 1.0
    DO i = 2, n
      factorial = factorial * i
    END DO
  END FUNCTION factorial

END MODULE generalHarmonicFunctions_m
