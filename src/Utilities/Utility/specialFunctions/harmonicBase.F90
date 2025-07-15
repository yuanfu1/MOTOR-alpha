!!--------------------------------------------------------------------------------------------------
! PROJECT         : MOTOR-DA.possionSolver_Kp.psMatrix
! AFFILIATION     : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring 
!                   Warning and Forecasting (GBA-MWF) Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Copilot on VScode initially on 2025-03-19
!   Modified by Yuanfu Xie, on 2025-03-20
!     for abstracting various harmonic base function with a polymorphism interface
!!--------------------------------------------------------------------------------------------------
!> @brief
!! This module contains a solver of streamfunction and velocity potential from vorticity and divergence
!! with u and v as boundary conditions.
!! @copyright (C) 2024 GBA-MWF, All rights reserved.
!! @note
!! @warning
!! @attention
MODULE harmonicBase_m
  USE, INTRINSIC :: iso_c_binding, ONLY: c_double, c_int
  ! USE generalHarmonicFunctions_m, ONLY: harmonicFunction !,sphericalHarmonic
  IMPLICIT NONE

  TYPE, ABSTRACT :: harmonics_t
  INTEGER(c_int) :: l, m
  REAL(c_double) :: value

  CONTAINS
    PROCEDURE(forward), DEFERRED :: harmonicFunction
    PROCEDURE(forward), DEFERRED :: dHarmonicDx
    PROCEDURE(forward), DEFERRED :: dHarmonicDy
  END TYPE harmonics_t

  ! Abstract interface for general harmonic functions
  ABSTRACT INTERFACE
  REAL(c_double) FUNCTION forward(this, l, m, x, y)
      IMPORT :: c_double, c_int, harmonics_t
      CLASS(harmonics_t) :: this
      INTEGER(c_int), INTENT(IN) :: l, m
      REAL(c_double), INTENT(IN) :: x, y
  END FUNCTION forward
  END INTERFACE

  TYPE, EXTENDS(harmonics_t) :: sphericalHarmonics_t
    CONTAINS
      PROCEDURE, PUBLIC :: harmonicFunction => sphericalHarmonic
      PROCEDURE, PUBLIC :: dHarmonicDx => dSphericalHarmonicDTheta
      PROCEDURE, PUBLIC :: dHarmonicDy => dSphericalHarmonicDPhi
  END TYPE sphericalHarmonics_t

  TYPE, EXTENDS(harmonics_t) :: planeHarmonics_t
    CONTAINS
      PROCEDURE, PUBLIC :: harmonicFunction => planeHarmonic
      PROCEDURE, PUBLIC :: dHarmonicDx => dPlaneHarmonicDx
      PROCEDURE, PUBLIC :: dHarmonicDy => dPlaneHarmonicDy
  END TYPE planeHarmonics_t

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
  REAL(c_double) FUNCTION sphericalHarmonic(this, l, m, x, y) ! replaced by x, y for theta, phi)
    IMPLICIT NONE
    CLASS(sphericalHarmonics_t) :: this
    INTEGER, INTENT(IN) :: l, m
    REAL(c_double), INTENT(IN) :: x, y ! replaced by x, y for theta, phi
    REAL(c_double) :: P_lm, normalization

    ! Compute the associated Legendre polynomial P_lm
    P_lm = legendreP(l, m, DCOS(x))

    ! Compute the normalization factor
    normalization = DSQRT((2.0_c_double * l + 1.0_c_double) / (4.0_c_double * 3.14159265358979323846_c_double) * &
                         factorial(l - m) / factorial(l + m))

    ! Compute the real-valued spherical harmonic
    IF (m > 0) THEN
      sphericalHarmonic = DSQRT(2.0_c_double) * normalization * P_lm * DCOS(m * y)
    ELSE IF (m < 0) THEN
      sphericalHarmonic = DSQRT(2.0_c_double) * normalization * P_lm * DSIN(ABS(m) * y)
    ELSE
      sphericalHarmonic = normalization * P_lm
    END IF
  END FUNCTION sphericalHarmonic

  !***********************************************************************
  ! Function: dSphericalHarmonicDTheta
  ! Description: Computes the derivative of spherical harmonic with respect to theta
  !***********************************************************************
  REAL(c_double) FUNCTION dSphericalHarmonicDTheta(this, l, m, x, y)
    IMPLICIT NONE
    CLASS(sphericalHarmonics_t) :: this
    INTEGER, INTENT(IN) :: l, m
    REAL(c_double), INTENT(IN) :: x, y
    REAL(c_double) :: P_lm, dP_lm, normalization

    ! Get the associated Legendre polynomial and its derivative
    P_lm = legendreP(l, m, DCOS(x))
    dP_lm = dLegendrePDTheta(l, m, x)
    
    ! Compute the normalization factor
    normalization = DSQRT((2.0_c_double * l + 1.0_c_double) / (4.0_c_double * 3.14159265358979323846_c_double) * &
                         factorial(l - m) / factorial(l + m))

    IF (m > 0) THEN
      dSphericalHarmonicDTheta = DSQRT(2.0_c_double) * normalization * dP_lm * DCOS(m * y)
    ELSE IF (m < 0) THEN
      dSphericalHarmonicDTheta = DSQRT(2.0_c_double) * normalization * dP_lm * DSIN(ABS(m) * y)
    ELSE
      dSphericalHarmonicDTheta = normalization * dP_lm
    END IF
  END FUNCTION dSphericalHarmonicDTheta

  !***********************************************************************
  ! Function: dSphericalHarmonicDPhi
  ! Description: Computes the derivative of spherical harmonic with respect to phi
  !***********************************************************************
  REAL(c_double) FUNCTION dSphericalHarmonicDPhi(this, l, m, x, y)
    IMPLICIT NONE
    CLASS(sphericalHarmonics_t) :: this
    INTEGER, INTENT(IN) :: l, m
    REAL(c_double), INTENT(IN) :: x, y
    REAL(c_double) :: P_lm, normalization

    P_lm = legendreP(l, m, DCOS(x))
    normalization = DSQRT((2.0_c_double * l + 1.0_c_double) / (4.0_c_double * 3.14159265358979323846_c_double) * &
                         factorial(l - m) / factorial(l + m))

    IF (m > 0) THEN
      dSphericalHarmonicDPhi = -m * DSQRT(2.0_c_double) * normalization * P_lm * DSIN(m * y)
    ELSE IF (m < 0) THEN
      dSphericalHarmonicDPhi = ABS(m) * DSQRT(2.0_c_double) * normalization * P_lm * DCOS(ABS(m) * y)
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
  REAL(c_double) FUNCTION planeHarmonic(this, l, m, x, y)
    IMPLICIT NONE
    CLASS(planeHarmonics_t) :: this
    INTEGER, INTENT(IN) :: l, m
    REAL(c_double), INTENT(IN) :: x, y
    REAL(c_double) :: r, theta

    ! Compute polar coordinates
    r = DSQRT(x*x + y*y)
    theta = DATAN2(y, x)

    ! Compute the real-valued plane harmonic
    ! Using r^l * DCOS(l*theta) as our harmonic function
    IF (r > 0.0_c_double) THEN
      planeHarmonic = r**l * DCOS(l*theta)
    ELSE
      planeHarmonic = 0.0_c_double
    END IF
  END FUNCTION planeHarmonic

  !***********************************************************************
  ! Function: dPlaneHarmonicDX
  ! Description: Computes the derivative of plane harmonic with respect to x
  !***********************************************************************
  REAL(c_double) FUNCTION dPlaneHarmonicDX(this, l, m, x, y)
    IMPLICIT NONE
    CLASS(planeHarmonics_t) :: this
    INTEGER, INTENT(IN) :: l, m
    REAL(c_double), INTENT(IN) :: x, y
    REAL(c_double) :: r, theta, r_pow

    r = DSQRT(x*x + y*y)
    theta = DATAN2(y, x)

    IF (r > 0.0_c_double) THEN
      r_pow = r**(l-1)
      dPlaneHarmonicDX = l*x/r * r_pow * DCOS(l*theta) - &
                         l*y/(r*r) * r_pow * DSIN(l*theta)
    ELSE
      dPlaneHarmonicDX = 0.0_c_double
    END IF
  END FUNCTION dPlaneHarmonicDX

  !***********************************************************************
  ! Function: dPlaneHarmonicDY
  ! Description: Computes the derivative of plane harmonic with respect to y
  !***********************************************************************
  REAL(c_double) FUNCTION dPlaneHarmonicDY(this, l, m, x, y)
    IMPLICIT NONE
    CLASS(planeHarmonics_t) :: this
    INTEGER, INTENT(IN) :: l, m
    REAL(c_double), INTENT(IN) :: x, y
    REAL(c_double) :: r, theta, r_pow

    r = DSQRT(x*x + y*y)
    theta = DATAN2(y, x)

    IF (r > 0.0_c_double) THEN
      r_pow = r**(l-1)
      dPlaneHarmonicDY = l*y/r * r_pow * DCOS(l*theta) + &
                         l*x/(r*r) * r_pow * DSIN(l*theta)
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
      legendreP = 1.0_c_double
      RETURN
    END IF

    pmm = 1.0_c_double
    IF (m > 0) THEN
      DO i = 1, m
        pmm = pmm * (-DSQRT(1.0_c_double - x * x)) * (2.0_c_double * i - 1.0_c_double) / i
      END DO
    END IF

    IF (l == m) THEN
      legendreP = pmm
      RETURN
    END IF

    pmmp1 = x * (2.0_c_double * m + 1.0_c_double) * pmm
    IF (l == m + 1) THEN
      legendreP = pmmp1
      RETURN
    END IF

    DO i = m + 2, l
      pll = ((2.0_c_double * i - 1.0_c_double) * x * pmmp1 - (i + m - 1.0_c_double) * pmm) / (i - m)
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
    cos_theta = DCOS(theta)
    
    ! The derivative with respect to theta is related to the derivative
    ! with respect to x = DCOS(theta) by the chain rule:
    ! dP/dtheta = -(DSIN theta) * dP/dx
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
END MODULE harmonicBase_m

