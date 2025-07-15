!***********************************************************************
! File: harmonicFunctions.F90
! Description: This file contains special harmonic functions.
! Author: [Yuanfu Xie]
! Date: [2025-03-18]
!***********************************************************************

module harmonicFunctions
  use, intrinsic :: iso_c_binding, only: c_double
  implicit none
  private
  public :: sphericalHarmonic, dYdTheta, dYdPhi

contains

  !***********************************************************************
  ! Function: sphericalHarmonic
  ! Description: Computes the real-valued spherical harmonic function Y_lm(theta, phi).
  ! Input: l - integer, the degree of the harmonic
  !        m - integer, the order of the harmonic
  !        theta - real, the polar angle in radians
  !        phi - real, the azimuthal angle in radians
  ! Output: result - real, the value of the spherical harmonic function
  !***********************************************************************
  real(c_double) function sphericalHarmonic(l, m, theta, phi)
    implicit none
    integer, intent(in) :: l, m
    real(c_double), intent(in) :: theta, phi
    real(c_double) :: P_lm, normalization

    ! Compute the associated Legendre polynomial P_lm
    P_lm = legendreP(l, m, cos(theta))

    ! Compute the normalization factor
    normalization = sqrt((2.0 * l + 1.0) / (4.0 * 3.14159265358979323846) * &
                         factorial(l - m) / factorial(l + m))

    ! Compute the real-valued spherical harmonic
    if (m > 0) then
      sphericalHarmonic = sqrt(2.0) * normalization * P_lm * cos(m * phi)
    else if (m < 0) then
      sphericalHarmonic = sqrt(2.0) * normalization * P_lm * sin(abs(m) * phi)
    else
      sphericalHarmonic = normalization * P_lm
    end if
  end function sphericalHarmonic

  !***********************************************************************
  ! Function: dYdTheta
  ! Description: Computes the derivative of the spherical harmonic function Y_lm(theta, phi) with respect to theta.
  ! Input: l - integer, the degree of the harmonic
  !        m - integer, the order of the harmonic
  !        theta - real, the polar angle in radians
  !        phi - real, the azimuthal angle in radians
  ! Output: result - real, the value of the derivative with respect to theta
  !***********************************************************************
  real(c_double) function dYdTheta(l, m, theta, phi)
    implicit none
    integer, intent(in) :: l, m
    real(c_double), intent(in) :: theta, phi
    real(c_double) :: P_lm, dP_lm_dTheta, normalization

    ! Compute the associated Legendre polynomial P_lm and its derivative
    P_lm = legendreP(l, m, cos(theta))
    dP_lm_dTheta = -sin(theta) * legendreP(l, m, cos(theta))

    ! Compute the normalization factor
    normalization = sqrt((2.0 * l + 1.0) / (4.0 * 3.14159265358979323846) * &
                         factorial(l - m) / factorial(l + m))

    ! Compute the derivative of the spherical harmonic with respect to theta
    if (m > 0) then
      dYdTheta = sqrt(2.0) * normalization * (dP_lm_dTheta * cos(m * phi) - m * P_lm * sin(m * phi) / tan(theta))
    else if (m < 0) then
      dYdTheta = sqrt(2.0) * normalization * (dP_lm_dTheta * sin(abs(m) * phi) + abs(m) * P_lm * cos(abs(m) * phi) / tan(theta))
    else
      dYdTheta = normalization * dP_lm_dTheta
    end if
  end function dYdTheta

  !***********************************************************************
  ! Function: dYdPhi
  ! Description: Computes the derivative of the spherical harmonic function Y_lm(theta, phi) with respect to phi.
  ! Input: l - integer, the degree of the harmonic
  !        m - integer, the order of the harmonic
  !        theta - real, the polar angle in radians
  !        phi - real, the azimuthal angle in radians
  ! Output: result - real, the value of the derivative with respect to phi
  !***********************************************************************
  real(c_double) function dYdPhi(l, m, theta, phi)
    implicit none
    integer, intent(in) :: l, m
    real(c_double), intent(in) :: theta, phi
    real(c_double) :: P_lm, normalization

    ! Compute the associated Legendre polynomial P_lm
    P_lm = legendreP(l, m, cos(theta))

    ! Compute the normalization factor
    normalization = sqrt((2.0 * l + 1.0) / (4.0 * 3.14159265358979323846) * &
                         factorial(l - m) / factorial(l + m))

    ! Compute the derivative of the spherical harmonic with respect to phi
    if (m > 0) then
      dYdPhi = -sqrt(2.0) * normalization * P_lm * m * sin(m * phi)
    else if (m < 0) then
      dYdPhi = sqrt(2.0) * normalization * P_lm * abs(m) * cos(abs(m) * phi)
    else
      dYdPhi = 0.0
    end if
  end function dYdPhi

  !***********************************************************************
  ! Function: legendreP
  ! Description: Computes the associated Legendre polynomial P_lm(x).
  ! Input: l - integer, the degree of the polynomial
  !        m - integer, the order of the polynomial
  !        x - real, the argument of the polynomial
  ! Output: result - real, the value of the associated Legendre polynomial
  !***********************************************************************
  real(c_double) function legendreP(l, m, x)
    implicit none
    integer, intent(in) :: l, m
    real(c_double), intent(in) :: x
    integer :: i
    real(c_double) :: pmm, pmmp1, pll

    if (m == 0 .and. l == 0) then
      legendreP = 1.0
      return
    end if

    pmm = 1.0
    if (m > 0) then
      do i = 1, m
        pmm = pmm * (-sqrt(1.0 - x * x)) * (2.0 * i - 1.0) / i
      end do
    end if

    if (l == m) then
      legendreP = pmm
      return
    end if

    pmmp1 = x * (2.0 * m + 1.0) * pmm
    if (l == m + 1) then
      legendreP = pmmp1
      return
    end if

    do i = m + 2, l
      pll = ((2.0 * i - 1.0) * x * pmmp1 - (i + m - 1.0) * pmm) / (i - m)
      pmm = pmmp1
      pmmp1 = pll
    end do

    legendreP = pll
  end function legendreP

  !***********************************************************************
  ! Function: factorial
  ! Description: Computes the factorial of a non-negative integer n.
  ! Input: n - integer, the number to compute the factorial of
  ! Output: result - real, the factorial of n
  !***********************************************************************
  real(c_double) function factorial(n)
    implicit none
    integer, intent(in) :: n
    integer :: i

    factorial = 1.0
    do i = 2, n
      factorial = factorial * i
    end do
  end function factorial

end module harmonicFunctions