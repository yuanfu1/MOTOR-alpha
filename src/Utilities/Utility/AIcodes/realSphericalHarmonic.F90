!!--------------------------------------------------------------------------------------------------
! PROJECT           : Application
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           : 2025-03-18   Created by Kimi AI chatBox.
! DESCRIPTION       : This is a simple Fortran program that computes the solution to the Laplace
!                     equation in spherical coordinates using spherical harmonics.
!                     The program defines a function to compute the associated Legendre polynomial
!                     P_l^m(x) and uses it to compute the solution for a given set of coefficients
!                     A(l, m) and B(l, m).
! USAGE             : Compile and run the program, then enter values for r, theta, and phi when prompted.
!                     The program will compute and print the solution.
! NOTES             : This program is for educational purposes and may need to be modified for specific
!                     applications.
! On 2025-03-18:
! Yuanfu Xie will pull the spherical harmonics part from this program to MOTOR as a component of a
! set of Harmonic functions.
!---------------------------------------------------------------------------------------------------
program laplace_solution
  implicit none
  integer, parameter :: max_l = 5, max_m = 5
  real, parameter :: pi = 3.141592653589793
  real :: r, theta, phi, solution
  real :: A(0:max_l, 0:max_m), B(0:max_l, 0:max_m)
  real :: dYdTheta
  integer :: l, m

  ! Initialize coefficients A and B (example values)
  do l = 0, max_l
      do m = 0, l
          A(l, m) = 1.0 / (l + m + 1)  ! Example values
          B(l, m) = 1.0 / (l + m + 2)  ! Example values
      end do
  end do

  ! Input values for r, theta, and phi
  print *, 'Enter r, theta (in radians), and phi (in radians):'
  read *, r, theta, phi

  ! Compute the solution
  solution = 0.0
  do l = 0, max_l
      do m = 0, l
          solution = solution + (A(l, m) * r**l + B(l, m) * r**(-(l+1))) * &
                     P_lm(l, m, cos(theta)) * cos(m * phi)
      end do
  end do

  ! Compute the derivative with respect to theta
  dYdTheta = 0.0
  do l = 0, max_l
      do m = 0, l
          dYdTheta = dYdTheta + (A(l, m) * r**l + B(l, m) * r**(-(l+1))) * &
                     dP_lm_dTheta(l, m, theta) * cos(m * phi)
      end do
  end do

  print *, 'Solution: ', solution
  print *, 'dYdTheta: ', dYdTheta

contains

  ! Function to compute associated Legendre polynomial P_l^m(x)
  real function P_lm(l, m, x)
      integer, intent(in) :: l, m
      real, intent(in) :: x
      real :: pmm, pmmp1, pll
      integer :: i

      if (m > l) then
          P_lm = 0.0
          return
      end if

      pmm = 1.0
      if (m > 0) then
          pmm = (1.0 - x**2)**(m/2.0)
          do i = 1, m
              pmm = -pmm * (2*i - 1) * sqrt(1.0 - x**2)
          end do
      end if

      pmmp1 = x * (2*m + 1) * pmm
      if (l == m) then
          P_lm = pmmp1
          return
      end if

      do i = m + 2, l
          pll = (x * (2*i - 1) * pmmp1 - (i + m - 1) * pmm) / (i - m)
          pmm = pmmp1
          pmmp1 = pll
      end do

      P_lm = pll
  end function P_lm

  ! Function to compute the derivative of the associated Legendre polynomial P_l^m(x) with respect to theta
  real function dP_lm_dTheta(l, m, theta)
      integer, intent(in) :: l, m
      real, intent(in) :: theta
      real :: x

      x = cos(theta)
      dP_lm_dTheta = -sin(theta) * P_lm(l, m, x)
  end function dP_lm_dTheta

end program laplace_solution