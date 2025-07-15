!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA.Analytic.F90
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring 
!                     Warning and Forecasting (GBA-MWF) Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : V 0.0
! Development docs  : "A Mixed Solver for Stream Function and Velocity Potential with Velocity Components
!                     as Boundary Conditions and Its Adjoint Solver.docs"
! HISTORY           :
!   Created by Yuanfu Xie on 2024-12-09
!   Updated by Yuanfu Xie on 2025-03-05
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This routine constructs spherical linear functions for Psi, Chi of a mixed solver of U and V
!! @copyright (C) 2020 GBA-MWF, All rights reserved.
! @note
! @warning
SUBROUTINE addLinearSolutions(psi, chi, u, v, sg)
  USE, INTRINSIC :: iso_c_binding, ONLY: c_double, c_int
  USE parameters_m, ONLY: EarthRadius
  USE singleGrid_m, ONLY: singleGrid_t
  USE State_m, ONLY: State_t
  USE State2NC_m
  USE YAMLRead_m

  IMPLICIT NONE

  TYPE(singleGrid_t), INTENT(IN) :: sg
  ! u and v are the velocity components to meet and 
  ! uwind and vwind are velocity components computed from the input psi and chi
  REAL(c_double), DIMENSION(sg%vlevel,sg%num_cell), INTENT(IN) :: u,v
  REAL(c_double), DIMENSION(sg%vlevel,sg%num_cell), INTENT(INOUT) :: psi,chi

  ! Local variables:
  INTEGER(c_int) :: i,j,nintegral
  REAL(c_double) :: coef(2,2),determinant,dll(2),matrix(2,2),r(2,2),cs,cs2,lon
  REAL(c_double), ALLOCATABLE :: uwind(:,:),vwind(:,:)
  REAL(c_double), ALLOCATABLE :: integral(:)

  ! For debugging plotting:
  CHARACTER(LEN=1024) :: yamlFile
  TYPE(state_t) :: states


  ! Plotting:----------------------------------------------------------------
  ! Get the yaml file:
  CALL getarg(1, yamlFile)
  ! Initialize the state:
  CALL states%initialize(yamlFile, sg)
  ! End plotting:------------------------------------------------------------

  ! Determining the general solution:

  ! Find an interior point to calculate the latlon increments for latlon grid:
  DO i=1,sg%num_cell
    IF (sg%cell_type(i) .EQ. 0) EXIT
  END DO
  dll(1) = sg%cell_cntr(1,sg%cell_adnb(3,i))-sg%cell_cntr(1,i)
  dll(2) = sg%cell_cntr(2,sg%cell_adnb(2,i))-sg%cell_cntr(2,i)

  ! Calculate the velocity at the cell center from input psi and chi:
  ALLOCATE(uwind(sg%vlevel,sg%num_cell),vwind(sg%vlevel,sg%num_cell))
  CALL uvVelocityOnInterior(psi, chi, uwind, vwind, sg)

  ALLOCATE(integral(sg%num_cell))
  integral = 0.0D0

  matrix = 0.0D0; r = 0.0D0
  DO i=1,sg%num_cell

    ! Integral of 1/cos(theta): THIS WORKS ONLY FOR Lat-lon grid
    nintegral = NINT((sg%cell_cntr(1,i)-sg%cell_cntr(1,1))/dll(1))
    integral(i) = 0.0D0
    DO j=1,nintegral
      ! Starting the interior cell south boundary:
      integral(i) = integral(i) + dll(1)/DCOS(sg%cell_cntr(1,1)+(DBLE(j)-0.5D0)*dll(1))

      ! Plotting the integral
      states%fields(1)%DATA(1,i,3) = integral(i)
    END DO

    ! For all interior point to find the least square solution determining the general solution:
    IF (sg%cell_type(i) .EQ. 1) THEN

      ! The corrent coefficient matrix:
      cs = DCOS(sg%cell_cntr(1,i))
      cs2 = cs**2
      lon = sg%cell_cntr(2,i)
      matrix(1,1) = matrix(1,1) + 1.0D0/cs2
      matrix(1,2) = matrix(1,2) + lon/cs2
      matrix(2,2) = matrix(2,2) + lon**2/cs2
      r(1,1) = r(1,1) +     (u(1,i)-uwind(1,i))/cs
      r(2,1) = r(2,1) - lon*(u(1,i)-uwind(1,i))/cs
      r(1,2) = r(1,2) +     (v(1,i)-vwind(1,i))/cs
      r(2,2) = r(2,2) + lon*(v(1,i)-vwind(1,i))/cs
    END IF
  END DO

  ! Correct form:
  determinant = matrix(1,1)*matrix(2,2) - matrix(1,2)*matrix(1,2)
  coef(1,1) = (matrix(2,2)*r(1,1)+matrix(1,2)*r(2,1))/determinant
  coef(1,2) = (matrix(1,2)*r(1,1)+matrix(1,1)*r(2,1))/determinant
  coef(2,1) = (matrix(2,2)*r(1,2)-matrix(1,2)*r(2,2))/determinant
  coef(2,2) =(-matrix(1,2)*r(1,2)+matrix(1,1)*r(2,2))/determinant

  ! Add the general solutions to both psi and chi:
  DO i=1,sg%num_cell
    psi(:,i) = psi(:,i) + EarthRadius*sg%cell_cntr(2,i)*(coef(2,1)+coef(1,2)*integral(i))
    chi(:,i) = chi(:,i) + EarthRadius*sg%cell_cntr(2,i)*(coef(1,1)+coef(2,2)*integral(i))
  END DO

  DEALLOCATE(uwind, vwind, integral)
END SUBROUTINE addLinearSolutions