
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
!   Updated by Yuanfu Xie on 2025-04-23 by defining the analytic function of Psi and Chi without sg
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This routine constructs analytic functions of Psi, Chi and voricity and divergence with U and V at
!  cell centers for testing a latitude and longitude grid.
!! @copyright (C) 2020 GBA-MWF, All rights reserved.
! @note
! @warning
SUBROUTINE LinearLatCase(npts,nlvl,psi,chi,vor,div,u,v,ll)
  USE, INTRINSIC :: iso_c_binding, ONLY: c_double, c_int
  USE parameters_m, ONLY: EarthRadius

  IMPLICIT NONE

  INTEGER(c_int), INTENT(IN) :: npts,nlvl
  REAL(c_double), INTENT(IN) :: ll(2,npts) ! longitude and latitude
  REAL(c_double), INTENT(OUT) :: psi(nlvl,npts),chi(nlvl,npts), &
    vor(nlvl,npts),div(nlvl,npts), &
    u(nlvl,npts), v(nlvl,npts)

  ! Local variables:
  INTEGER(c_int) :: i,info
  REAL(c_double) :: dlambda,ddtheta ! longitude and latitude increments

  ! Refer to https://en.wikipedia.org/wiki/Del_in_cylindrical_and_spherical_coordinates
  ! First to test a stream function as quadratic function of longitude only:
  psi = 0.0D0; chi = 0.0D0; vor = 0.0D0; div = 0.0D0; u = 0.0D0; v = 0.0D0
  ! Note: the analytic function of u and v are velocity in longitude and latitute directions at edges:
  ! Not the normal and tangential directions after the most recent modification:
  DO i=1,npts
    psi(:,i) = 0.5D0*EarthRadius*ll(2,i)**2
    vor(:,i) = 1.0D0/DCOS(ll(1,i))**2/EarthRadius
      v(:,i) = ll(2,i)/DCOS(ll(1,i))
  END DO
END SUBROUTINE LinearLatCase

! ! Check if the operators from Psi and Chi to U V are correct:
! SUBROUTINE checkPsiChi2UV_s(psi, chi, u, v, sg, errMax)
!   USE, INTRINSIC :: iso_c_binding, ONLY: c_double, c_int
!   USE SingleGrid_m, ONLY: SingleGrid_t

!   IMPLICIT NONE
!   REAL(r_kind), INTENT(IN) :: psi(nlvl,npts),chi(nlvl,npts)
!   REAL(r_kind), INTENT(IN) :: u(nlvl,npts),v(nlvl,npts)
!   REAL(r_kind), INTENT(OUT) :: errMax

!   ! Local variables:
!   INTEGER(i_kind) :: imx,i
!   REAL(r_kind), ALLOCATABLE :: uwind(:,:), vwind(:,:)

!   ! Allocate memory for computed wind:
!   ALLOCATE(uwind(nlvl,npts),vwind(nlvl,npts))
!   ! Convert to U and V:
!   CALL uvVelocityOnInterior(psi,chi,uwind,vwind,sg)

!   ! Find the max error in V component only (extension is needed for other test cases):
!   imx = 0
!   errMax = 0.0D0
!   DO i=1,npts
!     IF (sg%cell_type(i) .LE. 2) THEN
!       IF (ABS(v(1,i)-vwind(1,i)) .GT. errMax) THEN
!         imx = i
!         errMax = ABS(v(1,i)-vwind(1,i))
!       END IF
!     END IF
!   END DO
!   WRITE(*,1) errMax,imx,v(1,imx),vwind(1,imx)
! 1     FORMAT('V computed - anaV: ',E14.6,' at:',I6,' true/computed: ',2E14.6)
!   PRINT*,'Calculated U component: ',MAXVAL(ABS(uwind(1,:)))
! END SUBROUTINE checkPsiChi2UV_s