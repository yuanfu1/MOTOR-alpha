! Description:
!> @file
!!   Calculate Legendre coefficients for a given phase function.
!
!> @brief
!!   Calculate Legendre coefficients for a given phase function.
!!
!! @details
!!   Calculate the first nmom Legendre coefficients for the given
!!   phase function using Gaussian quadrature.
!!
!!   The default quadrature size is 1000 points. This can be changed
!!   by supplying the ngauss argument. For greater efficiency when
!!   making multiple calls to this subroutine it is possible to
!!   calculate the Gaussian quadrature points and weights outside of
!!   this subroutine and pass them in. This can be done using the
!!   gauss_quad subroutine contained in rttov_scattering_mod as follows:
!!
!!   CALL gauss_quad(-1._jprb, 1._jprb, q(1:ngauss), w(1:ngauss))
!!
!!   On output the REAL arrays q and w contain the quadrature points and
!!   weights respectively.
!!
!! @param[out]  err       status on exit
!! @param[in]   pha       phase function values at angles in phangle array
!! @param[in]   phangle   array of angles on which phase function is defined
!! @param[in]   nmom      number of Legendre coefficients to calculate
!! @param[out]  legcoef   calculated Legendre coefficients
!! @param[in]   ngauss    size of Gaussian quadrature to use, optional, default 1000
!! @param[in]   q         input precomputed Gaussian quadrature points, optional
!! @param[in]   w         input precomputed Gaussian quadrature weights, optional
!
! Copyright:
!    This software was developed within the context of
!    the EUMETSAT Satellite Application Facility on
!    Numerical Weather Prediction (NWP SAF), under the
!    Cooperation Agreement dated 7 December 2016, between
!    EUMETSAT and the Met Office, UK, by one or more partners
!    within the NWP SAF. The partners in the NWP SAF are
!    the Met Office, ECMWF, DWD and MeteoFrance.
!
!    Copyright 2016, EUMETSAT, All Rights Reserved.
!
SUBROUTINE rttov_legcoef_calc(err, pha, phangle, nmom, legcoef, ngauss, q, w)
!INTF_OFF
#include "throw.h"
!INTF_ON
  USE parkind1, ONLY : jpim, jprb

!INTF_OFF
  USE rttov_const, ONLY : deg2rad
  USE rttov_scattering_mod, ONLY : &
    gauss_quad,                    &
    spline_interp,                 &
    normalise,                     &
    calc_legendre_coef_gauss
!INTF_ON

  IMPLICIT NONE

  INTEGER(jpim), INTENT(OUT)            :: err
  REAL(jprb),    INTENT(IN)             :: pha(:)
  REAL(jprb),    INTENT(IN)             :: phangle(SIZE(pha))
  INTEGER(jpim), INTENT(IN)             :: nmom
  REAL(jprb),    INTENT(INOUT)          :: legcoef(:)
  INTEGER(jpim), INTENT(IN),   OPTIONAL :: ngauss
  REAL(jprb),    INTENT(IN),   OPTIONAL :: q(:), w(:)
!INTF_END

#include "rttov_errorreport.interface"

  INTEGER(jpim)           :: ngauss1, nphangle, nmom1
  REAL(jprb), ALLOCATABLE :: q1(:), w1(:), pha_interp(:)

  TRY

  nphangle = SIZE(phangle)

  IF (PRESENT(q)) THEN
    ! If quadrature is supplied check arguments and make use of it
    IF (.NOT. PRESENT(w)) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0, 'w argument is mandatory if q argument is supplied')
    ENDIF
    IF (SIZE(q) /= SIZE(w)) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0, 'q and w arguments must be the same size')
    ENDIF

    ngauss1 = SIZE(q)
    ALLOCATE(q1(ngauss1), w1(ngauss1), pha_interp(ngauss1), stat=err)
    THROWM(err.NE.0, 'allocation failure')

    q1 = q
    w1 = w
  ELSE
    ! Otherwise set size of Gaussian quadrature and calculate it
    ngauss1 = 1000
    IF (PRESENT(ngauss)) ngauss1 = ngauss

    ALLOCATE(q1(ngauss1), w1(ngauss1), pha_interp(ngauss1), stat=err)
    THROWM(err.NE.0, 'allocation failure')
    CALL gauss_quad(-1._jprb, 1._jprb, q1, w1)
  ENDIF

  CALL spline_interp(nphangle, COS(phangle(nphangle:1:-1) * deg2rad), &
                     pha(nphangle:1:-1), ngauss1, q1, pha_interp)
  CALL normalise(ngauss1, w1, pha_interp)
  CALL calc_legendre_coef_gauss(q1, w1, pha_interp, nmom, nmom, nmom1, legcoef)

  DEALLOCATE(q1, w1, pha_interp, stat=err)
  THROWM(err.NE.0, 'deallocation failure')

  CATCH
END SUBROUTINE rttov_legcoef_calc
