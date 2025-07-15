!
SUBROUTINE RTTOV_SUBLAYER_TL( &
            & z1,    &
            & z1_tl, &
            & z2,    &
            & z2_tl, &
            & x1,    &
            & x1_tl, &
            & x2,    &
            & x2_tl, &
            & w1,    &
            & w1_tl, &
            & w2,    &
            & w2_tl)
!
!-----------------------------------------------------------------------
!
!     Description:
!
!     Tangent linear code (by Niels Bormann) for the segment of the
!     piecewise weighted averaging interpolator described in the next
!     paragraph.
!
!     The module RTTOV_SUBLAYER determines weight coefficient contributions
!     w1 and w2 to assign to input domain (e.g. NWP model) variables at
!     x1 and x2. Weights are determined from integration over the
!     intersecting segment (y1,y2) of the ranges (x1,x2) for the input
!     domain and (z1,z2) for the output domain. Integrals are approximated
!     via the trapezoidal rule:
!
!         integral of f(x)=w(x)*t(x) from y1 to y2
!
!                                      = (f(y1)+f(y2))/2*(y2-y1)
!                                      = w(y1)*t(y1)+w(y2)*t(y2)
!                                      = w1*t(x1)+w2*t(x2)
!                                      = w1*t1+w2*t2
!
!     This is synonomous to having an integrand linear in x.
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
!    Copyright 2010, EUMETSAT, All Rights Reserved.
!
!     Method:
!     - Piecewise weighted interpolation in ln(P).
!     - Journal reference:
!       Rochon, Y.J., L. Garand, D.S. Turner, and S. Polavarapu.
!       Q. J. of the Royal Met. Soc., vol 133, 1547-1558, 2007.
!       (www.interscience.wiley.com) DOI:10.1002/qj.117
!
!     History:
!
!     The original version of this routine (sublayer.f90) was developed
!     and tested by Yves J. Rochon with contributions from L. Garand,
!     D.S. Turner, and S. Polavarapu of the Atmospheric Science and
!     Technology Directorate, Environment Canada (EC).
!
!     Copyright 2006, Environment Canada, All Rights Reserved.
!
!     Please provide appropriate acknowledgement through referencing
!     of the corresponding journal article (see above), and of the
!     RTTOV-9 Science and Validation Report (in preparation), available
!     from the NWP SAF website at
!     http://www.metoffice.gov.uk/research/interproj/nwpsaf/rtm/rtm_rttov9.html
!
!     Version   Date        Comment
!     -------   ----        -------
!     1         10/2005     Original F77 code by Y.J. Rochon with
!                           contributions from S. Polavarapu. (EC)
!     2         10/2006     Adaptation for F90 and for improved consistency
!                           with RTTOV code format (Y.J. Rochon, EC).
!                           Argument names of original code were not changed.
!     ----------------------------------------------------------------------
!                           Application within RTTOV (for NWP SAF)
!
!     3         12/2006     Some restructuring. (Peter Rayer)
!     4         01/2008     Tangent linear version for TL of input/output
!                           levels (lgradp option), Niels Bormann
!     5         02/2008     Updating of header and comments segments.
!                           (Y.J. Rochon, EC)
!     6         04/2008     Commenting out of tests for zero values
!                           done for speed optimisation (lines preceeded
!                           by !!$). (D. Salmond and M. Hamrud, ECMWF)
!
!     Arguments: (reflects changes under item 4 above)
!
!     Input
!
!     z1.......................Outer boundary of output domain half-layer
!                              (above or below z2)
!     z2.......................Inner boundary of output domain half-layer
!                              (position of reference output domain level)
!     z3.......................Calling routine only - here use z2
!
!     x1.......................Upper boundary of input domain layer (x1<x2)
!     x2.......................Lower boundary from input domain layer
!     z1_tl, z2_tl, x1_tl, x2_tl ....... Tangent linear of above
!
!     Output
!
!     w1.......................Weight assigned to variable at upper input level
!     w2.......................Weight assigned to variable at lower input level
!     w1_tl, w2_tl ............Tangent linear of above
!
!     Other
!
!     y1......................Upper boundary of integral range (y1<y2)
!     y2......................Lower boundary of integral range
!     errorstatus.............Returns error status to allow soft failure
!
!     External functions and subroutines: None
!
!     Related calling routine: RTTOV_LAYERAVG_TL
!
!     Assumptions:
!
!     1) x1<x2
!
!     2) The ranges (z1,z2) and (x1,x2) overlap. The overlap region
!     will be identified as (y1,y2) with y1<y2.
!
!     Comments: (related to above description section)
!
!     1) w(y1) and w(y2) are obtained by linear interpolation of the linear
!     weighting function w with w(z1)=0 and w(z2)= 1.
!
!     2) The w1 and w2 above are determined by expanding t(y1) and t(y2)
!     as linear functions of t(x1)=t1 and t(x2)=t2.
!
!     3) The factor of 1/2 in
!
!                          (f(y1)+f(y2))/2*(y2-y1)
!                           = w(y1)*t(y1)+w(y2)*t(y2)
!
!     is omitted as normalization is performed in the calling routine
!     RTTOV_LAYERAVG.
!
!     4) See RTTOV_SUBLAYER for a detailled description of steps reflected
!     in the condensed representation of the integration.
!
!--------------------------------------------------------------------
!
!
  USE parkind1, ONLY : jprb
!INTF_OFF
  USE parkind1, ONLY : jpim
!INTF_ON
  IMPLICIT NONE
!
!     --- Subroutine arguments ---------------------
!
  REAL(KIND=jprb), INTENT(IN)  :: z1   , z2   , z1_tl, z2_tl
  REAL(KIND=jprb), INTENT(IN)  :: x1   , x1_tl, x2   , x2_tl
  REAL(KIND=jprb), INTENT(OUT) :: w1   , w2
  REAL(KIND=jprb), INTENT(OUT) :: w1_tl, w2_tl
!INTF_END
!
!     --- Local scalars and arrays -----------------
!
  REAL   (KIND=jprb) :: y1    , y2    , d, w10 , w20  , dz , dz_tl , dx, dy, dzd, dzd_tl
  REAL   (KIND=jprb) :: y1_tl , y2_tl , dx_tl, d_tl, dy_tl, dxd, dxd_tl
  REAL   (KIND=jprb) :: w10_tl, w20_tl
  INTEGER(KIND=jpim) :: ibot  , itop
!
!     -------------------------------------------------
!     1. Identify and set upper and lower boundaries of
!        integration/averaging layers (y1 and y2) and
!        related perturbations.
!     -------------------------------------------------
!
  itop = 0
  ibot = 0
  IF (z1 .LT. z2) THEN
    y1    = z1
    y1_tl = z1_tl
    IF (x1 .GT. z1) THEN
      y1    = x1
      y1_tl = x1_tl
      itop  = 1
    ENDIF
    y2    = z2
    y2_tl = z2_tl
    IF (x2 .LT. z2) THEN
      y2    = x2
      y2_tl = x2_tl
      ibot  = 1
    ENDIF
  ELSE
    y1    = z2
    y1_tl = z2_tl
    IF (x1 .GT. z2) THEN
      y1    = x1
      y1_tl = x1_tl
      itop  = 1
    ENDIF
    y2    = z1
    y2_tl = z1_tl
    IF (x2 .LT. z1) THEN
      y2    = x2
      y2_tl = x2_tl
      ibot  = 1
    ENDIF
  ENDIF
!
!     ----------------------
!     2. Set weights for TLM
!     ----------------------
!
  dy     = y2 - y1
  dy_tl  = y2_tl - y1_tl
  dz     = z1 - z2
  dz_tl  = z1_tl - z2_tl
! Thickness of output layers all non-zero - divisor dz is positive
  dzd    = 1.0_JPRB / dz
  dzd_tl =  - 1. * dzd ** 2 * dz_tl
  w10    = (z1 - y1) * dzd * dy
  w10_tl = (z1_tl - y1_tl) * dzd * dy + (z1 - y1) * dzd_tl * dy + (z1 - y1) * dzd * dy_tl
  w20    = (z1 - y2) * dzd * dy
  w20_tl = (z1_tl - y2_tl) * dzd * dy + (z1 - y2) * dzd_tl * dy + (z1 - y2) * dzd * dy_tl
  dx     = (x2 - x1)
  dx_tl  = (x2_tl - x1_tl)
! Thickness of input layers all non-zero - divisor dx is positive
  dxd    = 1.0_JPRB / dx
  dxd_tl =  - dxd ** 2 * dx_tl
  d = (x2 - z2) * dxd
  d_tl   = (x2_tl - z2_tl) * dxd + (x2 - z2) * dxd_tl
  IF (z1 .LT. z2 .AND. ibot .EQ. 0) THEN
    w1    = w10 + w20 * d
    w1_tl = w10_tl + w20_tl * d + w20 * d_tl
    w2    = w20 * (1.0_JPRB - d)
    w2_tl = w20_tl * (1.0_JPRB - d) - w20 * d_tl
  ELSE IF (z1 .GT. z2 .AND. itop .EQ. 0) THEN
    w2    = w20 + w10 * (1.0_JPRB - d)
    w2_tl = w20_tl + w10_tl * (1.0_JPRB - d) - w10 * d_tl
    w1    = w10 * d
    w1_tl = w10_tl * d + w10 * d_tl
  ELSE
    w1    = w10
    w1_tl = w10_tl
    w2    = w20
    w2_tl = w20_tl
  ENDIF
!
!
END SUBROUTINE RTTOV_SUBLAYER_TL
