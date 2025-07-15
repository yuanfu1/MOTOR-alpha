!
SUBROUTINE RTTOV_SUBLAYER( &
            & z1, &
            & z2, &
            & x1, &
            & x2, &
            & w1, &
            & w2)
!
!-----------------------------------------------------------------------
!
!     Description:
!
!     Determine weight coefficient contributions w1 and w2 to assign
!     to input domain (e.g. NWP model) variables at x1 and x2. Weights
!     are determined from integration over the intersecting segment (y1,y2)
!     of the ranges (x1,x2) for the input domain and (z1,z2) for the
!     output domain. Integrals are approximated via the trapezoidal rule:
!
!         integral of f(x)=w(x)*t(x) from y1 to y2
!
!                                      = (f(y1)+f(y2))/2*abs(y2-y1)
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
!     4         01/2008     Removal of parts related to the calc of
!                           gradient contributions from sensitivity
!                           to changes in the vertical coordinate.
!                           This is now done in other routines.
!                           (Alan Geer and Neils Bormann, ECMWF)
!     5         02/2008     Updating of header and comments segments and
!                           removal of unused declared variables.
!                           Addition of expanded code version
!                           in comment section (item 4) which describes
!                           code used (in condensed form).
!                           (Y.J. Rochon, EC)
!     6         04/2008     Commenting out of tests for zero values
!                           done for speed optimization (lines preceeded
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
!
!     Output
!
!     w1.......................Weight assigned to variable at upper input level
!     w2.......................Weight assigned to variable at lower input level
!
!     Other
!
!     y1......................Upper boundary of integral range (y1<y2)
!     y2......................Lower boundary of integral range
!     errorstatus.............Returns error status to allow soft failure
!
!     External functions and subroutines: None
!
!     Related calling routine: RTTOV_LAYERAVG
!
!     Assumptions:
!
!     1) x1<x2
!
!     2) z1<z2 or z1>z2
!
!     3) The ranges (z1,z2) and (x1,x2) overlap. The overlap region
!     will be identified as (y1,y2) with y1<y2.
!
!     4) Weights w at (z1,z2) are (0,1).
!
!     Comments: (related to above description section)
!
!     1) w(y1) and w(y2) are obtained by linear interpolation of the linear
!     weighting function w with w(z1) and w(z2).
!
!     2) The w1 and w2 above are determined by expanding t(y1) and t(y2)
!     as linear functions of t(x1)=t1 and t(x2)=t2.
!
!     3) The factor of 1/2 in
!
!                          (f(y1)+f(y2))/2*abs(y2-y1)
!                           = w(y1)*t(y1)+w(y2)*t(y2)
!
!     is omitted as normalization is performed in the calling routine
!     RTTOV_LAYERAVG.
!
!     4) This version of the code was prepared considering the following
!     imposed conditions:
!
!        (w(z1),w(z2))=(0,1) implying  d1=(y1-x1)=0 from y1=x1 or
!                                      w(y1)=0 from w(x1)=0 and y1=z1.
!
!     5) A detailled description of steps reflected in the code but
!     for general (w(z1),w(z2)) pairs is provided below.
!
!     The lines beginning with !! denote the code in uncondensed form.
!
!     For this following description, if (z1>z2), switch z1 and z2
!     so that z1<z2.
!
!     Given y1, y2, dzd=1/(z2-z1), dxd=1/(x2-x1),
!
!     i) First determine averaging/interpolation weights w(y1) and w(y2) of
!
!                      f = w(y1)*t(y1)+w(y2)*t(y2)
!                        = wy1*t(y1)+wy2*t(y2)
!
!     for averaging/interpolation from t(y1) and t(y2) using the following
!     steps.
!
!     Weights are obtained by linear interpolation of the linear
!     weighting function w with w(z1)=wgt1 and w(z2)=wgt2.
!
!!      dzd=dzd*(wgt2-wgt1)
!!      wx1=wgt1+(y1-z1)*dzd
!!      wx2=wx1+(y2-y1)*dzd    !  wx2=wgt1+(y2-z1)*dzd
!
!!      dy=y2-y1
!!      wy1=dy*wx1
!!      wy2=dy*wx2
!
!     ii) Determine w1 and w2 of f=w1*t(x1)+w2*t(x2).
!
!     The w1 and w2 above are determined by expanding t(y1) and t(y2)
!     in f = w(y1)*t(y1)+w(y2)*t(y2) as linear functions of t(x1)=t1
!     and t(x2)=t2.
!
!           t(y1) = t1+(y1-x1)*(t2-t1)/(x2-x1) = t1+(y1-x1)*dxd*(t2-t1)
!           t(y2) = t2+(y2-x2)*(t2-t1)/(x2-x1) = t2+(y2-x2)*dxd*(t2-t1)
!
!     Therefore,
!
!           f = wy1*[t1+(y1-x1)*dxd*(t2-t1)]
!              +wy2*[t2+(y2-x2)*dxd*(t2-t1)]
!
!             = t1*[wy1-wy1*(y1-x1)*dxd-wy2*(y2-x2)*dxd]
!              +t2*[wy1*(y1-x1)*dxd+wy2+wy2*(y2-x2)*dxd]
!
!
!     Aside: w1+w2 = wy1+wy2 = 2 * w[(y2+y1)/2]
!
!!      dxy1=y1-x1
!!      d1=dxy1*dxd
!
!!      dxy2=y2-x2
!!      d2=dxy2*dxd
!
!!      w1=wy1*(1.0-d1)-wy2*d2
!!      w2=wy1+wy2-w1         ! w2=wy1*d1+wy2*(1.0+d2)
!
!--------------------------------------------------------------------
!
  USE parkind1, ONLY : jprb
!INTF_OFF
  USE parkind1, ONLY : jpim
!INTF_ON
  IMPLICIT NONE
!
!     --- Subroutine arguments ---------------------------
!
  REAL(KIND=jprb), INTENT(IN)  :: z1, z2, x1, x2
  REAL(KIND=jprb), INTENT(OUT) :: w1, w2
!INTF_END
!
!     --- Local scalars and arrays -----------------------
!
  REAL   (KIND=jprb) :: y1  , y2  , d, w10, w20, dz, dx, dy, dzd, dxd, dzddy
  INTEGER(KIND=jpim) :: ibot, itop
!
!     -------------------------------------------------
!     1. Identify and set upper and lower boundaries of
!        integration/averaging layers (y1 and y2)
!     -------------------------------------------------
!
  itop = 0
  ibot = 0
  IF (z1 .LT. z2) THEN
    y1 = z1
    IF (x1 .GT. z1) THEN
      y1   = x1
      itop = 1
    ENDIF
    y2 = z2
    IF (x2 .LT. z2) THEN
      y2   = x2
      ibot = 1
    ENDIF
  ELSE
    y1 = z2
    IF (x1 .GT. z2) THEN
      y1   = x1
      itop = 1
    ENDIF
    y2 = z1
    IF (x2 .LT. z1) THEN
      y2   = x2
      ibot = 1
    ENDIF
  ENDIF
!
!     ---------------------------------------
!     2. Set weights for forward interpolator
!     ---------------------------------------
!
  dy    = y2 - y1
  dz    = z1 - z2
! Thickness of output layers all non-zero - divisor dz is positive
  dzd   = 1.0_JPRB / dz
  dzddy = dzd * dy
  w10   = (z1 - y1) * dzddy
  w20   = (z1 - y2) * dzddy
  dx    = (x2 - x1)
! Thickness of input layers all non-zero - divisor dx is positive
  dxd   = 1.0_JPRB / dx
  d = (x2 - z2) * dxd
  IF (z1 .LT. z2 .AND. ibot .EQ. 0) THEN
    w1 = w10 + w20 * d
    w2 = w20 * (1.0_JPRB - d)
  ELSE IF (z1 .GT. z2 .AND. itop .EQ. 0) THEN
    w2 = w20 + w10 * (1.0_JPRB - d)
    w1 = w10 * d
  ELSE
    w1 = w10
    w2 = w20
  ENDIF
!
!
END SUBROUTINE RTTOV_SUBLAYER
