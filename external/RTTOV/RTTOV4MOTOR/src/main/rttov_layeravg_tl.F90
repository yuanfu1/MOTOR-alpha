SUBROUTINE RTTOV_LAYERAVG_TL( &
            & PX1,    &
            & PX1_TL, &
            & PX2,    &
            & PX2_TL, &
            & KN1,    &
            & KN2,    &
            & PZ,     &
            & PZ_TL,  &
            & kstart, &
            & kend,   &
            & interp_mode)
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
! Current Code Owner: SAF NWP
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".
!
!-----------------------------------------------------------------------
!
!     Description:
!
!     Tangent linear code (by Niels Bormann) for the piecewise
!     weighted averaging interpolator described in the next
!     paragraph.
!
!     The module RTTOV_LAYERAVG determines profile interpolation
!     weights by considering integrations over of series of
!     segments [PX1(KI-1),PX1(KI+1)] using  piecewise linear
!     weighting having weights of zero at  KI-1 and KI+1 and
!     max weight at KI. KI ranges from 1 to KN1.
!
!     Method:
!     - Piecewise weighted interpolation in ln(P).
!     - Journal reference:
!       Rochon, Y.J., L. Garand, D.S. Turner, and S. Polavarapu.
!       Jacobian mapping between coordinate systems in data assimilation,
!       Q. J. of the Royal Met. Soc., vol 133, 1547-1558, 2007.
!       (www.interscience.wiley.com) DOI:10.1002/qj.117
!
!     History:
!
!     The original version of this routine (layeravg.f90) was developed
!     and tested by Yves J. Rochon with contributions from L. Garand,
!     D.S. Turner, and S. Polavarapu of the Atmospheric Science and
!     Technology Directorate, Environment Canada (EC).
!
!     Copyright 2006, Environment Canada, All Rights Reserved.
!
!     The RTTOV-9 development team made significant improvements and
!     changes to the code of the forward interpolator and TLM adjoint
!     pair (i.e. the 'rttov_layeravg*' and 'rttov_sublayer*' routines).
!     Contributors consist of Niels Bormann, Alan Geer, Deborah Salmond,
!     Mats Hamrud, Roger Saunders, Peter Rayer, and Pascal Brunel.
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
!     4         11/2007     Assuming extrapolation at constant amount, better
!                           deal with output layers intersecting the
!                           highest or lowest input domain levels.
!                           In the process, significantly changed the
!                           preparation setup for calling rttov_sublayer.
!                           Added comments / tidied. See comment 2 below.
!                           (A. Geer, ECMWF)
!     5         01/2008     Cleaning. Removal of variables and code
!                           related to calculation towards TL and AD
!                           contributions. (Niels Bormann)
!                           See comment 1 below.
!     6         02/2008     Updating of header section and
!                           additions/modifications of comments.
!                           (Y.J. Rochon, EC)
!     7         04/2008     Optimisations: (some done earlier)
!                           - Re-write to reduce number of IF tests
!                           - Inlining of SUBLAYER
!                           - KSTART & KEND to remove multiplies by zero
!                           - Levels loop over KN1 moved into this routine
!                           (D. Salmond and M. Hamrud, ECMWF)
!
!     Arguments:
!
!     Input
!
!     PX1(KN1).................Levels of output domain (e.g. lnP; in increasing values)
!     PX2(KN2).................Levels of input domain (e.g. lnP; in increasing values)
!     PX1_TL(KN1)..............Perturbations on PX1
!     PX2_TL(KN2)..............Perturbations on PX2
!     KN1......................Dimension of PX1.
!     KN2......................Dimension of other arrays.
!
!     Output
!
!     PZ(KN2,KN1)..............Resultant accumulated weighting factors for
!                              interpolation from input to output domains.
!     PZ_TL(KN2,KN1)...........TL of PZ
!     kstart(KN1)..............Start index for relevant PZ row array segment.
!     kend(KN1)................End index for relevant PZ row array segment.
!     Errorstatus..............Returns error status to allow a soft failure.
!
!     External functions and subroutines: RTTOV_SUBLAYER
!
!     Assumptions:
!
!     1) PX1(i)<PX1(i+1) & PX2(i)<PX2(i+1)
!
!     Comments:
!
!     1) The calculations related to the TL and AD contributions
!     (and related arguments) were removed as this is now done
!     using separate routines.
!
!     2) The behaviour has been changed from the original Rochon code
!     in dealing with weighting integrals when extending outside the input
!     level range. The input profile is now extrapolated at constant value.
!
!     The impact is of most practical significance for instruments where
!     the weighting function peaks at or near the surface such as SSMI.
!
!     The new approach increases the weights of contribution from
!     the lowest and highest input domain levels for output
!     layers intersecting these input domain boundaries.
!     It consists of applying contant value extrapolation by
!     introducing 'fake' layers. For the surface of the input domain,
!     as example, this implies creating a virtual surface layer which
!     extends to the lower boundary of the output domain layer which
!     contains the input domain surface. This increases the
!     contributing weight of the surface which would be otherwise
!     underestimated in the original code due to the interpolator
!     actually doing piecewise weighted averaging.
!
!--------------------------------------------------------------------
!
  USE parkind1, ONLY : jpim, jprb
!INTF_OFF
  USE rttov_const, ONLY : &
    interp_loglinear
!INTF_ON
  IMPLICIT NONE
!
! --- Subroutine arguments ---------------------
!
  INTEGER(KIND=jpim), INTENT(IN)  :: KN1, KN2
  REAL   (KIND=jprb), INTENT(IN)  :: PX1   (KN1     ), PX1_TL(KN1     )
  REAL   (KIND=jprb), INTENT(IN)  :: PX2   (KN2     ), PX2_TL(KN2     )
  REAL   (KIND=jprb), INTENT(OUT) :: PZ    (KN2, KN1), PZ_TL (KN2, KN1)
  INTEGER(KIND=jpim), INTENT(OUT) :: kstart(KN1     ), kend  (KN1     )
  INTEGER(KIND=jpim), INTENT(IN)  :: interp_mode
!INTF_END
#include "rttov_sublayer_tl.interface"
!
! --- Local scalars and arrays ----------------
!
  INTEGER(KIND=jpim) :: J, KI   , istartn, istartnp1 , istart    , iend
  REAL   (KIND=jprb) :: Z1   , Z2   , Z3     , ZW_THIS   , ZW_NEXT   , ZSUM   , ZSUM_R
  REAL   (KIND=jprb) :: Z1_tl, Z2_tl, Z3_tl  , ZW_THIS_tl, ZW_NEXT_tl, ZSUM_tl, ZSUM_R_tl
  REAL   (KIND=jprb) :: zz   (kn2)
  REAL   (KIND=jprb) :: zz_tl(kn2)
  REAL   (KIND=jprb) ::      &
    & y1   , y2   , d, w10       , w20       , dz12   , dz32     , dx   , dy   , dzd12   , dzd32   , dxd   (kn2), dzddy
  REAL   (KIND=jprb) :: y1_tl, y2_tl, d_tl   , w10_tl    , w20_tl    , dz12_tl, dz32_tl  , dx_tl, dy_tl, dzd12_tl,      &
    & dzd32_tl, dxd_tl(kn2), dzddy_tl
  INTEGER(KIND=jpim) :: ibot , itop
  INTEGER(KIND=jpim) :: jthis

IF (interp_mode == interp_loglinear) THEN

!-----------------------------------------------------------------------------
! Log-linear interpolation
!-----------------------------------------------------------------------------

  jthis = 1

  DO KI = 1, KN1

    DO WHILE (px2(jthis) < px1(ki))
      IF (jthis == kn2) EXIT
      jthis = jthis + 1
    ENDDO

    IF (jthis == 1 .AND. px1(ki) <= px2(jthis)) THEN
      kstart(ki) = jthis
      kend(ki) = jthis
      pz(jthis, ki) = 1._jprb
      pz_tl(jthis, ki) = 0._jprb
    ELSE IF (jthis == kn2 .AND. px1(ki) >= px2(jthis)) THEN
      kstart(ki) = jthis
      kend(ki) = jthis
      pz(jthis, ki) = 1._jprb
      pz_tl(jthis, ki) = 0._jprb
    ELSE
      kstart(ki) = jthis - 1
      kend(ki) = jthis
      pz(jthis-1, ki) = (px2(jthis) - px1(ki)) / (px2(jthis) - px2(jthis-1))
      pz(jthis, ki) = 1._jprb - pz(jthis-1, ki)
      pz_tl(jthis-1, ki) = (px2_tl(jthis) - px1_tl(ki)) / (px2(jthis) - px2(jthis-1)) - &
          (px2(jthis) - px1(ki)) * (px2_tl(jthis) - px2_tl(jthis-1)) / (px2(jthis) - px2(jthis-1))**2_jpim
      pz_tl(jthis, ki) = - pz_tl(jthis-1, ki)
    ENDIF

  ENDDO

ELSE

!-----------------------------------------------------------------------------
! Rochon interpolation
!-----------------------------------------------------------------------------

! -----------------
! 1. Initialization
! -----------------
  istartn = 0
  DO j = 1, kn2 - 1
    dx = (px2(j + 1) - px2(j))
    dxd(j)    = 1.0_JPRB / dx
    dx_tl     = (px2_tl(j + 1) - px2_tl(j))
    dxd_tl(j) =  - dx_tl / (dx ** 2)
  ENDDO
! ----------------------------------------------
! 2. Loop through output ('reference') levels
!    for linear interpolator contribution to TLM
! ----------------------------------------------
  DO KI = 1, KN1
! ----------------------------------------------------------
! 2.1 Set integration/averaging range boundaries and related
!     perturbations for output domain
!     Ranges of integration are z1 to z3: (z1,z3)
! ----------------------------------------------------------
! --- Central output level (ki) ----------------------
    z2    = px1(ki)
    z2_tl = px1_tl(ki)
! --- Output level ki-1 (upper bound of weighting intergral) ---
    IF (ki .EQ. 1) THEN
      z1    = 2.0_JPRB * z2 - px1(ki + 1)
      z1_tl = 2.0_JPRB * z2_tl - px1_tl(ki + 1)
    ELSE
      z1    = px1(ki - 1)
      z1_tl = px1_tl(ki - 1)
    ENDIF
! --- Output level ki+1 (lower bound of weighting integral) ---
    IF (ki .EQ. kn1) THEN
      z3    = 2.0_JPRB * z2 - z1
      z3_tl = 2.0_JPRB * z2_tl - z1_tl
    ELSE
      z3    = px1(ki + 1)
      z3_tl = px1_tl(ki + 1)
    ENDIF
    dz12     = z1 - z2
    dz32     = z3 - z2
    dzd12    = 1.0_JPRB / dz12
    dzd32    = 1.0_JPRB / dz32
    dz12_tl  = z1_tl - z2_tl
    dz32_tl  = z3_tl - z2_tl
    dzd12_tl =  - dz12_tl / (dz12 ** 2)
    dzd32_tl =  - dz32_tl / (dz32 ** 2)
    istart   = kn2
    iend     = 1
! ----------------------------------------------------------------
! 2.2 Loop through input ('available') levels, starting initially
!     from a fake layer 0, used for extrapolation beyond the upper
!     bound of the input levels.
! ----------------------------------------------------------------
    IF (istartn == 0) THEN
      zz(1)    = 0.0_JPRB
      zz_tl(1) = 0.0_JPRB
      j = 0
! ----------------------------------------------------------------
! 2.2.1 Loop through upper and lower integral for fake extra level
!       above highest input level = extrapolation
!       of input variable at constant values.
! ----------------------------------------------------------------
      IF (px2(j + 1) .GT. z1) THEN
! --- Input layer is involved in upper integral (ki-1 to ki) ---
        CALL rttov_sublayer_tl( &
              & z1,            &
              & z1_tl,         &
              & z2,            &
              & z2_tl,         &
              & z1,            &
              & z1_tl,         &
              & px2(j + 1),    &
              & px2_tl(j + 1), &
              & zw_this,       &
              & zw_this_tl,    &
              & zw_next,       &
              & zw_next_tl)
! --- Add weight contributions from the
!     fake layer to the top input layer. ---------------------
        zz(1)        = zz(1) + zw_this
        zz(j + 1)    = zz(j + 1) + zw_next
        zz_tl(1)     = zz_tl(1) + zw_this_tl
        zz_tl(j + 1) = zz_tl(j + 1) + zw_next_tl
        istart       = 1
      ENDIF
      IF (px2(j + 1) .GE. z2) THEN
! --- Input layer is involved in lower integral (ki to ki+1) ---
        CALL rttov_sublayer_tl( &
              & z3,            &
              & z3_tl,         &
              & z2,            &
              & z2_tl,         &
              & z1,            &
              & z1_tl,         &
              & px2(j + 1),    &
              & px2_tl(j + 1), &
              & zw_this,       &
              & zw_this_tl,    &
              & zw_next,       &
              & zw_next_tl)
! --- Add weight contributions from the
!     fake layer to the top input layer. ---------------------
        zz(1)        = zz(1) + zw_this
        zz(j + 1)    = zz(j + 1) + zw_next
        zz_tl(1)     = zz_tl(1) + zw_this_tl
        zz_tl(j + 1) = zz_tl(j + 1) + zw_next_tl
        istart       = 1
      ENDIF
      istartnp1 = 1
    ELSE
      istartnp1        = istartn
      zz(istartnp1)    = 0.0_JPRB
      zz_tl(istartnp1) = 0.0_JPRB
    ENDIF
! --------------------------------------------------------
! 2.2.2 Loop through upper and lower integral accumulating
!       weight contributions.
! --------------------------------------------------------
    DO j = istartnp1, kn2 - 1
      zz(j + 1)    = 0.0_JPRB
      zz_tl(j + 1) = 0.0_JPRB
      IF (px2(j) .LT. z3) THEN
        IF ((px2(j) .LE. z2) .AND. (px2(j + 1) .GT. z1)) THEN
! --- Input layer is involved in upper integral (ki-1 to ki) ---
! --- Following code from rttov_sublayer_tl inlined for optimisation --
!
!     call rttov_sublayer_tl(z1, z1_tl, z2,  z2_tl, &
!           & px2(j), px2_tl(j), px2(j+1), px2_tl(j+1), &
!           & zw_this, zw_this_tl, zw_next, zw_next_tl)
          IF (px2(j) .GT. z1) THEN
            y1    = px2(j)
            y1_tl = px2_tl(j)
          ELSE
            y1    = z1
            y1_tl = z1_tl
          ENDIF
          IF (px2(j + 1) .LT. z2) THEN
            y2    = px2(j + 1)
            y2_tl = px2_tl(j + 1)
            ibot  = 1
          ELSE
            y2    = z2
            y2_tl = z2_tl
            ibot  = 0
          ENDIF
! --- Set weight contributions ----------------------
          dy       = y2 - y1
          dzddy    = dzd12 * dy
          w10      = (z1 - y1) * dzddy
          w20      = (z1 - y2) * dzddy
          dy_tl    = y2_tl - y1_tl
          dzddy_tl = dzd12 * dy_tl + dzd12_tl * dy
          w10_tl   = (z1 - y1) * dzddy_tl + (z1_tl - y1_tl) * dzddy
          w20_tl   = (z1 - y2) * dzddy_tl + (z1_tl - y2_tl) * dzddy
          IF (ibot .EQ. 0) THEN
            d = (px2(j + 1) - z2) * dxd(j)
            zw_this    = w10 + w20 * d
            zw_next    = w20 * (1.0_JPRB - d)
            d_tl       = (px2(j + 1) - z2) * dxd_tl(j) + (px2_tl(j + 1) - z2_tl) * dxd(j)
            zw_this_tl = w10_tl + w20 * d_tl + w20_tl * d
            zw_next_tl =  - w20 * d_tl + w20_tl * (1.0_JPRB - d)
          ELSE
            zw_this    = w10
            zw_next    = w20
            zw_this_tl = w10_tl
            zw_next_tl = w20_tl
          ENDIF
! --- End of inlined rttov_sublayer_tl -------------
! --- Add weight contributions ----------------------
          zz(j)        = zz(j) + zw_this
          zz(j + 1)    = zz(j + 1) + zw_next
          zz_tl(j)     = zz_tl(j) + zw_this_tl
          zz_tl(j + 1) = zz_tl(j + 1) + zw_next_tl
          istart       = min(j, istart)
          iend         = j + 1
        ENDIF
        IF ((px2(j + 1) .GE. z2)) THEN
! --- Input layer is involved in lower integral (ki to ki+1) ---
! --- Following code from rttov_sublayer_tl inlined for optimisation --
!
!     call rttov_sublayer_tl(z3, z3_tl, z2, z2_tl,&
!        & px2(j), px2_tl(j), px2(j+1), px2_tl(j+1), &
!        & zw_this, zw_this_tl, zw_next, zw_next_tl)
          IF (px2(j) .GT. z2) THEN
            y1    = px2(j)
            y1_tl = px2_tl(j)
            itop  = 1
          ELSE
            y1    = z2
            y1_tl = z2_tl
            itop  = 0
          ENDIF
          IF (px2(j + 1) .LT. z3) THEN
            y2    = px2(j + 1)
            y2_tl = px2_tl(j + 1)
          ELSE
            y2    = z3
            y2_tl = z3_tl
          ENDIF
! --- Set weight contributions ----------------------
          dy       = y2 - y1
          dzddy    = dzd32 * dy
          w10      = (z3 - y1) * dzddy
          w20      = (z3 - y2) * dzddy
          dy_tl    = y2_tl - y1_tl
          dzddy_tl = dzd32 * dy_tl + dzd32_tl * dy
          w10_tl   = (z3 - y1) * dzddy_tl + (z3_tl - y1_tl) * dzddy
          w20_tl   = (z3 - y2) * dzddy_tl + (z3_tl - y2_tl) * dzddy
          IF (itop .EQ. 0) THEN
            d = (px2(j + 1) - z2) * dxd(j)
            zw_next    = w20 + w10 * (1.0_JPRB - d)
            zw_this    = w10 * d
            d_tl       = (px2(j + 1) - z2) * dxd_tl(j) + (px2_tl(j + 1) - z2_tl) * dxd(j)
            zw_next_tl = w20_tl - w10 * d_tl + w10_tl * (1.0_JPRB - d)
            zw_this_tl = w10 * d_tl + w10_tl * d
          ELSE
            zw_this    = w10
            zw_next    = w20
            zw_this_tl = w10_tl
            zw_next_tl = w20_tl
          ENDIF
! --- End of inlined rttov_sublayer_tl -------------
! --- Add weight contributions ----------------------
          zz(j)        = zz(j) + zw_this
          zz(j + 1)    = zz(j + 1) + zw_next
          zz_tl(j)     = zz_tl(j) + zw_this_tl
          zz_tl(j + 1) = zz_tl(j + 1) + zw_next_tl
          IF (j < istart) istart = j
          iend = j + 1
        ENDIF
      ELSE
        EXIT
      ENDIF
    ENDDO
    j = kn2
    IF ((px2(j) .LT. z3)) THEN
! ----------------------------------------------------------
! 2.2.3 Loop through upper and lower integral for fake extra
!       level below lowest input level = extrapolation
!       of input variable at constant value.
! ----------------------------------------------------------
      IF ((px2(j) .LE. z2)) THEN
! --- Input layer is involved in upper integral (ki-1 to ki) ---
        CALL rttov_sublayer_tl( &
              & z1,         &
              & z1_tl,      &
              & z2,         &
              & z2_tl,      &
              & px2(j),     &
              & px2_tl(j),  &
              & z3,         &
              & z3_tl,      &
              & zw_this,    &
              & zw_this_tl, &
              & zw_next,    &
              & zw_next_tl)
! --- Add weight contributions from the
!     fake layer to the bottom input layer. ----------------
        zz(j)      = zz(j) + zw_this
        zz(kn2)    = zz(kn2) + zw_next
        zz_tl(j)   = zz_tl(j) + zw_this_tl
        zz_tl(kn2) = zz_tl(kn2) + zw_next_tl
        iend       = kn2
      ENDIF
      IF ((px2(j) .LT. z3)) THEN
! --- Input layer is involved in lower integral (ki to ki+1) ---
        CALL rttov_sublayer_tl( &
              & z3,         &
              & z3_tl,      &
              & z2,         &
              & z2_tl,      &
              & px2(j),     &
              & px2_tl(j),  &
              & z3,         &
              & z3_tl,      &
              & zw_this,    &
              & zw_this_tl, &
              & zw_next,    &
              & zw_next_tl)
! --- Add weight contributions from the
!     fake layer to the bottom input layer. ------------------
        zz(j)      = zz(j) + zw_this
        zz(kn2)    = zz(kn2) + zw_next
        zz_tl(j)   = zz_tl(j) + zw_this_tl
        zz_tl(kn2) = zz_tl(kn2) + zw_next_tl
        iend       = kn2
      ENDIF
    ENDIF
! -------------------------------------------------------------
! 2.3 Normalize sum of weights to unity (instead of calculating
!     and dividing by weighting denominator)
! -------------------------------------------------------------
    zsum = sum(zz(istart:iend))
    zsum_tl                = sum(zz_tl(istart:iend))
    zsum_r                 = 1.0_JPRB / zsum
    zsum_r_tl              =  - zsum_tl / (zsum * zsum)
    pz(istart:iend, ki)    = zz(istart:iend) * zsum_r
    pz_tl(istart:iend, ki) = zz_tl(istart:iend) * zsum_r + zz(istart:iend) * zsum_r_tl
! --- Save time in the input layer loop by ignoring irrelevant layers ---
    IF (istart > 1) istartn = istart
    kstart(ki) = istart
    kend(ki)   = iend
  ENDDO

ENDIF ! Interpolation mode

END SUBROUTINE RTTOV_LAYERAVG_TL
