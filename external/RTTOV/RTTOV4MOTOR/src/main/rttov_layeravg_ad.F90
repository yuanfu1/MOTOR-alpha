SUBROUTINE RTTOV_LAYERAVG_AD( &
            & PX1,    &
            & PX1_AD, &
            & PX2,    &
            & PX2_AD, &
            & KN1,    &
            & KN2,    &
            & PZ,     &
            & PZ_AD,  &
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
!     Adjoint of TL code (by Niels Bormann) for the piecewise
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
!     5         01/2008     Adjoint version for AD of input/output
!                           levels (lgradp option), Niels Bormann
!     6         02/2008     Updating of header section and
!                           additions/modifications of comments.
!                           (Y.J. Rochon, EC)
!     7         04/2008     Optimisations: (some done earlier)
!                           - Re-write to reduce number of IF tests
!                           - Inlining of SUBLAYER
!                           - KSTART & KEND to remove multiplies by zero
!                           - Levels loop over KN1 moved into this routine
!                           (D. Salmond and M. Hamrud, ECMWF)
!     8         06/2009     Corrected bug in indexing (T. Wilhelmsson)
!     Arguments:
!
!     Input
!
!     PX1(KN1).................Levels of output domain (e.g. lnP; in increasing values)
!     PZ_AD(KN2)...............Perturbations on PZ
!     PX2(KN2).................Levels of input domain (e.g. lnP; in increasing values)
!     KN1......................Dimension of PX1.
!     KN2......................Dimension of other arrays.
!     KI.......................Identifies region of interest:
!                              PX1(KI-1) to PX1(KI+1)
!
!     Output
!     PX1_AD(KN1)..............AD of output levels
!     PX2_AD(KN2)..............AD of input levels
!     PZ(KN2,KN1)..............Resultant accumulated weighting factors for
!                              interpolation from input to output domains.
!     kstart(KN1)..............Start index for relevant PZ row array segment.
!     kend(KN1)................End index for relevant PZ row array segment.
!     Errorstatus..............Returns error status to allow a soft failure.
!
!     External functions and subroutines: RTTOV_SUBLAYER_AD
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
  INTEGER(KIND=jpim), INTENT(IN)    :: KN1, KN2
  REAL   (KIND=jprb), INTENT(IN)    :: PX1   (KN1     ), PX2   (KN2     )
  REAL   (KIND=jprb), INTENT(INOUT) :: PX1_AD(KN1     ), PX2_AD(KN2     )
  REAL   (KIND=jprb), INTENT(INOUT) :: PZ    (KN2, KN1), PZ_AD (KN2, KN1)
  INTEGER(KIND=jpim), INTENT(OUT)   :: kstart(KN1     ), kend  (KN1     )
  INTEGER(KIND=jpim), INTENT(IN)    :: interp_mode
!INTF_END
#include "rttov_sublayer.interface"
#include "rttov_sublayer_ad.interface"
!
! --- Local scalars and arrays ----------------
!
  INTEGER(KIND=jpim) :: J, KI, istart(kn1), istartn, istartnp1
  REAL   (KIND=jprb) :: Z1, Z2, Z3, Z1_AD  , Z2_AD    , Z3_AD
  REAL   (KIND=jprb) :: ZW_THIS    , ZW_NEXT   , ZSUM(KN1), ZSUM_AD
  REAL   (KIND=jprb) :: ZW_THIS_AD , ZW_NEXT_AD
  REAL   (KIND=jprb) :: px2_this_ad
  REAL   (KIND=jprb) :: px2_next_ad
  REAL   (KIND=jprb) :: zz(KN2, KN1)
  REAL   (KIND=jprb) ::      &
    & y1   , y2, d   , w10    , w20      , dz12   , dz32   , dx   , dy   , dzd12   , dzd32   , dxd(kn2), dzddy
  REAL   (KIND=jprb) ::      &
    & y1_ad, y2_ad     , d_ad, w10_ad , w20_ad   , dz12_ad, dz32_ad, dx_ad, dy_ad, dzd12_ad, dzd32_ad, dxd_ad, dzddy_ad
  INTEGER(KIND=jpim) :: ibot , itop
  INTEGER(KIND=jpim) :: jthis

IF (interp_mode == interp_loglinear) THEN

!-----------------------------------------------------------------------------
! Log-linear interpolation
!-----------------------------------------------------------------------------

  jthis = 1

  DO KI = 1, KN1  ! No need to reverse loop for AD

    DO WHILE (px2(jthis) < px1(ki))
      IF (jthis == kn2) EXIT
      jthis = jthis + 1
    ENDDO

    IF (jthis == 1 .AND. px1(ki) <= px2(jthis)) THEN
      kstart(ki) = jthis
      kend(ki) = jthis
      pz(jthis, ki) = 1._jprb
      pz_ad(jthis, ki) = 0._jprb
    ELSE IF (jthis == kn2 .AND. px1(ki) >= px2(jthis)) THEN
      kstart(ki) = jthis
      kend(ki) = jthis
      pz(jthis, ki) = 1._jprb
      pz_ad(jthis, ki) = 0._jprb
    ELSE
      kstart(ki) = jthis - 1
      kend(ki) = jthis
      pz(jthis-1, ki) = (px2(jthis) - px1(ki)) / (px2(jthis) - px2(jthis-1))
      pz(jthis, ki) = 1._jprb - pz(jthis-1, ki)

      pz_ad(jthis-1, ki) = pz_ad(jthis-1, ki) - pz_ad(jthis, ki)
      px2_ad(jthis) = px2_ad(jthis) + pz_ad(jthis-1, ki) * (1._jprb - pz(jthis-1, ki)) / (px2(jthis) - px2(jthis-1))
      px1_ad(ki) = px1_ad(ki) - pz_ad(jthis-1, ki) / (px2(jthis) - px2(jthis-1))
      px2_ad(jthis-1) = px2_ad(jthis-1) + pz_ad(jthis-1, ki) * pz(jthis-1, ki) / (px2(jthis) - px2(jthis-1))
    ENDIF

  ENDDO

ELSE

!-----------------------------------------------------------------------------
! Rochon interpolation
!-----------------------------------------------------------------------------

! -----------------
! 1. Initialization
! -----------------
  istart = 0
  px2_ad = 0.0_jprb
  kstart = 1
  kend   = 1
  DO j = 1, kn2 - 1
    dx     = (px2(j + 1) - px2(j))
    dxd(j) = 1.0_JPRB / dx
  ENDDO
! ----------------------------------------------
! 2. Loop through output ('reference') levels
!    for linear interpolator contribution to TLM
! -----------------------------------------------
  DO KI = 1, KN1
! -----------------------------------------------------------------
! 2.1 Set integration/averaging range boundaries for output domain.
!     Ranges of integration are z1 to z3: (z1,z3)
! -----------------------------------------------------------------
! --- Central output level (ki) ----------------------
    z2 = px1(ki)
! --- Output level ki-1 (upper bound of weighting intergral) ---
    IF (ki .EQ. 1) THEN
      z1 = 2.0 * z2 - px1(ki + 1)
    ELSE
      z1 = px1(ki - 1)
    ENDIF
! --- Output level ki+1 (lower bound of weighting integral) ---
    IF (ki .EQ. kn1) THEN
      z3 = 2.0 * z2 - z1
    ELSE
      z3 = px1(ki + 1)
    ENDIF
    dz12          = z1 - z2
    dzd12         = 1.0_JPRB / dz12
    dz32          = z3 - z2
    dzd32         = 1.0_JPRB / dz32
    zz(1:kn2, ki) = 0.0_jprb
    kstart(ki)    = kn2
    kend(ki)      = 1
! ----------------------------------------------------------------
! 2.2 Loop through input ('available') levels, starting initially
!     from a fake layer 0, used for extrapolation beyond the upper
!     bound of the input levels.
! ----------------------------------------------------------------
    istartn       = istart(ki)
    IF (istartn == 0) THEN
      j = 0
! ----------------------------------------------------------------
! 2.2.1 Loop through upper and lower integral for fake extra level
!       above highest input level = extrapolation
!       of input variable at constant values.
! ----------------------------------------------------------------
      IF ((px2(j + 1) .GT. z1)) THEN
! --- Input layer is involved in upper integral (ki-1 to ki) ---
        CALL rttov_sublayer( &
              & z1,         &
              & z2,         &
              & z1,         &
              & px2(j + 1), &
              & zw_this,    &
              & zw_next)
! --- Add weight contributions from the
!     fake layer to the top input layer. ---------------------
        zz(1, ki)     = zz(1, ki) + zw_this
        zz(j + 1, ki) = zz(j + 1, ki) + zw_next
        kstart(ki)    = 1
      ENDIF
      IF ((px2(j + 1) .GE. z2)) THEN
! --- Input layer is involved in lower integral (ki to ki+1) ---
        CALL rttov_sublayer( &
              & z3,         &
              & z2,         &
              & z1,         &
              & px2(j + 1), &
              & zw_this,    &
              & zw_next)
! --- Add weight contributions from the
!     fake layer to the top input layer. ---------------------
        zz(1, ki)     = zz(1, ki) + zw_this
        zz(j + 1, ki) = zz(j + 1, ki) + zw_next
        kstart(ki)    = 1
      ENDIF
    ENDIF
    IF (istartn == 0) THEN
      istartnp1 = 1
    ELSE
      istartnp1 = istartn
    ENDIF
! --------------------------------------------------------
! 2.2.2 Loop through upper and lower integral accumulating
!       weight contributions.
! --------------------------------------------------------
    DO j = istartnp1, kn2 - 1
      IF ((px2(j) .LT. z3)) THEN
        IF ((px2(j) .LE. z2) .AND. (px2(j + 1) .GT. z1)) THEN
! --- Input layer is involved in upper integral (ki-1 to ki) ---
! --- Following code from rttov_sublayer inlined for optimisation --
!
!     call rttov_sublayer(z3, z2, &
!                & px2(j), px2(j+1), &
!                & zw_this, zw_next)
          IF (px2(j) .GT. z1) THEN
            y1 = px2(j)
          ELSE
            y1 = z1
          ENDIF
          IF (px2(j + 1) .LT. z2) THEN
            y2   = px2(j + 1)
            ibot = 1
          ELSE
            ibot = 0
            y2   = z2
          ENDIF
! --- Set weight contributions ----------------------
          dy    = y2 - y1
          dzddy = dzd12 * dy
          w10   = (z1 - y1) * dzddy
          w20   = (z1 - y2) * dzddy
          IF (ibot .EQ. 0) THEN
            d = (px2(j + 1) - z2) * dxd(j)
            zw_this = w10 + w20 * d
            zw_next = w20 * (1.0_JPRB - d)
          ELSE
            zw_this = w10
            zw_next = w20
          ENDIF
! --- End of inlined rttov_sublayer  ----------------
! --- Add weight contributions ----------------------
          zz(j, ki)     = zz(j, ki) + zw_this
          zz(j + 1, ki) = zz(j + 1, ki) + zw_next
          kstart(ki)    = min(j, kstart(ki))
          kend(ki)      = j + 1
        ENDIF
        IF ((px2(j + 1) .GE. z2)) THEN
! --- Input layer is involved in lower integral (ki to ki+1) ---
! --- Following code from rttov_sublayer_tl inlined for optimisation --
!
!     call rttov_sublayer(z1, z2, &
!        & px2(j), px2(j+1), &
!        & zw_this, zw_next)
          IF (px2(j) .GT. z2) THEN
            y1   = px2(j)
            itop = 1
          ELSE
            y1   = z2
            itop = 0
          ENDIF
          IF (px2(j + 1) .LT. z3) THEN
            y2 = px2(j + 1)
          ELSE
            y2 = z3
          ENDIF
! --- Set weight contributions ----------------------
          dy    = y2 - y1
          dzddy = dzd32 * dy
          w10   = (z3 - y1) * dzddy
          w20   = (z3 - y2) * dzddy
          IF (itop .EQ. 0) THEN
            d = (px2(j + 1) - z2) * dxd(j)
            zw_next = w20 + w10 * (1.0_JPRB - d)
            zw_this = w10 * d
          ELSE
            zw_this = w10
            zw_next = w20
          ENDIF
! --- End of inlined rttov_sublayer -----------------
! --- Add weight contributions ----------------------
          zz(j, ki)     = zz(j, ki) + zw_this
          zz(j + 1, ki) = zz(j + 1, ki) + zw_next
          IF (j < kstart(ki)) kstart(ki) = j
          kend(ki) = j + 1
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
        CALL rttov_sublayer( &
              & z1,      &
              & z2,      &
              & px2(j),  &
              & z3,      &
              & zw_this, &
              & zw_next)
! --- Add weight contributions from the
!     fake layer to the bottom input layer. ----------------
        zz(j, ki)   = zz(j, ki) + zw_this
        zz(kn2, ki) = zz(kn2, ki) + zw_next
        kend(ki)    = kn2
      ENDIF
      IF ((px2(j) .LT. z3)) THEN
! --- Input layer is involved in lower integral (ki to ki+1) ---
        CALL rttov_sublayer( &
              & z3,      &
              & z2,      &
              & px2(j),  &
              & z3,      &
              & zw_this, &
              & zw_next)
! --- Add weight contributions from the
!     fake layer to the bottom input layer. ------------------
        zz(j, ki)   = zz(j, ki) + zw_this
        zz(kn2, ki) = zz(kn2, ki) + zw_next
        kend(ki)    = kn2
      ENDIF
    ENDIF
! -------------------------------------------------------------
! 2.3 Normalize sum of weights to unity (instead of calculating
!     and dividing by weighting denominator)
! -------------------------------------------------------------
    zsum(ki)                    = sum(zz(kstart(ki):kend(ki), ki))
    pz(kstart(ki):kend(ki), ki) = zz(kstart(ki):kend(ki), ki) / zsum(ki)
! --- Save time in the input layer loop by ignoring irrelevant layers ---
    IF (kstart(ki) > 1 .AND. ki < kn1) istart(ki + 1) = kstart(ki)
  ENDDO


! ------------------------------------------------------
! 3. AD Part (for linearized contribution to non-linear
!             interpolator component)
! ------------------------------------------------------
  DO KI = KN1, 1,  - 1
! --- Initialization --------------------------------
    z1_ad    = 0.0_JPRB
    z2_ad    = 0.0_JPRB
    z3_ad    = 0.0_JPRB
    dzddy_ad = 0.0_JPRB
! -----------------------------------------------------------------
! 3.1 Set integration/averaging range boundaries for output domain.
!     Ranges of integration are z1 to z3: (z1,z3)
! -----------------------------------------------------------------
! --- Central output level (ki) ----------------------
    z2       = px1(ki)
! --- Output level ki-1 (upper bound of weighting intergral) ---
    IF (ki .EQ. 1) THEN
      z1 = 2.0 * z2 - px1(ki + 1)
    ELSE
      z1 = px1(ki - 1)
    ENDIF
! --- Output level ki+1 (lower bound of weighting integral) ---
    IF (ki .EQ. kn1) THEN
      z3 = 2.0 * z2 - z1
    ELSE
      z3 = px1(ki + 1)
    ENDIF
    dz12 = z1 - z2
    dzd12 = 1.0_JPRB / dz12
    dz32 = z3 - z2
    dzd32 = 1.0_JPRB / dz32
! -----------------
! 3.2 Normalization
! -----------------
    zsum_ad = Sum( - pz_ad(kstart(ki):kend(ki), ki) * zz(kstart(ki):kend(ki), ki)) / (zsum(ki) ** 2)
    pz_ad(kstart(ki):kend(ki), ki) = pz_ad(kstart(ki):kend(ki), ki) / zsum(ki)
    pz_ad(kstart(ki):kend(ki), ki) = pz_ad(kstart(ki):kend(ki), ki) + zsum_ad
    zsum_ad = 0.0_JPRB
! ----------------------------------------------------------------
! 3.3 Loop through input ('available') levels, starting initially
!     from a fake layer KN2, used for extrapolation beyond the upper
!     bound of the input levels.
! ----------------------------------------------------------------
    j = kn2
    IF ((px2(j) .LT. z3)) THEN
      px2_next_ad = 0.0_JPRB
! ----------------------------------------------------------
! 3.3.1 Loop through upper and lower integral for fake extra
!       level below lowest input level = extrapolation
!       of input variable at constant value.
! ----------------------------------------------------------
      IF ((px2(j) .LT. z3)) THEN
! --- Input layer is involved in lower integral (ki to ki+1) ---
        zw_next_ad = pz_ad(kn2, ki)
        zw_this_ad = pz_ad(j, ki)
        CALL rttov_sublayer_ad( &
              & z3,          &
              & z3_ad,       &
              & z2,          &
              & z2_ad,       &
              & px2(j),      &
              & px2_ad(j),   &
              & z3,          &
              & px2_next_ad, &
              & zw_this,     &
              & zw_this_ad,  &
              & zw_next,     &
              & zw_next_ad)
      ENDIF
      IF ((px2(j) .LE. z2)) THEN
! --- Input layer is involved in upper integral (ki-1 to ki) ---
        zw_next_ad = pz_ad(kn2, ki)
        zw_this_ad = pz_ad(j, ki)
        CALL rttov_sublayer_ad( &
              & z1,          &
              & z1_ad,       &
              & z2,          &
              & z2_ad,       &
              & px2(j),      &
              & px2_ad(j),   &
              & z3,          &
              & px2_next_ad, &
              & zw_this,     &
              & zw_this_ad,  &
              & zw_next,     &
              & zw_next_ad)
      ENDIF
! --- Add weight contributions from the
!     fake layer to the bottom input layer. ---------------------
      z3_ad       = z3_ad + px2_next_ad
      px2_next_ad = 0.0_JPRB
    ENDIF
! --------------------------------------------------------
! 3.3.2 Loop through upper and lower integral accumulating
!       weight contributions.
! --------------------------------------------------------
    istartn = istart(ki)
    IF (istartn == 0) THEN
      istartnp1 = 1
    ELSE
      istartnp1 = istartn
    ENDIF
    DO j = kn2 - 1, istartnp1,  - 1
      IF ((px2(j) .LT. z3)) THEN
        IF ((px2(j + 1) .GE. z2)) THEN
! --- Input layer is involved in lower integral (ki to ki+1) ---
          zw_next_ad = pz_ad(j + 1, ki)
          zw_this_ad = pz_ad(j, ki)
! --- Following code from rttov_sublayer_ad inlined for optimisation --
!
!     call rttov_sublayer_ad(z3, z3_ad, &
!              & z2, z2_ad,&
!              & px2(j), px2_ad(j), px2(j+1), px2_ad(j+1), &
!              & zw_this, zw_this_ad, zw_next, zw_next_ad)
          IF (px2(j) .GT. z2) THEN
            y1   = px2(j)
            itop = 1
          ELSE
            y1   = z2
            itop = 0
          ENDIF
          IF (px2(j + 1) .LT. z3) THEN
            y2 = px2(j + 1)
          ELSE
            y2 = z3
          ENDIF
! --- Set weight contributions ----------------------
          dy    = y2 - y1
          dzddy = dzd32 * dy
          w10   = (z3 - y1) * dzddy
          w20   = (z3 - y2) * dzddy
          dx    = (px2(j + 1) - px2(j))
          IF (itop .EQ. 0) THEN
            d = (px2(j + 1) - z2) * dxd(j)
            zw_next = w20 + w10 * (1.0_JPRB - d)
            zw_this = w10 * d
          ELSE
            zw_this = w10
            zw_next = w20
          ENDIF
! --- AD part ---
          IF (itop .EQ. 0) THEN
            w10_ad        = d * zw_this_ad
            d_ad          = w10 * zw_this_ad
            w20_ad        = zw_next_ad
            w10_ad        = w10_ad + (1.0_JPRB - d) * zw_next_ad
            d_ad          = d_ad - w10 * zw_next_ad
            px2_ad(j + 1) = px2_ad(j + 1) + dxd(j) * d_ad
            z2_ad         = z2_ad - dxd(j) * d_ad
            dxd_ad        = (px2(j + 1) - z2) * d_ad
          ELSE
            w10_ad = zw_this_ad
            w20_ad = zw_next_ad
            dxd_ad = 0.0_JPRB
          ENDIF
          zw_this_ad = 0.0_JPRB
          zw_next_ad = 0.0_JPRB
          z3_ad      = z3_ad + dzddy * w20_ad
          y2_ad      =  - dzddy * w20_ad
          dzddy_ad   = (z3 - y2) * w20_ad
          z3_ad      = z3_ad + dzddy * w10_ad
          y1_ad      =  - dzddy * w10_ad
          dzddy_ad   = dzddy_ad + (z3 - y1) * w10_ad
          dy_ad      = dzd32 * dzddy_ad
          dzd32_ad   = dy * dzddy_ad
          y2_ad      = y2_ad + dy_ad
          y1_ad      = y1_ad - dy_ad
          w10_ad     = 0.0_JPRB
          w20_ad     = 0.0_JPRB
          IF (px2(j + 1) .LT. z3) THEN
            px2_ad(j + 1) = px2_ad(j + 1) + y2_ad
          ELSE
            z3_ad = z3_ad + y2_ad
          ENDIF
          IF (px2(j) .GT. z2) THEN
            px2_ad(j) = px2_ad(j) + y1_ad
          ELSE
            z2_ad = z2_ad + y1_ad
          ENDIF
          dz32_ad       =  - dzd32_ad / (dz32 ** 2)
          z3_ad         = z3_ad + dz32_ad
          z2_ad         = z2_ad - dz32_ad
          dx_ad         =  - dxd_ad / (dx ** 2)
! --- Add weight contributions for input levels ----------
          px2_ad(j + 1) = px2_ad(j + 1) + dx_ad
          px2_ad(j)     = px2_ad(j) - dx_ad
! --- End of inlined rttov_sublayer_ad -----------------
        ENDIF
        IF ((px2(j) .LE. z2) .AND. (px2(j + 1) .GT. z1)) THEN
! --- Input layer is involved in upper integral (ki-1 to ki) ---
          zw_next_ad = pz_ad(j + 1, ki)
          zw_this_ad = pz_ad(j, ki)
! --- Following code from rttov_sublayer_ad inlined for optimisation --
!
!     call rttov_sublayer_ad(z1, z1_ad, &
!              & z2, z2_ad,&
!              & px2(j), px2_ad(j), px2(j+1), px2_ad(j+1), &
!              & zw_this, zw_this_ad, zw_next, zw_next_ad)
          IF (px2(j) .GT. z1) THEN
            y1 = px2(j)
          ELSE
            y1 = z1
          ENDIF
          IF (px2(j + 1) .LT. z2) THEN
            y2   = px2(j + 1)
            ibot = 1
          ELSE
            y2   = z2
            ibot = 0
          ENDIF
! --- Set weights contributions -------------------------
          dy    = y2 - y1
          dzddy = dzd12 * dy
          w10   = (z1 - y1) * dzddy
          w20   = (z1 - y2) * dzddy
          dx    = (px2(j + 1) - px2(j))
          IF (ibot .EQ. 0) THEN
            d = (px2(j + 1) - z2) * dxd(j)
            zw_this = w10 + w20 * d
            zw_next = w20 * (1.0_JPRB - d)
          ELSE
            zw_this = w10
            zw_next = w20
          ENDIF
! --- AD part ---
          IF (ibot .EQ. 0) THEN
            w20_ad        = (1.0_JPRB - d) * zw_next_ad
            d_ad          =  - w20 * zw_next_ad
            w10_ad        = zw_this_ad
            w20_ad        = w20_ad + d * zw_this_ad
            d_ad          = d_ad + w20 * zw_this_ad
            px2_ad(j + 1) = px2_ad(j + 1) + dxd(j) * d_ad
            z2_ad         = z2_ad - dxd(j) * d_ad
            dxd_ad        = (px2(j + 1) - z2) * d_ad
          ELSE
            w10_ad = zw_this_ad
            w20_ad = zw_next_ad
            dxd_ad = 0.0_JPRB
          ENDIF
          zw_this_ad = 0.0_JPRB
          zw_next_ad = 0.0_JPRB
          z1_ad      = z1_ad + dzddy * w20_ad
          y2_ad      =  - dzddy * w20_ad
          dzddy_ad   = (z1 - y2) * w20_ad
          z1_ad      = z1_ad + dzddy * w10_ad
          y1_ad      =  - dzddy * w10_ad
          dzddy_ad   = dzddy_ad + (z1 - y1) * w10_ad
          dy_ad      = dzd12 * dzddy_ad
          dzd12_ad   = dy * dzddy_ad
          y2_ad      = y2_ad + dy_ad
          y1_ad      = y1_ad - dy_ad
          w10_ad     = 0.0_JPRB
          w20_ad     = 0.0_JPRB
          IF (px2(j + 1) .LT. z2) THEN
            px2_ad(j + 1) = px2_ad(j + 1) + y2_ad
          ELSE
            z2_ad = z2_ad + y2_ad
          ENDIF
          IF (px2(j) .GT. z1) THEN
            px2_ad(j) = px2_ad(j) + y1_ad
          ELSE
            z1_ad = z1_ad + y1_ad
          ENDIF
          dz12_ad       =  - dzd12_ad / (dz12 ** 2)
          z1_ad         = z1_ad + dz12_ad
          z2_ad         = z2_ad - dz12_ad
          dx_ad         =  - dxd_ad / (dx ** 2)
! --- Add weight contributions for input levels -------
          px2_ad(j + 1) = px2_ad(j + 1) + dx_ad
          px2_ad(j)     = px2_ad(j) - dx_ad
! --- End of inlined rttov_sublayer_ad -------------
        ENDIF
      ENDIF
    ENDDO
    IF (istartn == 0) THEN
      j = 0
! ----------------------------------------------------------
! 3.3.3 Loop through upper and lower integral for fake extra
!       level above highest input level = extrapolation
!       of input variable at constant value.
! ----------------------------------------------------------
      px2_this_ad = 0.0_JPRB
      IF ((px2(j + 1) .GE. z2)) THEN
! --- Input layer is involved in lower integral (ki to ki+1) ---
        zw_next_ad = pz_ad(j + 1, ki)
        zw_this_ad = pz_ad(1, ki)
        CALL rttov_sublayer_ad( &
              & z3,            &
              & z3_ad,         &
              & z2,            &
              & z2_ad,         &
              & z1,            &
              & px2_this_ad,   &
              & px2(j + 1),    &
              & px2_ad(j + 1), &
              & zw_this,       &
              & zw_this_ad,    &
              & zw_next,       &
              & zw_next_ad)
      ENDIF
      IF ((px2(j + 1) .GT. z1)) THEN
! --- Input layer is involved in upper integral (ki-1 to ki) ---
        zw_next_ad = pz_ad(j + 1, ki)
        zw_this_ad = pz_ad(1, ki)
        CALL rttov_sublayer_ad( &
              & z1,            &
              & z1_ad,         &
              & z2,            &
              & z2_ad,         &
              & z1,            &
              & px2_this_ad,   &
              & px2(j + 1),    &
              & px2_ad(j + 1), &
              & zw_this,       &
              & zw_this_ad,    &
              & zw_next,       &
              & zw_next_ad)
      ENDIF
! --- Add weight contributions from the
!     fake layer to the top input layer. ----------------
      z1_ad       = z1_ad + px2_this_ad
      px2_this_ad = 0.0_JPRB
    ENDIF
! ----------------------------------------------------------
! 3.4 Add accumulated weight contributions for output levels
! ----------------------------------------------------------
! --- Output level ki+1 (lower bound of weighting integral) ---
    IF (ki .EQ. kn1) THEN
      z2_ad = z2_ad + 2.0 * z3_ad
      z1_ad = z1_ad - z3_ad
    ELSE
      px1_ad(ki + 1) = px1_ad(ki + 1) + z3_ad
    ENDIF
    z3_ad = 0.0_JPRB
! --- Output level ki-1 (upper bound of weighting intergral) ---
    IF (ki .EQ. 1) THEN
      z2_ad          = z2_ad + 2.0 * z1_ad
      px1_ad(ki + 1) = px1_ad(ki + 1) - z1_ad
    ELSE
      px1_ad(ki - 1) = px1_ad(ki - 1) + z1_ad
    ENDIF
    z1_ad      = 0.0_JPRB
! --- Output level ki ----
    px1_ad(ki) = px1_ad(ki) + z2_ad
    z2_ad      = 0.0_JPRB
  ENDDO

ENDIF ! Interpolation mode

END SUBROUTINE RTTOV_LAYERAVG_AD
