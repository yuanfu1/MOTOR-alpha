! Description:
!> @file
!!   Compute optical parameters for ice cloud from Baran 2014 parameterisation
!
!> @brief
!!   Compute optical parameters for ice cloud from Baran 2014 parameterisation
!!
!! @param[in]     optp             Baran optical property coefficient structure
!! @param[in]     ichn             channel index
!! @param[in]     t_in             layer temperature (K)
!! @param[in]     iwc_in           layer ice water content (g/cm3)
!! @param[out]    abso             computed absorption coefficient
!! @param[out]    sca              computed scattering coefficient
!! @param[out]    bpr              computed backscattering parameter
!! @param[out]    asym             computed asymmetry parameter
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
SUBROUTINE rttov_baran2014_calc_optpar (optp, ichn, t_in, iwc_in, abso, sca, bpr, asym)

  USE rttov_types, ONLY : rttov_optp_baran
  USE parkind1, ONLY : jpim, jprb
!INTF_OFF
  USE mod_rttov_baran2014_icldata
!INTF_ON
  IMPLICIT NONE

  TYPE(rttov_optp_baran),INTENT(IN)  :: optp
  INTEGER(KIND=jpim)    ,INTENT(IN)  :: ichn
  REAL(KIND=jprb)       ,INTENT(IN)  :: t_in
  REAL(KIND=jprb)       ,INTENT(IN)  :: iwc_in
  REAL(KIND=jprb)       ,INTENT(OUT) :: abso
  REAL(KIND=jprb)       ,INTENT(OUT) :: sca
  REAL(KIND=jprb)       ,INTENT(OUT) :: bpr
  REAL(KIND=jprb)       ,INTENT(OUT) :: asym
!INTF_END

  REAL(KIND=jprb)  :: absi, absj, bpri, asymi
  REAL(KIND=jprb)  :: scai, scaj, bprj, asymj

  REAL(KIND=jprb)     :: dx_dwn
  REAL(KIND=jprb)     :: liwc
  INTEGER(KIND=jpim)  :: iwn, jwn
  REAL(KIND=jprb)     :: iwc,t

!- End of header --------------------------------------------------------

! Do not need to test null values for IWC because
! this is done in the calling surbroutine opdpscattir
!    abso = 0.0_jprb
!    sca  = 0.0_jprb
!    bpr  = 0.0_jprb
!    asym = 0.0_jprb
!
!    IF( IWC .LE. 0._jprb ) THEN
!      RETURN
!    ENDIF

    if (iwc_in < baran2014_iwc_min) then
      iwc = baran2014_iwc_min
    else if (iwc_in > baran2014_iwc_max) then
      iwc = baran2014_iwc_max
    else
      iwc = iwc_in
    endif

    if (T_in < baran2014_temp_min) then
      T = baran2014_temp_min
    else if (T_in > baran2014_temp_max) then
      T = baran2014_temp_max
    else
      T = T_in
    endif

    LIWC = LOG10(iwc)
    iwn = optp%iwn(ichn)
    jwn = optp%jwn(ichn)
    dx_dwn = optp%dx_dwn(ichn)

    absi = baran2014_regcoef_abs(1_jpim,iwn)    + &
         & baran2014_regcoef_abs(2_jpim,iwn)*T + &
         & baran2014_regcoef_abs(3_jpim,iwn)*LIWC + &
         & baran2014_regcoef_abs(4_jpim,iwn)*T*T + &
         & baran2014_regcoef_abs(5_jpim,iwn)*LIWC*LIWC + &
         & baran2014_regcoef_abs(6_jpim,iwn)*T*LIWC

    absj = baran2014_regcoef_abs(1_jpim,jwn)    + &
         & baran2014_regcoef_abs(2_jpim,jwn)*T + &
         & baran2014_regcoef_abs(3_jpim,jwn)*LIWC + &
         & baran2014_regcoef_abs(4_jpim,jwn)*T*T + &
         & baran2014_regcoef_abs(5_jpim,jwn)*LIWC*LIWC + &
         & baran2014_regcoef_abs(6_jpim,jwn)*T*LIWC

    absi = 10**absi
    absj = 10**absj

    abso = absi + (absj - absi) * dx_dwn

    scai = baran2014_regcoef_sca(1_jpim,iwn)    + &
         & baran2014_regcoef_sca(2_jpim,iwn)*T + &
         & baran2014_regcoef_sca(3_jpim,iwn)*LIWC + &
         & baran2014_regcoef_sca(4_jpim,iwn)*T*T + &
         & baran2014_regcoef_sca(5_jpim,iwn)*LIWC*LIWC + &
         & baran2014_regcoef_sca(6_jpim,iwn)*T*LIWC

    scaj = baran2014_regcoef_sca(1_jpim,jwn)    + &
         & baran2014_regcoef_sca(2_jpim,jwn)*T + &
         & baran2014_regcoef_sca(3_jpim,jwn)*LIWC + &
         & baran2014_regcoef_sca(4_jpim,jwn)*T*T + &
         & baran2014_regcoef_sca(5_jpim,jwn)*LIWC*LIWC + &
         & baran2014_regcoef_sca(6_jpim,jwn)*T*LIWC

    scai = 10**scai
    scaj = 10**scaj

    sca = scai + (scaj - scai) * dx_dwn

    bpri = baran2014_regcoef_bpr(1_jpim,iwn)    + &
         & baran2014_regcoef_bpr(2_jpim,iwn)*T + &
         & baran2014_regcoef_bpr(3_jpim,iwn)*LIWC
    bprj = baran2014_regcoef_bpr(1_jpim,jwn)    + &
         & baran2014_regcoef_bpr(2_jpim,jwn)*T + &
         & baran2014_regcoef_bpr(3_jpim,jwn)*LIWC

    bpr = bpri + (bprj - bpri) * dx_dwn

    asymi = baran2014_regcoef_asym(1_jpim,iwn)    + &
          & baran2014_regcoef_asym(2_jpim,iwn)*T + &
          & baran2014_regcoef_asym(3_jpim,iwn)*LIWC
    asymj = baran2014_regcoef_asym(1_jpim,jwn)    + &
          & baran2014_regcoef_asym(2_jpim,jwn)*T + &
          & baran2014_regcoef_asym(3_jpim,jwn)*LIWC

    asym = asymi + (asymj - asymi) * dx_dwn
    if( asym .GT. 1.0_jprb) asym = 1.0_jprb

END SUBROUTINE
