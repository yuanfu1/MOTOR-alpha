! Description:
!> @file
!!   TL of optical parameters for ice cloud from Baran 2018 parameterisation
!
!> @brief
!!   TL of optical parameters for ice cloud from Baran 2018 parameterisation
!!
!! @details
!!   This also computes the direct model values.
!!
!! @param[in]     optp             Baran optical property coefficient structure
!! @param[in]     ichn             channel index
!! @param[in]     t_in             layer temperature (K)
!! @param[in]     iwc_in           layer ice water content (g/cm3)
!! @param[in]     t_in_tl          layer temperature perturbation
!! @param[in]     iwc_in_tl        layer ice water content perturbation
!! @param[out]    abso             computed absorption coefficient
!! @param[out]    sca              computed scattering coefficient
!! @param[out]    bpr              computed backscattering parameter
!! @param[out]    asym             computed asymmetry parameter
!! @param[out]    abso_tl          absorption coefficient perturbation
!! @param[out]    sca_tl           scattering coefficient perturbation
!! @param[out]    bpr_tl           backscattering parameter perturbation
!! @param[out]    asym_tl          asymmetry parameter perturbation
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
!    Copyright 2019, EUMETSAT, All Rights Reserved.
!
SUBROUTINE rttov_baran2018_calc_optpar_tl (optp, ichn, &
                       & T_in, IWC_in, T_in_tl, IWC_in_tl, &
                       & abso, sca, bpr, asym ,&
                       & abso_tl, sca_tl, bpr_tl, asym_tl)

  USE rttov_types, ONLY : rttov_optp_baran
  USE parkind1, ONLY : jpim, jprb
!INTF_OFF
  USE mod_rttov_baran2018_icldata
!INTF_ON
  IMPLICIT NONE

  TYPE(rttov_optp_baran) ,INTENT(IN)  :: optp
  INTEGER(KIND=jpim)     ,INTENT(IN)  :: ichn
  REAL(KIND=jprb)        ,INTENT(IN)  :: t_in, t_in_tl
  REAL(KIND=jprb)        ,INTENT(IN)  :: iwc_in, iwc_in_tl
  REAL(KIND=jprb)        ,INTENT(OUT) :: abso, abso_tl
  REAL(KIND=jprb)        ,INTENT(OUT) :: sca, sca_tl
  REAL(KIND=jprb)        ,INTENT(OUT) :: bpr, bpr_tl
  REAL(KIND=jprb)        ,INTENT(OUT) :: asym, asym_tl
!INTF_END

  REAL(KIND=jprb)  :: absi, absj, bpri, asymi
  REAL(KIND=jprb)  :: scai, scaj, bprj, asymj
  REAL(KIND=jprb)  :: absi_tl, absj_tl, bpri_tl, asymi_tl
  REAL(KIND=jprb)  :: scai_tl, scaj_tl, bprj_tl, asymj_tl

  REAL(KIND=jprb)     :: dx_dwn
  REAL(KIND=jprb)     :: liwc, liwc_tl
  INTEGER(KIND=jpim)  :: iwn, jwn
  REAL(KIND=jprb)     :: iwc,t,iwc_tl,t_tl

!- End of header --------------------------------------------------------

    if (iwc_in < baran2018_iwc_min) then
      iwc = baran2018_iwc_min
      iwc_tl = 0._jprb
    else if (iwc_in > baran2018_iwc_max) then
      iwc = baran2018_iwc_max
      iwc_tl = 0._jprb
    else
      iwc = iwc_in
      iwc_tl = iwc_in_tl
    endif

    if (T_in < baran2018_temp_min) then
      T = baran2018_temp_min
      T_tl = 0._jprb
    else if (T_in > baran2018_temp_max) then
      T = baran2018_temp_max
      T_tl = 0._jprb
    else
      T = T_in
      T_tl = T_in_tl
    endif

    LIWC = LOG10(iwc)
    LIWC_TL = iwc_tl / ( iwc * LOG(10.0_jprb) )

    iwn = optp%iwn(ichn)
    jwn = optp%jwn(ichn)
    dx_dwn = optp%dx_dwn(ichn)

    absi = baran2018_regcoef_abs(1_jpim,iwn)    + &
         & baran2018_regcoef_abs(2_jpim,iwn)*T + &
         & baran2018_regcoef_abs(3_jpim,iwn)*LIWC + &
         & baran2018_regcoef_abs(4_jpim,iwn)*T*T + &
         & baran2018_regcoef_abs(5_jpim,iwn)*LIWC*LIWC + &
         & baran2018_regcoef_abs(6_jpim,iwn)*T*LIWC
    absi_tl =  &
         & baran2018_regcoef_abs(2_jpim,iwn)*T_tl + &
         & baran2018_regcoef_abs(3_jpim,iwn)*LIWC_tl + &
         & baran2018_regcoef_abs(4_jpim,iwn)*2.0_jprb*T*T_tl + &
         & baran2018_regcoef_abs(5_jpim,iwn)*2.0_jprb*LIWC*LIWC_tl + &
         & baran2018_regcoef_abs(6_jpim,iwn)*(T*LIWC_tl + T_tl*LIWC )

    absj = baran2018_regcoef_abs(1_jpim,jwn)    + &
         & baran2018_regcoef_abs(2_jpim,jwn)*T + &
         & baran2018_regcoef_abs(3_jpim,jwn)*LIWC + &
         & baran2018_regcoef_abs(4_jpim,jwn)*T*T + &
         & baran2018_regcoef_abs(5_jpim,jwn)*LIWC*LIWC + &
         & baran2018_regcoef_abs(6_jpim,jwn)*T*LIWC
    absj_tl =  &
         & baran2018_regcoef_abs(2_jpim,jwn)*T_tl + &
         & baran2018_regcoef_abs(3_jpim,jwn)*LIWC_tl + &
         & baran2018_regcoef_abs(4_jpim,jwn)*2.0_jprb*T*T_tl + &
         & baran2018_regcoef_abs(5_jpim,jwn)*2.0_jprb*LIWC*LIWC_tl + &
         & baran2018_regcoef_abs(6_jpim,jwn)*(T*LIWC_tl + T_tl*LIWC )

    absi = 10**absi
    absi_tl = absi * absi_tl * LOG(10.0_jprb)
    absj = 10**absj
    absj_tl = absj * absj_tl * LOG(10.0_jprb)

    abso = absi + (absj - absi) * dx_dwn
    abso_tl = absi_tl + (absj_tl - absi_tl) * dx_dwn


    scai = baran2018_regcoef_sca(1_jpim,iwn)    + &
         & baran2018_regcoef_sca(2_jpim,iwn)*T + &
         & baran2018_regcoef_sca(3_jpim,iwn)*LIWC + &
         & baran2018_regcoef_sca(4_jpim,iwn)*T*T + &
         & baran2018_regcoef_sca(5_jpim,iwn)*LIWC*LIWC + &
         & baran2018_regcoef_sca(6_jpim,iwn)*T*LIWC
    scai_tl =  &
         & baran2018_regcoef_sca(2_jpim,iwn)*T_tl + &
         & baran2018_regcoef_sca(3_jpim,iwn)*LIWC_tl + &
         & baran2018_regcoef_sca(4_jpim,iwn)*2.0_jprb*T*T_tl + &
         & baran2018_regcoef_sca(5_jpim,iwn)*2.0_jprb*LIWC*LIWC_tl + &
         & baran2018_regcoef_sca(6_jpim,iwn)*(T*LIWC_tl + T_tl*LIWC )

    scaj = baran2018_regcoef_sca(1_jpim,jwn)    + &
         & baran2018_regcoef_sca(2_jpim,jwn)*T + &
         & baran2018_regcoef_sca(3_jpim,jwn)*LIWC + &
         & baran2018_regcoef_sca(4_jpim,jwn)*T*T + &
         & baran2018_regcoef_sca(5_jpim,jwn)*LIWC*LIWC + &
         & baran2018_regcoef_sca(6_jpim,jwn)*T*LIWC
    scaj_tl =  &
         & baran2018_regcoef_sca(2_jpim,jwn)*T_tl + &
         & baran2018_regcoef_sca(3_jpim,jwn)*LIWC_tl + &
         & baran2018_regcoef_sca(4_jpim,jwn)*2.0_jprb*T*T_tl + &
         & baran2018_regcoef_sca(5_jpim,jwn)*2.0_jprb*LIWC*LIWC_tl + &
         & baran2018_regcoef_sca(6_jpim,jwn)*(T*LIWC_tl + T_tl*LIWC )

    scai = 10**scai
    scai_tl = scai * scai_tl * LOG(10.0_jprb)
    scaj = 10**scaj
    scaj_tl = scaj * scaj_tl * LOG(10.0_jprb)

    sca = scai + (scaj - scai) * dx_dwn
    sca_tl = scai_tl + (scaj_tl - scai_tl) * dx_dwn


    bpri = baran2018_regcoef_bpr(1_jpim,iwn)    + &
         & baran2018_regcoef_bpr(2_jpim,iwn)*T + &
         & baran2018_regcoef_bpr(3_jpim,iwn)*LIWC
    bpri_tl =  &
         & baran2018_regcoef_bpr(2_jpim,iwn)*T_tl + &
         & baran2018_regcoef_bpr(3_jpim,iwn)*LIWC_tl
    bprj = baran2018_regcoef_bpr(1_jpim,jwn)    + &
         & baran2018_regcoef_bpr(2_jpim,jwn)*T + &
         & baran2018_regcoef_bpr(3_jpim,jwn)*LIWC
    bprj_tl =  &
         & baran2018_regcoef_bpr(2_jpim,jwn)*T_tl + &
         & baran2018_regcoef_bpr(3_jpim,jwn)*LIWC_tl

    bpr = bpri + (bprj - bpri) * dx_dwn
    bpr_tl = bpri_tl + (bprj_tl - bpri_tl) * dx_dwn


    asymi = baran2018_regcoef_asym(1_jpim,iwn)    + &
         & baran2018_regcoef_asym(2_jpim,iwn)*T + &
         & baran2018_regcoef_asym(3_jpim,iwn)*LIWC
    asymi_tl =  &
         & baran2018_regcoef_asym(2_jpim,iwn)*T_tl + &
         & baran2018_regcoef_asym(3_jpim,iwn)*LIWC_tl
    asymj = baran2018_regcoef_asym(1_jpim,jwn)    + &
         & baran2018_regcoef_asym(2_jpim,jwn)*T + &
         & baran2018_regcoef_asym(3_jpim,jwn)*LIWC
    asymj_tl =  &
         & baran2018_regcoef_asym(2_jpim,jwn)*T_tl + &
         & baran2018_regcoef_asym(3_jpim,jwn)*LIWC_tl

    asym = asymi + (asymj - asymi) * dx_dwn
    asym_tl = asymi_tl + (asymj_tl - asymi_tl) * dx_dwn

    if( asym .GT. 1.0_jprb) then
      asym = 1.0_jprb
      asym_tl = 0.0_jprb
    endif

END SUBROUTINE
