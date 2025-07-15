! Description:
!> @file
!!   AD of optical parameters for ice cloud from Baran 2018 parameterisation
!
!> @brief
!!   AD of optical parameters for ice cloud from Baran 2018 parameterisation
!!
!! @param[in]     optp             Baran optical property coefficient structure
!! @param[in]     ichn             channel index
!! @param[in]     t_in             layer temperature (K)
!! @param[in]     iwc_in           layer ice water content (g/cm3)
!! @param[in,out] t_in_ad          layer temperature increment
!! @param[in,out] iwc_in_ad        layer ice water content increment
!! @param[in,out] abs_ad           absorption coefficient increment
!! @param[in,out] sca_ad           scattering coefficient increment
!! @param[in,out] bpr_ad           backscattering parameter increment
!! @param[in,out] asym_ad          asymmetry parameter increment
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
SUBROUTINE rttov_baran2018_calc_optpar_ad (optp, ichn, &
                       & T_in, iwc_in, T_in_ad, iwc_in_ad, &
                       & abs_ad, sca_ad, bpr_ad, asym_ad)

  USE rttov_types, ONLY : rttov_optp_baran
  USE parkind1, ONLY : jpim, jprb
!INTF_OFF 
  USE mod_rttov_baran2018_icldata
!INTF_ON
  IMPLICIT NONE

  TYPE(rttov_optp_baran) ,INTENT(IN)  :: optp
  INTEGER(KIND=jpim)     ,INTENT(IN)  :: ichn
  REAL(KIND=jprb)        ,INTENT(IN)  :: t_in
  REAL(KIND=jprb)        ,INTENT(IN)  :: iwc_in
  REAL(KIND=jprb)        ,INTENT(INOUT):: t_in_ad
  REAL(KIND=jprb)        ,INTENT(INOUT):: iwc_in_ad

  REAL(KIND=jprb)        ,INTENT(INOUT) :: abs_ad
  REAL(KIND=jprb)        ,INTENT(INOUT) :: sca_ad
  REAL(KIND=jprb)        ,INTENT(INOUT) :: bpr_ad
  REAL(KIND=jprb)        ,INTENT(INOUT) :: asym_ad
!INTF_END

  REAL(KIND=jprb)  :: absi, absj!, brp, bpri, bprj
  REAL(KIND=jprb)  :: scai, scaj, asym, asymi, asymj
  REAL(KIND=jprb)  :: absi_ad, absj_ad, bpri_ad, asymi_ad
  REAL(KIND=jprb)  :: scai_ad, scaj_ad, bprj_ad, asymj_ad

  REAL(KIND=jprb)     :: dx_dwn
  REAL(KIND=jprb)     :: liwc, liwc_ad
  INTEGER(KIND=jpim)  :: iwn, jwn
  REAL(KIND=jprb)     :: iwc,t
  REAL(KIND=jprb)     :: iwc_ad,t_ad

  LOGICAL :: test5, test6, test7, test8

!- End of header --------------------------------------------------------

    test5 =.FALSE.
    test6 =.FALSE.
    test7 =.FALSE.
    test8 =.FALSE.

    if (iwc_in < baran2018_iwc_min) then
      iwc = baran2018_iwc_min
      test5=.TRUE.
    else if (iwc_in > baran2018_iwc_max) then
      iwc = baran2018_iwc_max
      test6=.TRUE.
    else
      iwc = iwc_in
    endif

    if (T_in < baran2018_temp_min) then
      T = baran2018_temp_min
      test7=.TRUE.
    else if (T_in > baran2018_temp_max) then
      T = baran2018_temp_max
      test8=.TRUE.
    else
      T = T_in
    endif

    LIWC = LOG10(iwc)

    iwn = optp%iwn(ichn)
    jwn = optp%jwn(ichn)
    dx_dwn = optp%dx_dwn(ichn)

    t_ad     = 0.0_jprb
    iwc_ad   = 0.0_jprb
    liwc_ad  = 0.0_jprb
    absi_ad  = 0.0_jprb
    absj_ad  = 0.0_jprb
    bpri_ad  = 0.0_jprb
    bprj_ad  = 0.0_jprb
    scai_ad  = 0.0_jprb
    scaj_ad  = 0.0_jprb
    asymi_ad = 0.0_jprb
    asymj_ad = 0.0_jprb


! ASYM
    asymi = baran2018_regcoef_asym(1_jpim,iwn)    + &
         &  baran2018_regcoef_asym(2_jpim,iwn)*T  + &
         &  baran2018_regcoef_asym(3_jpim,iwn)*LIWC
    asymj = baran2018_regcoef_asym(1_jpim,jwn)    + &
         &  baran2018_regcoef_asym(2_jpim,jwn)*T  + &
         &  baran2018_regcoef_asym(3_jpim,jwn)*LIWC

    asym = asymi + (asymj - asymi) * dx_dwn

    if( asym .GT. 1.0_jprb) then
      asym_ad = 0.0_jprb
    endif

    asymi_ad = asymi_ad + (1.0_jprb - dx_dwn)*asym_ad
    asymj_ad = asymj_ad + asym_ad * dx_dwn

    T_ad =     T_ad    + baran2018_regcoef_asym(2_jpim,iwn)*asymi_ad
    LIWC_ad =  LIWC_ad + baran2018_regcoef_asym(3_jpim,iwn)*asymi_ad
    asymi_ad = 0.0_jprb
    T_ad =     T_ad    + baran2018_regcoef_asym(2_jpim,jwn)*asymj_ad
    LIWC_ad =  LIWC_ad + baran2018_regcoef_asym(3_jpim,jwn)*asymj_ad
    asymj_ad = 0.0_jprb



! BPR
!     bpri = baran2018_regcoef_bpr(1_jpim,iwn)    + &
!          & baran2018_regcoef_bpr(2_jpim,iwn)*T  + &
!          & baran2018_regcoef_bpr(3_jpim,iwn)*LIWC
!     bprj = baran2018_regcoef_bpr(1_jpim,jwn)    + &
!          & baran2018_regcoef_bpr(2_jpim,jwn)*T  + &
!          & baran2018_regcoef_bpr(3_jpim,jwn)*LIWC

!     bpr = bpri + (bprj - bpri) * dx_dwn

    bpri_ad = bpri_ad + (1.0_jprb - dx_dwn)*bpr_ad
    bprj_ad = bprj_ad + bpr_ad * dx_dwn


    T_ad =     T_ad    + baran2018_regcoef_bpr(2_jpim,iwn)*bpri_ad
    LIWC_ad =  LIWC_ad + baran2018_regcoef_bpr(3_jpim,iwn)*bpri_ad
    bpri_ad = 0.0_jprb
    T_ad =     T_ad    + baran2018_regcoef_bpr(2_jpim,jwn)*bprj_ad
    LIWC_ad =  LIWC_ad + baran2018_regcoef_bpr(3_jpim,jwn)*bprj_ad
    bprj_ad = 0.0_jprb


! SCA
    scai = baran2018_regcoef_sca(1_jpim,iwn)    + &
         & baran2018_regcoef_sca(2_jpim,iwn)*T + &
         & baran2018_regcoef_sca(3_jpim,iwn)*LIWC + &
         & baran2018_regcoef_sca(4_jpim,iwn)*T*T + &
         & baran2018_regcoef_sca(5_jpim,iwn)*LIWC*LIWC + &
         & baran2018_regcoef_sca(6_jpim,iwn)*T*LIWC
    scaj = baran2018_regcoef_sca(1_jpim,jwn)    + &
         & baran2018_regcoef_sca(2_jpim,jwn)*T + &
         & baran2018_regcoef_sca(3_jpim,jwn)*LIWC + &
         & baran2018_regcoef_sca(4_jpim,jwn)*T*T + &
         & baran2018_regcoef_sca(5_jpim,jwn)*LIWC*LIWC + &
         & baran2018_regcoef_sca(6_jpim,jwn)*T*LIWC
    scai = 10**scai
    scaj = 10**scaj

    scai_ad = scai_ad + (1.0_jprb - dx_dwn)*sca_ad
    scaj_ad = scaj_ad + sca_ad * dx_dwn
    scai_ad = scai * scai_ad * LOG(10.0_jprb)
    scaj_ad = scaj * scaj_ad * LOG(10.0_jprb)

    T_ad =     T_ad    + (baran2018_regcoef_sca(2_jpim,iwn) + 2.0_jprb*baran2018_regcoef_sca(4_jpim,iwn)*T + &
                          baran2018_regcoef_sca(6_jpim,iwn)*LIWC)*scai_ad
    LIWC_ad =  LIWC_ad + (baran2018_regcoef_sca(3_jpim,iwn) + 2.0_jprb*baran2018_regcoef_sca(5_jpim,iwn)*LIWC + &
                          baran2018_regcoef_sca(6_jpim,iwn)*T)*scai_ad
    scai_ad = 0.0_jprb
    T_ad =     T_ad    + (baran2018_regcoef_sca(2_jpim,jwn) + 2.0_jprb*baran2018_regcoef_sca(4_jpim,jwn)*T + &
                          baran2018_regcoef_sca(6_jpim,jwn)*LIWC)*scaj_ad
    LIWC_ad =  LIWC_ad + (baran2018_regcoef_sca(3_jpim,jwn) + 2.0_jprb*baran2018_regcoef_sca(5_jpim,jwn)*LIWC + &
                          baran2018_regcoef_sca(6_jpim,jwn)*T)*scaj_ad
    scaj_ad = 0.0_jprb


! ABS
    absi = baran2018_regcoef_abs(1_jpim,iwn)    + &
         & baran2018_regcoef_abs(2_jpim,iwn)*T + &
         & baran2018_regcoef_abs(3_jpim,iwn)*LIWC + &
         & baran2018_regcoef_abs(4_jpim,iwn)*T*T + &
         & baran2018_regcoef_abs(5_jpim,iwn)*LIWC*LIWC + &
         & baran2018_regcoef_abs(6_jpim,iwn)*T*LIWC
    absj = baran2018_regcoef_abs(1_jpim,jwn)    + &
         & baran2018_regcoef_abs(2_jpim,jwn)*T + &
         & baran2018_regcoef_abs(3_jpim,jwn)*LIWC + &
         & baran2018_regcoef_abs(4_jpim,jwn)*T*T + &
         & baran2018_regcoef_abs(5_jpim,jwn)*LIWC*LIWC + &
         & baran2018_regcoef_abs(6_jpim,jwn)*T*LIWC
    absi = 10**absi
    absj = 10**absj


    absi_ad = absi_ad + (1.0_jprb - dx_dwn)*abs_ad
    absj_ad = absj_ad + abs_ad * dx_dwn
    absi_ad = absi * absi_ad * LOG(10.0_jprb)
    absj_ad = absj * absj_ad * LOG(10.0_jprb)

    T_ad =     T_ad    + (baran2018_regcoef_abs(2_jpim,iwn) + 2.0_jprb*baran2018_regcoef_abs(4_jpim,iwn)*T + &
                          baran2018_regcoef_abs(6_jpim,iwn)*LIWC)*absi_ad
    LIWC_ad =  LIWC_ad + (baran2018_regcoef_abs(3_jpim,iwn) + 2.0_jprb*baran2018_regcoef_abs(5_jpim,iwn)*LIWC + &
                          baran2018_regcoef_abs(6_jpim,iwn)*T)*absi_ad
    absi_ad = 0.0_jprb
    T_ad =     T_ad    + (baran2018_regcoef_abs(2_jpim,jwn) + 2.0_jprb*baran2018_regcoef_abs(4_jpim,jwn)*T + &
                          baran2018_regcoef_abs(6_jpim,jwn)*LIWC)*absj_ad
    LIWC_ad =  LIWC_ad + (baran2018_regcoef_abs(3_jpim,jwn) + 2.0_jprb*baran2018_regcoef_abs(5_jpim,jwn)*LIWC + &
                          baran2018_regcoef_abs(6_jpim,jwn)*T)*absj_ad
    absj_ad = 0.0_jprb


    iwc_ad = iwc_ad + liwc_ad / ( iwc * LOG(10.0_jprb) )
    liwc_ad = 0.0_jprb

    if (.NOT. (test7 .OR. test8)) then
      T_in_ad = T_in_ad + T_ad
    endif
    if (.NOT. (test5 .OR. test6)) then
      iwc_in_ad = iwc_in_ad + iwc_ad
    endif
    iwc_ad=0.0_jprb
    T_ad=0.0_jprb

END SUBROUTINE
