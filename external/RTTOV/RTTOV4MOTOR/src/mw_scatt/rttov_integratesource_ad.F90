!      
Subroutine rttov_integratesource_ad (&        
     & ccthres,       &! in
     & nlevels,       &! in
     & nchannels,     &! in
     & nprofiles,     &! in
     & nprofilesad,   &! in
     & lprofiles,     &! in
     & angles,        &! in
     & scatt_aux,     &! in
     & scatt_aux_ad,  &! inout
     & dp,            &! in
     & dp_ad,         &! inout
     & dm,            &! in
     & dm_ad,         &! inout
     & j_do,          &! inout
     & j_do_ad,       &! inout
     & j_up,          &! inout
     & j_up_ad)        ! inout 

  ! Description:
  ! integrate source in Eddington
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
  !    Copyright 2002, EUMETSAT, All Rights Reserved.
  !
  ! Method:
  ! - Bauer, P., 2002: Microwave radiative transfer modeling in clouds and precipitation.
  !     Part I: Model description.
  !     NWP SAF Report No. NWPSAF-EC-TR-005, 27 pp.
  ! - Chevallier, F., and P. Bauer, 2003:
  !     Model rain and clouds over oceans: comparison with SSM/I observations.
  !     Mon. Wea. Rev., 131, 1240-1255.
  ! - Moreau, E., P. Bauer and F. Chevallier, 2002: Microwave radiative transfer modeling in clouds and precipitation.
  !     Part II: Model evaluation.
  !     NWP SAF Report No. NWPSAF-EC-TR-006, 27 pp.
  !
  ! Current Code Owner: SAF NWP
  !
  ! History:
  ! Version   Date     Comment
  ! -------   ----     -------
  !  1.0       09/2002   Initial version     (E. Moreau)
  !  1.1       05/2003   RTTOV7.3 compatible (F. Chevallier)
  !  1.2       03/2004   Added polarimetry   (R. Saunders)
  !  1.3       08/2004   Polarimetry fixes   (U. O'Keefe)
  !  1.4       11/2004   Clean-up            (P. Bauer)
  !  1.5       11/2007   RTTOV9 version      (A. Geer)
  !  1.6       07/2008   Clear sky speed-ups (A. Geer)
  !  1.7       03/2010   Use mclayer rather than min_ssa (A. Geer)
  !  1.8       15/2014   Cure numerical instability (A. Geer)
  !  1.10      07/2020   Bugfix downward aa/cm/cp (H. Xie / A. Geer)
  !
  ! Code Description:
  !   Language:           Fortran 90.
  !   Software Standards: "European Standards for Writing and
  !     Documenting Exchangeable Fortran 90 Code".
  !
  ! Declarations:
  ! Modules used:
  ! Imported Type Definitions:

  Use rttov_types, Only :    &
       & rttov_profile_scatt_aux    ,&
       & rttov_geometry 

  Use parkind1, Only : jpim     ,jprb
!INTF_OFF
  Use rttov_const, Only: adk_adjoint, adk_k

  USE YOMHOOK, ONLY: LHOOK , DR_HOOK
!INTF_ON
  Implicit None

!* Subroutine arguments:
  Real    (Kind=jprb), Intent (in) :: ccthres
  Integer (Kind=jpim), Intent (in) :: nlevels      ! Number of levels
  Integer (Kind=jpim), Intent (in) :: nprofiles    ! Number of profiles
  Integer (Kind=jpim), Intent (in) :: nchannels    ! Number of channels*profiles=radiances
  Integer (Kind=jpim), Intent (in) :: nprofilesad  ! Number of profiles in adjoint
  Integer (Kind=jpim), Intent (in) :: lprofiles (nchannels) ! Profile indices

  Type (rttov_profile_scatt_aux), Intent (in)    :: scatt_aux          ! Auxiliary profile variables for RTTOV_SCATT
  Type (rttov_profile_scatt_aux), Intent (inout) :: scatt_aux_ad       ! Auxiliary profile variables for RTTOV_SCATT
  Type (rttov_geometry),          Intent (in)    :: angles (nprofiles) ! Zenith angles

  Real (Kind=jprb), Intent (in)   , dimension (nchannels,nlevels) :: dp       ! D+ for boundary conditions
  Real (Kind=jprb), Intent (in)   , dimension (nchannels,nlevels) :: dm       ! D- for boundary conditions
  Real (Kind=jprb), Intent (inout), dimension (nchannels,nlevels) :: j_do     ! Downward source terms 
  Real (Kind=jprb), Intent (inout), dimension (nchannels,nlevels) :: j_up     ! Upward source terms
  Real (Kind=jprb), Intent (inout), dimension (nchannels,nlevels) :: dp_ad    ! D+ for boundary conditions
  Real (Kind=jprb), Intent (inout), dimension (nchannels,nlevels) :: dm_ad    ! D- for boundary conditions
  Real (Kind=jprb), Intent (inout), dimension (nchannels,nlevels) :: j_do_ad  ! Downward source terms 
  Real (Kind=jprb), Intent (inout), dimension (nchannels,nlevels) :: j_up_ad  ! Upward source terms

!INTF_END

!* Local variables
  Real    (Kind=jprb) :: ja1, jb1, jc1, jd1, aa_up, aa_do, apm, bb, cp_do, cm_do, cp_up, cm_up, cpm, ztmp
  Real    (Kind=jprb) :: ja2, jb2, jc2, jd2
  Real    (Kind=jprb) :: ja_ad, jb_ad, jc_ad, jd_ad, aa_up_ad, aa_do_ad, apm_ad, bb_ad
  Real    (Kind=jprb) :: cp_do_ad, cm_do_ad, cp_up_ad, cm_up_ad, cpm_ad, ztmp_ad
  Integer (Kind=jpim) :: iprof, ichan
  Integer (Kind=jpim) :: iprofad, adk, ii
  Logical             :: lstable
    
  REAL(KIND=JPRB) :: ZHOOK_HANDLE

  !- End of header --------------------------------------------------------

  IF (LHOOK) CALL DR_HOOK('RTTOV_INTEGRATESOURCE_AD',0_jpim,ZHOOK_HANDLE)

  if (nprofilesad == nprofiles) then 
    adk = adk_adjoint   ! Adjoint mode
  else if (nprofilesad == nchannels) then
    adk = adk_k         ! K mode
  endif 

  !* Channels * Profiles
  do ii=1,nlevels
    do ichan = 1, nchannels
      iprof = lprofiles (ichan)
      if (adk == adk_adjoint) then
        iprofad = iprof
      else if (adk == adk_k) then
        iprofad = ichan
      endif

      if (ii >= scatt_aux % mclayer(ichan) .and. scatt_aux % cfrac (iprof) > ccthres ) then

        !* Reset
        aa_do_ad = 0.0_JPRB
        aa_up_ad = 0.0_JPRB
        apm_ad   = 0.0_JPRB
        bb_ad    = 0.0_JPRB
        cp_do_ad = 0.0_JPRB
        cm_do_ad = 0.0_JPRB
        cp_up_ad = 0.0_JPRB
        cm_up_ad = 0.0_JPRB
        cpm_ad   = 0.0_JPRB
        ja_ad  = 0.0_JPRB
        jb_ad  = 0.0_JPRB
        jc_ad  = 0.0_JPRB
        jd_ad  = 0.0_JPRB

        !* Coefficients
        apm   = 1.5_JPRB * scatt_aux % asm (ichan,ii) * scatt_aux % ssa (ichan,ii) &
             & * angles (iprof) % coszen * scatt_aux % b1 (ichan,ii) / scatt_aux % h (ichan,ii)
        aa_up = scatt_aux % b0 (ichan,ii) - apm
        aa_do = scatt_aux % b0 (ichan,ii) + apm
        bb  = scatt_aux % b1 (ichan,ii)
        cpm = 1.5_JPRB * scatt_aux % asm (ichan,ii) &
          & * angles (iprof) % coszen * scatt_aux % lambda (ichan,ii) / scatt_aux % h (ichan,ii)
        cp_up = dp (ichan,ii) * scatt_aux % ssa (ichan,ii) * (1.0_JPRB - cpm)
        cm_up = dm (ichan,ii) * scatt_aux % ssa (ichan,ii) * (1.0_JPRB + cpm)
        cp_do = dp (ichan,ii) * scatt_aux % ssa (ichan,ii) * (1.0_JPRB + cpm)
        cm_do = dm (ichan,ii) * scatt_aux % ssa (ichan,ii) * (1.0_JPRB - cpm)

        ja1  = 0.0_JPRB
        jb1  = 0.0_JPRB
        jc1  = 0.0_JPRB
        jd1  = 0.0_JPRB

        ja2  = 0.0_JPRB
        jb2  = 0.0_JPRB
        jc2  = 0.0_JPRB
        jd2  = 0.0_JPRB

        !* Downward radiance source terms
        !* FORWARD PART
        ja1   = 1.0_JPRB - scatt_aux % tau (ichan,ii)
        jb1   = angles (iprof) % coszen / scatt_aux % ext (ichan,ii) * (1.0_JPRB - scatt_aux % tau (ichan,ii)) &
               & - scatt_aux % tau (ichan,ii) * scatt_aux % dz (iprof,ii) 

        lstable = abs(scatt_aux % ext (ichan,ii) - scatt_aux % lambda (ichan,ii) * angles (iprof) % coszen) > 1E-7_JPRB
        if (lstable) then
          ztmp  = exp (scatt_aux % dz (iprof,ii) * (scatt_aux % lambda (ichan,ii) - scatt_aux % ext (ichan,ii) / &
            &angles (iprof) % coszen))
          jc1   = scatt_aux % ext (ichan,ii) &
            & / (scatt_aux % lambda (ichan,ii) * angles (iprof) % coszen - scatt_aux % ext (ichan,ii)) *&
            & (ztmp - 1.0_JPRB)
        else
          ! Numerically unstable case needs an alternative formulation, valid only for very small dz*(lambda-ext/coszen)
          jc1   = scatt_aux % ext (ichan,ii) * scatt_aux % dz (iprof,ii) / angles (iprof) % coszen
        endif

        ztmp  = exp (scatt_aux % dz (iprof,ii) * (scatt_aux % lambda (ichan,ii) + scatt_aux % ext (ichan,ii) / &
         &angles (iprof) % coszen))
        jd1   = scatt_aux % ext (ichan,ii) &
              & / (scatt_aux % lambda (ichan,ii) * angles (iprof) % coszen + scatt_aux % ext (ichan,ii)) * &
              &(1.0_JPRB - 1.0_JPRB / ztmp) 

        j_do (ichan,ii) = ja1 * aa_do + jb1 * bb + jc1 * cp_do + jd1 * cm_do

        !* Upward radiance source terms
        ja2   = 1.0_JPRB - scatt_aux % tau (ichan,ii)
        jb2   = scatt_aux % dz (iprof,ii) - angles (iprof) % coszen / scatt_aux % ext (ichan,ii) &
              & * (1.0_JPRB - scatt_aux % tau (ichan,ii)) 

        ztmp  = exp (scatt_aux % dz (iprof,ii) * scatt_aux % lambda (ichan,ii))
        jc2   = scatt_aux % ext (ichan,ii) / (scatt_aux % ext (ichan,ii) + scatt_aux % lambda (ichan,ii) &
               & * angles (iprof) % coszen) * (ztmp  - scatt_aux % tau (ichan,ii)) 
        if (lstable) then
          jd2 = scatt_aux % ext (ichan,ii) / (scatt_aux % ext (ichan,ii) - scatt_aux % lambda (ichan,ii) &
               & * angles (iprof) % coszen) * (1.0_JPRB / ztmp  - scatt_aux % tau (ichan,ii)) 
        else
          ! Numerically unstable case needs an alternative formulation, valid only for very small dz*(lambda-ext/coszen)
          jd2 = scatt_aux % ext (ichan,ii) * scatt_aux % dz (iprof,ii) / (ztmp * angles (iprof) % coszen) 
        endif

        j_up (ichan,ii) = ja2 * aa_up + jb2 * bb + jc2 * cp_up + jd2 * cm_up

        !* ADJOINT PART
        !* Upward radiance source terms

        ja_ad    = ja_ad    + j_up_ad (ichan,ii) * aa_up
        aa_up_ad = aa_up_ad + j_up_ad (ichan,ii) * ja2
        jb_ad    = jb_ad    + j_up_ad (ichan,ii) * bb
        bb_ad    = bb_ad    + j_up_ad (ichan,ii) * jb2
        jc_ad    = jc_ad    + j_up_ad (ichan,ii) * cp_up
        cp_up_ad = cp_up_ad + j_up_ad (ichan,ii) * jc2
        jd_ad    = jd_ad    + j_up_ad (ichan,ii) * cm_up
        cm_up_ad = cm_up_ad + j_up_ad (ichan,ii) * jd2
        
        j_up_ad (ichan,ii) = 0.0_JPRB

        if (lstable) then
          scatt_aux_ad % ext (ichan,ii) = scatt_aux_ad % ext (ichan,ii) &
               & - jd_ad * scatt_aux % lambda (ichan,ii) * angles (iprof) % coszen &
               & / ( (scatt_aux % ext (ichan,ii) - scatt_aux % lambda (ichan,ii) * angles (iprof) % coszen) ** 2) &
               & * (1.0_JPRB/ ztmp  - scatt_aux % tau (ichan,ii)) 
          scatt_aux_ad % lambda (ichan,ii) = scatt_aux_ad % lambda (ichan,ii) &
               & + jd_ad  * angles (iprof) % coszen * scatt_aux % ext (ichan,ii) &
               & / ( (scatt_aux % ext (ichan,ii) - scatt_aux % lambda (ichan,ii) * angles (iprof) % coszen) ** 2) &
               & * (1.0_JPRB/ ztmp  - scatt_aux % tau (ichan,ii)) 
          ztmp_ad  = -1.0_JPRB * jd_ad  / ztmp  / ztmp  * scatt_aux % ext (ichan,ii) &
               & / (scatt_aux % ext (ichan,ii) - scatt_aux % lambda (ichan,ii) * angles (iprof) % coszen) 
          scatt_aux_ad % tau (ichan,ii) = scatt_aux_ad % tau (ichan,ii) &
              & - jd_ad  * scatt_aux % ext (ichan,ii) &
               & / (scatt_aux % ext (ichan,ii) - scatt_aux % lambda (ichan,ii) * angles (iprof) % coszen) 
        else
          scatt_aux_ad % ext (ichan,ii) = scatt_aux_ad % ext (ichan,ii) &
               & + jd_ad * scatt_aux % dz (iprof,ii) / ztmp / angles (iprof) % coszen
          ztmp_ad  = -1.0_JPRB * jd_ad * scatt_aux % ext (ichan,ii) * scatt_aux % dz (iprof,ii) &
               & / ztmp  / ztmp / angles (iprof) % coszen
          scatt_aux_ad % dz (iprofad,ii) = scatt_aux_ad % dz (iprofad,ii) &
               & + jd_ad * scatt_aux % ext (ichan,ii) / ztmp / angles (iprof) % coszen
        endif
        jd_ad  = 0.0_JPRB
    
        scatt_aux_ad % ext (ichan,ii) = scatt_aux_ad % ext (ichan,ii) &
             & + jc_ad * (ztmp  - scatt_aux % tau (ichan,ii)) &
             & * scatt_aux % lambda (ichan,ii) * angles (iprof) % coszen & 
             & / ( (scatt_aux % ext (ichan,ii) + scatt_aux % lambda (ichan,ii) * angles (iprof) % coszen) ** 2)
        scatt_aux_ad % lambda (ichan,ii) = scatt_aux_ad % lambda (ichan,ii) &
             & - jc_ad * (ztmp  - scatt_aux % tau (ichan,ii)) &
             & * scatt_aux % ext (ichan,ii) * angles (iprof) % coszen  &
             & / ( (scatt_aux % ext (ichan,ii) + scatt_aux % lambda (ichan,ii) * angles (iprof) % coszen) ** 2) 
        ztmp_ad  = ztmp_ad  + jc_ad  * scatt_aux % ext (ichan,ii) &
             & / (scatt_aux % ext (ichan,ii) + scatt_aux % lambda (ichan,ii) * angles (iprof) % coszen) 
        scatt_aux_ad % tau (ichan,ii) = scatt_aux_ad % tau (ichan,ii) - jc_ad  * scatt_aux % ext (ichan,ii) &
             & / (scatt_aux % ext (ichan,ii) + scatt_aux % lambda (ichan,ii) * angles (iprof) % coszen) 
        jc_ad    = 0.0_JPRB
    
        scatt_aux_ad % dz     (iprofad,ii) = scatt_aux_ad % dz     (iprofad,ii) + ztmp_ad  * &
         &scatt_aux % lambda (ichan,ii) * ztmp 
        scatt_aux_ad % lambda (ichan,ii) = scatt_aux_ad % lambda (ichan,ii) + ztmp_ad  * &
         &scatt_aux % dz     (iprof,ii) * ztmp 
        ztmp_ad  = 0.0_JPRB

        scatt_aux_ad % dz  (iprofad,ii) = scatt_aux_ad % dz  (iprofad,ii) + jb_ad 
        scatt_aux_ad % ext (ichan,ii) = scatt_aux_ad % ext (ichan,ii) + jb_ad  / &
         &scatt_aux % ext (ichan,ii) / scatt_aux % ext (ichan,ii) &
                                   & * angles (iprof) % coszen * (1.0_JPRB - scatt_aux % tau (ichan,ii)) 
        scatt_aux_ad % tau (ichan,ii) = scatt_aux_ad % tau (ichan,ii) + jb_ad  * &
         &angles (iprof) % coszen / scatt_aux % ext (ichan,ii)
        jb_ad  = 0.0_JPRB

        scatt_aux_ad % tau (ichan,ii)  = scatt_aux_ad % tau (ichan,ii)  - ja_ad 
        ja_ad  = 0.0_JPRB

        !* Downward radiance source terms

        ja_ad    = ja_ad    + j_do_ad (ichan,ii) * aa_do
        aa_do_ad = aa_do_ad + j_do_ad (ichan,ii) * ja1
        jb_ad    = jb_ad    + j_do_ad (ichan,ii) * bb
        bb_ad    = bb_ad    + j_do_ad (ichan,ii) * jb1
        jc_ad    = jc_ad    + j_do_ad (ichan,ii) * cp_do
        cp_do_ad = cp_do_ad + j_do_ad (ichan,ii) * jc1
        jd_ad    = jd_ad    + j_do_ad (ichan,ii) * cm_do
        cm_do_ad = cm_do_ad + j_do_ad (ichan,ii) * jd1
        j_do_ad (ichan,ii) = 0.0_JPRB

        ztmp  = exp (scatt_aux % dz (iprof,ii) * (scatt_aux % lambda (ichan,ii) + scatt_aux % ext (ichan,ii) / &
         &angles (iprof) % coszen))

        scatt_aux_ad % ext (ichan,ii) = scatt_aux_ad % ext (ichan,ii) &
             & + jd_ad * (1.0_JPRB - 1.0_JPRB / ztmp ) * scatt_aux % lambda (ichan,ii) * angles (iprof) % coszen &
             & / ( (scatt_aux % lambda (ichan,ii) * angles (iprof) % coszen + scatt_aux % ext (ichan,ii)) ** 2)
        scatt_aux_ad % lambda (ichan,ii) = scatt_aux_ad % lambda (ichan,ii) &
             & - jd_ad * angles (iprof) % coszen * scatt_aux % ext (ichan,ii) * (1.0_JPRB - 1.0_JPRB / ztmp ) &
             & / ( (scatt_aux % lambda (ichan,ii) * angles (iprof) % coszen + scatt_aux % ext (ichan,ii)) ** 2)
        ztmp_ad  = jd_ad  * scatt_aux % ext (ichan,ii) / ztmp / ztmp &
             & / (scatt_aux % lambda (ichan,ii) * angles (iprof) % coszen + scatt_aux % ext (ichan,ii)) 
        jd_ad  = 0.0_JPRB
        
        scatt_aux_ad % dz (iprofad,ii) = scatt_aux_ad % dz (iprofad,ii) &
             & + ztmp_ad * (scatt_aux % lambda (ichan,ii) + scatt_aux % ext (ichan,ii) / &
             & angles (iprof) % coszen) * ztmp  
        scatt_aux_ad % lambda (ichan,ii) = scatt_aux_ad % lambda (ichan,ii) &
             & + ztmp_ad  * scatt_aux % dz (iprof,ii) * ztmp  
        scatt_aux_ad % ext    (ichan,ii) = scatt_aux_ad % ext    (ichan,ii) &
             & + ztmp_ad  * scatt_aux % dz (iprof,ii) / angles (iprof) % coszen * ztmp  
        ztmp_ad  = 0.0_JPRB

        if (lstable) then 
          ztmp = exp (scatt_aux % dz (iprof,ii) * (scatt_aux % lambda (ichan,ii) - scatt_aux % ext (ichan,ii) / &
              &angles (iprof) % coszen))
         
          scatt_aux_ad % ext    (ichan,ii) = scatt_aux_ad % ext    (ichan,ii) &
              & + jc_ad  / (scatt_aux % lambda (ichan,ii) * angles (iprof) % coszen - &
              &scatt_aux % ext (ichan,ii)) * (ztmp  - 1.0_JPRB) &
              & * (1.0_JPRB + scatt_aux % ext (ichan,ii) / (scatt_aux % lambda (ichan,ii) * &
              &angles (iprof) % coszen - scatt_aux % ext (ichan,ii)) ) 
          scatt_aux_ad % lambda (ichan,ii) = scatt_aux_ad % lambda (ichan,ii) &
              & - jc_ad  * angles (iprof) % coszen * scatt_aux % ext (ichan,ii) * (ztmp - 1.0_JPRB) &
              & / (scatt_aux % lambda (ichan,ii) * angles (iprof) % coszen - scatt_aux % ext (ichan,ii)) &
              & / (scatt_aux % lambda (ichan,ii) * angles (iprof) % coszen - scatt_aux % ext (ichan,ii))  
          ztmp_ad  = jc_ad  * scatt_aux % ext (ichan,ii) &
              & / (scatt_aux % lambda (ichan,ii) * angles (iprof) % coszen - scatt_aux % ext (ichan,ii))  
          jc_ad  = 0.0_JPRB

          scatt_aux_ad % dz (iprofad,ii) = scatt_aux_ad % dz (iprofad,ii) &
              & + ztmp_ad  * (scatt_aux % lambda (ichan,ii) - scatt_aux % ext (ichan,ii) / &
              &angles (iprof) % coszen) * ztmp  
          scatt_aux_ad % lambda(ichan,ii) = scatt_aux_ad % lambda( ichan,ii) &
              & + ztmp_ad  *  scatt_aux % dz     (iprof,ii) * ztmp  
          scatt_aux_ad % ext(ichan,ii) = scatt_aux_ad % ext(ichan,ii) &
              & - ztmp_ad  *  scatt_aux % dz     (iprof,ii) / angles (iprof) % coszen * ztmp  
          ztmp_ad  = 0.0_JPRB
        else
          scatt_aux_ad % ext (ichan,ii) = scatt_aux_ad % ext (ichan,ii) &
              & + jc_ad * scatt_aux % dz (iprof,ii)/ angles (iprof) % coszen
          scatt_aux_ad % dz (iprofad,ii) = scatt_aux_ad % dz (iprofad,ii) &
              & + jc_ad * scatt_aux % ext (ichan,ii) / angles (iprof) % coszen
          jc_ad  = 0.0_JPRB      
        endif
        
        scatt_aux_ad % ext (ichan,ii) = scatt_aux_ad % ext (ichan,ii) &
            & - jb_ad  / scatt_aux % ext (ichan,ii) / scatt_aux % ext (ichan,ii) * &
            &angles (iprof) % coszen * (1.0_JPRB - scatt_aux % tau (ichan,ii))
        scatt_aux_ad % tau(ichan,ii) = scatt_aux_ad % tau(ichan,ii)  &
            & - jb_ad  * (angles (iprof) % coszen / scatt_aux % ext (ichan,ii) + scatt_aux % dz (iprof,ii))
        scatt_aux_ad % dz (iprofad,ii) = scatt_aux_ad % dz (iprofad,ii) - jb_ad  * scatt_aux % tau (ichan,ii)
        jb_ad  = 0.0_JPRB

        scatt_aux_ad % tau (ichan,ii)  = scatt_aux_ad % tau (ichan,ii) - ja_ad
        ja_ad  = 0.0_JPRB

        dm_ad (ichan,ii) = dm_ad (ichan,ii) + cm_do_ad * scatt_aux % ssa (ichan,ii) * (1.0_JPRB - cpm)
        scatt_aux_ad % ssa (ichan,ii) = scatt_aux_ad % ssa (ichan,ii) + cm_do_ad * dm (ichan,ii) * (1.0_JPRB - cpm)
        cpm_ad = cpm_ad - cm_do_ad * dm (ichan,ii) * scatt_aux % ssa (ichan,ii)
        cm_do_ad = 0.0_JPRB

        dp_ad (ichan,ii) = dp_ad (ichan,ii) + cp_do_ad * scatt_aux % ssa (ichan,ii) * (1.0_JPRB + cpm)
        scatt_aux_ad % ssa (ichan,ii) = scatt_aux_ad % ssa (ichan,ii) + cp_do_ad * dp (ichan,ii) * (1.0_JPRB + cpm)
        cpm_ad = cpm_ad + cp_do_ad  * dp (ichan,ii) * scatt_aux % ssa (ichan,ii)
        cp_do_ad  = 0.0_JPRB

        dm_ad (ichan,ii) = dm_ad (ichan,ii) + cm_up_ad * scatt_aux % ssa (ichan,ii) * (1.0_JPRB + cpm)
        scatt_aux_ad % ssa (ichan,ii) = scatt_aux_ad % ssa (ichan,ii) + cm_up_ad * dm (ichan,ii) * (1.0_JPRB + cpm)
        cpm_ad = cpm_ad + cm_up_ad * dm (ichan,ii) * scatt_aux % ssa (ichan,ii)
        cm_up_ad = 0.0_JPRB

        dp_ad (ichan,ii) = dp_ad (ichan,ii) + cp_up_ad * scatt_aux % ssa (ichan,ii) * (1.0_JPRB - cpm)
        scatt_aux_ad % ssa (ichan,ii) = scatt_aux_ad % ssa (ichan,ii) + cp_up_ad * dp (ichan,ii) * (1.0_JPRB - cpm)
        cpm_ad = cpm_ad - cp_up_ad  * dp (ichan,ii) * scatt_aux % ssa (ichan,ii)
        cp_up_ad  = 0.0_JPRB

        scatt_aux_ad % asm    (ichan,ii) = scatt_aux_ad % asm    (ichan,ii) + cpm_ad * 1.5_JPRB &
                                       & * angles (iprof) % coszen * scatt_aux % lambda (ichan,ii) / scatt_aux % h (ichan,ii)
        scatt_aux_ad % lambda (ichan,ii) = scatt_aux_ad % lambda (ichan,ii) + cpm_ad * 1.5_JPRB &
                                       & * angles (iprof) % coszen * scatt_aux % asm    (ichan,ii) / scatt_aux % h (ichan,ii)
        scatt_aux_ad % h      (ichan,ii) = scatt_aux_ad % h      (ichan,ii) - cpm_ad * 1.5_JPRB &
                                       & * angles (iprof) % coszen * scatt_aux % asm (ichan,ii) &
                                       & * scatt_aux % lambda (ichan,ii) / scatt_aux % h (ichan,ii) / scatt_aux % h (ichan,ii)
        cpm_ad = 0.0_JPRB

        scatt_aux_ad % b1 (ichan,ii) = scatt_aux_ad % b1 (ichan,ii) + bb_ad
        bb_ad  = 0.0_JPRB

        scatt_aux_ad % b0 (ichan,ii) = scatt_aux_ad % b0 (ichan,ii) + aa_do_ad + aa_up_ad
        apm_ad                       = apm_ad                       + aa_do_ad - aa_up_ad
        aa_do_ad = 0._JPRB
        aa_up_ad = 0._JPRB

        scatt_aux_ad % asm (ichan,ii) = scatt_aux_ad % asm (ichan,ii) + apm_ad  * 1.5_JPRB * scatt_aux % ssa (ichan,ii) &
                                   & * angles (iprof) % coszen * scatt_aux % b1 (ichan,ii) / scatt_aux % h (ichan,ii)
        scatt_aux_ad % ssa (ichan,ii) = scatt_aux_ad % ssa (ichan,ii) + apm_ad  * 1.5_JPRB * scatt_aux % asm (ichan,ii) &
                                   & * angles (iprof) % coszen * scatt_aux % b1 (ichan,ii) / scatt_aux % h (ichan,ii)
        scatt_aux_ad % b1  (ichan,ii) = scatt_aux_ad % b1  (ichan,ii) + apm_ad  * 1.5_JPRB * scatt_aux % asm (ichan,ii) &
                                   & * scatt_aux % ssa (ichan,ii) * angles (iprof) % coszen / scatt_aux % h  (ichan,ii)
        scatt_aux_ad % h   (ichan,ii) = scatt_aux_ad % h   (ichan,ii) - apm_ad  * 1.5_JPRB * scatt_aux % asm (ichan,ii) &
                                   & * scatt_aux % ssa (ichan,ii) * angles (iprof) % coszen * scatt_aux % b1 (ichan,ii) &
                                   & / scatt_aux % h( ichan,ii) / scatt_aux % h (ichan,ii)
        apm_ad = 0._JPRB
      endif
    end do
  end do

  IF (LHOOK) CALL DR_HOOK('RTTOV_INTEGRATESOURCE_AD',1_jpim,ZHOOK_HANDLE)

End subroutine rttov_integratesource_ad
