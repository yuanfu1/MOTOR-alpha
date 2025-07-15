!      
Subroutine rttov_integratesource_tl (&        
     & ccthres,       &! in
     & nlevels,       &! in
     & nchannels,     &! in
     & nprofiles,     &! in
     & lprofiles,     &! in
     & angles,        &! in
     & scatt_aux,     &! in
     & scatt_aux_tl,  &! in
     & dp,            &! in
     & dp_tl,         &! in
     & dm,            &! in
     & dm_tl,         &! in
     & j_do,          &! inout
     & j_do_tl,       &! inout
     & j_up,          &! inout
     & j_up_tl)        ! inout 

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
  !  1.6       07/2008   Speed-ups / tidied  (A. Geer)
  !  1.7       03/2010   Use mclayer rather than min_ssa (A. Geer)
  !  1.8       11/2012   Trajectory multiplication order changed so as to reproduce AD exactly (A. Geer)
  !  1.9       01/2014   Numerical instability fixed (A. Geer)
  !  1.10      07/2020   Bugfix downward aa/cm/cp (H. Xie / A. Geer)
  !
  ! Code Description:
  !   Language:           Fortran 90.
  !   Software Standards: "European Standards for Writing and
  !   Documenting Exchangeable Fortran 90 Code".
  !
  ! Declarations:
  ! Modules used:
  ! Imported Type Definitions:

  Use rttov_types, Only :    &
       & rttov_profile_scatt_aux    ,&
       & rttov_geometry 

  Use parkind1, Only : jpim     ,jprb
!INTF_OFF
  USE YOMHOOK, ONLY: LHOOK , DR_HOOK
!INTF_ON
  Implicit None

!* Subroutine arguments:
  Real    (Kind=jprb), Intent (in) :: ccthres
  Integer (Kind=jpim), Intent (in) :: nlevels               ! Number of levels
  Integer (Kind=jpim), Intent (in) :: nprofiles             ! Number of profiles
  Integer (Kind=jpim), Intent (in) :: nchannels             ! Number of channels*profiles=radiances
  Integer (Kind=jpim), Intent (in) :: lprofiles (nchannels) ! Profile indices

  Type (rttov_profile_scatt_aux), Intent(in) :: scatt_aux         ! Auxiliary profile variables for RTTOV_SCATT
  Type (rttov_profile_scatt_aux), Intent(in) :: scatt_aux_tl      ! Auxiliary profile variables for RTTOV_SCATT
  Type (rttov_geometry),          Intent(in) :: angles (nprofiles)! Zenith angles

  Real (Kind=jprb), Intent (in)   , dimension (nchannels,nlevels) :: dp       ! D+ for boundary conditions
  Real (Kind=jprb), Intent (in)   , dimension (nchannels,nlevels) :: dm       ! D- for boundary conditions
  Real (Kind=jprb), Intent (inout), dimension (nchannels,nlevels) :: j_do     ! Downward source terms 
  Real (Kind=jprb), Intent (inout), dimension (nchannels,nlevels) :: j_up     ! Upward source terms
  Real (Kind=jprb), Intent (in)   , dimension (nchannels,nlevels) :: dp_tl    ! D+ for boundary conditions
  Real (Kind=jprb), Intent (in)   , dimension (nchannels,nlevels) :: dm_tl    ! D- for boundary conditions
  Real (Kind=jprb), Intent (inout), dimension (nchannels,nlevels) :: j_do_tl  ! Downward source terms 
  Real (Kind=jprb), Intent (inout), dimension (nchannels,nlevels) :: j_up_tl  ! Upward source terms

!INTF_END

!* Local variables
  Real    (Kind=jprb) :: ja   , jb   , jc   , jd   , aa_up   , aa_do   , apm   , bb
  Real    (Kind=jprb) :: ja_tl, jb_tl, jc_tl, jd_tl, aa_up_tl, aa_do_tl, apm_tl, bb_tl
  Real    (Kind=jprb) :: cp_do   , cm_do   , cp_up   , cm_up    , cpm    , ztmp
  Real    (Kind=jprb) :: cp_do_tl, cm_do_tl, cp_up_tl, cm_up_tl , cpm_tl , ztmp_tl
  Integer (Kind=jpim) :: iprof, ichan, ilayer
  Logical             :: lstable
  
  REAL(KIND=JPRB) :: ZHOOK_HANDLE

  !- End of header --------------------------------------------------------

  IF (LHOOK) CALL DR_HOOK('RTTOV_INTEGRATESOURCE_TL',0_jpim,ZHOOK_HANDLE)

  !* Channels * Profiles
  do ilayer=1,nlevels
    do ichan = 1, nchannels
      iprof = lprofiles (ichan)
    
      if (ilayer >= scatt_aux % mclayer(ichan) .and. scatt_aux % cfrac (iprof) > ccthres ) then

        !* Coefficients
        apm_tl   = 1.5_JPRB * angles( iprof) % coszen &
                & * (scatt_aux_tl % asm (ichan,ilayer) * scatt_aux    % ssa (ichan,ilayer) * scatt_aux    % b1 (ichan,ilayer)  &
                & +  scatt_aux    % asm (ichan,ilayer) * scatt_aux_tl % ssa (ichan,ilayer) * scatt_aux    % b1 (ichan,ilayer)  &
                & +  scatt_aux    % asm (ichan,ilayer) * scatt_aux    % ssa (ichan,ilayer) * scatt_aux_tl % b1 (ichan,ilayer)) &
                & /  scatt_aux % h (ichan,ilayer) &
                & - 1.5_JPRB * scatt_aux % asm (ichan,ilayer) * scatt_aux % ssa (ichan,ilayer) &
                & *  angles (iprof) % coszen * scatt_aux % b1 (ichan,ilayer) * scatt_aux_tl % h (ichan,ilayer) &
                & / (scatt_aux % h (ichan,ilayer) * scatt_aux % h (ichan,ilayer))
        aa_up_tl = scatt_aux_tl % b0  (ichan,ilayer) - apm_tl
        aa_do_tl = scatt_aux_tl % b0  (ichan,ilayer) + apm_tl

        apm   = 1.5_JPRB * scatt_aux % asm (ichan,ilayer) &
             & * scatt_aux % ssa (ichan,ilayer) * angles (iprof) % coszen &
             & * scatt_aux % b1 (ichan,ilayer) / scatt_aux % h (ichan,ilayer)
        aa_up = scatt_aux % b0 (ichan,ilayer) - apm
        aa_do = scatt_aux % b0 (ichan,ilayer) + apm

        bb_tl  = scatt_aux_tl % b1 (ichan,ilayer)
        bb     = scatt_aux    % b1 (ichan,ilayer)

        cpm_tl = 1.5_JPRB * angles (iprof) % coszen  &
                & * (scatt_aux_tl % asm (ichan,ilayer) * scatt_aux    % lambda (ichan,ilayer) / scatt_aux    % h( ichan,ilayer) &
                & +  scatt_aux    % asm (ichan,ilayer) * scatt_aux_tl % lambda (ichan,ilayer) / scatt_aux    % h (ichan,ilayer) &
                & -  scatt_aux    % asm (ichan,ilayer) * scatt_aux    % lambda (ichan,ilayer) * scatt_aux_tl % h (ichan,ilayer) &
                & / (scatt_aux    % h   (ichan,ilayer) * scatt_aux    % h      (ichan,ilayer)))
        cpm    = 1.5_JPRB * scatt_aux % asm (ichan,ilayer) * angles (iprof) % coszen &
                & * scatt_aux % lambda (ichan,ilayer) / scatt_aux % h (ichan,ilayer)

        cp_up_tl = (dp_tl (ichan,ilayer) * scatt_aux % ssa (ichan,ilayer) + dp (ichan,ilayer) * scatt_aux_tl % ssa (ichan,ilayer)) &
               & * (1.0_JPRB - cpm) &
               & - dp (ichan,ilayer) * scatt_aux % ssa (ichan,ilayer) * cpm_tl
        cp_up    = dp (ichan,ilayer) * scatt_aux % ssa (ichan,ilayer) * (1.0_JPRB - cpm)

        cm_up_tl = (dm_tl (ichan,ilayer) * scatt_aux % ssa (ichan,ilayer) + dm (ichan,ilayer) * scatt_aux_tl % ssa (ichan,ilayer)) &
               & * (1.0_JPRB + cpm)&
               & + dm (ichan,ilayer) * scatt_aux % ssa (ichan,ilayer) * cpm_tl
        cm_up    = dm (ichan,ilayer) * scatt_aux % ssa (ichan,ilayer) * (1.0_JPRB + cpm)

        cp_do_tl = (dp_tl (ichan,ilayer) * scatt_aux % ssa (ichan,ilayer) + dp (ichan,ilayer) * scatt_aux_tl % ssa (ichan,ilayer)) &
               & * (1.0_JPRB + cpm) &
               & + dp (ichan,ilayer) * scatt_aux % ssa (ichan,ilayer) * cpm_tl
        cp_do    = dp (ichan,ilayer) * scatt_aux % ssa (ichan,ilayer) * (1.0_JPRB + cpm)

        cm_do_tl = (dm_tl (ichan,ilayer) * scatt_aux % ssa (ichan,ilayer) + dm (ichan,ilayer) * scatt_aux_tl % ssa (ichan,ilayer)) &
               & * (1.0_JPRB - cpm)&
               & - dm (ichan,ilayer) * scatt_aux % ssa (ichan,ilayer) * cpm_tl
        cm_do    = dm (ichan,ilayer) * scatt_aux % ssa (ichan,ilayer) * (1.0_JPRB - cpm)

        !* Downward radiance source terms
        ja_tl  = -1.0_JPRB * scatt_aux_tl % tau (ichan,ilayer)
        ja     =  1.0_JPRB - scatt_aux    % tau (ichan,ilayer)
 
        jb_tl  = -1.0_JPRB * angles(iprof) % coszen &
                & * (scatt_aux_tl % ext (ichan,ilayer) / (scatt_aux % ext (ichan,ilayer) * scatt_aux % ext (ichan,ilayer))*&
                &(1.0_JPRB - scatt_aux % tau (ichan,ilayer)) &
                & +  scatt_aux_tl % tau (ichan,ilayer) /  scatt_aux % ext (ichan,ilayer)) &
                & -  scatt_aux_tl % tau (ichan,ilayer) *  scatt_aux % dz  (iprof,ilayer) - &
                & scatt_aux % tau (ichan,ilayer) * scatt_aux_tl % dz (iprof,ilayer) 
        jb     = angles (iprof) % coszen / scatt_aux % ext (ichan,ilayer) * (1.0_JPRB - scatt_aux % tau (ichan,ilayer)) &
                & - scatt_aux % tau (ichan,ilayer) * scatt_aux % dz (iprof,ilayer) 

        lstable = abs(scatt_aux % ext (ichan,ilayer) - scatt_aux % lambda (ichan,ilayer) * angles (iprof) % coszen) > 1E-7_JPRB
        if (lstable) then

          ztmp     = exp (scatt_aux % dz (iprof,ilayer) * (scatt_aux % lambda (ichan,ilayer) - scatt_aux % ext (ichan,ilayer) / &
           &angles (iprof) % coszen))
          ztmp_tl  = ztmp * (scatt_aux_tl % dz (iprof,ilayer) * (scatt_aux    % lambda (ichan,ilayer) - &
           &scatt_aux    % ext (ichan,ilayer) / angles (iprof) % coszen) &
                  &              + scatt_aux    % dz (iprof,ilayer) * (scatt_aux_tl % lambda (ichan,ilayer) - &
                  &scatt_aux_tl % ext (ichan,ilayer) / angles (iprof) % coszen)) 

          jc_tl = (scatt_aux_tl % ext (ichan,ilayer)    * scatt_aux % lambda (ichan,ilayer) &
              & -  scatt_aux_tl % lambda (ichan,ilayer) * scatt_aux % ext (ichan,ilayer)) &
              & * ( ztmp - 1.0_JPRB ) * angles (iprof) % coszen &
              & / ((scatt_aux % lambda (ichan,ilayer) * angles (iprof) % coszen - scatt_aux % ext (ichan,ilayer)) **2 )&
              & + ztmp_tl * scatt_aux % ext (ichan,ilayer) &
              & / (scatt_aux % lambda (ichan,ilayer) * angles (iprof) % coszen - scatt_aux % ext (ichan,ilayer))     
          jc    =  scatt_aux % ext (ichan,ilayer) &
              & / (scatt_aux % lambda (ichan,ilayer) * angles (iprof) % coszen - scatt_aux % ext (ichan,ilayer)) &
              & * (ztmp - 1.0_JPRB)
 
        else
        
          ! Numerically unstable case needs an alternative formulation, valid only for very small dz*(lambda-ext/coszen)
          jc_tl = ( scatt_aux_tl % ext (ichan,ilayer) * scatt_aux % dz (iprof,ilayer) &
                & + scatt_aux % ext (ichan,ilayer) * scatt_aux_tl % dz (iprof,ilayer) ) / angles (iprof) % coszen
          jc    = scatt_aux % ext (ichan,ilayer) * scatt_aux % dz (iprof,ilayer) / angles (iprof) % coszen 

        endif

        ztmp    = exp (scatt_aux % dz (iprof,ilayer) * (scatt_aux % lambda (ichan,ilayer) + scatt_aux % ext (ichan,ilayer) / &
                & angles (iprof) % coszen))
        ztmp_tl = ztmp * ( scatt_aux_tl % dz (iprof,ilayer) & 
              & * (scatt_aux % lambda (ichan,ilayer) + scatt_aux % ext (ichan,ilayer) / angles (iprof) % coszen) &
              & + scatt_aux % dz (iprof,ilayer) * (scatt_aux_tl % lambda (ichan,ilayer) &
              & + scatt_aux_tl % ext (ichan,ilayer) / angles (iprof) % coszen)) 

        jd_tl = (scatt_aux_tl % ext (ichan,ilayer)    * scatt_aux % lambda (ichan,ilayer) &
            & -  scatt_aux_tl % lambda (ichan,ilayer) * scatt_aux % ext (ichan,ilayer) )  &
            & * (1.0_JPRB - 1.0_JPRB / ztmp) * angles (iprof) % coszen &
            & / ((scatt_aux % lambda (ichan,ilayer) * angles (iprof) % coszen + scatt_aux % ext (ichan,ilayer)) ** 2) &
            & +  scatt_aux % ext    (ichan,ilayer) * ztmp_tl  / ztmp  / ztmp &
            & / (scatt_aux % lambda (ichan,ilayer) * angles (iprof) % coszen + scatt_aux % ext (ichan,ilayer))    
        jd    =    scatt_aux % ext  (ichan,ilayer) &
            & / (scatt_aux % lambda (ichan,ilayer) * angles (iprof) % coszen + scatt_aux % ext (ichan,ilayer)) &
            & * (1.0_JPRB - 1.0_JPRB / ztmp )
  
        j_do_tl (ichan,ilayer) = ja_tl * aa_do + ja * aa_do_tl + jb_tl * bb    + jb * bb_tl  &
                             & + jc_tl * cp_do + jc * cp_do_tl + jd_tl * cm_do + jd * cm_do_tl
        j_do    (ichan,ilayer) = ja * aa_do + jb * bb + jc * cp_do + jd * cm_do

        !* Upward radiance source terms

        ja_tl  = -1.0_JPRB * scatt_aux_tl % tau (ichan,ilayer)
        ja     =  1.0_JPRB - scatt_aux    % tau (ichan,ilayer)
       
        jb_tl  = angles (iprof) % coszen  &
                & * (scatt_aux_tl % ext (ichan,ilayer) / (scatt_aux % ext (ichan,ilayer) * scatt_aux    % ext (ichan,ilayer)) * &
                &(1.0_JPRB - scatt_aux % tau (ichan,ilayer)) &
                & +  scatt_aux_tl % tau (ichan,ilayer) /  scatt_aux % ext (ichan,ilayer)) + scatt_aux_tl % dz  (iprof,ilayer) 
        jb     =  scatt_aux    % dz  (iprof,ilayer) - angles (iprof) % coszen / scatt_aux % ext (ichan,ilayer) * &
         &(1.0_JPRB - scatt_aux % tau (ichan,ilayer)) 

        ztmp     = exp (scatt_aux % dz (iprof,ilayer) * scatt_aux % lambda (ichan,ilayer))
        ztmp_tl  = (scatt_aux_tl % dz (iprof,ilayer) * scatt_aux    % lambda (ichan,ilayer) &
             &   +  scatt_aux    % dz (iprof,ilayer) * scatt_aux_tl % lambda (ichan,ilayer)) * ztmp 

        jc_tl  =  (  scatt_aux_tl % ext   (ichan,ilayer) * scatt_aux % lambda (ichan,ilayer)   &
                & -  scatt_aux_tl % lambda(ichan,ilayer) * scatt_aux % ext    (ichan,ilayer) ) &
                & * ( ztmp - scatt_aux % tau (ichan,ilayer)) * angles(iprof) % coszen &
                & / ( (scatt_aux % ext (ichan,ilayer) + scatt_aux % lambda (ichan,ilayer) * angles (iprof) % coszen) ** 2) &
                & + (ztmp_tl - scatt_aux_tl % tau (ichan,ilayer)) * scatt_aux % ext (ichan,ilayer) &
                & / (scatt_aux % ext (ichan,ilayer) + scatt_aux % lambda (ichan,ilayer) * angles(iprof) % coszen ) 
        jc     = scatt_aux % ext (ichan,ilayer) / (scatt_aux % ext (ichan,ilayer) + scatt_aux % lambda (ichan,ilayer) &
                & * angles (iprof) % coszen) * (ztmp - scatt_aux % tau (ichan,ilayer)) 


        if(lstable) then

          jd_tl =  ( (scatt_aux_tl % ext (ichan,ilayer) * (-1.0_JPRB) * scatt_aux % lambda (ichan,ilayer) * angles (iprof) % coszen &
                & + scatt_aux_tl % lambda (ichan,ilayer) * scatt_aux % ext (ichan,ilayer) * angles (iprof) % coszen ) &
                & * (1.0_JPRB/ztmp - scatt_aux % tau (ichan,ilayer) ) &
                & / (scatt_aux % ext (ichan,ilayer) - scatt_aux % lambda (ichan,ilayer) * angles (iprof) % coszen) &
                & +  scatt_aux % ext (ichan,ilayer) * (-1.0_JPRB * ztmp_tl / ztmp / ztmp  - scatt_aux_tl % tau (ichan,ilayer))) &
                & / (scatt_aux % ext (ichan,ilayer) - scatt_aux % lambda (ichan,ilayer) * angles (iprof) % coszen) 
          jd    = scatt_aux    % ext (ichan,ilayer) / (scatt_aux % ext (ichan,ilayer) - scatt_aux % lambda (ichan,ilayer) &
                & * angles (iprof) % coszen) * (1.0_JPRB / ztmp  - scatt_aux % tau (ichan,ilayer)) 

        else

          ! Numerically unstable case needs an alternative formulation, valid only for very small dz*(lambda-ext/coszen)
          jd_tl = ( scatt_aux % dz (iprof,ilayer) * scatt_aux_tl % ext (ichan,ilayer) &
                & - scatt_aux % dz (iprof,ilayer) * scatt_aux % ext (ichan,ilayer) * ztmp_tl / ztmp &
                & + scatt_aux % ext (ichan,ilayer) * scatt_aux_tl % dz (iprof,ilayer) ) &
                & / (ztmp * angles (iprof) % coszen)
          jd    = scatt_aux % ext (ichan,ilayer) * scatt_aux % dz (iprof,ilayer) / (ztmp * angles (iprof) % coszen)
  
        endif

        j_up_tl (ichan,ilayer) = ja_tl * aa_up + ja * aa_up_tl + jb_tl * bb    + jb * bb_tl &
                             & + jc_tl * cp_up + jc * cp_up_tl + jd_tl * cm_up + jd * cm_up_tl
        j_up    (ichan,ilayer) = ja * aa_up + jb * bb + jc * cp_up + jd * cm_up

      end if
    end do
  end do
  
  IF (LHOOK) CALL DR_HOOK('RTTOV_INTEGRATESOURCE_TL',1_jpim,ZHOOK_HANDLE)

End subroutine rttov_integratesource_tl
