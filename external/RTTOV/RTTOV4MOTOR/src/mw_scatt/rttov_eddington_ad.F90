!
Subroutine rttov_eddington_ad ( & 
     & ccthres,           &! in
     & nlevels,           &! in
     & nchannels,         &! in
     & nprofiles,         &! in
     & nprofilesad,       &! in
     & chanprof,          &! in
     & angles,            &! in
     & scatt_aux,         &! in
     & scatt_aux_ad,      &! inout
     & rad_cld,           &! out
     & rad_cld_ad)         ! inout 

  ! Description:
  ! AD of routine 
  ! to compute Eddington approximation to RT
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
  ! - Bauer, P., 2002: Microwave radiative transfer modeling in clouds and 
  !     precipitation. Part I: Model description.
  !     NWP SAF Report No. NWPSAF-EC-TR-005, 27 pp.
  ! - Chevallier, F., and P. Bauer, 2003:
  !     Model rain and clouds over oceans: comparison with SSM/I observations.
  !     Mon. Wea. Rev., 131, 1240-1255.
  ! - Moreau, E., P. Bauer and F. Chevallier, 2002: Microwave radiative transfer 
  !     modeling in clouds and precipitation. Part II: Model evaluation.
  !     NWP SAF Report No. NWPSAF-EC-TR-006, 27 pp.
  !
  ! Current Code Owner: SAF NWP
  !
  ! History:
  ! Version   Date     Comment
  ! -------   ----     -------
  !  1.0       09/2002   Initial version     (P. Bauer, E. Moreau)
  !  1.1       05/2003   RTTOV7.3 compatible (F. Chevallier)
  !  1.2       03/2004   Added polarimetry   (R. Saunders)
  !  1.3       08/2004   Polarimetry fixes   (U. O'Keefe)
  !  1.4       11/2004   Clean-up            (P. Bauer)
  !  1.5       11/2005   Limit lines to 132 characters (J. Cameron)
  !  1.6       11/2007   RTTOV9 version      (A. Geer)
  !  1.7       07/2008   Speed up and make consistent minimum SSA (A. Geer)
  !  1.8       07/2008   Clear sky speed-ups (A. Geer)
  !  1.9       06/2009   Bug fixed in index (T. Wilhelmsson)
  !  1.10      03/2010   Use mclayer rather than min_ssa (A. Geer)
  !  1.11      11/2017   R/T now done with radiances, not Tb (A. Geer)
  !
  ! Code Description:
  !   Language:           Fortran 90.
  !   Software Standards: "European Standards for Writing and
  !   Documenting Exchangeable Fortran 90 Code".
  !
  ! Declarations:
  ! Modules used:
  ! Imported Type Definitions:

  Use rttov_types, Only :      &
       & rttov_chanprof       ,&
       & rttov_geometry       ,&
       & rttov_profile_scatt_aux

  Use parkind1, Only : jpim     ,jprb
!INTF_OFF
  ! Use rttov_const, Only: adk_adjoint, adk_k

  USE YOMHOOK, ONLY: LHOOK , DR_HOOK
!INTF_ON
  Implicit None


!* Subroutine arguments:
  Real    (Kind=jprb), Intent (in) :: ccthres
  Integer (Kind=jpim), Intent (in) :: nlevels      ! Number of levels
  Integer (Kind=jpim), Intent (in) :: nprofiles    ! Number of profiles
  Integer (Kind=jpim), Intent (in) :: nchannels    ! Number of channels*profiles=radiances
  Integer (Kind=jpim), Intent (in) :: nprofilesad  ! Number of profiles in adjoint variables
  type(rttov_chanprof), intent (in) :: chanprof(nchannels)  ! Channel and profile indices

  Type (rttov_geometry),          Intent (in)    :: angles       (nprofiles) ! Zenith angles
  Type (rttov_profile_scatt_aux), Intent (in)    :: scatt_aux                ! Auxiliary profile variables for RTTOV_SCATT
  Type (rttov_profile_scatt_aux), Intent (inout) :: scatt_aux_ad             ! Auxiliary profile variables for RTTOV_SCATT
  Real (Kind=jprb),               Intent (out)   :: rad_cld       (nchannels) ! Radiances
  Real (Kind=jprb),               Intent (inout) :: rad_cld_ad    (nchannels) ! Radiances

!INTF_END

#include "rttov_boundaryconditions.interface"
#include "rttov_integratesource.interface"
#include "rttov_boundaryconditions_ad.interface"
#include "rttov_integratesource_ad.interface"

!* Local variables
  Real (Kind=jprb), Dimension (nchannels,nlevels) :: dp      ! D+ for boundary conditions
  Real (Kind=jprb), Dimension (nchannels,nlevels) :: dm      ! D- for boundary conditions
  Real (Kind=jprb), Dimension (nchannels,nlevels) :: j_up    ! Upward radiance source terms
  Real (Kind=jprb), Dimension (nchannels,nlevels) :: j_do    ! Downward radiance source terms
  Real (Kind=jprb), Dimension (nchannels,nlevels) :: dp_ad   ! D+ for boundary conditions
  Real (Kind=jprb), Dimension (nchannels,nlevels) :: dm_ad   ! D- for boundary conditions
  Real (Kind=jprb), Dimension (nchannels,nlevels) :: j_up_ad ! Upward radiance source terms
  Real (Kind=jprb), Dimension (nchannels,nlevels) :: j_do_ad ! Downward radiance source terms


  Real (Kind=jprb), Dimension (nchannels,0:nlevels) :: irad_do1       ! Downward radiances 
  Real (Kind=jprb), Dimension (nchannels,0:nlevels) :: irad_do2       ! Downward radiances
  Real (Kind=jprb), Dimension (nchannels,nlevels+1) :: irad_up        ! Downward radiances 
  Real (Kind=jprb), Dimension (nchannels)              :: irad_sfc       ! Inward radiances at surface
  Real (Kind=jprb), Dimension (nchannels,nlevels+1) :: tau_t          ! Transmittances integrated over all levels
  Real (Kind=jprb), Dimension (nchannels)              :: ftop, ftop_ad  ! Downward radiances 
  Real (Kind=jprb), Dimension (nchannels,0:nlevels) :: irad_do2_ad    ! Downward radiances 
  Real (Kind=jprb), Dimension (nchannels,0:nlevels) :: irad_do1_ad    ! Downward radiances 
  Real (Kind=jprb), Dimension (nchannels,nlevels+1) :: irad_up_ad     ! Upward radiances 
  Real (Kind=jprb), Dimension (nchannels)              :: irad_sfc_ad    ! Inward radiances at surface
  Real (Kind=jprb), Dimension (nchannels,nlevels+1) :: tau_t_ad       ! Transmittances integrated over all levels

  Real (Kind=jprb), Dimension (nchannels,nlevels) :: j_part  ! Part of the linear in tau
  Real (Kind=jprb), Dimension (nchannels,nlevels) :: j_part_ad  

  Integer (Kind=jpim) :: ilayer, jlayer, iprof, ichan
  ! Integer (Kind=jpim) :: adk, iprofad

  REAL(KIND=JPRB) :: ZHOOK_HANDLE

  !- End of header --------------------------------------------------------

  IF (LHOOK) CALL DR_HOOK('RTTOV_EDDINGTON_AD',0_jpim,ZHOOK_HANDLE)        

  ! if (nprofilesad == nprofiles) then 
  !   adk = adk_adjoint   ! Adjoint mode
  ! else if (nprofilesad == nchannels) then
  !   adk = adk_k         ! K mode
  ! endif 

  j_up (:,:) = 0.0_JPRB
  j_do (:,:) = 0.0_JPRB
 
  !* Channels * Profiles      
  do ichan = 1, nchannels
    iprof = chanprof(ichan)%prof

    !* Top/bottom
    irad_sfc   (ichan) = scatt_aux % ems_cld (ichan) * scatt_aux % bsfc (ichan)

    if (scatt_aux % cfrac(iprof) > ccthres ) then

      !* Clear-sky source terms (but only actually used if SSA is neglible)
      !* Linear in tau for consistency with RTTOV9
      do ilayer = 1, nlevels
        if (ilayer < scatt_aux % mclayer(ichan)) then
         if ( scatt_aux % delta (ichan,ilayer) > 1E-8 ) then

           j_part (ichan,ilayer) = &
             &   (scatt_aux % b0 (ichan,ilayer) - scatt_aux % bn (ichan,ilayer)) &
             & * (1.0_JPRB  - scatt_aux % tau (ichan,ilayer)) &
             & / scatt_aux % delta (ichan,ilayer)

           j_up (ichan,ilayer)   = scatt_aux % bn (ichan,ilayer) &
             & - scatt_aux % b0 (ichan,ilayer) * scatt_aux % tau (ichan,ilayer) &
             & + j_part(ichan,ilayer)

           j_do (ichan,ilayer)   = scatt_aux % b0 (ichan,ilayer) &
             & - scatt_aux % bn (ichan,ilayer) * scatt_aux % tau (ichan,ilayer) &
             & - j_part(ichan,ilayer)

         else

           ! Linear in tau formula is computationally unstable for small delta - avoid
           j_up (ichan,ilayer) = scatt_aux % b0 (ichan,ilayer) &
            &                  * (1.0_JPRB - scatt_aux % tau (ichan,ilayer)) 
           j_do (ichan,ilayer) = scatt_aux % b0 (ichan,ilayer) &
            &                  * (1.0_JPRB - scatt_aux % tau (ichan,ilayer)) 

         endif
       endif
      enddo


      !* Downward radiance at cloud top
      irad_do1 (ichan,0) = scatt_aux % btop (ichan)

      Do ilayer = 1, scatt_aux % mclayer (ichan) - 1
        irad_do1 (ichan,ilayer) = irad_do1 (ichan,ilayer-1) &
          &                     * scatt_aux % tau (ichan,ilayer) + j_do (ichan,ilayer)
      End do
  
      ftop (ichan) = irad_do1 (ichan,scatt_aux % mclayer (ichan) - 1)
    endif
  End do

  !* Get D+, D- from boundary conditions

  Call rttov_boundaryconditions (&
    & ccthres,               &! in
    & nlevels,               &! in
    & nchannels,             &! in
    & chanprof%prof,         &! in
    & scatt_aux,             &! in
    & ftop,                  &! in
    & dp,                    &! out
    & dm)                     ! out 

  !* Integrate radiance source terms
  Call rttov_integratesource (&
    & ccthres,               &! in
    & nlevels,               &! in
    & nchannels,             &! in
    & nprofiles,             &! in
    & chanprof%prof,         &! in
    & angles,                &! in
    & scatt_aux,             &! in
    & dp,                    &! in
    & dm,                    &! in
    & j_do,                  &! inout
    & j_up)                   ! inout 

!* Integrate downward radiances/transmittance
  irad_do2 (:,0)         = scatt_aux % btop (:)
  irad_up  (:,nlevels+1) = irad_sfc   (:)
  tau_t    (:,nlevels+1) = 1.0_JPRB
  
  Do ilayer = 1, nlevels
     jlayer = nlevels + 1 - ilayer
  
     irad_do2 (:,ilayer) = irad_do2 (:,ilayer-1) * scatt_aux % tau (:,ilayer) + j_do (:,ilayer)
     irad_up  (:,jlayer) = irad_up  (:,jlayer+1) * scatt_aux % tau (:,jlayer) + j_up (:,jlayer)
     
     tau_t    (:,jlayer) = tau_t    (:,jlayer+1) * scatt_aux % tau (:,jlayer)
  Enddo

  rad_cld (:) = irad_up (:,1) + scatt_aux % ref_cld (:) * irad_do2 (:,nlevels) * tau_t (:,1)

!* ADJOINT PART
  irad_up_ad    (:,:) = 0.0_JPRB
  irad_do1_ad   (:,:) = 0.0_JPRB
  irad_do2_ad   (:,:) = 0.0_JPRB
  tau_t_ad      (:,:) = 0.0_JPRB
  j_up_ad       (:,:) = 0.0_JPRB
  j_do_ad       (:,:) = 0.0_JPRB 
  irad_sfc_ad   (:)   = 0.0_JPRB
  dp_ad         (:,:) = 0.0_JPRB
  dm_ad         (:,:) = 0.0_JPRB 
  ftop_ad       (:)   = 0.0_JPRB 

!* Integrate downward radiances/transmittance
  irad_up_ad (:,1) = irad_up_ad (:,1) + rad_cld_ad (:) 
  scatt_aux_ad % ref_cld (:) = scatt_aux_ad % ref_cld (:) &
    & + irad_do2 (:,nlevels) * tau_t (:,1) * rad_cld_ad (:) 
  irad_do2_ad (:,nlevels) = irad_do2_ad (:,nlevels) &
    & + scatt_aux % ref_cld (:) * tau_t (:,1) * rad_cld_ad (:)
  tau_t_ad (:,1)             = tau_t_ad (:,1) &
    & + scatt_aux % ref_cld (:) * irad_do2 (:,nlevels) * rad_cld_ad (:)
  rad_cld_ad (:)   = 0.0_JPRB
  
  Do ilayer = nlevels, 1, -1
     jlayer = nlevels + 1 - ilayer
  
     tau_t_ad           (:,jlayer+1) = tau_t_ad (:,jlayer+1) + scatt_aux % tau (:,jlayer) &
                                   & * tau_t_ad (:,jlayer)
     scatt_aux_ad % tau (:,jlayer)   = scatt_aux_ad % tau (:,jlayer) + tau_t (:,jlayer+1) &
                                   & * tau_t_ad (:,jlayer)
     tau_t_ad           (:,jlayer)   = 0.0_JPRB
     
     irad_up_ad         (:,jlayer+1) = irad_up_ad (:,jlayer+1) + scatt_aux % tau (:,jlayer) &
                                   & * irad_up_ad (:,jlayer)
     scatt_aux_ad % tau (:,jlayer)   = scatt_aux_ad % tau (:,jlayer) + irad_up  (:,jlayer+1) &
                                   & * irad_up_ad (:,jlayer)
     j_up_ad            (:,jlayer)   = j_up_ad (:,jlayer) + irad_up_ad (:,jlayer)
     irad_up_ad         (:,jlayer)   = 0.0_JPRB
     
     irad_do2_ad        (:,ilayer-1) = irad_do2_ad (:,ilayer-1) + scatt_aux % tau (:,ilayer) &
                                   & * irad_do2_ad (:,ilayer)
     scatt_aux_ad % tau (:,ilayer)   = scatt_aux_ad % tau (:,ilayer) + irad_do2 (:,ilayer-1) &
                                   & * irad_do2_ad (:,ilayer)
     j_do_ad            (:,ilayer)   = j_do_ad (:,ilayer) + irad_do2_ad (:,ilayer)
     irad_do2_ad        (:,ilayer)   = 0.0_JPRB     
  Enddo

  tau_t_ad (:,nlevels+1) = 0.0_JPRB

  irad_sfc_ad (:)              = irad_sfc_ad (:) + irad_up_ad (:,nlevels+1)
  irad_up_ad  (:,nlevels+1) = 0.0_JPRB 
  
  scatt_aux_ad % btop (:) = scatt_aux_ad % btop (:) + irad_do2_ad (:,0) 
  irad_do2_ad   (:,0) = 0.0_JPRB

  !* Get D+, D- from boundary conditions

  Call rttov_integratesource_ad (&
    & ccthres,               &! in
    & nlevels,               &! in
    & nchannels,             &! in
    & nprofiles,             &! in
    & nprofilesad,           &! in
    & chanprof%prof,         &! in
    & angles,                &! in
    & scatt_aux,             &! in
    & scatt_aux_ad,          &! in
    & dp,                    &! in
    & dp_ad,                 &! in
    & dm,                    &! in
    & dm_ad,                 &! in
    & j_do,                  &! inout
    & j_do_ad,               &! inout
    & j_up,                  &! inout
    & j_up_ad)                ! inout 

  Call rttov_boundaryconditions_ad (&
    & ccthres,               &! in
    & nlevels,               &! in
    & nchannels,             &! in
    & nprofiles,             &! in
    & nprofilesad,           &! in
    & chanprof%prof,         &! in
    & scatt_aux,             &! in
    & scatt_aux_ad,          &! in
    & ftop,                  &! in
    & ftop_ad,               &! in
    & dp,                    &! out
    & dp_ad,                 &! out
    & dm,                    &! out
    & dm_ad)                  ! out 

!* Channels * Profiles      
  do ichan = 1, nchannels
    iprof = chanprof(ichan)%prof

    ! if (adk == adk_adjoint) then
    !   iprofad = iprof  
    ! else if (adk == adk_k) then
    !   iprofad = ichan  
    ! endif

    if (scatt_aux % cfrac(iprof) > ccthres) then
    
      !* Downward radiance at cloud top
      irad_do1_ad (ichan,scatt_aux % mclayer (ichan) - 1) = &
      & irad_do1_ad (ichan,scatt_aux % mclayer (ichan) - 1) + ftop_ad (ichan)
      ftop_ad (ichan) = 0.0_JPRB
     
      Do ilayer = scatt_aux % mclayer (ichan) - 1, 1, -1
        irad_do1_ad        (ichan,ilayer-1) = irad_do1_ad (ichan,ilayer-1) + &
         & scatt_aux % tau (ichan,ilayer) * irad_do1_ad (ichan,ilayer)
        scatt_aux_ad % tau (ichan,ilayer)   = scatt_aux_ad % tau (ichan,ilayer) + &
         & irad_do1 (ichan,ilayer-1) * irad_do1_ad (ichan,ilayer)
        j_do_ad            (ichan,ilayer)   = j_do_ad (ichan,ilayer) + irad_do1_ad (ichan,ilayer)
        irad_do1_ad (ichan,ilayer)   = 0.0_JPRB
      End do
     
      scatt_aux_ad % btop (ichan) = scatt_aux_ad % btop (ichan) + irad_do1_ad (ichan,0)
      irad_do1_ad (ichan,0) = 0.0_JPRB         
     
      !* Clear-sky source terms

      j_part_ad(ichan,:) = 0.0_JPRB
      do ilayer = 1, nlevels
        if (ilayer < scatt_aux % mclayer(ichan)) then
          if ( scatt_aux % delta (ichan,ilayer) > 1E-8 ) then

            scatt_aux_ad % b0 (ichan,ilayer)  = scatt_aux_ad % b0 (ichan,ilayer) &
              & + j_do_ad (ichan,ilayer)
            scatt_aux_ad % bn (ichan,ilayer)  = scatt_aux_ad % bn (ichan,ilayer) &
              & - j_do_ad (ichan,ilayer) * scatt_aux % tau (ichan,ilayer)
            scatt_aux_ad % tau (ichan,ilayer)   = scatt_aux_ad % tau (ichan,ilayer)  &
              & - j_do_ad (ichan,ilayer) * scatt_aux % bn (ichan,ilayer)
            j_part_ad (ichan,ilayer) = j_part_ad(ichan,ilayer) - j_do_ad (ichan,ilayer)
            j_do_ad (ichan,ilayer)   = 0.0_JPRB

            scatt_aux_ad % bn (ichan,ilayer)  = scatt_aux_ad % bn (ichan,ilayer) &
             & + j_up_ad (ichan,ilayer)
            scatt_aux_ad % b0 (ichan,ilayer)  = scatt_aux_ad % b0 (ichan,ilayer) &
             & - j_up_ad (ichan,ilayer) * scatt_aux % tau (ichan,ilayer)
            scatt_aux_ad % tau (ichan,ilayer) = scatt_aux_ad % tau (ichan,ilayer)    &
              & - j_up_ad (ichan,ilayer) * scatt_aux % b0 (ichan,ilayer)
            j_part_ad (ichan,ilayer) = j_part_ad(ichan,ilayer) + j_up_ad (ichan,ilayer)
            j_up_ad (ichan,ilayer)   = 0.0_JPRB

            scatt_aux_ad % b0 (ichan,ilayer) = scatt_aux_ad % b0 (ichan,ilayer)      &
             & + j_part_ad (ichan,ilayer) * (1.0_JPRB  - scatt_aux % tau (ichan,ilayer)) &
             & / scatt_aux % delta (ichan,ilayer)
            scatt_aux_ad % bn (ichan,ilayer) = scatt_aux_ad % bn (ichan,ilayer)      &
             & - j_part_ad (ichan,ilayer) * (1.0_JPRB  - scatt_aux % tau (ichan,ilayer)) &
             & / scatt_aux % delta (ichan,ilayer)
            scatt_aux_ad % tau (ichan,ilayer)  = scatt_aux_ad % tau (ichan,ilayer)       &
             & - j_part_ad (ichan,ilayer)                                                &
             & * (scatt_aux % b0 (ichan,ilayer) - scatt_aux % bn (ichan,ilayer))         &
             & / scatt_aux % delta (ichan,ilayer)       
            scatt_aux_ad % delta (ichan,ilayer) = scatt_aux_ad % delta (ichan,ilayer)    &
             & - j_part_ad (ichan,ilayer) * j_part (ichan,ilayer)                        &
             & / scatt_aux % delta (ichan,ilayer)
            j_part_ad(ichan,ilayer) = 0.0_JPRB

          else
            ! Linear in tau formula is computationally unstable for small delta - avoid
            j_up_ad (ichan,ilayer) = j_up_ad (ichan,ilayer) + j_do_ad (ichan,ilayer)
            j_do_ad (ichan,ilayer) = 0.0_JPRB

            scatt_aux_ad % b0  (ichan,ilayer) = scatt_aux_ad % b0 (ichan,ilayer) &
             & + (1.0_JPRB - scatt_aux % tau (ichan,ilayer)) * j_up_ad (ichan,ilayer)
            scatt_aux_ad % tau (ichan,ilayer) = scatt_aux_ad % tau (ichan,ilayer)    &
             & - scatt_aux % b0 (ichan,ilayer) * j_up_ad (ichan,ilayer)
            j_up_ad (ichan,ilayer) = 0.0_JPRB
          endif
        endif
      enddo      
    endif 
             
    !* Top/bottom
    scatt_aux_ad % ems_cld (ichan) = scatt_aux_ad % ems_cld (ichan) &
      & + scatt_aux % bsfc (ichan) * irad_sfc_ad (ichan)
    scatt_aux_ad % bsfc (ichan) = scatt_aux_ad % bsfc (ichan) &
      & + scatt_aux % ems_cld (ichan) * irad_sfc_ad (ichan)
    irad_sfc_ad (ichan) = 0.0_JPRB

  End do

  if (lhook) call dr_hook('RTTOV_EDDINGTON_AD',1_jpim,zhook_handle)

End Subroutine rttov_eddington_ad
