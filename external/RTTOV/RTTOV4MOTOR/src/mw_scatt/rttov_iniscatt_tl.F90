!
Subroutine rttov_iniscatt_tl (&
      & errorstatus,       &! out
      & opts,              &! in
      & opts_scatt,        &! in
      & lreflectivity,     &! in
      & nlevels,           &! in
      & nchannels,         &! in
      & nprofiles,         &! in
      & chanprof,          &! in
      & frequencies,       &! in
      & profiles,          &! in  
      & profiles_tl,       &! in  
      & cld_profiles,      &! in 
      & cld_profiles_tl,   &! in 
      & coef_rttov,        &! in
      & coef_scatt,        &! in
      & transmission,      &! in
      & transmission_tl,   &! in
      & calcemiss,         &! in
      & usercfrac,         &! in
      & angles,            &! out
      & scatt_aux,         &! inout
      & scatt_aux_tl)       ! inout 

  !
  ! Description:
  ! Calculates some variables related to the input precipitation profile
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
  !
  ! Current Code Owner: SAF NWP
  !
  ! History:
  ! Version   Date     Comment
  ! -------   ----     -------
  !   1.0    09/2002      Initial version     (F. Chevallier)
  !   1.1    05/2003      RTTOV7.3 compatible (F. Chevallier)
  !   1.2    03/2004      Added polarimetry   (R. Saunders)
  !   1.3    08/2004      Polarimetry fixes   (U. O'Keeffe)
  !   1.4    11/2004      Clean-up            (P. Bauer)
  !   1.5    10/2005      Fixes for rttov8 indexing   (U. O'Keeffe)
  !   1.6    11/2005      Add errorstatus to arguments (J. Cameron)   
  !   1.7    09/2006      Use zccmax_tl instead of iccmax index (A. Doherty)
  !   1.8    11/2007      RTTOV9 / cleanup    (A. Geer) 
  !   1.9    03/2008      Revised cloud partitioning (A. Geer)  
  !   1.10   03/2009      Safety check on cloud fraction (A. Geer)
  !   1.11   11/2009      User may supply average cloud fraction (A. Geer)
  !   1.12   11/2009      RTTOV transmittances / optical depths come on levels (A Geer)
  !   1.13   11/2012      Remove "old cloud" option; new subroutine for hydrometeors (A Geer)
  !   1.14   11/2017      R/T now done with radiances, not Tb (A. Geer)
  !
  ! Code Description:
  !   Language:           Fortran 90.
  !   Software Standards: "European Standards for Writing and
  !   Documenting Exchangeable Fortran 90 Code".
  !
  ! Declaratiochannelsns:
  ! Modules used:
  ! Imported Type Definitions:

  Use rttov_types, Only :           &
       & rttov_coef                ,&
       & rttov_scatt_coef          ,&
       & rttov_transmission        ,&
       & rttov_geometry            ,&
       & rttov_profile_scatt_aux   ,&
       & rttov_profile             ,&
       & rttov_profile_cloud       ,&
       & rttov_chanprof            ,&
       & rttov_options             ,&
       & rttov_options_scatt

  Use parkind1, Only : jpim, jplm
!INTF_OFF
  Use rttov_types, Only : rttov_transmission_aux

  Use rttov_const, Only:      &
       & errorstatus_success, &
       & errorstatus_fatal,   &
       & gravity,             &
       & pressure_top,        &
       & rgc,                 &
       & mair,                &
       & tcosmic,             &
       & max_scatt_optical_depth, &
       & min_reflectivity

  Use parkind1, Only : jprb

  USE YOMHOOK, ONLY: LHOOK , DR_HOOK
!INTF_ON
  Implicit None

!* Subroutine arguments:
  Type(rttov_options), Intent (in) :: opts       ! RTTOV options
  Type(rttov_options_scatt), Intent (in) :: opts_scatt ! RTTOV_SCATT options
  Logical (Kind=jplm), Intent (in) :: lreflectivity    ! Computation of radar reflectivity, not TB
  Integer (Kind=jpim), Intent (in) :: nlevels    ! Number of levels
  Integer (Kind=jpim), Intent (in) :: nprofiles  ! Number of profiles
  Integer (Kind=jpim), Intent (in) :: nchannels  ! Number of channels*profiles=radiances
  Integer (Kind=jpim), Intent(out) :: errorstatus             ! Error return code
  Integer (Kind=jpim), Intent (in) :: frequencies (nchannels) ! Frequency indices
  Type(rttov_chanprof), Intent(in) :: chanprof    (nchannels) ! Channel and profile indices
  Logical (Kind=jplm), Intent (in) :: calcemiss   (nchannels) ! Emissivity flags
  Logical (Kind=jplm), Intent (in) :: usercfrac               ! User has supplied cloud fraction

  Type (rttov_profile),             Intent (in)    :: profiles        (nprofiles)   ! Atmospheric profiles
  Type (rttov_profile),             Intent (in)    :: profiles_tl     (nprofiles)   ! Atmospheric profiles
  Type (rttov_coef),                Intent (in)    :: coef_rttov                    ! RTTOV Coefficients
  Type (rttov_scatt_coef),          Intent (in)    :: coef_scatt                    ! RTTOV_SCATT Coefficients
  Type (rttov_profile_cloud),       Intent (in)    :: cld_profiles    (nprofiles)   ! Cloud profiles
  Type (rttov_profile_cloud),       Intent (in)    :: cld_profiles_tl (nprofiles)   ! Cloud profiles
  Type (rttov_transmission),        Intent (in)    :: transmission                  ! Transmittances and optical depths
  Type (rttov_transmission),        Intent (in)    :: transmission_tl               ! Transmittances and optical depths
  Type (rttov_geometry),            Intent (out)   :: angles          (nprofiles)   ! Zenith angles
  Type (rttov_profile_scatt_aux),   Intent (inout) :: scatt_aux                     ! Auxiliary profile variables
  Type (rttov_profile_scatt_aux),   Intent (inout) :: scatt_aux_tl                  ! Auxiliary profile variables

!INTF_END

!* Local variables
  Integer (Kind=jpim) :: ilayer, iprof, ichan, ichanid
  Real    (Kind=jprb) :: p1, p2, pm, p1_tl, p2_tl, pm_tl, dp2dz

  Real (Kind=jprb), Dimension (nprofiles,nlevels)   :: presf  ! Pressure levels [hPa]
  Real (Kind=jprb), Dimension (nprofiles,nlevels+1) :: presfh ! Half-level pressure levels [hPa]
  Real (Kind=jprb), Dimension (nprofiles,nlevels)   :: dztop  ! Thickness of top half of level [m]
  Real (Kind=jprb), Dimension (nprofiles,nlevels)   :: dzbot  ! Thickness of bottom half of level [m]  
  Real (Kind=jprb), Dimension (2:nlevels)           :: dzr    ! Thickness of RTTOV level [m]
  Real (Kind=jprb), Dimension (nchannels,nlevels)   :: od     ! Single layer optical depths on our levels
  Real (Kind=jprb), Dimension (nlevels+1)           :: od_rttov ! Single layer optical depths on RTTOV levels
  Real (Kind=jprb), Dimension (nprofiles,nlevels)   :: dztop_tl    
  Real (Kind=jprb), Dimension (nprofiles,nlevels)   :: dzbot_tl  
  Real (Kind=jprb), Dimension (2:nlevels)           :: dzr_tl  
  Real (Kind=jprb), Dimension (nchannels,nlevels)   :: od_tl       
  Real (Kind=jprb), Dimension (nlevels+1)           :: od_rttov_tl 
  Real (Kind=jprb), Dimension (nchannels)           :: zod_up_cld   ! Optical depth from top of the atmosphere 
  Real (Kind=jprb), Dimension (nprofiles,nlevels)   :: presf_tl     ! Pressure levels [hPa]
  Real (Kind=jprb), Dimension (nprofiles,nlevels+1) :: presfh_tl    ! Half-level pressure levels [hPa]
  Real (Kind=jprb), Dimension (nchannels)           :: zod_up_cld_tl    ! Optical depth from top of the atmosphere 
  Real (Kind=jprb), Dimension (nprofiles)           :: tbd_prof, tbd_prof_tl ! Boundary temperature, by profile
  
  Type (rttov_transmission_aux) :: transmissioncld, transmissioncld_tl             ! Clear+cloud transmittances with cloud

  Character (len=80) :: errMessage
  Character (len=18) :: NameOfRoutine = 'rttov_iniscatt_tl '
  
  REAL(KIND=JPRB) :: ZHOOK_HANDLE

#include "rttov_mieproc_tl.interface"
#include "rttov_iniedd_tl.interface"
#include "rttov_calcemis_mw.interface"
#include "rttov_calcemis_mw_tl.interface"
#include "rttov_hydro_tl.interface"
#include "rttov_setgeometry.interface"
#include "rttov_errorreport.interface"
  !- End of header --------------------------------------------------------

  if (lhook) call dr_hook('RTTOV_INISCATT_TL',0_jpim,zhook_handle)

  errorstatus = errorstatus_success

  allocate (transmissioncld    % thermal_path1)
  allocate (transmissioncld    % thermal_path1 % tau_surf (0:0,nchannels))
  allocate (transmissioncld_tl % thermal_path1)
  allocate (transmissioncld_tl % thermal_path1 % tau_surf (0:0,nchannels))

  dp2dz = - rgc / gravity / mair

  scatt_aux    % ext (:,:) = 0.0_JPRB
  scatt_aux    % ssa (:,:) = 0.0_JPRB
  scatt_aux    % asm (:,:) = 0.0_JPRB
  scatt_aux    % zef (:,:) = min_reflectivity
  scatt_aux_tl % ext (:,:) = 0.0_JPRB
  scatt_aux_tl % ssa (:,:) = 0.0_JPRB
  scatt_aux_tl % asm (:,:) = 0.0_JPRB
  scatt_aux_tl % zef (:,:) = 0.0_JPRB

  !* Security on user-defined pressures
  Do iprof = 1, nprofiles
     Do ilayer = 1, nlevels
        If (profiles (iprof) % p (ilayer) >= pressure_top) Then
            presf_tl (iprof,ilayer) = profiles_tl (iprof) % p (ilayer)
            presf    (iprof,ilayer) = profiles    (iprof) % p (ilayer)
        else
            presf_tl (iprof,ilayer) = 0.0_JPRB
            presf    (iprof,ilayer) = pressure_top
        Endif
     Enddo
     Do ilayer = 1, nlevels + 1
        If (cld_profiles(iprof) % ph (ilayer) >= pressure_top) Then
            presfh_tl (iprof,ilayer) = cld_profiles_tl (iprof) % ph (ilayer)
            presfh    (iprof,ilayer) = cld_profiles    (iprof) % ph (ilayer)
        else
            presfh_tl (iprof,ilayer) = 0.0_JPRB
            presfh    (iprof,ilayer) = pressure_top
        Endif
     Enddo
  Enddo

  !* Geometric variables
  Call rttov_setgeometry ( &
    & .FALSE._jplm,    & ! in (plane-parallel)
    & profiles,        & ! in
    & coef_rttov,      & ! in
    & angles)            ! out

  !* Temperature at layer boundaries (K) 
  do ichan = 1, nchannels
    iprof = chanprof(ichan) % prof
    scatt_aux_tl % tbd (ichan,nlevels+1) = profiles_tl (iprof) % s2m % t
    scatt_aux    % tbd (ichan,nlevels+1) = profiles (iprof) % s2m % t
    scatt_aux_tl % tbd (ichan,1)         = profiles_tl (iprof) % t (1)
    scatt_aux    % tbd (ichan,1)         = profiles (iprof) % t (1)
    scatt_aux_tl % tsfc (ichan)          = profiles_tl (iprof) % skin % t
    scatt_aux    % tsfc (ichan)          = profiles (iprof) % skin % t
    scatt_aux_tl % tcosmic (ichan)       = 0.0_JPRB
    scatt_aux    % tcosmic (ichan)       = tcosmic
  enddo
  
  do ilayer = 1, nlevels-1
    do iprof = 1, nprofiles
      p1_tl = presf_tl  (iprof,ilayer+1)
      p1    = presf     (iprof,ilayer+1)
      p2_tl = presf_tl  (iprof,ilayer  )
      p2    = presf     (iprof,ilayer  )
      pm_tl = presfh_tl (iprof,ilayer+1)
      pm    = presfh    (iprof,ilayer+1)

      tbd_prof_tl(iprof) = profiles_tl (iprof) % t (ilayer+1) + (profiles_tl (iprof) % t (ilayer)     &
                       & - profiles_tl (iprof) % t (ilayer+1)) / log(p2/p1) * log(pm/p1)              &
                       & + (profiles   (iprof) % t (ilayer)   -  profiles    (iprof) % t (ilayer+1))  &
                       & / (-1.0_JPRB *  log(p2/p1) * log(p2/p1) ) * (p2_tl / p2 - p1_tl / p1) * log(pm/p1) & 
                       & + (profiles   (iprof) % t (ilayer)   -  profiles    (iprof) % t  (ilayer+1)) &
                       & / log(p2/p1) * (pm_tl / pm - p1_tl / p1) 
      tbd_prof(iprof) = profiles (iprof) % t (ilayer+1) + (profiles (iprof) % t (ilayer) &
                    & - profiles (iprof) % t (ilayer+1)) / log(p2/p1) * log(pm/p1) 
    enddo
    do ichan = 1, nchannels
      iprof = chanprof(ichan) % prof
      scatt_aux_tl % tbd (ichan,ilayer+1) = tbd_prof_tl (iprof)
      scatt_aux    % tbd (ichan,ilayer+1) = tbd_prof    (iprof)
    enddo
  enddo

  ! Convert temperatures to "effective" - i.e. apply band corrections
  if(coef_rttov%ff_val_bc) then
    do ichan = 1, nchannels
      ichanid = chanprof(ichan) % chan
      scatt_aux_tl%tbd(ichan,:)   = coef_rttov%ff_bcs(ichanid) * scatt_aux_tl%tbd(ichan,:)
      scatt_aux_tl%tsfc(ichan)    = coef_rttov%ff_bcs(ichanid) * scatt_aux_tl%tsfc(ichan)    
      scatt_aux_tl%tcosmic(ichan) = coef_rttov%ff_bcs(ichanid) * scatt_aux_tl%tcosmic(ichan)    
      scatt_aux%tbd(ichan,:)   = coef_rttov%ff_bco(ichanid) + coef_rttov%ff_bcs(ichanid) * scatt_aux%tbd(ichan,:)
      scatt_aux%tsfc(ichan)    = coef_rttov%ff_bco(ichanid) + coef_rttov%ff_bcs(ichanid) * scatt_aux%tsfc(ichan)    
      scatt_aux%tcosmic(ichan) = coef_rttov%ff_bco(ichanid) + coef_rttov%ff_bcs(ichanid) * scatt_aux%tcosmic(ichan)           
    enddo
  endif

  !* Nadir heights (km)
  Do ilayer = nlevels, 1, -1
     Do iprof = 1, nprofiles
        p1_tl = presfh_tl (iprof,ilayer+1)
        p1    = presfh    (iprof,ilayer+1)
        p2_tl = presfh_tl (iprof,ilayer  )
        p2    = presfh    (iprof,ilayer  )
        pm_tl = presf_tl  (iprof,ilayer  )
        pm    = presf     (iprof,ilayer  )

        If (p1 <= p2) then
           errorstatus = errorstatus_fatal
           Write( errMessage, '( "iniscatt : problem with user-defined pressure layering")' )
           Call rttov_errorreport (errorstatus, errMessage, NameOfRoutine)
           if (lhook) call dr_hook('RTTOV_INISCATT',1_jpim,zhook_handle)
           Return
        End If

        scatt_aux_tl % dz (iprof,ilayer) = dp2dz * (Log(p2/p1) * profiles_tl (iprof) % t (ilayer) &
                         & + (p2_tl / p2 - p1_tl / p1) * profiles (iprof) % t (ilayer))        
        scatt_aux    % dz (iprof,ilayer) = dp2dz * Log(p2/p1) * profiles (iprof) % t (ilayer) 

        dzbot_tl (iprof,ilayer) = dp2dz * Log(pm/p1) * profiles_tl (iprof) % t (ilayer) &
                              & + dp2dz * (pm_tl / pm - p1_tl / p1) * profiles (iprof) % t (ilayer)
        dztop_tl (iprof,ilayer) = dp2dz * Log(p2/pm) * profiles_tl (iprof) % t (ilayer) &
                              & + dp2dz * (p2_tl / p2 - pm_tl / pm) * profiles (iprof) % t (ilayer)

        dzbot (iprof,ilayer) = dp2dz * Log(pm/p1) * profiles (iprof) % t (ilayer)
        dztop (iprof,ilayer) = dp2dz * Log(p2/pm) * profiles (iprof) % t (ilayer)

     Enddo
  Enddo
 
  !* Get single layer optical depths (at nadir and in hPa-1) and put onto model half levels
  Do ichan = 1, nchannels
     iprof = chanprof(ichan) % prof
        
     ! Top RTTOV level to space   
     od_rttov_tl (1)      = -1.0_jprb / transmission % tau_levels (1,ichan)  &
                        & * transmission_tl % tau_levels (1,ichan) 
     od_rttov (1)         = -1.0_jprb * log( transmission % tau_levels (1,ichan) )     
     Do ilayer = 2, nlevels 
        od_rttov_tl (ilayer) = transmission_tl % tau_levels (ilayer-1,ichan) &
                        & / transmission % tau_levels (ilayer-1,ichan)    &
                        & - transmission_tl % tau_levels (ilayer,ichan)   &
                        & / transmission % tau_levels (ilayer,ichan) 
        od_rttov (ilayer) = log( transmission % tau_levels (ilayer-1,ichan) ) &
                        & - log( transmission % tau_levels (ilayer,ichan) )
     Enddo
     ! Surface to bottom RTTOV (full pressure) level
     od_rttov_tl (nlevels+1) = transmission_tl % tau_levels (nlevels,ichan) &
                        & / transmission % tau_levels (nlevels,ichan)    &
                        & - transmission_tl % tau_total (ichan)   &
                        & / transmission % tau_total (ichan) 
     od_rttov (nlevels+1) = log( transmission % tau_levels (nlevels,ichan) ) &
                        & - log( transmission % tau_total (ichan) )

     ! RTTOV layer thickness
     Do ilayer = 2, nlevels  
       dzr(ilayer)    = dzbot(iprof,ilayer-1)    + dztop(iprof,ilayer)
       dzr_tl(ilayer) = dzbot_tl(iprof,ilayer-1) + dztop_tl(iprof,ilayer)
     Enddo
 
     ! Re-allocate optical depths between half pressure levels       
     od_tl (ichan,1)      = od_rttov_tl(1) & 
                        & + od_rttov_tl(2) * dzbot(iprof,1) / dzr(2) &
                        & + dzbot_tl(iprof,1) * od_rttov(2) / dzr(2) &
                        & - dzr_tl(2) * od_rttov(2) * dzbot(iprof,1) / dzr(2)**2
     od (ichan,1)         = od_rttov(1) &
                        & + od_rttov(2) * dzbot(iprof,1) / dzr(2)     
     Do ilayer = 2, nlevels - 1  
        od_tl (ichan,ilayer) = od_rttov_tl(ilayer) * dztop(iprof,ilayer) / dzr(ilayer)  &
                        & + dztop_tl(iprof,ilayer) * od_rttov(ilayer)    / dzr(ilayer)  &
                        & - dzr_tl(ilayer) * od_rttov(ilayer) * dztop(iprof,ilayer)     &
                        & / dzr(ilayer)**2                                              & 
                        & + od_rttov_tl(ilayer+1) * dzbot(iprof,ilayer) / dzr(ilayer+1) &
                        & + dzbot_tl(iprof,ilayer) * od_rttov(ilayer+1) / dzr(ilayer+1) &
                        & - dzr_tl(ilayer+1) * od_rttov(ilayer+1) * dzbot(iprof,ilayer) &
                        & / dzr(ilayer+1)**2
        od (ichan,ilayer) = od_rttov(ilayer)   * dztop(iprof,ilayer) / dzr(ilayer)   &
                        & + od_rttov(ilayer+1) * dzbot(iprof,ilayer) / dzr(ilayer+1) 
     Enddo
     od_tl (ichan,nlevels) = od_rttov_tl(nlevels)  * dztop(iprof,nlevels) / dzr(nlevels) &
                        & + dztop_tl(iprof,nlevels) * od_rttov(nlevels)   / dzr(nlevels) &
                        & - dzr_tl(nlevels) * od_rttov(nlevels) * dztop(iprof,nlevels)   &
                        & / dzr(nlevels)**2                                              &
                        & + od_rttov_tl(nlevels+1)
     od (ichan,nlevels)   = od_rttov(nlevels)  * dztop(iprof,nlevels) / dzr(nlevels) &
                        & + od_rttov(nlevels+1) 


  Enddo

  !* Optical depths in km-1 and at nadir
  Do ilayer = 1,nlevels
     Do ichan = 1, nchannels
        iprof = chanprof(ichan) % prof
     
        scatt_aux_tl % ext (ichan,ilayer) = od_tl (ichan,ilayer) &
                                        & / scatt_aux % dz (iprof,ilayer) * angles (iprof) % coszen &
                                        & - od    (ichan,ilayer) &
                                        & * scatt_aux_tl % dz(iprof,ilayer) / (scatt_aux % dz (iprof,ilayer) &
                                        & * scatt_aux % dz (iprof,ilayer)) &
                                        & * angles (iprof) % coszen 
        scatt_aux    % ext (ichan,ilayer) = od (ichan,ilayer)  &
                                       & / scatt_aux % dz (iprof,ilayer) * angles (iprof) % coszen 

        If (scatt_aux    % ext (ichan,ilayer) < 1.0e-10_JPRB) Then
            scatt_aux_tl % ext (ichan,ilayer) = 0.0_JPRB
            scatt_aux    % ext (ichan,ilayer) = 1.0e-10_JPRB
        Endif
     Enddo 
  Enddo

  !* Creation of the hydrometeor profiles
  Call rttov_hydro_tl (           &
       & opts_scatt%cc_threshold, &! in
       & nlevels,                 &! in
       & nprofiles,               &! in   
       & usercfrac,               &! in
       & lreflectivity,           &! in
       & presf,                   &! in 
       & presf_tl,                &! in 
       & profiles,                &! in
       & profiles_tl,             &! in
       & cld_profiles,            &! in
       & cld_profiles_tl,         &! in
       & opts_scatt,              &! in
       & coef_scatt,              &! in
       & scatt_aux,               &! inout 
       & scatt_aux_tl)             ! inout 

  !* Cloud/rain absorption/scattering parameters
  Call rttov_mieproc_tl (         &
       & nlevels,                 &! in
       & nchannels,               &! in
       & nprofiles,               &! in
       & frequencies,             &! in
       & chanprof%prof,           &! in
       & lreflectivity,           &! in
       & profiles,                &! in
       & cld_profiles,            &! in
       & cld_profiles_tl,         &! in
       & opts_scatt,              &! in
       & coef_scatt,              &! in
       & coef_rttov,              &! in
       & chanprof,                &! in
       & scatt_aux,               &! inout
       & scatt_aux_tl)             ! inout 
       
  !* Scattering parameters for Eddington RT
  Call rttov_iniedd_tl(           &
       & lreflectivity,           &! in
       & opts_scatt%cc_threshold, &! in
       & nlevels,                 &! in
       & nchannels ,              &! in
       & nprofiles ,              &! in
       & chanprof,                &! in
       & angles    ,              &! in
       & coef_rttov,              &! in
       & scatt_aux ,              &! inout
       & scatt_aux_tl)             ! inout 
       
  !* Surface emissivities
  zod_up_cld_tl (:) = 0.0_JPRB
  zod_up_cld    (:) = 0.0_JPRB
  
  Do ichan = 1, nchannels
     iprof = chanprof(ichan) % prof
     
     Do ilayer = 1, nlevels     
        zod_up_cld_tl (ichan) = zod_up_cld_tl (ichan) &
                            & + scatt_aux_tl % ext (ichan,ilayer) * scatt_aux    % dz (iprof,ilayer)  &
                            & + scatt_aux    % ext (ichan,ilayer) * scatt_aux_tl % dz (iprof,ilayer)
        zod_up_cld    (ichan) = zod_up_cld    (ichan) &
                            & + scatt_aux    % ext (ichan,ilayer) * scatt_aux    % dz (iprof,ilayer)  
     Enddo
     if (zod_up_cld (ichan) >= max_scatt_optical_depth) then
         zod_up_cld    (ichan) = max_scatt_optical_depth
         zod_up_cld_tl (ichan) =  0.0_JPRB              
     endif
     
     transmissioncld    % thermal_path1 % tau_surf (0,ichan) = &
                                         & Exp(-1.0_JPRB * zod_up_cld (ichan) / angles (iprof) % coszen)
     transmissioncld_tl % thermal_path1 % tau_surf (0,ichan) = &
                                         & -1.0_JPRB * zod_up_cld_tl (ichan)  / angles (iprof) % coszen &
                                         & * transmissioncld % thermal_path1 % tau_surf (0,ichan)
  Enddo

  Call rttov_calcemis_mw(         &
       & opts,                    &! in
       & profiles,                &! in
       & angles,                  &! in
       & coef_rttov,              &! in
       & chanprof,               &! in
       & transmissioncld,         &! in
       & calcemiss,               &! in
       & scatt_aux % ems_cld,     &! inout
       & scatt_aux % ref_cld,     &! out
       & errorstatus          )    ! inout 
       
  Call rttov_calcemis_mw_tl(      &
       & opts,                    &! in
       & profiles,                &! in
       & profiles_tl,             &! in
       & angles,                  &! in
       & coef_rttov,              &! in
       & chanprof,               &! in
       & transmissioncld,         &! in
       & transmissioncld_tl,      &! in
       & calcemiss,               &! in
       & scatt_aux_tl % ems_cld,  &! inout
       & scatt_aux_tl % ref_cld)   ! out 

  !* Hemispheric emissivity (= Fastem's effective emissivity)
  scatt_aux_tl % ems_bnd (:) = scatt_aux_tl % ems_cld (:)
  scatt_aux    % ems_bnd (:) = scatt_aux    % ems_cld (:)
  scatt_aux_tl % ref_bnd (:) = scatt_aux_tl % ref_cld (:)
  scatt_aux    % ref_bnd (:) = scatt_aux    % ref_cld (:)

  !* Deallocate
  Deallocate (transmissioncld    % thermal_path1 % tau_surf)
  Deallocate (transmissioncld    % thermal_path1)
  Deallocate (transmissioncld_tl % thermal_path1 % tau_surf)
  Deallocate (transmissioncld_tl % thermal_path1)

  if (lhook) call dr_hook('RTTOV_INISCATT_TL',1_jpim,zhook_handle)

End Subroutine rttov_iniscatt_tl
