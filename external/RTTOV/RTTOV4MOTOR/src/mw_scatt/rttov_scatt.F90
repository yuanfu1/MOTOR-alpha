! Description:
!> @file
!!   Runs RTTOV-SCATT direct model
!
!> @brief
!!   Runs RTTOV-SCATT direct model
!!
!! @details
!!   Computes microwave multi-channel radiances and brightness
!!   temperatures for many profiles per call in a cloudy and/or rainy sky.
!!   If the optional reflectivity argument is supplied, radar reflectivities
!!   are calculated.
!!   
!!   The same profiles structure as standard RTTOV is used to input the
!!   pressure, temperature, water vapour and surface variables. The additional
!!   cloud and hydrometeor profiles are input via the cld_profiles structure:
!!   these profiles specify the the constituent amounts on the "full" pressure
!!   levels defined in the profiles structure and they apply to a domain
!!   bounded by "half" pressure levels defined in the cld_profiles%ph(:) array.
!!   The half-pressure levels lie between the full pressure levels with the
!!   top-most half-pressure typically being zero (space) and the last half-
!!   pressure being the surface pressure, giving (nlevels+1) half-pressures
!!   in total.
!!
!!   RTTOV-SCATT allows a limited number of options to be set in the opts_scatt
!!   argument which configure the simulations.
!!
!!   The chanprof and frequencies arguments should be populated by calling the
!!   rttov_scatt_setupindex subroutine. Both should be allocated with size equal
!!   to the total number of channels being simulated.
!!
!!   The optional cfrac argument can be used to return the cloud fractions used
!!   by RTTOV-SCATT. These are calculated internally unless opts_scatt%lusercfrac
!!   is true in which case cfrac contains the input cloud fractions on exit.
!!
!!   If the emis_retrieval_terms argument is supplied the structure will be
!!   filled with data required for all-sky emissivity retrievals. See the
!!   rttov_scatt_emis_retrieval subroutine which performs these calculations.
!!
!!
!!   The RTTOV-SCATT methodology and validation is described in the following:
!!
!!   Bauer, P., 2002: Microwave radiative transfer modeling in clouds and precipitation.
!!     Part I: Model description.
!!     NWP SAF Report No. NWPSAF-EC-TR-005, 21 pp.
!!
!!   Moreau, E., P. Bauer and F. Chevallier, 2002: Microwave radiative transfer modeling in clouds and precipitation.
!!     Part II: Model evaluation.
!!     NWP SAF Report No. NWPSAF-EC-TR-006, 27 pp.
!!
!!   Chevallier, F., and P. Bauer, 2003:
!!     Model rain and clouds over oceans:comparison with SSM/I observations. Mon. Wea. Rev., 131, 1240-1255.
!!
!!   Smith, E. A., P. Bauer, F. S. Marzano, C. D. Kummerow, D. McKague, A. Mugnai, G. Panegrossi, 2002:
!!     Intercomparison of microwave radiative transfer models for precipitating clouds.
!!     IEEE Trans. Geosci. Remote Sens. 40, 541-549.
!!
!!   Bauer, P., Moreau, E., Chevallier, F., O'Keeffe, U., 2006:
!!     Multiple-scattering microwave radiative transfer for data assimilation applications
!!     Quart. J. R. Meteorol. Soc. 132, 1259-1281
!!
!!   Geer, A.J., Bauer, P. and O'Dell, C.W., 2009:
!!     A Revised Cloud Overlap Scheme for Fast Microwave Radiative Transfer in Rain and Cloud.
!!     Journal of Applied Meteorology and Climatology. 48, 2257-2270
!!
!! @param[out]    errorstatus             status on exit
!! @param[in]     opts_scatt              RTTOV-SCATT options to configure the simulations
!! @param[in]     nlevels                 number of levels in profiles structure
!! @param[in]     chanprof                specifies channels and profiles to simulate
!! @param[in]     frequencies             frequency indices for each channel
!! @param[in]     profiles                input atmospheric profiles and surface variables
!! @param[in]     cld_profiles            input cloud and hydrometeor profiles
!! @param[in]     coef_rttov              RTTOV coefficients structure
!! @param[in]     coef_scatt              RTTOV-SCATT Hydrotable structure
!! @param[in]     calcemis                flags for internal RTTOV surface emissivity calculation
!! @param[in,out] emissivity              input/output surface emissivities
!! @param[in,out] radiance                output radiances and corresponding BTs
!! @param[out]    cfrac                   cloud fraction used by RTTOV-SCATT, optional
!! @param[in,out] emis_retrieval_terms    output data for emissivity retrievals, optional
!! @param[in,out] reflectivity            output data for radar simulator, optional
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
!    Copyright 2015, EUMETSAT, All Rights Reserved.
!
Subroutine rttov_scatt(   &
      & errorstatus,       &! out
      & opts_scatt,        &! in
      & nlevels,           &! in
      & chanprof,          &! in
      & frequencies,       &! in
      & profiles,          &! in
      & cld_profiles,      &! in
      & coef_rttov,        &! in
      & coef_scatt,        &! in
      & calcemis,          &! in
      & emissivity,        &! inout
      & radiance,          &! inout 
      & cfrac,             &! out, optional
      & emis_retrieval_terms, &! inout, optional
      & reflectivity)       ! inout, optional

  ! History:
  ! Version   Date     Comment
  ! -------   ----     -------
  !   1.0    09/2002      Initial version     (F. Chevallier)
  !   1.1    05/2003      RTTOV7.3 compatible (F. Chevallier)
  !   1.2    03/2004      Added polarimetry   (R. Saunders)
  !   1.3    08/2004      Polarimetry fixes   (U. O'Keeffe)
  !   1.4    11/2004      Clean-up            (P. Bauer)
  !   1.5    07/2005      Polarimetry fixes after re-write (U O'Keeffe)
  !   1.6    11/2005      Add errorstatus to iniscatt arguments and use a temporary
  !                       radiance type for the calcpolarisation call (J Cameron)
  !   1.7    11/2007      RTTOV9 version      (A. Geer)
  !   1.8    07/2008      Clear sky speed-ups (A. Geer)
  !   1.9    10/2008      Revised cloud partitioning (A. Geer)
  !   1.10   11/2009      User may supply the average cloud fraction (A. Geer)
  !   1.11   11/2009      Adapted for RTTOV10 (A. Geer)
  !   1.12   04/2010      Tidied up after code cleaning (A. Geer)
  !   1.13   02/2010      Revised cloud partitioning is now the default (A. Geer)
  !   1.14   08/2010      Fix for polarimetric channels until they can be 
  !                       hanled properly by RTTOV_SCATT (W. Bell)
  !          12/2010      Return transmittances (A. Geer)
  !          11/2012      RTTOV11; remove "old cloud fraction" approach (A. Geer)
  !          02/2013      Emissivity estimation from satellite observation + rttov emis out (F. Baordo)
  !          06/2015      Allow foam_fraction input (L-F Meunier)
  !          11/2017      R/T can now be done with radiances, not Tb (A. Geer)
  !          02/2018      Emissivity retrieval terms returned through interface (A. Geer)
  !   2.0    10/2018      Flexible hydrometeors (A Geer)
  !
  Use rttov_types, Only :    &
       & rttov_coefs          ,&
       & rttov_scatt_coef     ,&
       & rttov_profile        ,&
       & rttov_profile_cloud  ,&
       & rttov_radiance       ,&
       & rttov_reflectivity   ,&
       & rttov_chanprof       ,&
       & rttov_emissivity     ,&
       & rttov_options_scatt  ,&
       & rttov_scatt_emis_retrieval_type

  Use parkind1, Only : jpim, jprb, jplm
!INTF_OFF
  Use rttov_types, Only :            &
       & rttov_geometry             ,&
       & rttov_profile_scatt_aux    ,&
       & rttov_transmission         ,&
       & rttov_radiance2            ,&
       & rttov_options,&
       & eddington_sfc_type

  Use rttov_const, Only :   &
       & errorstatus_success ,&
       & errorstatus_fatal, &
       & sensor_id_po, &
       & min_reflectivity, &
       & min_radiance_radar

  USE YOMHOOK, ONLY: LHOOK , DR_HOOK
!INTF_ON
  Implicit None


!* Subroutine arguments:
  Type(rttov_options_scatt), Intent(in) :: opts_scatt  ! RTTOV-SCATT options
  Integer (Kind=jpim), Intent (in)      :: nlevels     ! Number of levels
  Type(rttov_chanprof), Intent (in)     :: chanprof(:) ! Channel and profile indices
  Type(rttov_profile),  Intent (in)     :: profiles(:)
  Integer (Kind=jpim), Intent (in)      :: frequencies (size(chanprof)) ! Frequency indices
  Integer (Kind=jpim), Intent (out)     :: errorstatus                  ! Error return flag

  Logical (Kind=jplm),    Intent (in)    :: calcemis   (size(chanprof))  ! Switch for emmissivity calculation
  Type(rttov_emissivity), Intent (inout) :: emissivity (size(chanprof))  ! Surface emissivity

  Type (rttov_coefs),         Intent (in)    :: coef_rttov               ! RTTOV Coefficients
  Type (rttov_scatt_coef),    Intent (in)    :: coef_scatt               ! RTTOV_SCATT Coefficients
  Type (rttov_profile_cloud), Intent (in)    :: cld_profiles (size(profiles)) ! Cloud profiles
  Type (rttov_radiance),      Intent (inout) :: radiance                 ! Radiances
 
  Real (Kind=jprb), optional, Intent (out)  :: cfrac (size(profiles))  ! Cloud fraction (diagnostic)

  Type (rttov_scatt_emis_retrieval_type), optional, Intent (inout) :: emis_retrieval_terms ! Optional for all-sky emis retrievals
  
  Type (rttov_reflectivity), optional, Intent (inout) :: reflectivity ! Optional for radar

!INTF_END

#include "rttov_direct.interface"
#include "rttov_iniscatt.interface"
#include "rttov_eddington.interface"
#include "rttov_errorreport.interface"
#include "rttov_calcbt.interface"
#include "rttov_scatt_emis_terms.interface"

  Integer (Kind=jpim), target :: sa__mclayer (size(chanprof))
 
  Real (Kind=jprb), target :: t__tau_total  (size(chanprof))
  Real (Kind=jprb), target :: t__tau_levels (nlevels,size(chanprof))

  Real (Kind=jprb), target :: sa__cfrac   (size(profiles))
  Real (Kind=jprb), target :: sa__ems_bnd (size(chanprof))
  Real (Kind=jprb), target :: sa__ref_bnd (size(chanprof))
  Real (Kind=jprb), target :: sa__ems_cld (size(chanprof))
  Real (Kind=jprb), target :: sa__ref_cld (size(chanprof))
  
  Real (Kind=jprb), target :: sa__tbd (size(chanprof),nlevels+1)
  Real (Kind=jprb), target :: sa__tsfc (size(chanprof))
  Real (Kind=jprb), target :: sa__tcosmic (size(chanprof))
  
  Real (Kind=jprb), target :: sa__delta  (size(chanprof),nlevels)
  Real (Kind=jprb), target :: sa__tau    (size(chanprof),nlevels)
  Real (Kind=jprb), target :: sa__int_tau(size(chanprof),nlevels)
  Real (Kind=jprb), target :: sa__ext    (size(chanprof),nlevels)
  Real (Kind=jprb), target :: sa__ssa    (size(chanprof),nlevels)
  Real (Kind=jprb), target :: sa__asm    (size(chanprof),nlevels)
  Real (Kind=jprb), target :: sa__zef    (size(chanprof),nlevels)
  Real (Kind=jprb), target :: sa__lambda (size(chanprof),nlevels)
  Real (Kind=jprb), target :: sa__h      (size(chanprof),nlevels)
  
  Real (Kind=jprb), target :: sa__b0     (size(chanprof),nlevels)
  Real (Kind=jprb), target :: sa__b1     (size(chanprof),nlevels)
  Real (Kind=jprb), target :: sa__bn     (size(chanprof),nlevels)
  Real (Kind=jprb), target :: sa__btop   (size(chanprof))
  Real (Kind=jprb), target :: sa__bsfc   (size(chanprof))
  Real (Kind=jprb), target :: sa__dz     (size(profiles),nlevels)
  Real (Kind=jprb), target :: sa__hydro  (size(profiles),nlevels,cld_profiles(1)%nhydro)

  Real (Kind=jprb), target :: sf__t_down (size(chanprof))
  Real (Kind=jprb), target :: sf__t_up   (size(chanprof))
  Real (Kind=jprb), target :: sf__tau    (size(chanprof))

  Real (Kind=jprb), target :: sr__upclear(size(chanprof))
  Real (Kind=jprb), target :: sr__dnclear(size(chanprof))
  Real (Kind=jprb), target :: sr__refldnclear(size(chanprof))
  Real (Kind=jprb), target :: sr__up  (nlevels-1, size(chanprof))
  Real (Kind=jprb), target :: sr__down(nlevels-1, size(chanprof))
  Real (Kind=jprb), target :: sr__surf(nlevels-1, size(chanprof))
  
!* Local variables:
  Integer (Kind=jpim) :: nprofiles, nchannels
  Logical (Kind=jplm) :: lpolarimetric(size(chanprof)), lthermal(size(chanprof))
  Integer (Kind=jpim) :: iprof, ichan, ilayer
  Integer (Kind=jpim) :: ichan_act
  Real    (Kind=jprb) :: rad_cld     (size(chanprof))
  Real    (Kind=jprb) :: zlayers(nlevels) !altitude from ground in km
  Logical (Kind=jplm) :: lreflectivity
  
  ! Variables for emissivity calculations
  Type (eddington_sfc_type) :: sfc_terms ! Downward radiance source terms, Upward radiance source terms, Total transmittancs
  Type (rttov_radiance2)    :: radiance2               
  Type (rttov_transmission)      :: transmission
  Type (rttov_geometry)          :: angles (size(profiles))
  Type (rttov_profile_scatt_aux) :: scatt_aux  
                         
  Type(rttov_options)    :: opts

  Character (len=80) :: errMessage
  Character (len=15) :: NameOfRoutine = 'rttov_scatt '
  
  REAL(KIND=JPRB) :: ZHOOK_HANDLE
  
  !- End of header --------------------------------------------------------

  if (lhook) call dr_hook('RTTOV_SCATT',0_jpim,zhook_handle)

  lreflectivity = present(reflectivity)

  nprofiles = size(profiles)
  nchannels = size(chanprof)

  errorstatus = errorstatus_success
  
  transmission % tau_total  => t__tau_total
  transmission % tau_levels => t__tau_levels
 
  scatt_aux % cfrac    => sa__cfrac
  scatt_aux % ems_bnd  => sa__ems_bnd
  scatt_aux % ref_bnd  => sa__ref_bnd
  scatt_aux % ems_cld  => sa__ems_cld
  scatt_aux % ref_cld  => sa__ref_cld
  scatt_aux % tbd      => sa__tbd
  scatt_aux % tsfc     => sa__tsfc
  scatt_aux % tcosmic  => sa__tcosmic
  scatt_aux % mclayer  => sa__mclayer
  scatt_aux % delta    => sa__delta
  scatt_aux % tau      => sa__tau
  scatt_aux % ext      => sa__ext
  scatt_aux % ssa      => sa__ssa
  scatt_aux % asm      => sa__asm
  scatt_aux % int_tau  => sa__int_tau
  scatt_aux % zef      => sa__zef
  scatt_aux % lambda   => sa__lambda
  scatt_aux % h        => sa__h
  scatt_aux % b0       => sa__b0
  scatt_aux % b1       => sa__b1
  scatt_aux % bn       => sa__bn
  scatt_aux % btop     => sa__btop
  scatt_aux % bsfc     => sa__bsfc
  scatt_aux % dz       => sa__dz
  scatt_aux % hydro    => sa__hydro

  sfc_terms % down   => sf__t_down
  sfc_terms % up     => sf__t_up
  sfc_terms % tau    => sf__tau

  radiance2 % upclear => sr__upclear
  radiance2 % dnclear => sr__dnclear
  radiance2 % refldnclear => sr__refldnclear
  radiance2 % up      => sr__up  
  radiance2 % down    => sr__down
  radiance2 % surf    => sr__surf
   
  ! Check inputs
  ! ------------
  Do iprof = 1, nprofiles
    if (  profiles(iprof) % s2m % p /= cld_profiles(iprof) % ph(nlevels+1)  ) then
      errorstatus = errorstatus_fatal
      Write( errMessage, '( "Surface pressure and lowest half level should be identical")' )
    endif
    if ( cld_profiles(iprof) % nhydro /= coef_scatt % nhydro ) then
      errorstatus = errorstatus_fatal
      Write( errMessage, '( "Number of hydrometeors differs between inputs and scattering coefficients ")' )
    endif
    if ( cld_profiles(iprof) % nhydro_frac /= coef_scatt % nhydro .and. &
         cld_profiles(iprof) % nhydro_frac /= 1_JPIM ) then
      errorstatus = errorstatus_fatal
      Write( errMessage, '( "Number of hydrometeor fractions should be 1 or nhydro ")' )
    endif
    if(errorstatus == errorstatus_fatal) then
      Call rttov_errorreport (errorstatus_fatal, errMessage, NameOfRoutine)
      IF (LHOOK) CALL DR_HOOK('RTTOV_SCATT',1_jpim,ZHOOK_HANDLE)
      Return
    endif
  End do

  ! Identify polarimetric channels for fix to use only the clear sky part of RT calculation 
  do ichan = 1, nchannels
    ichan_act = chanprof(ichan)%chan
    lpolarimetric(ichan) = ( (coef_rttov % coef % id_sensor == sensor_id_po) &
                   & .and.   (coef_rttov % coef % fastem_polar(ichan_act) + 1_jpim >= 6_jpim) )
  enddo
      
  !* 1.   Gas absorption

  ! Profiles will be interpolated from model/RTTOV-SCATT levels to 
  ! RTTOV coefficient levels within RTTOV itself.
  opts%interpolation%addinterp   = .true.
  opts%rt_ir%addclouds           = .false.
  opts%rt_ir%addsolar            = .false.
  opts%rt_ir%addaerosl           = .false.
  opts%rt_ir%pc%addpc            = .false.
  opts%rt_ir%pc%addradrec        = .false.
  opts%rt_mw%clw_data            = .false.

  opts%rt_mw%fastem_version           = opts_scatt%fastem_version
  opts%rt_mw%supply_foam_fraction     = opts_scatt%supply_foam_fraction
  opts%rt_all%ozone_data              = opts_scatt%ozone_data
  opts%rt_all%use_t2m_opdep           = opts_scatt%use_t2m_opdep
  opts%rt_all%use_q2m                 = opts_scatt%use_q2m
  opts%rt_all%addrefrac               = opts_scatt%addrefrac
  opts%rt_all%rad_down_lin_tau        = opts_scatt%rad_down_lin_tau
  opts%rt_all%dtau_test               = opts_scatt%dtau_test
  opts%config                         = opts_scatt%config
  opts%interpolation%interp_mode      = opts_scatt%interp_mode
  opts%interpolation%reg_limit_extrap = opts_scatt%reg_limit_extrap
  opts%interpolation%lgradp           = opts_scatt%lgradp

  Call rttov_direct(               &!
         & errorstatus,            &! out
         & chanprof,               &! in
         & opts,                   &! in
         & profiles,               &! in
         & coef_rttov,             &! in
         & transmission,           &! inout
         & radiance,               &! inout
         & radiance2 = radiance2, &! inout 
         & calcemis   = calcemis,  &! in
         & emissivity = emissivity) ! inout

  If ( errorstatus == errorstatus_fatal ) Then
     Write( errMessage, '( "error in rttov_direct")' )
     Call rttov_errorreport (errorstatus_fatal, errMessage, NameOfRoutine)
     IF (LHOOK) CALL DR_HOOK('RTTOV_SCATT',1_jpim,ZHOOK_HANDLE)
     Return
  End If

  scatt_aux % ems_cld (:) = emissivity (:) % emis_in
  scatt_aux % ref_cld (:) = 1.0_JPRB - emissivity (:) % emis_in

  !* 2.   Initialisations for Eddington
  Call rttov_iniscatt(            &
         & errorstatus,           &! out
         & opts,                  &! in
         & opts_scatt,            &! in
         & lreflectivity,         &! in
         & nlevels,               &! in
         & nchannels,             &! in
         & nprofiles,             &! in
         & chanprof,              &! in
         & frequencies,           &! in
         & profiles,              &! in
         & cld_profiles,          &! in
         & coef_rttov%coef,       &! in
         & coef_scatt,            &! in
         & transmission,          &! in
         & calcemis,              &! in
         & opts_scatt%lusercfrac, &! in
         & angles,                &! out
         & scatt_aux)              ! inout

  If ( errorstatus == errorstatus_fatal ) Then
     Write( errMessage, '( "error in rttov_iniscatt")' )
     Call rttov_errorreport (errorstatus_fatal, errMessage, NameOfRoutine)
     IF (LHOOK) CALL DR_HOOK('RTTOV_SCATT',1_jpim,ZHOOK_HANDLE)
     Return
  End If

  if (.not. lreflectivity) then

    !* 3. Eddington (in temperature space)
    Call rttov_eddington(  &
          & opts_scatt%cc_threshold, &! in
          & nlevels,                 &! in
          & nchannels,               &! in
          & nprofiles,               &! in
          & chanprof,                &! in
          & angles,                  &! in
          & scatt_aux,               &! in
          & rad_cld,                 &! out
          & sfc_terms = sfc_terms)    ! inout, optional, Upward and downward radiance source terms, Total transmittances

    ! Emissivity retrieval terms
    if(present(emis_retrieval_terms)) then
      call rttov_scatt_emis_terms( &
          & opts_scatt%cc_threshold, &! in
          & chanprof,                &! in
          & coef_rttov,              &! in
          & scatt_aux,               &! in
          & emissivity,              &! in
          & transmission,            &! in
          & radiance2,               &! in
          & sfc_terms,               &! in
          & emis_retrieval_terms)     ! inout
    endif

    !* 4.  Combine clear and cloudy parts
    do ichan = 1, nchannels
      iprof = chanprof(ichan)%prof

      if (scatt_aux % cfrac (iprof) > opts_scatt % cc_threshold .and. .not. lpolarimetric(ichan)) then

        radiance % total (ichan) = rad_cld (ichan) * scatt_aux % cfrac (iprof) &
          & + radiance % clear (ichan) * (1.0_JPRB - scatt_aux % cfrac (iprof))

      else

        radiance % total (ichan) = radiance % clear (ichan)

      endif
    enddo

    ! Return the cloud fraction actually used - this is diagnostic output
    ! only provided by the forward model.
    if(present(cfrac)) then
      cfrac(:) = scatt_aux % cfrac (:)
    endif

  else

    reflectivity % zef(:,:)     = 0._JPRB
    reflectivity % azef(:,:)    = 0._JPRB
    radiance % geometric_height(:,:) = 0._JPRB
    radiance % total (:)        = min_radiance_radar
    radiance % clear (:)        = min_radiance_radar

    do ichan = 1, nchannels

      iprof            = chanprof(ichan)%prof
      zlayers(nlevels) = profiles(iprof) % elevation

      do ilayer = nlevels, 1, -1

        if(ilayer < nlevels) zlayers(ilayer) = zlayers(ilayer+1) + scatt_aux%dz(iprof,ilayer+1)

        if (scatt_aux % zef (ichan,ilayer) > min_reflectivity) then

          reflectivity % zef (ilayer,ichan)  = scatt_aux % zef (ichan,ilayer)

          reflectivity % azef (ilayer,ichan) = scatt_aux % zef (ichan,ilayer) + &
                & 2*10*log10( scatt_aux % int_tau (ichan,ilayer) )

        else
          reflectivity % zef (ilayer,ichan)  = min_reflectivity
          reflectivity % azef (ilayer,ichan) = min_reflectivity
        endif

        ! Approximate altitude of the middle of the RTTOV-SCATT layer (corresponding to the IFS full level), back from km to m
        radiance % geometric_height (ilayer,ichan) = 1000.0_JPRB * (zlayers(ilayer) + scatt_aux%dz(iprof,ilayer)/2.)

      enddo
    enddo

  endif

  lthermal=.true.
  call rttov_calcbt(chanprof, coef_rttov%coef, lthermal, radiance)

  if (lhook) call dr_hook('RTTOV_SCATT',1_jpim,zhook_handle)

End Subroutine rttov_scatt
