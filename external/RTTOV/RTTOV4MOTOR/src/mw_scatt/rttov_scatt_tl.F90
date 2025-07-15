! Description:
!> @file
!!   Runs RTTOV-SCATT tangent linear (TL) model
!
!> @brief
!!   Runs RTTOV-SCATT tangent linear (TL) model
!!
!! @details
!!   Given input profile and cloud/hydrometeor profile structures and
!!   corresponding profile perturbations, computes the resulting radiance
!!   perturbation from the tangent linear (TL) to the direct model
!!   evaluated at the input profile.
!!
!!   If the optional reflectivity and reflectivity_tl arguments are supplied
!!   then the perturbations in radar reflectivities are calculated.
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
!!
!!   The methodology and validation is described in the following:
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
!! @param[out]    errorstatus      status on exit
!! @param[in]     opts_scatt       RTTOV-SCATT options to configure the simulations
!! @param[in]     nlevels          number of levels in profiles structure
!! @param[in]     chanprof         specifies channels and profiles to simulate
!! @param[in]     frequencies      frequency indices for each channel
!! @param[in]     profiles         input atmospheric profiles and surface variables
!! @param[in]     cld_profiles     input cloud and hydrometeor profiles
!! @param[in]     coef_rttov       RTTOV coefficients structure
!! @param[in]     coef_scatt       RTTOV-SCATT Hydrotable structure
!! @param[in]     calcemis         flags for internal RTTOV surface emissivity calculation
!! @param[in,out] emissivity       input/output surface emissivities
!! @param[in]     profiles_tl      input atmospheric profile and surface variable perturbations
!! @param[in]     cld_profiles_tl  input cloud and hydrometeor profile perturbations
!! @param[in,out] emissivity_tl    input/output surface emissivity perturbations
!! @param[in,out] radiance         output radiances and corresponding BTs
!! @param[in,out] radiance_tl      output radiance and BT perturbations
!! @param[in,out] reflectivity     output data for radar simulator, optional
!! @param[in,out] reflectivity_tl  output radar reflectivity perturbations, optional
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
Subroutine rttov_scatt_tl( &
     & errorstatus,        &! out
     & opts_scatt,         &! in
     & nlevels,            &! in
     & chanprof,           &! in
     & frequencies,        &! in
     & profiles,           &! in  
     & cld_profiles,       &! in
     & coef_rttov,         &! in
     & coef_scatt,         &! in
     & calcemis,           &! in
     & emissivity,         &! inout
     & profiles_tl,        &! in
     & cld_profiles_tl,    &! in
     & emissivity_tl,      &! inout
     & radiance,           &! inout
     & radiance_tl,        &! inout
     & reflectivity,       &! inout, optional
     & reflectivity_tl)     ! inout, optional

  ! History:
  ! Version   Date     Comment
  ! -------   ----     -------
  !   1.0    09/2002      Initial version     (F. Chevallier)
  !   1.1    05/2003      RTTOV7.3 compatible (F. Chevallier)
  !   1.2    03/2004      Added polarimetry   (R. Saunders)
  !   1.3    08/2004      Polarimetry fixes   (U. O'Keeffe)
  !   1.4    11/2004      Clean-up            (P. Bauer)
  !   1.5    07/2005      Polarimetry fixes   (U. O'Keeffe)
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
  !   1.15   06/2015      Allow foam_fraction input (L-F Meunier)
  !   1.16   11/2017      R/T can now be done with radiances, not Tb (A. Geer) 
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
       & rttov_options_scatt

  Use parkind1, Only : jpim, jplm
!INTF_OFF
  Use parkind1, Only : jprb

  Use rttov_types, Only :            &
       & rttov_geometry             ,&
       & rttov_profile_scatt_aux    ,&
       & rttov_transmission         ,&
       & rttov_options

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
  Type(rttov_options_scatt), Intent(in) :: opts_scatt              ! RTTOV-SCATT options
  Integer (Kind=jpim), Intent (in)      :: nlevels                 ! Number of levels
  Type(rttov_chanprof),Intent (in)      :: chanprof(:)             ! Indices
  Type(rttov_profile), Intent (in)      :: profiles(:)             ! Atmospheric profiles
  Integer (Kind=jpim), Intent (in)      :: frequencies (size(chanprof)) ! Frequency indices
  Integer (Kind=jpim), Intent (out)     :: errorstatus                  ! Error return flag

  Logical (Kind=jplm),    Intent (in)    :: calcemis      (size(chanprof))  ! Switch for emissivity calculation
  Type(rttov_emissivity), Intent (inout) :: emissivity    (size(chanprof))  ! Surface emissivity
  Type(rttov_emissivity), Intent (inout) :: emissivity_tl (size(chanprof))  ! Surface emissivity
  
  Type (rttov_profile),       Intent (in)    :: profiles_tl     (size(profiles))
  Type (rttov_coefs),         Intent (in)    :: coef_rttov                  ! RTTOV Coefficients
  Type (rttov_scatt_coef),    Intent (in)    :: coef_scatt                  ! RTTOV_SCATT Coefficients
  Type (rttov_profile_cloud), Intent (in)    :: cld_profiles    (size(profiles)) ! Cloud profiles
  Type (rttov_profile_cloud), Intent (in)    :: cld_profiles_tl (size(profiles))
  Type (rttov_radiance),      Intent (inout) :: radiance                    ! Radiances
  Type (rttov_radiance),      Intent (inout) :: radiance_tl

  Type (rttov_reflectivity), optional, Intent (inout) :: reflectivity    ! Optional for radar
  Type (rttov_reflectivity), optional, Intent (inout) :: reflectivity_tl ! Optional for radar

!INTF_END

#include "rttov_tl.interface"
#include "rttov_iniscatt_tl.interface"
#include "rttov_eddington_tl.interface"
#include "rttov_errorreport.interface"
#include "rttov_calcbt.interface"
#include "rttov_calcbt_tl.interface"

  Integer (Kind=jpim), target :: sa__mclayer    (size(chanprof))
  Integer (Kind=jpim), target :: sa_tl__mclayer (size(chanprof))
   
  Real (Kind=jprb), target :: t__tau_total     (size(chanprof))
  Real (Kind=jprb), target :: t__tau_levels    (nlevels,size(chanprof))
  Real (Kind=jprb), target :: t_tl__tau_total  (size(chanprof))
  Real (Kind=jprb), target :: t_tl__tau_levels (nlevels,size(chanprof))
  
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

  Real (Kind=jprb), target :: sa_tl__cfrac   (size(profiles))
  Real (Kind=jprb), target :: sa_tl__ems_bnd (size(chanprof))
  Real (Kind=jprb), target :: sa_tl__ref_bnd (size(chanprof))
  Real (Kind=jprb), target :: sa_tl__ems_cld (size(chanprof))
  Real (Kind=jprb), target :: sa_tl__ref_cld (size(chanprof)) 
  Real (Kind=jprb), target :: sa_tl__tbd (size(chanprof),nlevels+1)
  Real (Kind=jprb), target :: sa_tl__tsfc (size(chanprof))
  Real (Kind=jprb), target :: sa_tl__tcosmic (size(chanprof))
  
  Real (Kind=jprb), target :: sa_tl__delta  (size(chanprof),nlevels)
  Real (Kind=jprb), target :: sa_tl__tau    (size(chanprof),nlevels)
  Real (Kind=jprb), target :: sa_tl__int_tau(size(chanprof),nlevels)
  Real (Kind=jprb), target :: sa_tl__ext    (size(chanprof),nlevels)
  Real (Kind=jprb), target :: sa_tl__ssa    (size(chanprof),nlevels)
  Real (Kind=jprb), target :: sa_tl__asm    (size(chanprof),nlevels)
  Real (Kind=jprb), target :: sa_tl__zef    (size(chanprof),nlevels)
  Real (Kind=jprb), target :: sa_tl__lambda (size(chanprof),nlevels)
  Real (Kind=jprb), target :: sa_tl__h      (size(chanprof),nlevels) 
  Real (Kind=jprb), target :: sa_tl__b0     (size(chanprof),nlevels)
  Real (Kind=jprb), target :: sa_tl__b1     (size(chanprof),nlevels)
  Real (Kind=jprb), target :: sa_tl__bn     (size(chanprof),nlevels)
  Real (Kind=jprb), target :: sa_tl__btop   (size(chanprof))
  Real (Kind=jprb), target :: sa_tl__bsfc   (size(chanprof))
  
  Real (Kind=jprb), target :: sa_tl__dz     (size(profiles),nlevels)
  Real (Kind=jprb), target :: sa_tl__hydro  (size(profiles),nlevels,cld_profiles(1)%nhydro)

!* Local variables:
  Integer (Kind=jpim) :: nprofiles, nchannels
  Logical (Kind=jplm) :: lpolarimetric(size(chanprof)), lthermal(size(chanprof))
  Integer (Kind=jpim) :: iprof, ichan, ichan_act, ilayer
  Real    (Kind=jprb) :: rad_cld        (size(chanprof))            
  Real    (Kind=jprb) :: rad_cld_tl     (size(chanprof))  
  Real    (Kind=jprb) :: zlayers(nlevels) !altitude from ground in km
  Logical (Kind=jplm) :: lreflectivity

  Type (rttov_transmission)      :: transmission, transmission_tl
  Type (rttov_geometry)          :: angles (size(profiles))
  Type (rttov_profile_scatt_aux) :: scatt_aux, scatt_aux_tl

  Type(rttov_options)    :: opts

  Character (len=80) :: errMessage
  Character (len=15) :: NameOfRoutine = 'rttov_scatt_tl '
  
  Real (KIND=JPRB) :: ZHOOK_HANDLE        

  !- End of header --------------------------------------------------------
  
  if (lhook) call dr_hook('RTTOV_SCATT_TL',0_jpim,zhook_handle)

  lreflectivity = present(reflectivity) .and. present(reflectivity_tl)

  nprofiles = size(profiles)
  nchannels = size(chanprof)

  errorstatus = errorstatus_success

  transmission % tau_total  => t__tau_total
  transmission % tau_levels => t__tau_levels

  transmission_tl % tau_total  => t_tl__tau_total
  transmission_tl % tau_levels => t_tl__tau_levels

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

  scatt_aux_tl % cfrac    => sa_tl__cfrac
  scatt_aux_tl % ems_bnd  => sa_tl__ems_bnd
  scatt_aux_tl % ref_bnd  => sa_tl__ref_bnd
  scatt_aux_tl % ems_cld  => sa_tl__ems_cld
  scatt_aux_tl % ref_cld  => sa_tl__ref_cld
  scatt_aux_tl % tbd      => sa_tl__tbd
  scatt_aux_tl % tsfc     => sa_tl__tsfc
  scatt_aux_tl % tcosmic  => sa_tl__tcosmic
  scatt_aux_tl % mclayer  => sa_tl__mclayer
  scatt_aux_tl % delta    => sa_tl__delta
  scatt_aux_tl % tau      => sa_tl__tau
  scatt_aux_tl % ext      => sa_tl__ext
  scatt_aux_tl % ssa      => sa_tl__ssa
  scatt_aux_tl % asm      => sa_tl__asm
  scatt_aux_tl % int_tau  => sa_tl__int_tau
  scatt_aux_tl % zef      => sa_tl__zef
  scatt_aux_tl % lambda   => sa_tl__lambda
  scatt_aux_tl % h        => sa_tl__h
  scatt_aux_tl % b0       => sa_tl__b0
  scatt_aux_tl % b1       => sa_tl__b1
  scatt_aux_tl % bn       => sa_tl__bn
  scatt_aux_tl % btop     => sa_tl__btop
  scatt_aux_tl % bsfc     => sa_tl__bsfc
  scatt_aux_tl % dz       => sa_tl__dz
  scatt_aux_tl % hydro    => sa_tl__hydro

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

!*         1.   Gas absorption

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

  Call rttov_tl(                      &
     & errorstatus,                   &! out
     & chanprof,                      &! in
     & opts,                          &! in
     & profiles,                      &! in
     & profiles_tl,                   &! in
     & coef_rttov,                    &! in
     & transmission,                  &! inout
     & transmission_tl,               &! inout
     & radiance,                      &! inout
     & radiance_tl,                   &! inout
     & calcemis      = calcemis,      &! in
     & emissivity    = emissivity,    &! inout
     & emissivity_tl = emissivity_tl)  ! inout

  If ( errorstatus == errorstatus_fatal ) Then
     Write( errMessage, '( "error in rttov_tl")' )
     Call rttov_errorreport (errorstatus_fatal, errMessage, NameOfRoutine)
     IF (LHOOK) CALL DR_HOOK('RTTOV_SCATT_TL',1_jpim,ZHOOK_HANDLE)
     Return
  End If
 
  scatt_aux_tl % ems_cld (:) = emissivity_tl (:) % emis_in
  scatt_aux    % ems_cld (:) = emissivity    (:) % emis_in
  scatt_aux_tl % ref_cld (:) = -1.0_JPRB * emissivity_tl (:) % emis_in
  scatt_aux    % ref_cld (:) =  1.0_JPRB - emissivity    (:) % emis_in

!*  2.   Initialisations for Eddington
  Call rttov_iniscatt_tl(        &
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
        & profiles_tl,           &! in
        & cld_profiles,          &! in
        & cld_profiles_tl,       &! in
        & coef_rttov%coef,       &! in
        & coef_scatt,            &! in
        & transmission,          &! in
        & transmission_tl,       &! in
        & calcemis,              &! in
        & opts_scatt%lusercfrac, &! in
        & angles,                &! out
        & scatt_aux,             &! inout
        & scatt_aux_tl)           ! inout

  If ( errorstatus == errorstatus_fatal ) Then
     Write( errMessage, '( "error in rttov_iniscatt_tl")' )
     Call rttov_errorreport (errorstatus_fatal, errMessage, NameOfRoutine)
     IF (LHOOK) CALL DR_HOOK('RTTOV_SCATT_TL',1_jpim,ZHOOK_HANDLE)
     Return
  End If

  if (.not. lreflectivity) then

    !* 3.   Eddington (in temperature space)
    Call rttov_eddington_tl(   &
          & opts_scatt%cc_threshold, &! in
          & nlevels,                 &! in
          & nchannels,               &! in
          & nprofiles,               &! in
          & chanprof,                &! in
          & angles,                  &! in
          & scatt_aux,               &! in
          & scatt_aux_tl,            &! in
          & rad_cld,                 &! out
          & rad_cld_tl)               ! out

    !* 4.   Combine clear and cloudy parts
    do ichan = 1, nchannels
      iprof = chanprof (ichan) % prof
      if (scatt_aux % cfrac (iprof) > opts_scatt % cc_threshold .and. .not. lpolarimetric(ichan)) then

        radiance_tl % total (ichan) = &
          &   scatt_aux % cfrac (iprof)                    * rad_cld_tl (ichan)  &
          & + (rad_cld (ichan) - radiance % clear (ichan)) * scatt_aux_tl % cfrac (iprof) &
          & + (1.0_JPRB - scatt_aux % cfrac (iprof))       * radiance_tl % clear (ichan)

        radiance % total (ichan) = rad_cld (ichan) * scatt_aux % cfrac (iprof)   &
          & + radiance % clear (ichan) * (1.0_JPRB - scatt_aux % cfrac (iprof))

      else
        radiance_tl % total (ichan) = radiance_tl % clear (ichan)
        radiance    % total (ichan) = radiance    % clear (ichan)
      endif
    enddo

  else

    reflectivity % zef(:,:)             = 0._JPRB
    reflectivity % azef(:,:)            = 0._JPRB
    radiance % geometric_height(:,:)    = 0._JPRB
    radiance % total (:)                = min_radiance_radar
    radiance % clear (:)                = min_radiance_radar
    reflectivity_tl % zef(:,:)          = 0._JPRB
    reflectivity_tl % azef(:,:)         = 0._JPRB
    radiance_tl % geometric_height(:,:) = 0._JPRB
    radiance_tl % total (:)             = 0._JPRB
    radiance_tl % clear (:)             = 0._JPRB

    do ichan = 1, nchannels

      iprof            = chanprof(ichan)%prof
      zlayers(nlevels) = profiles(iprof) % elevation

      do ilayer = nlevels, 1, -1

        if(ilayer < nlevels) zlayers(ilayer) = zlayers(ilayer+1) + scatt_aux%dz(iprof,ilayer+1)

        if (scatt_aux % zef (ichan,ilayer) > min_reflectivity) then

          reflectivity    % zef (ilayer,ichan)  = scatt_aux    % zef (ichan,ilayer)
          reflectivity_tl % zef (ilayer,ichan)  = scatt_aux_tl % zef (ichan,ilayer)

          reflectivity    % azef (ilayer,ichan) = scatt_aux    % zef (ichan,ilayer) + &
                & 2.0_JPRB*10.0_JPRB*log10( scatt_aux % int_tau (ichan,ilayer) )
          reflectivity_tl % azef (ilayer,ichan) = scatt_aux_tl % zef (ichan,ilayer) + &
                & 2.0_JPRB*10.0_JPRB/log(10.0_JPRB)/scatt_aux % int_tau (ichan,ilayer) * scatt_aux_tl % int_tau (ichan,ilayer)

        else
          reflectivity    % zef  (ilayer,ichan) = min_reflectivity
          reflectivity_tl % zef  (ilayer,ichan) = 0.0_JPRB
          reflectivity    % azef (ilayer,ichan) = min_reflectivity
          reflectivity_tl % azef (ilayer,ichan) = 0.0_JPRB
        endif

        ! Approximate altitude of the middle of the RTTOV-SCATT layer (corresponding to the IFS full level), back from km to m
        radiance % geometric_height (ilayer,ichan) = 1000.0_JPRB * (zlayers(ilayer) + scatt_aux%dz(iprof,ilayer)/2.)

      enddo
    enddo
  endif

  lthermal=.true.
  call rttov_calcbt(chanprof, coef_rttov%coef, lthermal, radiance)
  call rttov_calcbt_tl(chanprof, coef_rttov%coef, lthermal, radiance, radiance_tl)

  if (lhook) call dr_hook('RTTOV_SCATT_TL',1_jpim,zhook_handle) 
  
End Subroutine rttov_scatt_tl
