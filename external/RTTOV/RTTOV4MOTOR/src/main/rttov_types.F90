! Description:
!> @file
!!   Defines all derived types for RTTOV
!
!> @brief
!!   Defines all derived types for RTTOV
!!
!! @details
!!   This contains types that users will make use of in their code
!!   as well as types that RTTOV only uses internally.
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
MODULE rttov_types

  USE rttov_const, ONLY : &
      fastem_sp,          &
      ncldtyp,            &
      interp_rochon,      &
      ir_scatt_chou,      &
      vis_scatt_dom

  USE parkind1, ONLY : jpim, jprb, jplm

  IMPLICIT NONE

  ! ---------------------------------------------------------------------------
  ! User-level RTTOV structures
  ! ---------------------------------------------------------------------------

  !> Specify channels and profiles to simulate, declare as array of size
  !! nchanprof which is the total number of channels to simulate
  TYPE rttov_chanprof
    INTEGER(jpim) :: chan              !< Channel index
    INTEGER(jpim) :: prof              !< Profile index
  END TYPE

  !> Input/output surface emissivities, declare as array of size nchanprof,
  !! also used to specify per-channel specularity parameter for use with do_lambertian option
  TYPE rttov_emissivity
    REAL(jprb)    :: emis_in     = 0._jprb !< Input emissivity (0-1)
    REAL(jprb)    :: emis_out    = 0._jprb !< Output emissivity, value used by RTTOV (0-1)
    REAL(jprb)    :: specularity = 0._jprb !< Mix of specular/Lambertian for do_lambertian option (0-1, 0=Lambertian/1=specular)
  END TYPE

  !> Input/output surface reflectances, declare as array of size nchanprof,
  !! can also be used to specify cloud top BRDF for simple cloud scheme in VIS/NIR channels
  TYPE rttov_reflectance
    REAL(jprb)    :: refl_in          = 0._jprb !< Input direct sun-surface-satellite BRDF (>=0)
    REAL(jprb)    :: refl_out         = 0._jprb !< Output BRDF, value used by RTTOV (>=0)
    REAL(jprb)    :: diffuse_refl_in  = 0._jprb !< Input diffuse reflectance for downward scattered radiation (>=0)
    REAL(jprb)    :: diffuse_refl_out = 0._jprb !< Output diffuse reflectance, value used by RTTOV (>=0)
    REAL(jprb)    :: refl_cloud_top   = 0._jprb !< Optional, cloud top BRDF for simple cloud
  END TYPE

  !> Surface skin variables
  TYPE rttov_skin
    INTEGER(jpim) :: surftype          !< Surface type: 0=land, 1=sea, 2=sea-ice
    INTEGER(jpim) :: watertype         !< Water type: 0=fresh water, 1=ocean water
    REAL(jprb)    :: t                 !< Radiative skin temperature (K)
    REAL(jprb)    :: salinity          !< Practical ocean salinity unit (\%o) - FASTEM-4/5/6 only
    REAL(jprb)    :: foam_fraction     !< ocean foam fraction (0-1; only used if supply_foam_fraction is true)
    REAL(jprb)    :: snow_fraction     !< Snow coverage fraction for IR emissivity atlas (0-1)
    REAL(jprb)    :: soil_moisture     !< Soil moisture (m^3/m^3) - not currently used
    REAL(jprb)    :: fastem(fastem_sp) !< FASTEM land/sea-ice surface parameters
  END TYPE rttov_skin

  !> Surface 2m variables
  TYPE rttov_s2m
    REAL(jprb) :: t                    !< Temperature (K)
    REAL(jprb) :: q                    !< Water vapour (ppmv or kg/kg)
    REAL(jprb) :: o                    !< Ozone (ppmv or kg/kg) - not currently used
    REAL(jprb) :: p                    !< Surface pressure (hPa)
    REAL(jprb) :: u                    !< U 10m wind component (m/s)
    REAL(jprb) :: v                    !< V 10m wind component (m/s)
    REAL(jprb) :: wfetc                !< Wind fetch (metres)
  END TYPE rttov_s2m

  !> Atmospheric profiles on model pressure levels
  TYPE rttov_profile
    CHARACTER(LEN=128) :: id                       !< Optional profile ID string
    INTEGER(jpim) :: date(3)                       !< Year, Month, Day
    INTEGER(jpim) :: time(3)                       !< Hour, Minute, Second

    INTEGER(jpim) :: nlevels                       !< Number of atmospheric levels
    INTEGER(jpim) :: nlayers                       !< Number of atmospheric layers

    INTEGER(jpim) :: gas_units                     !< Units of gas profiles (0 or less => ppmv over dry air,
                                                   !! 1 => kg/kg over moist air, 2 => ppmv over moist air)

    REAL(jprb), POINTER :: p(:)          => NULL() !< Pressure (hPa)
    REAL(jprb), POINTER :: t(:)          => NULL() !< Temperature (K)
    REAL(jprb), POINTER :: q(:)          => NULL() !< Water vapour (ppmv or kg/kg)
    REAL(jprb), POINTER :: o3(:)         => NULL() !< O3 (ppmv or kg/kg)
    REAL(jprb), POINTER :: co2(:)        => NULL() !< CO2 (ppmv or kg/kg)
    REAL(jprb), POINTER :: n2o(:)        => NULL() !< N2O (ppmv or kg/kg)
    REAL(jprb), POINTER :: co(:)         => NULL() !< CO (ppmv or kg/kg)
    REAL(jprb), POINTER :: ch4(:)        => NULL() !< CH4 (ppmv or kg/kg)
    REAL(jprb), POINTER :: so2(:)        => NULL() !< SO2 (ppmv or kg/kg)
    REAL(jprb), POINTER :: clw(:)        => NULL() !< Cloud liquid water absorption only (kg/kg)

    LOGICAL(jplm)       :: mmr_cldaer              !< Cloud/aerosol units: False => num density cm-3 (aer), g.m-3 (cld);
                                                   !!                      True => kg/kg (cld, aer)
    REAL(jprb), POINTER :: aerosols(:,:) => NULL() !< Aerosol layer concentrations (units: see mmr_cldaer)
    REAL(jprb), POINTER :: cloud(:,:)    => NULL() !< Cloud liquid/ice water layer concentrations (units: see mmr_cldaer)
    REAL(jprb), POINTER :: cfrac(:)      => NULL() !< Layer cloud fraction (0-1)
    REAL(jprb), POINTER :: clwde(:)      => NULL() !< Cloud liquid water particle effective diameter (microns)
    INTEGER(jpim)       :: clwde_param             !< Cloud liquid water effective diameter parameterisation (1)
    INTEGER(jpim)       :: clw_scheme              !< Select liquid water cloud scheme (1-2)
    REAL(jprb), POINTER :: icede(:)      => NULL() !< Ice particle effective diameter (microns)
    INTEGER(jpim)       :: icede_param             !< Ice particle effective diameter parameterisation (1-4)
    INTEGER(jpim)       :: ice_scheme              !< Select ice cloud scheme (1-2)

    TYPE(rttov_skin) :: skin                       !< Surface skin variables
    TYPE(rttov_s2m)  :: s2m                        !< Surface 2m variables

    REAL(jprb) :: zenangle                         !< Satellite zenith angle (degrees)
    REAL(jprb) :: azangle                          !< Satellite azimuth angle (degrees)
    REAL(jprb) :: sunzenangle                      !< Solar azimuth angle (degrees)
    REAL(jprb) :: sunazangle                       !< Solar azimuth angle (degrees)
    REAL(jprb) :: elevation                        !< Surface elevation (km)
    REAL(jprb) :: latitude                         !< Latitude (degrees)
    REAL(jprb) :: longitude                        !< Longitude (degrees)

    REAL(jprb) :: Be                               !< Earth magnetic field strength (Gauss)
    REAL(jprb) :: cosbk                            !< Cosine of the angle between the Earth magnetic
                                                   !!   field and wave propagation direction

    REAL(jprb) :: ctp                              !< Black body (simple) cloud top pressure (hPa)
    REAL(jprb) :: cfraction                        !< Black body (simple) cloud fraction (0-1)
  END TYPE rttov_profile

  !> Additional atmospheric cloud/hydrometeor profile input for RTTOV-SCATT
  TYPE rttov_profile_cloud
    INTEGER(jpim) :: nlevels                         !< Number of atmospheric levels (same as in rttov_profile)
    INTEGER(jpim) :: nhydro                          !< Number of hydrometeor types
    INTEGER(jpim) :: nhydro_frac                     !< Number of hydrometeor fractions (should be 1 or nhydro)
    INTEGER(jpim), POINTER :: flux_conversion(:) => NULL() !< 0 => no flux conv, input units are kg/kg,
                                                           !< 1,2 => input kg/m2/s rain,snow
    REAL(jprb)    :: cfrac                           !< Average cloud fraction (only used if lusercfrac = TRUE, 0-1)
    REAL(jprb), POINTER :: ph(:)           => NULL() !< nlevels+1 of half-level model pressures (hPa)
    REAL(jprb), POINTER :: hydro(:,:)      => NULL() !< nlevels by ntypes of hydrometeor (kg/kg or (deprecated) kg/m2/s)
    REAL(jprb), POINTER :: hydro_frac(:,:) => NULL() !< nlevels of hydrometeor fraction (cloud cover) (0-1)
  END TYPE rttov_profile_cloud

  !> @internal Data for interpolating phase functions
  TYPE rttov_phasefn_int
    REAL(jprb)             :: zminphadiff
    REAL(jprb),    POINTER :: cosphangle(:) => NULL()
    INTEGER(jpim), POINTER :: iphangle(:)   => NULL()
  END TYPE rttov_phasefn_int

  !> Explicit optical parameters for IR scattering
  TYPE rttov_opt_param
    REAL(jprb),    POINTER :: abs(:,:)       => NULL() !< Absorption coef (nlayers,nchannels) (km-1)
    REAL(jprb),    POINTER :: sca(:,:)       => NULL() !< Scattering coef (nlayers,nchannels) (km-1)
    REAL(jprb),    POINTER :: bpr(:,:)       => NULL() !< b parameter (nlayers,nchannels) (no units)
    INTEGER(jpim)          :: nmom                     !< Number of Leg. coefs provided for each phase fn
    REAL(jprb),    POINTER :: legcoef(:,:,:) => NULL() !< Phase fn Leg. coefs (1:nmom+1,nlayers,nchannels)
                                                       !! - required for DOM scattering calculations
    REAL(jprb),    POINTER :: pha(:,:,:)     => NULL() !< Phase function (nphangle,nlayers,nchannels), should be
                                                       !! normalised such that integral over all angles is 4*pi
                                                       !! - required for solar scattering calculations
    REAL(jprb),    POINTER :: phangle(:)     => NULL() !< Angles over which phase fns defined (nphangle) (degrees)
                                                       !! - required for solar scattering calculations

    ! The following are for RTTOV internal purposes only
    TYPE(rttov_phasefn_int) :: phasefn_int             !< Data for interpolating phase functions - for internal use only
  END TYPE rttov_opt_param

  !> Output transmittances
  TYPE rttov_transmission
    REAL(jprb), POINTER  :: tau_total(:)             => NULL() !< Surface-satellite transmittance (channels)
    REAL(jprb), POINTER  :: tau_levels(:,:)          => NULL() !< Level-satellite transmittance (levels,channels)
    REAL(jprb), POINTER  :: tausun_total_path2(:)    => NULL() !< Sun-surface-satellite solar transmittance
    REAL(jprb), POINTER  :: tausun_levels_path2(:,:) => NULL() !< Sun-level-satellite solar transmittance for each level
    REAL(jprb), POINTER  :: tausun_total_path1(:)    => NULL() !< Surface-satellite solar transmittance
    REAL(jprb), POINTER  :: tausun_levels_path1(:,:) => NULL() !< Level-satellite solar transmittance for each level
    REAL(jprb), POINTER  :: tau_total_cld(:)         => NULL() !< Surface-satellite cloud-only transmittance (channels)
    REAL(jprb), POINTER  :: tau_levels_cld(:,:)      => NULL() !< Level-satellite cloud-only transmittance (levels,channels)
  END TYPE rttov_transmission

  !> Output radiances and corresponding brightness temperatures and reflectances (BRFs)
  !! Radiance unit: mW/m2/sr/cm-1; BT unit: K; BRFs are unitless.
  !! Array sizes are (nchannels) or (nlayers,nchannels)
  TYPE rttov_radiance
    REAL(jprb),    POINTER :: clear(:)              => NULL() !< Clear sky radiance
    REAL(jprb),    POINTER :: total(:)              => NULL() !< Cloudy radiance for given cloud
    REAL(jprb),    POINTER :: bt_clear(:)           => NULL() !< Brightness temp equivalent to clear radiance
    REAL(jprb),    POINTER :: bt(:)                 => NULL() !< Brightness temp equivalent to total radiance
    REAL(jprb),    POINTER :: refl_clear(:)         => NULL() !< Reflectance calculated from clear radiance
    REAL(jprb),    POINTER :: refl(:)               => NULL() !< Reflectance calculated from total radiance
    REAL(jprb),    POINTER :: overcast(:,:)         => NULL() !< Overcast radiance for opaque cloud at level bounding
                                                              !!   bottom of each layer
    REAL(jprb),    POINTER :: cloudy(:)             => NULL() !< 100% cloudy radiance for given cloud (simple cloud scheme)
                                                              !!   or same as total (addclouds/addaerosl true)
    INTEGER(jpim), POINTER :: quality(:)            => NULL() !< Flag to indicate whether they may be accuracy issues with
                                                              !!   simulated radiances
    LOGICAL(jplm)          :: plane_parallel                  !< Flag to indicate if strict plane-parallel was enforced
                                                              !!   e.g. by user-specified option or for DOM
    REAL(jprb),    POINTER :: geometric_height(:,:) => NULL() !< Geometric height of each level, direct model output
                                                              !!   only (nlevels,nchannels), units: m (NB the size is
                                                              !!   nchannels for consistency with other outputs, values
                                                              !!   are identical among channels for a given profile)
  END TYPE rttov_radiance

  !> Output radar reflectivities. Get the corresponding altitudes from the geometric_height in rttov_radiance
  !! Array sizes are (nchannels) or (nlevels,nchannels)
  TYPE rttov_reflectivity
    REAL(jprb),    POINTER  :: zef(:,:)    => NULL() !< Radar reflectivity
    REAL(jprb),    POINTER  :: azef(:,:)   => NULL() !< Attenuated radar reflectivity
  END TYPE rttov_reflectivity

  !> Secondary radiances optionally calculated in direct model only for clear-sky with no solar contribution
  !! Radiance unit: mW/m2/sr/cm-1. Array sizes are (nchannels) or (nlayers,nchannels)
  TYPE rttov_radiance2
    REAL(jprb), POINTER  :: upclear(:)     => NULL() !< Clear sky upwelling radiance without reflection term
    REAL(jprb), POINTER  :: dnclear(:)     => NULL() !< Clear sky downwelling radiance
    REAL(jprb), POINTER  :: refldnclear(:) => NULL() !< Reflected clear sky downwelling radiance
    REAL(jprb), POINTER  :: up(:,:)        => NULL() !< Sum( B * dT ) above cloud upwelling radiance from each layer
    REAL(jprb), POINTER  :: down(:,:)      => NULL() !< Sum( B / T**2 dT ) above cloud downwelling radiance from
                                                     !!   each layer
    REAL(jprb), POINTER  :: surf(:,:)      => NULL() !< Radiance at surface emitted from a black cloud
  END TYPE rttov_radiance2

  !> Output variables from RTTOV-SCATT enabling emissivity retrievals
  !! All arrays are size(chanprof)
  TYPE rttov_scatt_emis_retrieval_type
    REAL(jprb), POINTER  :: cfrac(:)    => NULL() !< RTTOV-SCATT effective cloud fraction
                                                  !!  (Tallsky = cfrac * Tcld + (1-cfrac) * Tclr
    REAL(jprb), POINTER  :: bsfc(:)     => NULL() !< Surface black-body radiance (i.e. Planck(Tsfc))
    REAL(jprb), POINTER  :: tau_cld(:)  => NULL() !< Along-path transmittance, surface to space (cloudy column)
    REAL(jprb), POINTER  :: up_cld(:)   => NULL() !< TOA upwelling radiance from atmosphere (not inc. surface
                                                  !!  emission or reflection)
    REAL(jprb), POINTER  :: down_cld(:) => NULL() !< SFC downwelling radiance (inc. cosmic term)
    REAL(jprb), POINTER  :: tau_clr(:)  => NULL() !< Along-path transmittance, surface to space (clear column)
    REAL(jprb), POINTER  :: up_clr(:)   => NULL() !< TOA upwelling radiance from atmosphere (not inc. surface
                                                  !!  emission or reflection)
    REAL(jprb), POINTER  :: down_clr(:) => NULL() !< SFC downwelling radiance (inc. cosmic term)
  END TYPE rttov_scatt_emis_retrieval_type

  !> Output PC scores and reconstructed radiances from PC-RTTOV or HTFRTC
  !! Radiance unit: mW/m2/sr/cm-1; BT unit: K.
  TYPE rttov_pccomp
    REAL(jprb), POINTER  :: clear_pcscores(:)      => NULL() !< Clear PC scores
    REAL(jprb), POINTER  :: total_pcscores(:)      => NULL() !< Cloudy PC scores for given cloud
    REAL(jprb), POINTER  :: overcast_pcscores(:,:) => NULL() !< Overcast PC scores for opaque cloud at level bounding
                                                             !<  bottom of each layerr
    REAL(jprb), POINTER  :: cloudy_pcscores(:)     => NULL() !< 100% cloudy PC scores for given cloud (simple cloud scheme)
                                                             !<  or same as total (addclouds true)
    REAL(jprb), POINTER  :: clear_pccomp(:)        => NULL() !< Clear radiances reconstructed using PCs
    REAL(jprb), POINTER  :: total_pccomp(:)        => NULL() !< Cloudy radiances reconstructed using PCs for given cloud
    REAL(jprb), POINTER  :: overcast_pccomp(:,:)   => NULL() !< Overcast radiances reconstructed using PCs for opaque cloud
                                                             !<  at level bounding bottom of each layer
    REAL(jprb), POINTER  :: cloudy_pccomp(:)       => NULL() !< 100% cloudy radiances reconstructed using PCs for given cloud
                                                             !< (simple cloud scheme) or as total (addclouds true)
    REAL(jprb), POINTER  :: bt_clear_pccomp(:)     => NULL() !< Brightness temp equivalent to clear radiances
                                                             !< reconstructed using PCs
    REAL(jprb), POINTER  :: bt_pccomp(:)           => NULL() !< Brightness temp equivalent to total radiances
                                                             !! reconstructed using PCs
  END TYPE rttov_pccomp

  !> Configuration options that apply to all flavours of RTTOV
  TYPE rttov_opts_config
    LOGICAL(jplm) :: apply_reg_limits = .FALSE.  !< Switch to restrict input profiles to coef training limits
    LOGICAL(jplm) :: verbose          = .TRUE.   !< Switch for verbose output
    LOGICAL(jplm) :: do_checkinput    = .TRUE.   !< Switch to apply internal profile checking

    ! Deprecated options: recommend default values
    LOGICAL(jplm) :: fix_hgpl         = .TRUE.   !< Switch to apply fix to match 2m p with elevation in geometry
                                                 !!   calculations (deprecated)
  END TYPE

  !> Options for PC-RTTOV
  TYPE rttov_opts_pc
    LOGICAL(jplm) :: addpc     = .FALSE.   !< Switch to enable PC-RTTOV
    INTEGER(jpim) :: ipcbnd    = 1         !< PC spectral band
    INTEGER(jpim) :: ipcreg    = -1        !< PC predictor channel set
    INTEGER(jpim) :: npcscores = -1        !< Number of PC scores to compute
    LOGICAL(jplm) :: addradrec = .FALSE.   !< Switch for calculation of reconstructed radiances
  END TYPE

  !> General radiative transfer options
  TYPE rttov_opts_rt_all
    LOGICAL(jplm) :: addrefrac              = .TRUE.   !< Switch to enable atmospheric refraction
    LOGICAL(jplm) :: switchrad              = .FALSE.  !< Switch for input units in AD/K models
    LOGICAL(jplm) :: use_t2m_opdep          = .TRUE.   !< Switch to enable use of 2m T variable in optical depth calculation
    LOGICAL(jplm) :: use_q2m                = .TRUE.   !< Switch to enable use of 2m q variable
    LOGICAL(jplm) :: do_lambertian          = .FALSE.  !< Switch for setting Lambertian reflection (IR and MW)
    LOGICAL(jplm) :: lambertian_fixed_angle = .TRUE.   !< Switch for fixed/parameterised effective angle for Lambertian option
    LOGICAL(jplm) :: plane_parallel         = .FALSE.  !< Switch to ignore atmospheric curvature
    LOGICAL(jplm) :: rad_down_lin_tau       = .TRUE.   !< Linear-in-tau or layer-mean for downwelling radiances

    LOGICAL(jplm) :: ozone_data             = .FALSE.  !< Switch to enable input of O3 profile
    LOGICAL(jplm) :: co2_data               = .FALSE.  !< Switch to enable input of CO2 profile
    LOGICAL(jplm) :: n2o_data               = .FALSE.  !< Switch to enable input of N2O profile
    LOGICAL(jplm) :: co_data                = .FALSE.  !< Switch to enable input of CO profile
    LOGICAL(jplm) :: ch4_data               = .FALSE.  !< Switch to enable input of CH4 profile
    LOGICAL(jplm) :: so2_data               = .FALSE.  !< Switch to enable input of SO2 profile

    ! Deprecated options: recommend default values
    LOGICAL(jplm) :: dtau_test              = .FALSE.  !< Switch to apply dtau test in transmit/integrate
                                                       !!   calculations (deprecated)
  END TYPE

  !> VIS/IR-only radiative transfer options
  TYPE rttov_opts_rt_ir
    TYPE(rttov_opts_pc) :: pc                            !< PC-RTTOV options

    INTEGER(jpim) :: solar_sea_brdf_model    = 2         !< Solar sea BRDF model (1-2)
    INTEGER(jpim) :: ir_sea_emis_model       = 2         !< IR sea emissivity model (1-2)
    LOGICAL(jplm) :: addsolar                = .FALSE.   !< Switch to enable solar simulations
    REAL(jprb)    :: rayleigh_max_wavelength = 2._jprb   !< Ignore Rayleigh scattering for channels at wavelengths
                                                         !!   above this (microns)
    REAL(jprb)    :: rayleigh_min_pressure   = 0._jprb   !< Ignore Rayleigh scattering at pressures below this (hPa)
    LOGICAL(jplm) :: rayleigh_single_scatt   = .TRUE.    !< Switch to enable Rayleigh single-scattering for VIS/NIR channels
    LOGICAL(jplm) :: do_nlte_correction      = .FALSE.   !< Switch to enable NLTE bias correction
    LOGICAL(jplm) :: addaerosl               = .FALSE.   !< Switch to enable IR aerosol calculations
    LOGICAL(jplm) :: user_aer_opt_param      = .FALSE.   !< Switch to supply aerosol optical properties explicitly per channel
    LOGICAL(jplm) :: addclouds               = .FALSE.   !< Switch to enable IR cloudy calculations
    LOGICAL(jplm) :: user_cld_opt_param      = .FALSE.   !< Switch to supply cloud optical properties explicitly per channel
    LOGICAL(jplm) :: grid_box_avg_cloud      = .TRUE.    !< Switch to supply grid-box average cloud concentration (true) or
                                                         !!   cloud concentration in cloudy fraction of each layer (false)
    REAL(jprb)    :: cldcol_threshold        = -1.0_jprb !< Ignore cloud columns with weights lower than this
    INTEGER(jpim) :: cloud_overlap           = 1         !< Cloud overlap scheme
    REAL(jprb)    :: cc_low_cloud_top        = 750._jprb !< Upper pressure limit for cloud_overlap_simple option (hPa)
    INTEGER(jpim) :: ir_scatt_model  = ir_scatt_chou     !< IR scattering model
    INTEGER(jpim) :: vis_scatt_model = vis_scatt_dom     !< VIS/NIR scattering model
    INTEGER(jpim) :: dom_nstreams            = 8         !< Number of DOM streams, must be even and not less than 2
    REAL(jprb)    :: dom_accuracy            = 0._jprb   !< Convergence criterion for termination of DOM azimuthal loop
    REAL(jprb)    :: dom_opdep_threshold     = 0._jprb   !< DOM ignores levels below this optical depth:
                                                         !!   10. is reasonable, not applied if <= 0
    LOGICAL(jplm) :: dom_rayleigh            = .FALSE.   !< DOM includes Rayleigh scattering for VIS/NIR channels
  END TYPE

  !> MW-only radiative transfer options
  TYPE rttov_opts_rt_mw
    INTEGER(jpim) :: fastem_version        = 6         !< FASTEM version (0-6); 0 => TESSEM2
    LOGICAL(jplm) :: supply_foam_fraction  = .FALSE.   !< Supply a foam fraction to FASTEM
    LOGICAL(jplm) :: clw_data              = .FALSE.   !< Switch to enable input of cloud liquid water profile
    INTEGER(jpim) :: clw_scheme            = 2         !< MW CLW scheme: 1 => Liebe, 2 => Rosenkranz, 3 => TKC
    REAL(jprb)    :: clw_cloud_top         = 322._jprb !< Lower pressure limit for MW CLW calculations (hPa)
  END TYPE

  !> Options for internal vertical interpolation and vertical grid setup
  TYPE rttov_opts_interp
    LOGICAL(jplm) :: addinterp        = .FALSE.        !< Switch to enable RTTOV interpolator
    INTEGER(jpim) :: interp_mode      = interp_rochon  !< Interpolation mode (1-5, see user guide)
    LOGICAL(jplm) :: lgradp           = .FALSE.        !< Switch to make pressure an active variable in TL/AD/K models

    ! Deprecated options: recommend default values
    LOGICAL(jplm) :: spacetop         = .TRUE.         !< Switch to assume space boundary at top-most input pressure
                                                       !!   level (deprecated)
    LOGICAL(jplm) :: reg_limit_extrap = .TRUE.         !< Switch to extrapolate input profiles using regression
                                                       !!   limits (deprecated)
  END TYPE

  !> HTFRTC options structure
  TYPE rttov_opts_htfrtc
    LOGICAL(jplm) :: htfrtc       = .FALSE. !< Switch to use htfrtc
    INTEGER(jpim) :: n_pc_in      = -1      !< Number of principal components to be used
    LOGICAL(jplm) :: reconstruct  = .FALSE. !< Switch to select reconstructed radiances
    LOGICAL(jplm) :: simple_cloud = .FALSE. !< Calculate simple cloud
    LOGICAL(jplm) :: overcast     = .FALSE. !< Calculate overcast cloud on all levels
  END TYPE rttov_opts_htfrtc

  !> Developer-only options structure
  TYPE rttov_opts_dev
    LOGICAL(jplm) :: do_opdep_calc      = .TRUE.  !< Enable/disable the RTTOV gas optical depth calculation
    LOGICAL(jplm) :: no_opt_param_tladk = .FALSE. !< Ignore cld/aer_opt_param_tl/ad/k even if passed as arguments
    ! MFASIS tuning parameters
    REAL(jprb)    :: od1_thresh = 0.5_jprb  ! Tuning parameters for determining which portion of ice cloud is
    REAL(jprb)    :: o_del1     = 0.15_jprb ! (potentially) part of mixed layer cloud
    REAL(jprb)    :: qq_1       = 1._jprb   ! tuning parameters for switching on/off mixed cloud correction
    REAL(jprb)    :: qq_2       = 0.3_jprb  !
    REAL(jprb)    :: x_del1     = 0.3_jprb  !
    REAL(jprb)    :: x_del2     = 0.2_jprb  !
  END TYPE rttov_opts_dev

  !> RTTOV options structure
  TYPE rttov_options
    TYPE(rttov_opts_config)  :: config          !< General configuration options
    TYPE(rttov_opts_rt_all)  :: rt_all          !< General RT options
    TYPE(rttov_opts_rt_ir)   :: rt_ir           !< VIS/IR RT options
    TYPE(rttov_opts_rt_mw)   :: rt_mw           !< MW RT options
    TYPE(rttov_opts_interp)  :: interpolation   !< Interpolation options
    TYPE(rttov_opts_htfrtc)  :: htfrtc_opts     !< HTFRTC options
    TYPE(rttov_opts_dev)     :: dev             !< Developer-only options
  END TYPE

  !> RTTOV-SCATT options structure: RTTOV-SCATT deliberately does
  !! not give user control over certain core RT options.
  TYPE rttov_options_scatt
    TYPE(rttov_opts_config) :: config                      !< General configuration options
    LOGICAL(jplm) :: lusercfrac            = .FALSE.       !< Switch to enable user-specified effective cloud fraction
    REAL(jprb)    :: cc_threshold          = 1.E-3_jprb    !< Minimum effective cloud fraction threshold to consider scattering
    REAL(jprb)    :: ice_polarisation      = 1.40_jprb     !< Polarised scattering factor for ice hydrometeors (<0 = no polarisation)
    LOGICAL(jplm) :: ozone_data            = .FALSE.       !< Switch to enable input of O3 profile
    LOGICAL(jplm) :: use_t2m_opdep         = .TRUE.        !< Switch to enable use of 2m T variable in optical depth calculation
    LOGICAL(jplm) :: use_q2m               = .TRUE.        !< Switch to enable use of 2m q variable
    LOGICAL(jplm) :: addrefrac             = .TRUE.        !< Switch to enable atmospheric refraction
    LOGICAL(jplm) :: rad_down_lin_tau      = .TRUE.        !< Linear-in-tau or layer-mean for downwelling radiances
    INTEGER(jpim) :: fastem_version        = 6             !< FASTEM version (1-6)
    LOGICAL(jplm) :: supply_foam_fraction  = .FALSE.       !< Supply a foam fraction to FASTEM
    INTEGER(jpim) :: interp_mode           = interp_rochon !< Interpolation mode (1-5, see user guide)
    LOGICAL(jplm) :: lgradp                = .FALSE.       !< Switch to make pressure an active variable in TL/AD/K models
    LOGICAL(jplm) :: hydro_cfrac_tlad      = .TRUE.        !< Switch for hydrometeor TL/AD sensitivity to effective cfrac
    LOGICAL(jplm) :: zero_hydro_tlad       = .FALSE.       !< Switch for hydrometeor TL/AD sensitivity in layers with zero
                                                           !!   hydrometeor concentration

    ! Deprecated options: recommend default values
    LOGICAL(jplm) :: reg_limit_extrap      = .TRUE.        !< Switch to extrapolate input profiles using regression
                                                           !!   limits (deprecated)
    LOGICAL(jplm) :: dtau_test             = .FALSE.       !< Switch to apply dtau test in transmit/integrate
                                                           !!   calculations (deprecated)
  END TYPE

  ! ---------------------------------------------------------------------------
  ! Internal RTTOV structures
  ! ---------------------------------------------------------------------------

  !> @internal Satellite and solar geometry
  TYPE rttov_geometry
    REAL(jprb) :: sinzen
    REAL(jprb) :: sinzen_sq
    REAL(jprb) :: coszen
    REAL(jprb) :: coszen_sq
    REAL(jprb) :: seczen
    REAL(jprb) :: seczen_sq
    REAL(jprb) :: sinview
    REAL(jprb) :: sinview_sq
    REAL(jprb) :: cosview_sq
    REAL(jprb) :: normzen
    REAL(jprb) :: coszen_sun
    REAL(jprb) :: sinzen_sun
    REAL(jprb) :: sinlat
    REAL(jprb) :: coslat
  END TYPE rttov_geometry

  !> @internal The actual predictor arrays
  TYPE rttov_path_pred
    REAL(jprb), POINTER :: mixedgas(:,:)    => NULL() ! (nmixed,  nlayers)
    REAL(jprb), POINTER :: watervapour(:,:) => NULL() ! (nwater,  nlayers)
    REAL(jprb), POINTER :: ozone(:,:)       => NULL() ! (nozone,  nlayers)
    REAL(jprb), POINTER :: wvcont(:,:)      => NULL() ! (nwvcont, nlayers)
    REAL(jprb), POINTER :: co2(:,:)         => NULL() ! (nco2,    nlayers)
    REAL(jprb), POINTER :: n2o(:,:)         => NULL() ! (nn2o,    nlayers)
    REAL(jprb), POINTER :: co(:,:)          => NULL() ! (nco,     nlayers)
    REAL(jprb), POINTER :: ch4(:,:)         => NULL() ! (nch4,    nlayers)
    REAL(jprb), POINTER :: so2(:,:)         => NULL() ! (nso2,    nlayers)
    REAL(jprb), POINTER :: pmc(:,:,:)       => NULL() ! pressure modulated cell (npmc,nlevels,nchannels)
  END TYPE rttov_path_pred

  !> @internal Predictors
  TYPE rttov_predictors
    ! the nxxxx could be set to 0 to indicate the abscence
    ! of the predictor, in that case there is no need to
    ! allocate the corresponding predictor
    INTEGER(jpim) :: nlevels   ! number of levels for predictors (all same)
    INTEGER(jpim) :: nmixed    ! number of variables for Mixed Gases
    INTEGER(jpim) :: nwater    ! number of variables for Water Vapour
    INTEGER(jpim) :: nozone    ! number of variables for Ozone
    INTEGER(jpim) :: nwvcont   ! number of variables for WV Continuum
    INTEGER(jpim) :: nco2      ! number of variables for CO2
    INTEGER(jpim) :: nn2o      ! number of variables for N2O
    INTEGER(jpim) :: nco       ! number of variables for CO
    INTEGER(jpim) :: nch4      ! number of variables for CH4
    INTEGER(jpim) :: nso2      ! number of variables for SO2
    INTEGER(jpim) :: npmc      ! number of variables for pressure modulated cell correction

    TYPE(rttov_path_pred), POINTER :: path1(:) => NULL() ! Predictors for surface-satellite path (always required)
    TYPE(rttov_path_pred), POINTER :: path2(:) => NULL() ! Predictors for sun-surface-satellite path (only required for solar)
  END TYPE rttov_predictors

  !> @internal Actual storage for gas fast coefficients 
  TYPE rttov_fast_coef_gas
    REAL(jprb), POINTER :: coef(:,:) => NULL()
  END TYPE rttov_fast_coef_gas

  !> @internal Fast coefficients
  !  Separate structure to allow for non-PW and PW if necessary
  TYPE rttov_fast_coef
    ! SIZE is fmv_gas, i.e. one element per gas.
    TYPE(rttov_fast_coef_gas), POINTER  :: gasarray(:) => NULL()

    ! SHAPE is (ncoef, levels): these point to coefs within gasarray for ease of use
    REAL(jprb), POINTER :: mixedgas(:,:)    => NULL()
    REAL(jprb), POINTER :: watervapour(:,:) => NULL()
    REAL(jprb), POINTER :: ozone(:,:)       => NULL()
    REAL(jprb), POINTER :: wvcont(:,:)      => NULL()
    REAL(jprb), POINTER :: co2(:,:)         => NULL()
    REAL(jprb), POINTER :: n2o(:,:)         => NULL()
    REAL(jprb), POINTER :: co(:,:)          => NULL()
    REAL(jprb), POINTER :: ch4(:,:)         => NULL()
    REAL(jprb), POINTER :: so2(:,:)         => NULL()
  END TYPE rttov_fast_coef

  !> @internal RTTOV-11.2 style fast coefficient structure to maintain backwards
  !! compatibility for later versions
  TYPE rttov_fast_coef_hdf_io
    ! SHAPE is (levels, channels, ncoef)
    REAL(jprb), POINTER :: mixedgas(:,:,:)    => NULL()
    REAL(jprb), POINTER :: watervapour(:,:,:) => NULL()
    REAL(jprb), POINTER :: ozone(:,:,:)       => NULL()
    REAL(jprb), POINTER :: wvcont(:,:,:)      => NULL()
    REAL(jprb), POINTER :: co2(:,:,:)         => NULL()
    REAL(jprb), POINTER :: n2o(:,:,:)         => NULL()
    REAL(jprb), POINTER :: co(:,:,:)          => NULL()
    REAL(jprb), POINTER :: ch4(:,:,:)         => NULL()
    REAL(jprb), POINTER :: so2(:,:,:)         => NULL()
  END TYPE rttov_fast_coef_hdf_io

  !> @internal NLTE coefficients
  TYPE rttov_nlte_coef
    REAL(jprb), POINTER :: coef(:,:,:,:)    => NULL() ! ncoef x nsat x nsol x nchan
    REAL(jprb), POINTER :: sol_zen_angle(:) => NULL()
    REAL(jprb), POINTER :: sat_zen_angle(:) => NULL()
    REAL(jprb), POINTER :: cos_sol(:)       => NULL()
    REAL(jprb), POINTER :: sec_sat(:)       => NULL()
    INTEGER(jpim)       :: ncoef, nsol, nsat, nchan
    INTEGER(jpim)       :: start_chan, end_chan
    REAL(jprb)          :: max_sat_angle
  END TYPE rttov_nlte_coef

  !> @internal Optical depth coefs
  TYPE rttov_coef
     ! Structure for the storage of RTTOV coefficients
     ! this may differ from what is stored in the coefficient files especially
     ! for the units (ie kg/kg to ppmv)
     ! Gases are separated in MxG WV O3
     ! Number of levels is the same for all gases (taken from MxG).
     !
    INTEGER(jpim)       :: id_platform         ! platform   (see user guide or rttov_const)
    INTEGER(jpim)       :: id_sat              ! satellite  (.....)
    INTEGER(jpim)       :: id_inst             ! instrument (.....)
    INTEGER(jpim)       :: id_sensor           ! 1 = Infrared, 2 = Microwave, 3 = High-resolution, 4 = Polarimeter
    INTEGER(jpim)       :: id_comp_lvl         ! RTTOV coefficient file version number
    INTEGER(jpim)       :: id_comp_pc          ! Principal component coefficient file version number
    INTEGER(jpim)       :: id_creation_date(3) ! YYYY MM DD
    CHARACTER (LEN=80)  :: id_creation         ! Creation comment
    CHARACTER (LEN=32)  :: id_common_name      ! usual name of the satellite
    CHARACTER (LEN=132) :: line_by_line(100)   ! readme LBL
    CHARACTER (LEN=132) :: readme_srf(100)     ! readme Spectral Response Function

       !FUNDAMENTAL_CONSTANTS section
    REAL(jprb)          :: fc_planck_c1        ! first radiation constant (mW/(m2*sr*cm-4))
    REAL(jprb)          :: fc_planck_c2        ! second radiation constant (cm*K)
    REAL(jprb)          :: fc_sat_height       ! satellite nominal altitude (km)

       !FAST_MODEL_VARIABLES section
    CHARACTER (LEN=32)     :: fmv_model_def  ! Predictors version (RTTOV7, RTTOV8, RTTOV9, RTTOV13)
    INTEGER(jpim)          :: fmv_model_ver  ! fast model version compatibility level
    INTEGER(jpim)          :: fmv_ori_nchn   ! number of channels in original file
    INTEGER(jpim)          :: fmv_chn        ! number of channels read from file into coef structure
    INTEGER(jpim)          :: fmv_gas        ! number of gases in file
    INTEGER(jpim), POINTER :: fmv_gas_id(:)  => NULL() ! gas id. number i gas_id list (fmv_gas)
    INTEGER(jpim), POINTER :: fmv_gas_pos(:) => NULL() ! respective position of each gas of gas_id list (ngases_max)
    INTEGER(jpim), POINTER :: fmv_var(:)     => NULL() ! number of variables/predictors by gas (fmv_gas)
    INTEGER(jpim), POINTER :: fmv_coe(:)     => NULL() ! number of coefficients by gas (fmv_gas)
    INTEGER(jpim), POINTER :: fmv_ncorr(:)   => NULL() ! number of coefs by gas for correction term (fmv_gas) (v13 only)
    INTEGER(jpim), POINTER :: fmv_lvl(:)     => NULL() ! number of levels(pres/absorber) by gas (fmv_gas) (v789 only)

    INTEGER(jpim)          :: nlevels             ! number of levels(pres/absorber) same for all gases
    INTEGER(jpim)          :: nlayers             ! number of layers(pres/absorber) nlevels-1
    LOGICAL(jplm)          :: inczeeman           ! Flag to include Zeeman effect for this sensor
    INTEGER(jpim)          :: nmixed              ! number of variables/predictors for Mixed Gases
    INTEGER(jpim)          :: nwater              ! number of variables/predictors for Water Vapour
    INTEGER(jpim)          :: nozone              ! number of variables/predictors for Ozone
    INTEGER(jpim)          :: nwvcont             ! number of variables/predictors for WV continuum
    INTEGER(jpim)          :: nco2                ! number of variables/predictors for CO2
    INTEGER(jpim)          :: nn2o                ! number of variables/predictors for N2O
    INTEGER(jpim)          :: nco                 ! number of variables/predictors for CO
    INTEGER(jpim)          :: nch4                ! number of variables/predictors for CH4
    INTEGER(jpim)          :: nso2                ! number of variables/predictors for SO2
    INTEGER(jpim)          :: ncmixed,  nccmixed  ! number of coefs for Mixed Gases, and correction term
    INTEGER(jpim)          :: ncwater,  nccwater  ! number of coefs for Water Vapour, and correction term
    INTEGER(jpim)          :: ncozone,  nccozone  ! number of coefs for Ozone, and correction term
    INTEGER(jpim)          :: ncwvcont, nccwvcont ! number of coefs for WV continuum, and correction term
    INTEGER(jpim)          :: ncco2,    nccco2    ! number of coefs for CO2, and correction term
    INTEGER(jpim)          :: ncn2o,    nccn2o    ! number of coefs for N2O, and correction term
    INTEGER(jpim)          :: ncco,     nccco     ! number of coefs for CO, and correction term
    INTEGER(jpim)          :: ncch4,    nccch4    ! number of coefs for CH4, and correction term
    INTEGER(jpim)          :: ncso2,    nccso2    ! number of coefs for SO2, and correction term

       !FILTER_FUNCTIONS section  array size is fmv_chn
    LOGICAL(jplm)          :: ff_val_bc                      ! are any band corrections to be applied?
    LOGICAL(jplm)          :: ff_val_gam                     ! any gamma corrections?
    INTEGER(jpim), POINTER :: ff_ori_chn(:)        => NULL() ! original chan number
    INTEGER(jpim), POINTER :: ff_val_chn(:)        => NULL() ! validity of the channel (1=OK)
    REAL(jprb),    POINTER :: ff_cwn(:)            => NULL() ! cental wave number (cm-1)
    REAL(jprb),    POINTER :: ff_bco(:)            => NULL() ! band correction offset (K)
    REAL(jprb),    POINTER :: ff_bcs(:)            => NULL() ! band correction slope (K/K)
    REAL(jprb),    POINTER :: ff_gam(:)            => NULL() ! gamma factor transmittance correction

       !TRANSMITTANCE_TRESHOLD section  array size is fmv_chn
    INTEGER(jpim), POINTER :: tt_val_chn(:)        => NULL()
    REAL(jprb),    POINTER :: tt_a0(:)             => NULL()
    REAL(jprb),    POINTER :: tt_a1(:)             => NULL()

       !PLANCK_WEIGHTED section array size if fmv_chn
    INTEGER(jpim), POINTER :: pw_val_chn(:)        => NULL() ! 0 => non-PW thermal coefs, 1 => PW thermal coefs

       !SOLAR_SPECTRUM section array size is fmv_chn
    INTEGER(jpim), POINTER :: ss_val_chn(:)        => NULL()
    REAL(jprb),    POINTER :: ss_solar_spectrum(:) => NULL() ! TOA solar irradiance (mW/m2/cm-1)
    REAL(jprb),    POINTER :: ss_rayleigh_ext(:)   => NULL() ! Channel-averaged Rayleigh extinction coefficients (km-1)
    REAL(jprb),    POINTER :: refl_visnir_ow(:)    => NULL() ! Vis/NIR ocean water refl, populated by init_coef
    REAL(jprb),    POINTER :: refl_visnir_fw(:)    => NULL() ! Vis/NIR fresh water refl, populated by init_coef

       !WATER_OPTICAL_CONSTANT section array size is fmv_chn
    COMPLEX(jprb), POINTER :: woc_waopc_ow(:)      => NULL() ! Refractive index of ocean water
    COMPLEX(jprb), POINTER :: woc_waopc_fw(:)      => NULL() ! Refractive index of fresh water

       !WAVE_SPECTRUM section array size is ws_nomega
       !Data used to compute the frequency spectrum of the JONSWAP
       !wave model surface wave.
    INTEGER(jpim)          :: ws_nomega
    REAL(jprb),    POINTER :: ws_npoint(:)         => NULL()
    REAL(jprb),    POINTER :: ws_k_omega(:)        => NULL()

       !FASTEM section
    INTEGER(jpim), POINTER :: fastem_polar(:)      => NULL() ! polarisation of each channel
    REAL(jprb),    POINTER :: pol_phi(:)           => NULL() ! angle for V/H mix for pol_id 7 (degrees)
    REAL(jprb),    POINTER :: pol_fac_v(:)         => NULL() ! factor for V-pol contribution
    REAL(jprb),    POINTER :: pol_fac_h(:)         => NULL() ! factor for H-pol contribution
       ! 0 = 0.5 V+H
       ! 1 = 90 - incident angle
       ! 2 = incident angle
       ! 3 = vertical
       ! 4 = horizontal
       ! 5 = 3rd element of Stokes vector
       ! 6 = 4th element of Stokes vector
       ! 7 = fixed mixture of V+H according to angle pol_phi:
       !       Tb = Tb_H * cos(pol_phi)^2 + Tb_V * sin(pol_phi)^2

       !SSIREM section     array size is fmv_chn
       ! ems =   ssirem_a0
       !       - ssirem_a1*(zen**ssirem_xzn1)
       !       - ssirem_a2*(zen**ssirem_xzn2)
       ! where zen is satellite zenith angle in degrees, divided by 60.
    INTEGER(jpim)          :: ssirem_ver                ! version number
    REAL(jprb),    POINTER :: ssirem_a0(:)    => NULL() ! constant coef
    REAL(jprb),    POINTER :: ssirem_a1(:)    => NULL() ! first coef
    REAL(jprb),    POINTER :: ssirem_a2(:)    => NULL() ! second coef
    REAL(jprb),    POINTER :: ssirem_xzn1(:)  => NULL() ! 1st exponent on zenith angle
    REAL(jprb),    POINTER :: ssirem_xzn2(:)  => NULL() ! 2nd exponent on zenith angle

       ! IREMIS - IR sea surface emissivity model section
    INTEGER(jpim)         :: iremis_version             ! version number
    REAL(jprb)            :: iremis_angle0              ! reference zenith angle (degrees)
    REAL(jprb)            :: iremis_tskin0              ! reference Tskin (K)
    INTEGER(jpim)         :: iremis_ncoef               ! version number
    REAL(jprb),   POINTER :: iremis_coef(:,:) => NULL() ! coefficients (iremis_ncoef,fmv_chn)

       !REFERENCE_PROFILE section  defined on mixed gases pressure levels
       ! gases are in the order of gas id codes
       ! unit for mr in coef file is ppmv over dry air
       ! unit for mr for optical depth calculations is ppmv over dry air
    REAL(jprb), POINTER   :: ref_prfl_p(:)    => NULL() ! pressure  (hPa)       (levels)
    REAL(jprb), POINTER   :: ref_prfl_t(:,:)  => NULL() ! temperature (K)       (levels, gases)
    REAL(jprb), POINTER   :: ref_prfl_mr(:,:) => NULL() ! mixing ratio (ppmv)   (levels, gases)
       !Background profiles to use when no optional input gas profile supplied
    REAL(jprb), POINTER   :: bkg_prfl_mr(:,:) => NULL() ! mixing ratio (ppmv)   (levels, gases)

       !PROFILE_ENVELOPE section
       ! gases are in the order of gas id codes
       ! unit for mr in coef file is ppmv over dry air
       ! unit for mr for optical depth calculations is ppmv over dry air
    REAL(jprb), POINTER   :: lim_prfl_p(:)      => NULL() ! pressure  (hPa)     (levels)
    REAL(jprb), POINTER   :: env_prfl_tmax(:)   => NULL() ! max temperature (K) (levels)
    REAL(jprb), POINTER   :: env_prfl_tmin(:)   => NULL() ! min temperature (K) (levels)
    REAL(jprb), POINTER   :: env_prfl_gmax(:,:) => NULL() ! max mixing r (ppmv) (levels, gases)
    REAL(jprb), POINTER   :: env_prfl_gmin(:,:) => NULL() ! min mixing r (ppmv) (levels, gases)
       ! Profile limits calculated from envelope
    REAL(jprb), POINTER   :: lim_prfl_tmax(:)   => NULL() ! max temperature (K) (levels)
    REAL(jprb), POINTER   :: lim_prfl_tmin(:)   => NULL() ! min temperature (K) (levels)
    REAL(jprb), POINTER   :: lim_prfl_gmax(:,:) => NULL() ! max mixing r (ppmv) (levels, gases)
    REAL(jprb), POINTER   :: lim_prfl_gmin(:,:) => NULL() ! min mixing r (ppmv) (levels, gases)

       !FAST_COEFFICIENTS/SOLAR_FAST_COEFFICIENTS section
       ! For non-PW instruments, "solar" will point to "thermal" coefs
       ! For instruments with solar-affected PW channels, both thermal and solar
       ! structures will be populated from coef file
       ! Fast coefficients now grouped per channel rather than per gas (RTTOV12)
       ! so thermal and solar are now arrays with size nchannels.
    TYPE(rttov_fast_coef), POINTER :: thermal(:)      => NULL() ! FAST_COEFFICIENTS for gas opdep prediction
    TYPE(rttov_fast_coef), POINTER :: thermal_corr(:) => NULL() ! Coefficients for correction term
    TYPE(rttov_fast_coef), POINTER :: solar(:)        => NULL() ! SOLAR_FAST_COEFFICIENTS when present
    TYPE(rttov_fast_coef), POINTER :: solar_corr(:)   => NULL() ! Coefficients for correction term
    LOGICAL(jplm)                  :: solarcoef                 ! .TRUE. if solar fast coefs present in file: need
                                                                ! to know whether solar points to thermal or is
                                                                ! allocated separately.
    LOGICAL(jplm)                  :: nltecoef = .FALSE.   ! .TRUE. if nlte coefs present
    TYPE(rttov_nlte_coef), POINTER :: nlte_coef  => NULL() ! nlte_coef

         !PRESSURE_MODULATED_CELL section
    LOGICAL(jplm)       :: pmc_shift = .FALSE.       ! .TRUE. if pmc shift coefs present
    REAL(jprb)          :: pmc_lengthcell            ! cell length (cm)
    REAL(jprb), POINTER :: pmc_pnominal(:) => NULL() ! nominal cell pressure (hPa) - nchannels
    REAL(jprb)          :: pmc_tempcell              ! cell temperature (K)
    REAL(jprb)          :: pmc_betaplus1             ! co2 band-average: self-HW/air-HW
    INTEGER(jpim)       :: pmc_nlay                  ! number of layers used
    INTEGER(jpim)       :: pmc_nvar                  ! number of variables used
    REAL(jprb), POINTER :: pmc_coef(:,:,:) => NULL() ! pressure moodulated cell corrections - nlevels, nchannels, nvariables
    REAL(jprb), POINTER :: pmc_ppmc(:)     => NULL() ! actual cell pressure (hPa) - nchannels

       ! Auxiliary variables
    REAL(jprb)               :: ratoe                      ! ratio (H+R)/R  H=sat height, R=Earth radius
    REAL(jprb), POINTER      :: planck1(:)       => NULL() ! C1 * Nu**3 (mW/(m2*sr*cm-1))
    REAL(jprb), POINTER      :: planck2(:)       => NULL() ! C2 * Nu (K)
    REAL(jprb), POINTER      :: frequency_ghz(:) => NULL() ! frequency in GHz

       ! other predictor variables see Science and Validation report
    REAL(jprb), POINTER      :: dp(:)      => NULL()                         ! interval between standard p levels (hPa)
    REAL(jprb), POINTER      :: dpp(:)     => NULL()                         ! pressure based variable (hPa**2)
    REAL(jprb), POINTER      :: tstar(:)   => NULL(), tstar_r(:)   => NULL() ! layer temp (K)
    REAL(jprb), POINTER      :: to3star(:) => NULL(), to3star_r(:) => NULL() ! layer temp for O3 calculations (K)
    REAL(jprb), POINTER      :: wstar(:)   => NULL(), wstar_r(:)   => NULL() ! layer WV  (ppmv)
    REAL(jprb), POINTER      :: ostar(:)   => NULL(), ostar_r(:)   => NULL() ! layer O3  (ppmv)
    REAL(jprb), POINTER      :: co2star(:) => NULL(), co2star_r(:) => NULL() ! layer co2 (ppmv)
    REAL(jprb), POINTER      :: n2ostar(:) => NULL(), n2ostar_r(:) => NULL() ! layer n2o (ppmv)
    REAL(jprb), POINTER      :: costar(:)  => NULL(), costar_r(:)  => NULL() ! layer co  (ppmv)
    REAL(jprb), POINTER      :: ch4star(:) => NULL(), ch4star_r(:) => NULL() ! layer ch4 (ppmv)
    REAL(jprb), POINTER      :: so2star(:) => NULL(), so2star_r(:) => NULL() ! layer so2 (ppmv)
    REAL(jprb), POINTER      :: tstar_wsum_r(:)   => NULL(), tstarmod_wsum_r(:) => NULL()
    REAL(jprb), POINTER      :: tstar_uwsum_r(:)  => NULL()
    REAL(jprb), POINTER      :: wstar_wsum_r(:)   => NULL(), wtstar_wsum_r(:)   => NULL()
    REAL(jprb), POINTER      :: ostar_wsum_r(:)   => NULL()
    REAL(jprb), POINTER      :: co2star_wsum_r(:) => NULL()
    REAL(jprb), POINTER      :: n2ostar_wsum_r(:) => NULL(), n2otstar_wsum_r(:) => NULL()
    REAL(jprb), POINTER      :: costar_wsum_r(:)  => NULL(), cotstar_wsum_r(:)  => NULL()
    REAL(jprb), POINTER      :: ch4star_wsum_r(:) => NULL(), ch4tstar_wsum_r(:) => NULL()
    REAL(jprb), POINTER      :: so2star_wsum_r(:) => NULL(), so2tstar_wsum_r(:) => NULL()

    INTEGER(jpim), POINTER   :: bounds(:,:,:,:) => NULL() 
  END TYPE rttov_coef

  !> @internal RTTOV-SCATT coefs
  TYPE rttov_scatt_coef
    ! Structure for the storage of RTTOV_SCATT coefficients (hydrotable files)
    INTEGER(jpim) :: rttov_version ! RTTOV compatibility version
    INTEGER(jpim) :: nhydro        ! Number of hydrometeors in computation
    INTEGER(jpim) :: mtype         ! Number of hydrometeors
    INTEGER(jpim) :: mfreqm        ! Number of frequencies
    INTEGER(jpim) :: mtemp         ! Number of temperature bins
    INTEGER(jpim) :: mwc           ! Number of water bins
    REAL(jprb), POINTER :: offset_temp(:) => NULL() ! temperature offset in table for each hydrometeor type (K)
    REAL(jprb)          :: offset_water             ! liquid/ice water offset in table
    REAL(jprb)          :: scale_water              ! log10(liquid/ice water) scaling factor in table
    REAL(jprb)          :: from_scale_water         ! 10**(1/scale_water)
    LOGICAL(jplm), POINTER :: is_frozen(:) => NULL() ! whether the hydrometeor type is frozen (T) or liquid (F)
    REAL(jprb), POINTER :: freq(:)      => NULL()   ! list of frequencies in hydrotable
    INTEGER(jpim), POINTER :: mpol(:)   => NULL()   ! list of polarisations in Hydro table (-1 = polarisation ignored)
    REAL(jprb), POINTER :: ext(:,:,:,:) => NULL()   ! extinction coefficent table
    REAL(jprb), POINTER :: ssa(:,:,:,:) => NULL()   ! single scattering albedo table
    REAL(jprb), POINTER :: asp(:,:,:,:) => NULL()   ! asymmetry parameter table
    REAL(jprb), POINTER :: zef(:,:,:,:) => NULL()   ! reflectivity (Z/Z0) parameter table (only present for radars)
  END TYPE rttov_scatt_coef

  !> @internal Surface and cloud fraction
  TYPE rttov_profile_aux_s
    INTEGER(jpim) :: nearestlev_surf ! nearest model level above surface
    REAL(jprb)    :: pfraction_surf  ! pressure fraction of surface in model layer (hPa)
    INTEGER(jpim) :: nearestlev_ctp  ! nearest model level above cloud top
    REAL(jprb)    :: pfraction_ctp   ! pressure fraction of cloud top pressure in layer (hPa)
    REAL(jprb)    :: cfraction       ! cloud fraction (0-1) 1 for 100% cloud cover
  END TYPE rttov_profile_aux_s

  !> @internal Auxiliary profile variables (user levels)
  TYPE rttov_profile_aux
    TYPE(rttov_profile_aux_s), POINTER :: s(:) => NULL() 

    ! MW CLW absorption
    REAL(jprb),    POINTER :: clw(:,:) => NULL() ! CLW * layer thickness (g.m^-3.km)

    ! Liebe MW CLW permittivity
    REAL(jprb),    POINTER :: debye_prof(:,:,:) => NULL() ! Debye terms

    ! Rosenkranz MW CLW permittivity
    REAL(jprb),    POINTER :: ros_eps_s(:,:)  => NULL() ! Static dielectric constant (epsilon(nu=0))
    REAL(jprb),    POINTER :: ros_dr(:,:)     => NULL() ! Rosenkranz parameter Delta_R
    REAL(jprb),    POINTER :: ros_gammar(:,:) => NULL() ! Rosenkranz parameter gamma_R
    REAL(jprb),    POINTER :: ros_db(:,:)     => NULL() ! Rosenkranz parameter Delta_B
    REAL(jprb),    POINTER :: ros_nu1(:,:)    => NULL() ! Rosenkranz parameter nu1
    COMPLEX(jprb), POINTER :: ros_z1(:,:)     => NULL() ! Rosenkranz parameter z1
    COMPLEX(jprb)          :: ros_z2                    ! Rosenkranz parameter z2 (this is a constant)
    REAL(jprb),    POINTER :: ros_log_abs_z1_sq(:,:) => NULL() ! LOG(ABS(z1)**2)
    COMPLEX(jprb), POINTER :: ros_log_z1star_z2(:,:) => NULL() ! LOG(CONJG(z1)*z2)
    COMPLEX(jprb), POINTER :: ros_div1(:,:)   => NULL() ! Precalculated complex division
    COMPLEX(jprb), POINTER :: ros_div2(:,:)   => NULL() ! Precalculated complex division

    ! TKC MW CLW permittivity
    REAL(jprb),    POINTER :: tkc_eps_s(:,:)   => NULL() ! Static dielectric constant (epsilon(nu=0))
    REAL(jprb),    POINTER :: tkc_delta_1(:,:) => NULL() ! TKC parameter delta_1
    REAL(jprb),    POINTER :: tkc_delta_2(:,:) => NULL() ! TKC parameter delta_2
    REAL(jprb),    POINTER :: tkc_tau_1(:,:)   => NULL() ! TKC parameter tau_1
    REAL(jprb),    POINTER :: tkc_tau_2(:,:)   => NULL() ! TKC parameter tau_2

    ! Variables used in water cloud effective diameter parameterisations (visible/IR scattering)
    REAL(jprb),    POINTER :: clw_dg(:,:)      => NULL() ! Generalized effective diameter (microns)
    REAL(jprb),    POINTER :: clw_dg_ref(:,:)  => NULL() ! Generalized effective diameter before clipping (microns)

    ! Variables used in ice effective diameter parameterisations (visible/IR scattering)
    REAL(jprb),    POINTER :: ice_dg(:,:)      => NULL() ! Generalized effective diameter (microns)
    REAL(jprb),    POINTER :: ice_dg_ref(:,:)  => NULL() ! Generalized effective diameter before clipping (microns)
    REAL(jprb),    POINTER :: fac1_ice_dg(:,:) => NULL() ! Intermediate variables used to compute the generalized diameter

    ! Variables used in aerosol scattering including rel. hum. calculation
    REAL(jprb),    POINTER :: relhum(:,:)   => NULL() ! Relative humidity (%)
    REAL(jprb),    POINTER :: tave(:,:)     => NULL() ! Layer average temperature (K)
    REAL(jprb),    POINTER :: wmixave(:,:)  => NULL() ! Layer average WV (ppmv)
    REAL(jprb),    POINTER :: xpresave(:,:) => NULL() ! Layer average pressure (hPa)
    REAL(jprb),    POINTER :: ppv(:,:)      => NULL() ! Partial pressure of WV (hPa)
    REAL(jprb),    POINTER :: esw(:,:)      => NULL() ! Saturated vapour pressure liquid (hPa)
    REAL(jprb),    POINTER :: esi(:,:)      => NULL() ! Saturated vapour pressure ice (hPa)
  END TYPE rttov_profile_aux

  !> @internal Auxiliary profile variables on coef levels
  TYPE rttov_profile_aux_coef
    TYPE(rttov_profile_aux_s), POINTER :: s(:) => NULL() 

    ! Variables used in predictor calculations
    REAL(jprb), POINTER :: t_layer(:,:)   => NULL()   ! avg layer temperature (K)
    REAL(jprb), POINTER :: w_layer(:,:)   => NULL()   ! avg layer humidity (ppmv)
    REAL(jprb), POINTER :: o3_layer(:,:)  => NULL()   ! avg layer ozone (ppmv)
    REAL(jprb), POINTER :: co2_layer(:,:) => NULL()
    REAL(jprb), POINTER :: n2o_layer(:,:) => NULL()
    REAL(jprb), POINTER :: co_layer(:,:)  => NULL()
    REAL(jprb), POINTER :: ch4_layer(:,:) => NULL()
    REAL(jprb), POINTER :: so2_layer(:,:) => NULL()
    REAL(jprb), POINTER :: dt(:,:)        => NULL(), dtabsdt(:,:)    => NULL() ! deviation from ref temperature prof (K)
    REAL(jprb), POINTER :: tr(:,:)        => NULL(), tr2(:,:)        => NULL() ! ratio t / ref_t
    REAL(jprb), POINTER :: tr_r(:,:)      => NULL(), tr_sqrt(:,:)    => NULL()
    REAL(jprb), POINTER :: tw(:,:)        => NULL(), tw_sqrt(:,:)    => NULL()
    REAL(jprb), POINTER :: tw_4rt(:,:)    => NULL(), tuw(:,:)        => NULL(), twr(:,:)       => NULL()
    REAL(jprb), POINTER :: wr(:,:)        => NULL(), wr_sqrt(:,:)    => NULL()
    REAL(jprb), POINTER :: wr_4rt(:,:)    => NULL(), wr_rsqrt(:,:)   => NULL()
    REAL(jprb), POINTER :: ww(:,:)        => NULL(), ww_sqrt(:,:)    => NULL()
    REAL(jprb), POINTER :: ww_4rt(:,:)    => NULL(), ww_rsqrt(:,:)   => NULL(), ww_r(:,:)      => NULL()
    REAL(jprb), POINTER :: wwr(:,:)       => NULL(), wwr_r(:,:)      => NULL()
    REAL(jprb), POINTER :: wrw_r(:,:)     => NULL(), wrwr_r(:,:)     => NULL()
    REAL(jprb), POINTER :: dto(:,:)       => NULL(), tro(:,:)        => NULL() ! deviation / ratio from ref ozone temp prof (K)
    REAL(jprb), POINTER :: or(:,:)        => NULL(), or_sqrt(:,:)    => NULL()
    REAL(jprb), POINTER :: ow(:,:)        => NULL(), ow_r(:,:)       => NULL()
    REAL(jprb), POINTER :: ow_sqrt(:,:)   => NULL(), ow_rsqrt(:,:)   => NULL(), ow_4rt(:,:)    => NULL()
    REAL(jprb), POINTER :: co2r(:,:)      => NULL(), co2r_sqrt(:,:)  => NULL(), co2w(:,:)      => NULL()
    REAL(jprb), POINTER :: n2or(:,:)      => NULL(), n2or_sqrt(:,:)  => NULL()
    REAL(jprb), POINTER :: n2or_4rt(:,:)  => NULL(), n2ow(:,:)       => NULL()
    REAL(jprb), POINTER :: n2owr(:,:)     => NULL(), n2ow_r(:,:)     => NULL()
    REAL(jprb), POINTER :: cor(:,:)       => NULL(), cor_sqrt(:,:)   => NULL(), cor_4rt(:,:)   => NULL()
    REAL(jprb), POINTER :: corw_r(:,:)    => NULL(), corw_rsqrt(:,:) => NULL(), corw_r4rt(:,:) => NULL()
    REAL(jprb), POINTER :: cow(:,:)       => NULL(), cow_rsqrt(:,:)  => NULL()
    REAL(jprb), POINTER :: cow_r4rt(:,:)  => NULL(), cowr(:,:)       => NULL()
    REAL(jprb), POINTER :: cowr_r(:,:)    => NULL(), cowr_4rt(:,:)   => NULL()
    REAL(jprb), POINTER :: ch4r(:,:)      => NULL(), ch4r_sqrt(:,:)  => NULL(), ch4r_4rt(:,:)  => NULL()
    REAL(jprb), POINTER :: ch4w(:,:)      => NULL(), ch4w_4rt(:,:)   => NULL()
    REAL(jprb), POINTER :: ch4w_r(:,:)    => NULL(), ch4wr(:,:)      => NULL(), ch4rw_r(:,:)   => NULL()
    REAL(jprb), POINTER :: so2r(:,:)      => NULL(), so2r_sqrt(:,:)  => NULL(), so2r_4rt(:,:)  => NULL()
    REAL(jprb), POINTER :: so2w(:,:)      => NULL(), so2w_sqrt(:,:)  => NULL()
    REAL(jprb), POINTER :: so2w_4rt(:,:)  => NULL(), so2w_r(:,:)     => NULL()
    REAL(jprb), POINTER :: so2wr(:,:)     => NULL(), so2wr_r(:,:)    => NULL()
    REAL(jprb), POINTER :: so2rw_r(:,:)   => NULL(), so2rwr_r(:,:)   => NULL()

    REAL(jprb), POINTER :: co2_cm(:) => NULL() ! (profiles) CO2 thickness at STP (for PMC shift)

    ! These are derived from pathsat and patheff in the rttov_raytracing structure
    ! but are only used for predictor calculations so are stored here
    REAL(jprb), POINTER :: pathsat_rsqrt(:,:) => NULL() ! (layers,profiles) Reciprocal of SQRT of pathsat
    REAL(jprb), POINTER :: pathsat_sqrt(:,:)  => NULL() ! (layers,profiles) SQRT of pathsat
    REAL(jprb), POINTER :: pathsat_4rt(:,:)   => NULL() ! (layers,profiles) 4thRT of pathsat
    REAL(jprb), POINTER :: patheff_rsqrt(:,:) => NULL() ! (layers,profiles) Reciprocal of SQRT of patheff
    REAL(jprb), POINTER :: patheff_sqrt(:,:)  => NULL() ! (layers,profiles) SQRT of patheff
    REAL(jprb), POINTER :: patheff_4rt(:,:)   => NULL() ! (layers,profiles) 4thRT of patheff
  END TYPE rttov_profile_aux_coef

  !> @internal Auxiliary profile variables for RTTOV_SCATT
  TYPE rttov_profile_scatt_aux
    REAL(jprb), POINTER :: cfrac(:)      => NULL() ! horizontal cloud fraction, one value used for all layers (0-1)
    REAL(jprb), POINTER :: ems_bnd(:)    => NULL() ! surface emissivity for boundary conditions
    REAL(jprb), POINTER :: ref_bnd(:)    => NULL() ! surface emissivity for boundary conditions
    REAL(jprb), POINTER :: ems_cld(:)    => NULL() ! surface emissivity taking into account cloud/rain impact on od
    REAL(jprb), POINTER :: ref_cld(:)    => NULL() ! surface reflectivity taking into account cloud/rain impact on od
    REAL(jprb), POINTER :: dz(:,:)       => NULL() ! layer depth   [km]
    REAL(jprb), POINTER :: tbd(:,:)      => NULL() ! (effective) temperature at layer boundary [K]
    REAL(jprb), POINTER :: tsfc(:)       => NULL() ! (effective) temperature at surface [K]
    REAL(jprb), POINTER :: tcosmic(:)    => NULL() ! (effective) temperature of cosmic background [K]
    REAL(jprb), POINTER :: hydro(:,:,:)  => NULL() ! hydrometeor water content (g/m3)
    INTEGER(jpim), POINTER :: mclayer(:) => NULL() ! upper level cloud layer
    REAL(jprb), POINTER :: delta(:,:)    => NULL() ! (= ext*dz/coszen)
    REAL(jprb), POINTER :: tau(:,:)      => NULL() ! optical depths (= exp(-delta))
    REAL(jprb), POINTER :: ext(:,:)      => NULL() ! extinction coefficient integreated over hydrometeor types
    REAL(jprb), POINTER :: ssa(:,:)      => NULL() ! single scattering albedo integreated over hydrometeor types
    REAL(jprb), POINTER :: asm(:,:)      => NULL() ! asymetry parameter integreated over hydrometeor types [-1,1]
    REAL(jprb), POINTER :: int_tau(:,:)  => NULL() ! integrated optical depths 
    REAL(jprb), POINTER :: zef(:,:)      => NULL() ! radar reflectivity 
    REAL(jprb), POINTER :: lambda(:,:)   => NULL() ! eddington approx. variable
                                                   ! (= sqrt( 3*ext*ext*(1-ssa)*(1-ssa*asm) )
    REAL(jprb), POINTER :: h (:,:)       => NULL() ! boundary condition variable (= 1.5_jprb*ext(1-ssa*asm))
    REAL(jprb), POINTER :: b0(:,:)       => NULL() ! lower level planck function
    REAL(jprb), POINTER :: b1(:,:)       => NULL() ! planck function gradient
    REAL(jprb), POINTER :: bn(:,:)       => NULL() ! upper level planck function
    REAL(jprb), POINTER :: btop(:)       => NULL() ! planck function entering atmosphere (cosmic radiation)
    REAL(jprb), POINTER :: bsfc(:)       => NULL() ! planck function at surface
  END TYPE rttov_profile_scatt_aux

  !> @internal Auxillary variables for RTTOV_SCATT (emissivity retrieval - internal)
  TYPE eddington_sfc_type
    REAL(jprb), POINTER  :: tau(:)  => NULL()
    REAL(jprb), POINTER  :: up(:)   => NULL() ! TOA upwelling radiance
    REAL(jprb), POINTER  :: down(:) => NULL() ! SFC downwelling radiance
  END TYPE eddington_sfc_type

  !> @internal Path optical depths as predicted or interpolated (unitless)
  TYPE rttov_opdp_path
    REAL(jprb), POINTER :: atm_level(:,:)       => NULL() ! neg optical depth for thermal radiation (levels to space),
                                                          ! size (levels, channels)
    REAL(jprb), POINTER :: sun_level_path2(:,:) => NULL() ! neg optical depth for solar radiation (levels to space) for
                                                          ! combined sun-surface-satellite path, size (levels, channels)
  END TYPE rttov_opdp_path

  !> @internal Transmissions and optical depths (unitless)
  TYPE rttov_path_transmission
    REAL(jprb), POINTER  :: tau_surf_p(:,:)         => NULL() ! Lambertian transmittance from surface (columns,channels)
    REAL(jprb), POINTER  :: tau_surf_p_r(:,:)       => NULL() ! reciprocal Lambertian trans. from surface (columns,channels)
    REAL(jprb), POINTER  :: tau_surf(:,:)           => NULL() ! transmittance from surface (columns,channels)
    REAL(jprb), POINTER  :: tau_surf_r(:,:)         => NULL() ! reciprocal transmittance from surface (columns,channels)
    REAL(jprb), POINTER  :: tau_level(:,:,:)        => NULL() ! transmittance from each standard pressure level
                                                              ! (levels,columns,channels)
    REAL(jprb), POINTER  :: tau_level_r(:,:,:)      => NULL() ! reciprocal transmittance from each standard pressure level
                                                              ! (levels,columns,channels)
    REAL(jprb), POINTER  :: tau_level_p(:,:,:)      => NULL() ! Lambertian transmittance from each standard pressure level
                                                              ! (levels,columns,channels)
    REAL(jprb), POINTER  :: tau_level_p_r(:,:,:)    => NULL() ! reciprocal Lambertian trans. from each standard pressure level
                                                              ! (levels,columns,channels)
    REAL(jprb), POINTER  :: od_singlelayer(:,:,:)   => NULL() ! single-layer optical depth
    REAL(jprb), POINTER  :: od_singlelayer_r(:,:,:) => NULL() ! reciprocal single-layer optical depth
    REAL(jprb), POINTER  :: od_sfrac(:,:)           => NULL() ! optical depth of partial layer above surface
    REAL(jprb), POINTER  :: od_sfrac_r(:,:)         => NULL() ! reciprocal optical depth of partial layer above surface
    REAL(jprb), POINTER  :: od_frac_ac(:,:)         => NULL() ! cloud/aerosol optical depth of partial layer above surface
    REAL(jprb), POINTER  :: tau_surf_ac(:,:)        => NULL() ! cloud/aerosol transmittance from surface

    REAL(jprb), POINTER  :: fac2(:,:,:)             => NULL() ! Mask for integration calculation: thermal and solar path1
  END TYPE

  !> @internal Auxiliary transmittance variables
  TYPE rttov_transmission_aux
    REAL(jprb), POINTER  :: fac1(:,:,:)   => NULL() ! Mask for integration calculation
    REAL(jprb), POINTER  :: surf_fac(:,:) => NULL() ! Mask for near-surface layer integration calculation

    TYPE(rttov_path_transmission), POINTER :: thermal_path1 => NULL() ! Thermal transmittances on surface-sensor path
    TYPE(rttov_path_transmission), POINTER :: solar_path2   => NULL() ! Solar transmittances on combined sun-surface-sensor path
    TYPE(rttov_path_transmission), POINTER :: solar_path1   => NULL() ! Solar transmittances on surface-sensor path
  END TYPE rttov_transmission_aux

  !> @internal Auxiliary radiance variables
  TYPE rttov_radiance_aux
    ! Auxiliary calculation arrays for RTE integration
    ! Direct model arrays need to be passed to TL AD and K codes
    ! array size is of (nchannels) or (nlevels, nchannels)
    REAL(jprb), POINTER :: air(:,:)        => NULL() ! Planck emission from atmospheric layers (mW/m2/sr/cm-1)
    REAL(jprb), POINTER :: surfair(:)      => NULL() ! Planck emission from near-surface layer (mW/m2/sr/cm-1)
    REAL(jprb), POINTER :: skin(:)         => NULL() ! Planck emission from surface skin (mW/m2/sr/cm-1)
    REAL(jprb), POINTER :: cosmic(:)       => NULL() ! Planck emission from CMBR (mW/m2/sr/cm-1)
    REAL(jprb), POINTER :: air_t_eff(:,:)  => NULL() ! Effective (band-corrected) temperature of atmospheric layers (K)
    REAL(jprb), POINTER :: surf_t_eff(:)   => NULL() ! Effective (band-corrected) temperature of near-surface layer (K)
    REAL(jprb), POINTER :: skin_t_eff(:)   => NULL() ! Effective (band-corrected) temperature of surface skin (K)
    REAL(jprb), POINTER :: cosmic_t_eff(:) => NULL() ! Effective (band-corrected) temperature of CMBR (K)

    ! Radiances (units: mW/m2/sr/cm-1)
    REAL(jprb), POINTER :: up(:,:,:)               => NULL() ! sum( B * dT )
    REAL(jprb), POINTER :: down(:,:,:)             => NULL() ! sum ( B / T**2 dT )
    REAL(jprb), POINTER :: down_p(:,:,:)           => NULL() ! sum ( B / T**2 dT ) for Lambertian downwelling radiance
    REAL(jprb), POINTER :: up_solar(:,:,:)         => NULL() ! sum( B * dT )
    REAL(jprb), POINTER :: down_solar(:,:,:)       => NULL() ! sum ( B / T**2 dT )
    REAL(jprb), POINTER :: meanrad_up(:,:)         => NULL()
    REAL(jprb), POINTER :: meanrad_down(:,:)       => NULL()
    REAL(jprb), POINTER :: meanrad_down_p(:,:)     => NULL()
    REAL(jprb), POINTER :: meanrad_up_solar(:,:)   => NULL()
    REAL(jprb), POINTER :: meanrad_down_solar(:,:) => NULL()
    REAL(jprb), POINTER :: down_ref(:,:,:)         => NULL()
    REAL(jprb), POINTER :: down_p_ref(:,:,:)       => NULL()
    REAL(jprb), POINTER :: down_ref_solar(:,:,:)   => NULL()
    REAL(jprb), POINTER :: cloudy(:,:)             => NULL()

    ! Quantities used for solar single-scattering
    REAL(jprb), POINTER :: fac1_2(:,:,:) => NULL()
    REAL(jprb), POINTER :: fac3_2(:,:)   => NULL()
    REAL(jprb), POINTER :: fac4_2(:,:,:) => NULL()
    REAL(jprb), POINTER :: fac5_2(:,:,:) => NULL()
    REAL(jprb), POINTER :: fac6_2(:,:,:) => NULL()
    REAL(jprb), POINTER :: fac7_2(:,:)   => NULL()
    REAL(jprb), POINTER :: fac4_3(:,:)   => NULL()
    REAL(jprb), POINTER :: fac5_3(:,:)   => NULL()
  END TYPE rttov_radiance_aux

  !> @internal Raytracing variables
  TYPE rttov_raytracing
    REAL(jprb), POINTER :: ltick(:,:)         => NULL() ! (levels,profiles) Layer thickness (km)
    REAL(jprb), POINTER :: hgpl(:,:)          => NULL() ! (levels,profiles) Level geometric height (km)
    REAL(jprb), POINTER :: dmair(:,:)         => NULL() ! (levels,profiles) Density of moist air (kg/m3)
    REAL(jprb), POINTER :: refractivity(:,:)  => NULL() ! (levels,profiles) Refractive index of air
    REAL(jprb), POINTER :: r(:,:)             => NULL() ! (levels,profiles)
    REAL(jprb), POINTER :: r_r(:,:)           => NULL() ! (levels,profiles)
    REAL(jprb), POINTER :: z_r(:,:)           => NULL() ! (layers,profiles)
    REAL(jprb), POINTER :: ratoesun(:,:)      => NULL() ! (layers,profiles)
    REAL(jprb), POINTER :: ratoesat(:,:)      => NULL() ! (layers,profiles)
    REAL(jprb), POINTER :: zasun(:,:)         => NULL() ! (layers,profiles) Sine of local angle of sun-surface path
    REAL(jprb), POINTER :: zasat(:,:)         => NULL() ! (layers,profiles) Sine of local angle of surface-satellite path
    REAL(jprb), POINTER :: int(:,:)           => NULL() ! (levels,profiles) Integrated layer values for dmair (hPa/(kg/m3))
    REAL(jprb), POINTER :: ztemp(:,:)         => NULL() ! (levels,profiles)
    REAL(jprb), POINTER :: ppw(:,:)           => NULL() ! (levels,profiles) Partial pressure of water vapour (hPa)
    REAL(jprb), POINTER :: dispco2(:,:)       => NULL() ! (levels,profiles) Correction factor for ref. index due to CO2
    REAL(jprb), POINTER :: pathsat(:,:)       => NULL() ! (layers,profiles) Secant of local angle of surface-satellite path
    REAL(jprb), POINTER :: pathsun(:,:)       => NULL() ! (layers,profiles) Secant of local angle of sun-surface path
    REAL(jprb), POINTER :: patheff(:,:)       => NULL() ! (layers,profiles) Sum of pathsat and pathsun
  END TYPE rttov_raytracing

  !> @internal Sea-surface solar BRDF model variables
  TYPE rttov_sunglint_s
    ! Profile-related quantities
    REAL(jprb) :: windsp
    REAL(jprb) :: wangl
    REAL(jprb) :: dazng
    REAL(jprb) :: zensat
    REAL(jprb) :: zensun

    ! Variables for Yoshimori wave facet model
    REAL(jprb) :: gamma_sq
    REAL(jprb) :: gamma_o
    REAL(jprb) :: gamma_p
    REAL(jprb) :: csi
    REAL(jprb) :: alfa
    REAL(jprb) :: omega
    REAL(jprb) :: gammax
    REAL(jprb) :: q_shad
    REAL(jprb) :: a_shad
    REAL(jprb) :: b_shad
    REAL(jprb) :: lambda_a
    REAL(jprb) :: lambda_b
    REAL(jprb) :: c_shad
    REAL(jprb) :: p_prime
    REAL(jprb) :: g_shad
    REAL(jprb) :: fac1
    REAL(jprb) :: pxy_gammaxy
    REAL(jprb) :: glint

    ! Variables for JONSWAP spectrum
    REAL(jprb) :: x_u
    REAL(jprb) :: alfa1
    REAL(jprb) :: omega_m

    ! Variables for Elfouhaily et al spectrum
    REAL(jprb) :: omega_c
  END TYPE rttov_sunglint_s

  !> @internal Sea-surface solar BRDF model variables
  TYPE rttov_sunglint
    TYPE(rttov_sunglint_s), POINTER :: s(:) => NULL() ! (nprofiles)

    ! Variables for JONSWAP spectrum
    REAL(jprb), POINTER :: beta(:,:)      => NULL() ! (nomega,nprofiles)
    REAL(jprb), POINTER :: psi(:,:)       => NULL() ! (nomega,nprofiles)

    ! Variables for Elfouhaily et al spectrum
    REAL(jprb), POINTER :: c(:,:)         => NULL() ! (nk,nprofiles)
    REAL(jprb), POINTER :: lpm(:,:)       => NULL() ! (nk,nprofiles)
    REAL(jprb), POINTER :: gamma_exp(:,:) => NULL() ! (nk,nprofiles)
    REAL(jprb), POINTER :: jp(:,:)        => NULL() ! (nk,nprofiles)
    REAL(jprb), POINTER :: fpexp(:,:)     => NULL() ! (nk,nprofiles)
    REAL(jprb), POINTER :: fm(:,:)        => NULL() ! (nk,nprofiles)
    REAL(jprb), POINTER :: dk(:,:)        => NULL() ! (nk,nprofiles)
    REAL(jprb), POINTER :: sk2(:,:)       => NULL() ! (0:nk,nprofiles)
  END TYPE rttov_sunglint

  !> @internal Phase function Legendre coefficients for DOM
  TYPE rttov_phasefn_lcoef
    REAL(jprb), POINTER :: legcoef(:,:,:) => NULL() ! Phase function Legendre coefficients (0:dom_nstr,0:1,nlay)
  END TYPE rttov_phasefn_lcoef

  !> @internal optical depths and other quantities for aerosols and clouds
  TYPE rttov_scatt_ir_aercld
    REAL(jprb), POINTER :: opdpsca(:,:)    => NULL() ! nadir scattering optical depth for all particle types (nlay,nchan)
    REAL(jprb), POINTER :: opdpabs(:,:)    => NULL() ! nadir absorption optical depth for all particle types (nlay,nchan)
    REAL(jprb), POINTER :: opdpscabpr(:,:) => NULL() ! nadir scattering op dep weighted by bpr for all part. types (nlay,nchan)
    REAL(jprb), POINTER :: opdp(:,:)       => NULL() ! nadir total layer optical depth (nlay,nchan)
    REAL(jprb), POINTER :: opdpsun(:,:)    => NULL() ! nadir total layer optical depth on solar effective path (nlay,nchan)
    REAL(jprb), POINTER :: phintup(:,:,:)  => NULL() ! interpolated phase functions for upward scattering (nparticles,nlay,nchan)
    REAL(jprb), POINTER :: phintdo(:,:,:)  => NULL() ! interpolated phase functions for downward scattering (nparticles,nlay,nchan)
    REAL(jprb), POINTER :: phtotup(:,:)    => NULL() ! combined interpolated phase function for upward scattering (nlay,nchan)
    REAL(jprb), POINTER :: phtotdo(:,:)    => NULL() ! combined interpolated phase function for downward scattering (nlay,nchan)
    REAL(jprb), POINTER :: partsca(:,:,:)  => NULL() ! scattering coefficient for each particle type (km-1) (nparticles,nlay,nchan)
    REAL(jprb), POINTER :: sca(:,:)        => NULL() ! total layer scattering coefficient (km-1) (nlay,nchan)
    REAL(jprb), POINTER :: partbpr(:,:,:)  => NULL() ! bpr value for each particle type (nparticles,nlay,nchan)
  END TYPE rttov_scatt_ir_aercld

  !> @internal Visible/IR scattering optical depths and related data
  TYPE rttov_transmission_scatt_ir
    TYPE(rttov_scatt_ir_aercld), POINTER :: aer => NULL()
    TYPE(rttov_scatt_ir_aercld), POINTER :: cld => NULL()
    REAL(jprb), POINTER :: ray_sca(:,:)           => NULL() ! Rayleigh nadir scattering extinction (nlay,nchan)
    REAL(jprb), POINTER :: phup(:,:,:)            => NULL() ! aer/aer+cld interp phasefn for upward scattering (0:1,nlay,nchan)
    REAL(jprb), POINTER :: phdo(:,:,:)            => NULL() ! aer/aer+cld interp phasefn for downward scattering (0:1,nlay,nchan)
    REAL(jprb), POINTER :: opdpext(:,:,:)         => NULL() ! aer/aer+cld extinction optical depth (0:1,nlay,nchan)
    REAL(jprb), POINTER :: opdpabs(:,:,:)         => NULL() ! nadir aer/aer+cld absorption optical depth (0:1,nlay,nchan)
    REAL(jprb), POINTER :: opdpsca(:,:,:)         => NULL() ! nadir aer/aer+cld scattering optical depth (0:1,nlay,nchan)
    REAL(jprb), POINTER :: opdpac(:,:,:)          => NULL() ! aer/aer+cld accumulated optical depth (nlev,0:ncolumn,nchan)
    REAL(jprb), POINTER :: opdpacl(:,:,:)         => NULL() ! aer/aer+cld layer optical depth (0:1,nlay,nchan)
    REAL(jprb), POINTER :: opdpacsun(:,:,:)       => NULL() ! aer/aer+cld accumulated op dep (solar path) (nlev,0:ncolumn,nchan)
    REAL(jprb), POINTER :: opdpaclsun(:,:,:)      => NULL() ! aer/aer+cld layer optical depth (solar path) (0:1,nlay,nchan)
    REAL(jprb), POINTER :: ssa_solar(:,:,:)       => NULL() ! aer/aer+cld single-scattering albedo (solar path) (0:1,nlay,nchan)
    REAL(jprb), POINTER :: ssa_thermal(:,:,:)     => NULL() ! aer/aer+cld single-scattering albedo (0:1,nlay,nchan)
    REAL(jprb), POINTER :: layerod_solar(:,:,:)   => NULL() ! nadir aer/aer+cld layer op deps w/ gas (solar path) (0:1,nlay,nchan)
    REAL(jprb), POINTER :: layerod_thermal(:,:,:) => NULL() ! nadir aer/aer+cld layer optical depths w/ gas (0:1,nlay,nchan)
    TYPE(rttov_phasefn_lcoef), POINTER :: phasefn(:) => NULL() ! phase function Legendre coefficients (nchan)
  END TYPE rttov_transmission_scatt_ir

  !> @internal PC coefs for each predictor channel set
  TYPE rttov_coef_pccomp1
    INTEGER(jpim)           :: fmv_pc_npred                  ! Number of predictors in the regression set
    INTEGER(jpim), POINTER  :: predictindex(:)     => NULL() ! Predictors channel indices
    REAL(jprb)   , POINTER  :: coefficients(:,:)   => NULL() ! Regression coefficients
    REAL(jprb)   , POINTER  :: coefficients_t(:,:) => NULL() ! Regression coefficients transposed
  END TYPE rttov_coef_pccomp1

  !> @internal PC eigenvectors
  TYPE rttov_coef_pccomp2
    REAL(jprb)   , POINTER  :: eigenvectors  (:,:) => NULL() ! Eigenvectors
    REAL(jprb)   , POINTER  :: eigenvectors_t(:,:) => NULL() ! Transposed Eigenvectors
  END TYPE rttov_coef_pccomp2

  !> @internal PC-RTTOV coefs
  TYPE rttov_coef_pccomp
    INTEGER(jpim)           :: fmv_pc_comp_pc                    ! PC file version number
    INTEGER(jpim)           :: fmv_pc_cld                        ! File trained for cloudy simulations
    INTEGER(jpim)           :: fmv_pc_aer                        ! File trained for aerosol simulations
    INTEGER(jpim)           :: fmv_pc_naer_types                 ! Number of aerosol types for aerosol simulations
    INTEGER(jpim)           :: fmv_pc_nlte                       ! File trained for NLTE simulations
    INTEGER(jpim)           :: fmv_pc_msets                      ! Maximum number of regression sets
    INTEGER(jpim)           :: fmv_pc_bands                      ! Number of bands
    INTEGER(jpim)           :: fmv_pc_mnum                       ! Maximum number of eigenvectors
    INTEGER(jpim)           :: fmv_pc_mchn                       ! Maximum number of channels
    INTEGER(jpim)           :: fmv_pc_nchn                       ! Number of channels
    INTEGER(jpim)           :: fmv_pc_nchn_noise                 ! Number of channels for which instrument noise is available
    INTEGER(jpim)           :: fmv_pc_nche                       ! Number of channels for which emissisity coefs are available
    INTEGER(jpim)           :: fmv_pc_gas                        ! Number of gases for which a reference profile is given
    INTEGER(jpim)           :: fmv_pc_gas_lim                    ! Number of gases for which min/max limit profiles are given
    INTEGER(jpim), POINTER  :: fmv_pc_sets   (:)       => NULL() ! Number of regression sets in each band
    INTEGER(jpim), POINTER  :: emiss_chn     (:)       => NULL() ! Number of channels for which emissivity coefficients are stored
    REAL   (jprb), POINTER  :: emiss_c1      (:)       => NULL() ! Emissivity coefficient
    REAL   (jprb), POINTER  :: emiss_c2      (:)       => NULL() ! Emissivity coefficient
    REAL   (jprb), POINTER  :: emiss_c3      (:)       => NULL() ! Emissivity coefficient
    REAL   (jprb), POINTER  :: emiss_c4      (:)       => NULL() ! Emissivity coefficient
    REAL   (jprb), POINTER  :: emiss_c5      (:)       => NULL() ! Emissivity coefficient
    REAL   (jprb), POINTER  :: emiss_c6      (:)       => NULL() ! Emissivity coefficient
    REAL   (jprb), POINTER  :: emiss_c7      (:)       => NULL() ! Emissivity coefficient
    REAL   (jprb), POINTER  :: emiss_c8      (:)       => NULL() ! Emissivity coefficient
    REAL   (jprb), POINTER  :: emiss_c9      (:)       => NULL() ! Emissivity coefficient
    INTEGER(jpim)           :: fmv_pc_nlev                       ! Number of reference profile levels
    REAL(jprb), POINTER     :: ref_pc_prfl_p (:)       => NULL() ! pressure  (hPa)       (levels)
    REAL(jprb), POINTER     :: ref_pc_prfl_mr(:,:)     => NULL() ! mixing ratio (ppmv)   (levels)
    REAL(jprb), POINTER     :: lim_pc_prfl_tmin(:)     => NULL() ! Profile limit :temperature (K)
    REAL(jprb), POINTER     :: lim_pc_prfl_tmax(:)     => NULL() ! Profile limit :temperature (K)
    REAL(jprb), POINTER     :: lim_pc_prfl_qmin(:)     => NULL() ! Profile limit :water vapour (ppmv)
    REAL(jprb), POINTER     :: lim_pc_prfl_qmax(:)     => NULL() ! Profile limit :water vapour (ppmv)
    REAL(jprb), POINTER     :: lim_pc_prfl_ozmin(:)    => NULL() ! Profile limit :ozone (ppmv)
    REAL(jprb), POINTER     :: lim_pc_prfl_ozmax(:)    => NULL() ! Profile limit :ozone (ppmv)
    REAL(jprb), POINTER     :: lim_pc_prfl_gasmin(:,:) => NULL() ! Profile limit :additional gases (ppmv)
    REAL(jprb), POINTER     :: lim_pc_prfl_gasmax(:,:) => NULL() ! Profile limit :additional gases (ppmv)
    REAL(jprb), POINTER     :: lim_pc_prfl_aermin(:,:) => NULL() ! Profile limit :aerosols (cm^-3) (layers)
    REAL(jprb), POINTER     :: lim_pc_prfl_aermax(:,:) => NULL() ! Profile limit :aerosols (cm^-3) (layers)
    REAL(jprb)              :: lim_pc_prfl_pmin                  ! Surface pressure (hPa)
    REAL(jprb)              :: lim_pc_prfl_pmax                  ! Surface pressure (hPa)
    REAL(jprb)              :: lim_pc_prfl_tsmin                 ! Surface temperature (K)
    REAL(jprb)              :: lim_pc_prfl_tsmax                 ! Surface temperature (K)
    REAL(jprb)              :: lim_pc_prfl_skmin                 ! Skin temperature (K)
    REAL(jprb)              :: lim_pc_prfl_skmax                 ! Skin temperature (K)
    REAL(jprb)              :: lim_pc_prfl_wsmin                 ! 10m wind speed (m/s)
    REAL(jprb)              :: lim_pc_prfl_wsmax                 ! 10m wind speed (m/s)
    REAL(jprb), POINTER     :: co2_pc_ref(:)           => NULL() ! Fixed co2 profile to be used in the computation of PC's (ppmv)
    REAL(jprb), POINTER     :: n2o_pc_ref(:)           => NULL() ! Fixed n2o profile to be used in the computation of PC's (ppmv)
    REAL(jprb), POINTER     :: co_pc_ref(:)            => NULL() ! Fixed co  profile to be used in the computation of PC's (ppmv)
    REAL(jprb), POINTER     :: ch4_pc_ref(:)           => NULL() ! Fixed ch4 profile to be used in the computation of PC's (ppmv)
    REAL(jprb), POINTER     :: co2_pc_min(:)           => NULL() ! Fixed co2 profile to be used in the computation of PC's (ppmv)
    REAL(jprb), POINTER     :: n2o_pc_min(:)           => NULL() ! Fixed n2o profile to be used in the computation of PC's (ppmv)
    REAL(jprb), POINTER     :: co_pc_min(:)            => NULL() ! Fixed co  profile to be used in the computation of PC's (ppmv)
    REAL(jprb), POINTER     :: ch4_pc_min(:)           => NULL() ! Fixed ch4 profile to be used in the computation of PC's (ppmv)
    REAL(jprb), POINTER     :: co2_pc_max(:)           => NULL() ! Fixed co2 profile to be used in the computation of PC's (ppmv)
    REAL(jprb), POINTER     :: n2o_pc_max(:)           => NULL() ! Fixed n2o profile to be used in the computation of PC's (ppmv)
    REAL(jprb), POINTER     :: co_pc_max(:)            => NULL() ! Fixed co  profile to be used in the computation of PC's (ppmv)
    REAL(jprb), POINTER     :: ch4_pc_max(:)           => NULL() ! Fixed ch4 profile to be used in the computation of PC's (ppmv)
    REAL(jprb), POINTER     :: noise_in(:)             => NULL() ! Noise values for the channels whose radiances are
                                                                 ! reconstructed using principal components
    REAL(jprb), POINTER     :: noise(:)                => NULL() ! Noise values for the channels whose radiances are
                                                                 ! used as predictors in the computation of principal components
    REAL(jprb), POINTER     :: noise_r(:)              => NULL() ! Reciprocal noise
    INTEGER(jpim), POINTER  :: ff_ori_chn_in(:)        => NULL()
    REAL(jprb),    POINTER  :: ff_cwn_in(:)            => NULL() ! central wave number of reconstructed radiances (cm-1)
    REAL(jprb),    POINTER  :: ff_bco_in (:)                     ! band correction offset (K)
    REAL(jprb),    POINTER  :: ff_bcs_in (:)           => NULL() ! band correction slope (K/K)
    REAL(jprb),    POINTER  :: planck1_in(:)           => NULL() ! C1 * Nu**3 (mW/(m2*sr*cm-1))
    REAL(jprb),    POINTER  :: planck2_in(:)           => NULL() ! C2 * Nu (K)
    TYPE(rttov_coef_pccomp1), POINTER:: pcreg(:,:)     => NULL()
    TYPE(rttov_coef_pccomp2), POINTER:: eigen(:)       => NULL()
  END TYPE rttov_coef_pccomp

  !> @internal Optical property data for one cld/aer particle type
  TYPE rttov_optp_data
    CHARACTER(LEN=4)       :: name                         ! Name of particle type
    REAL(jprb)             :: confac                       ! Unit conversion factor ([particle.cm^-3]/[g.m^-3])
    INTEGER(jpim)          :: nrelhum                      ! Number of rel. hum. values
    REAL(jprb),    POINTER :: relhum(:)          => NULL() ! Relative humidity values in % (nrelhum)
    INTEGER(jpim)          :: ndeff                        ! Number of Deff values
    REAL(jprb),    POINTER :: deff(:)            => NULL() ! Deff values in microns (ndeff)
    REAL(jprb),    POINTER :: abs(:,:,:)         => NULL() ! Absorption coefficients (nrelhum, ndeff, nchan)
    REAL(jprb),    POINTER :: sca(:,:,:)         => NULL() ! Scattering coefficients (nrelhum, ndeff, nchan)
    REAL(jprb),    POINTER :: bpr(:,:,:)         => NULL() ! bpr parameters for Chou-scaling (nrelhum, ndeff, nchan)
    REAL(jprb),    POINTER :: pha(:,:,:,:)       => NULL() ! Phase functions (nphangle, nrelhum, ndeff, nchan_pha)
    INTEGER(jpim), POINTER :: nmom(:,:)          => NULL() ! Number of Legendre moments (nrelhum, nchan)
    REAL(jprb),    POINTER :: legcoef(:,:,:,:)   => NULL() ! Legendre coefficients (maxnmom, nrelhum, ndeff, nchan)
  END TYPE rttov_optp_data

  !> @internal Top-level structure for one set of cld/aer optical properties
  TYPE rttov_optp
    INTEGER(jpim)           :: version                     ! File format version number
    INTEGER(jpim)           :: id                          ! 1/2/0 for OPAC/CAMS/user aerosol files,
                                                           ! 0/1 for CLW using OPAC/Segelstein (1984) ref index,
                                                           ! 0 for Baum ice properties
    INTEGER(jpim)           :: nchan                       ! Number of channels
    INTEGER(jpim)           :: nchan_pha                   ! Number of solar-affected channels
    INTEGER(jpim), POINTER  :: chan_pha(:)       => NULL() ! List of solar-affected channel numbers (nchan_pha)
    INTEGER(jpim), POINTER  :: chan_pha_index(:) => NULL() ! Indices of channels in solar-affected channel list (nchan)
    INTEGER(jpim)           :: maxnmom                     ! Maximum number of Legendre moments
    INTEGER(jpim)           :: nphangle                    ! Number of phase function angles
    REAL(jprb),    POINTER  :: phangle(:)        => NULL() ! Phase function angles (nphangle)
    TYPE(rttov_phasefn_int) :: phfn_int                    ! Contains info for interpolating phase functions
    INTEGER(jpim)           :: ntypes                      ! Number of particle types
    TYPE(rttov_optp_data), POINTER :: data(:)    => NULL() ! Optical property data (ntypes)
  END TYPE rttov_optp

  !> @internal Baran scheme variables
  TYPE rttov_optp_baran
    ! interpolation factors for frequency
    INTEGER(jpim), POINTER  :: iwn(:)            => NULL()
    INTEGER(jpim), POINTER  :: jwn(:)            => NULL()
    REAL(jprb),    POINTER  :: dx_dwn(:)         => NULL()

    REAL(jprb),    POINTER  :: q(:)              => NULL() ! Gaussian quadrature to calculate phase fn Leg. coefs
    REAL(jprb),    POINTER  :: w(:)              => NULL() !
    TYPE(rttov_phasefn_int) :: phfn_int                    ! Phase fn interpolation data
  END TYPE rttov_optp_baran

  !> @internal Visible/IR scattering coefs
  TYPE rttov_coef_scatt
    TYPE(rttov_optp)       :: optp_aer
    TYPE(rttov_optp)       :: optp_wcl_opac
    TYPE(rttov_optp)       :: optp_wcl_deff
    TYPE(rttov_optp)       :: optp_icl_baum
    TYPE(rttov_optp_baran) :: optp_icl_baran2014
    TYPE(rttov_optp_baran) :: optp_icl_baran2018
  END TYPE rttov_coef_scatt

  !> @internal MFASIS LUT data for one channel
  TYPE rttov_mfasis_lut
    INTEGER(jpim)       :: nluts               ! number of LUTs for this channel
    REAL(jprb), POINTER :: qint(:,:) => NULL() ! integrated water vapour values for each LUT (2,nluts)
    REAL(jprb), POINTER :: data(:,:) => NULL() ! LUT(s) for a single channel (all dimensions concatenated per LUT) (:,nluts)
  END TYPE rttov_mfasis_lut

  !> @internal MFASIS LUT axis definition
  TYPE rttov_mfasis_axis
    CHARACTER(LEN=32)   :: name        ! axis variable name
    INTEGER(jpim)       :: dim_type    ! type of dimension (interp. method):
                                       !   0 : theta+ Fourier index (k)
                                       !   1 : theta- Fourier index (l)
                                       !   2 : Albedo (no interp., use Jonkheid 2012)
                                       !   3 : tau-like (interpolate in logarithm)
                                       !   4 : Reff-like (interpolate linearly)
                                       !   5 : scattering angle (constr. lin. interp.)
    INTEGER(jpim)       :: nvalues
    REAL(jprb), POINTER :: values(:) => NULL() ! axis values (nvalues)
  END TYPE rttov_mfasis_axis

  !> @internal MFASIS cloud or aerosol LUT structure
  TYPE rttov_coef_mfasis
    INTEGER(jpim) :: file_type                                ! 1=>clouds, 2=>aerosols
    INTEGER(jpim) :: version                                  ! version number, may be useful for handling
                                                              !   updates to LUT format later and checking
                                                              !   compatibility with code
    CHARACTER(LEN=132) :: readme_lut(100)                     ! description of LUT creation, etc

    INTEGER(jpim)          :: ndims                           ! number of LUT dimensions
    TYPE(rttov_mfasis_axis), POINTER :: lut_axes(:) => NULL() ! table axis values (ndims)

    INTEGER(jpim)          :: nparticles                      ! number of particle types in the table (for
                                                              !   clouds this is 2, water and ice cloud; for
                                                              !   aerosol this has yet to be decided)
    INTEGER(jpim), POINTER :: aer_types(:)          => NULL() ! the indices of the RTTOV aerosol types used
                                                              !   in training aerosol LUTs (nparticles)

    INTEGER(jpim)          :: clw_scheme                      ! liquid cloud scheme used for training LUT
    INTEGER(jpim)          :: ice_scheme                      ! ice cloud scheme used for training LUT
    REAL(jprb)             :: maxzenangle                     ! maximum satellite and solar zenith angle used for training LUT

    INTEGER(jpim)          :: nchannels                       ! number of supported channels
    INTEGER(jpim), POINTER :: channel_list(:)       => NULL() ! supported channel numbers (nchannels)
    INTEGER(jpim)          :: nchannels_coef                  ! number of channels in the associated rtcoef file
    INTEGER(jpim), POINTER :: channel_lut_index(:)  => NULL() ! index of each supported channel in lut(:)  (nchannels_coef)
    TYPE(rttov_mfasis_lut), POINTER :: lut(:)       => NULL() ! LUT(s) for each supported channel (nchannels)
  END TYPE rttov_coef_mfasis

  !> @internal HT-FRTC scheme coefficients (Optical properties + PC related)
  TYPE rttov_coef_htfrtc
    INTEGER(jpim)          :: n_f
    REAL(jprb), POINTER    :: freq(:)             => NULL()
    INTEGER(jpim)          :: n_gas_l
    INTEGER(jpim), POINTER :: gasid_l(:)          => NULL()
    INTEGER(jpim)          :: n_p
    REAL(jprb), POINTER    :: p(:)                => NULL()
    INTEGER(jpim)          :: n_val_l
    INTEGER(jpim)          :: n_b
    REAL(jprb), POINTER    :: val_b(:)            => NULL()
    INTEGER(jpim)          :: n_lt
    REAL(jprb), POINTER    :: val_lt(:)           => NULL()
    REAL(jprb), POINTER    :: coef_l(:,:,:,:)     => NULL()
    INTEGER(jpim)          :: n_cont
    REAL(jprb), POINTER    :: coef_ct(:,:)        => NULL()
    REAL(jprb), POINTER    :: coef_ctt(:,:,:)     => NULL()
    REAL(jprb), POINTER    :: coef_b(:,:)         => NULL()
    REAL(jprb), POINTER    :: coef_lt(:)          => NULL()
    INTEGER(jpim)          :: n_ssemp
    REAL(jprb), POINTER    :: coef_ssemp(:,:)     => NULL()
    INTEGER(jpim)          :: n_iremis
    REAL(jprb), POINTER    :: coef_iremis(:,:)    => NULL()
    INTEGER(jpim)          :: n_pc
    INTEGER(jpim)          :: n_pc_oc
    REAL(jprb), POINTER    :: coef_pdt(:,:)       => NULL()
    REAL(jprb), POINTER    :: val_mean(:)         => NULL()
    REAL(jprb), POINTER    :: val_norm(:)         => NULL()
    INTEGER(jpim)          :: n_ch
    REAL(jprb), POINTER    :: sensor_freq(:)      => NULL()
    REAL(jprb), POINTER    :: ch_mean(:)          => NULL()
    REAL(jprb), POINTER    :: pc(:,:)             => NULL()
    REAL(jprb), POINTER    :: mixed_ref_frac(:,:) => NULL()
    INTEGER(jpim)          :: n_mftlb
    REAL(jprb), POINTER    :: mftlb(:)            => NULL()
    INTEGER(jpim), POINTER :: addf(:,:)           => NULL()
    INTEGER(jpim), POINTER :: addch(:,:)          => NULL()
  END TYPE rttov_coef_htfrtc

  !> @internal RTTOV coefs
  TYPE rttov_coefs
    LOGICAL(jplm)           :: initialised = .FALSE.
    TYPE(rttov_coef)        :: coef
    TYPE(rttov_coef_scatt)  :: coef_scatt
    TYPE(rttov_coef_pccomp) :: coef_pccomp
    TYPE(rttov_coef_mfasis) :: coef_mfasis_cld
    TYPE(rttov_coef_mfasis) :: coef_mfasis_aer
    TYPE(rttov_coef_htfrtc) :: coef_htfrtc
  END TYPE rttov_coefs

  !> @internal IR cloud column variables
  TYPE rttov_ircld
    INTEGER(jpim), POINTER  :: ncolumn(:)        => NULL()
    INTEGER(jpim), POINTER  :: ncolumnref(:)     => NULL()
    INTEGER(jpim), POINTER  :: iloop(:)          => NULL()
    INTEGER(jpim), POINTER  :: icount(:)         => NULL()
    INTEGER(jpim), POINTER  :: icouncol(:)       => NULL()
    INTEGER(jpim), POINTER  :: icount1(:)        => NULL()
    REAL(jprb)   , POINTER  :: xcolclr(:)        => NULL()
    INTEGER(jpim), POINTER  :: icldarr   (:,:,:) => NULL()
    REAL(jprb)   , POINTER  :: xcolref1  (:,:,:) => NULL()
    REAL(jprb)   , POINTER  :: xcolref2  (:,:,:) => NULL()
    INTEGER(jpim), POINTER  :: indexcol  (:,:)   => NULL()
    INTEGER(jpim), POINTER  :: icount1ref(:,:)   => NULL()
    INTEGER(jpim), POINTER  :: iloopin   (:,:)   => NULL()
    INTEGER(jpim), POINTER  :: iflag     (:,:)   => NULL()
    REAL(jprb)   , POINTER  :: xcol      (:,:)   => NULL()
    REAL(jprb)   , POINTER  :: xcolminref(:,:)   => NULL()
    REAL(jprb)   , POINTER  :: xcolref   (:,:)   => NULL()
    REAL(jprb)   , POINTER  :: cldcfr    (:,:)   => NULL()
    REAL(jprb)   , POINTER  :: maxcov    (:,:)   => NULL()
    REAL(jprb)   , POINTER  :: xcolmax   (:,:)   => NULL()
    REAL(jprb)   , POINTER  :: xcolmin   (:,:)   => NULL()
    REAL(jprb)   , POINTER  :: a         (:,:)   => NULL()
    REAL(jprb)   , POINTER  :: ntotref   (:,:)   => NULL()
    LOGICAL(jplm), POINTER  :: flag      (:,:)   => NULL()
  END TYPE rttov_ircld

  !> @internal Profiles input to Discrete Ordinates algorithm
  TYPE rttov_profile_dom
    INTEGER(jpim)          :: nlayers              ! Number of layers in DOM profile
    LOGICAL(jplm)          :: surface              ! Flag to indicate if surface is visible
    REAL(jprb),    POINTER :: layerod(:) => NULL() ! Layer total optical depth (nlayers)
    INTEGER(jpim), POINTER :: laymap(:)  => NULL() ! Mapping from DOM layer to user layer
  END TYPE

  !> @internal Direct model internal state of Discrete Ordinates algorithm
  TYPE rttov_dom_state
    ! Anything of size 0:naz is 0:0 for thermal channels, 0:nstr-1 for solar

    ! Solar only
    INTEGER(jpim), POINTER :: nazloops(:)      => NULL() ! 0:ncolumns
    REAL(jprb),    POINTER :: z(:,:,:,:)       => NULL() ! nstr,0:profiles(1)%nlayers,0:1,0:nstr-1

    ! Thermal only
    REAL(jprb),    POINTER :: y0(:,:,:)        => NULL() ! nstr,nlayers,0:ncolumns
    REAL(jprb),    POINTER :: y1(:,:,:)        => NULL() ! nstr,nlayers,0:ncolumns

    ! Thermal and solar
    REAL(jprb),    POINTER :: thisrad(:)       => NULL() ! 0:ncolumns
    REAL(jprb),    POINTER :: radsurfup(:)     => NULL() ! 0:ncolumns
    REAL(jprb),    POINTER :: radsurfup_sum(:) => NULL() ! 0:ncolumns
    REAL(jprb),    POINTER :: x(:,:,:)         => NULL() ! nstr*maxnlayers, 0:naz, 0:ncolumns
    REAL(jprb),    POINTER :: xp(:,:,:,:,:)    => NULL() ! nstr/2,nstr/2,0:profiles(1)%nlayers,0:1,0:nstr-1
    REAL(jprb),    POINTER :: eval(:,:,:,:)    => NULL() ! nstr/2,profiles(1)%nlayers,0:1,0:naz
  END TYPE

  !> @internal Direct model internal state of MFASIS
  TYPE rttov_mfasis_refl
    REAL(jprb)             :: refl                       ! Reflectance computed by MFASIS
    REAL(jprb)             :: refl_wv(3)                 ! Reflectance computed by MFASIS for different water vapour profiles
    REAL(jprb),    POINTER :: refl_lin_coef(:) => NULL() ! Sensitivity of reflectance (matrix elements of linear computations)
  END TYPE rttov_mfasis_refl

  !> @internal RTTOV internal state
  TYPE rttov_traj
!
! Hold RTTOV trajectory; these variables can have counterparts in TL, AD, K,
! and/or their dimensions are known before running RTTOV (nlevels, nprofiles, nchannels)
! it is possible to allocate these variables from outside RTTOV
!
    TYPE(rttov_profile),      POINTER :: profiles_coef(:) => NULL()
    TYPE(rttov_profile),      POINTER :: profiles_int(:)  => NULL()
    TYPE(rttov_predictors)            :: predictors
    TYPE(rttov_profile_aux)           :: aux_prof
    TYPE(rttov_profile_aux_coef)      :: aux_prof_coef
    TYPE(rttov_raytracing)            :: raytracing
    TYPE(rttov_raytracing)            :: raytracing_coef
    TYPE(rttov_opdp_path)             :: opdp_path
    TYPE(rttov_opdp_path)             :: opdp_path_coef

    REAL(jprb),               POINTER :: diffuse_refl(:)  => NULL() ! Surface refl for downwelling radiation (nchanprof)
    REAL(jprb),               POINTER :: fresnrefl(:)     => NULL() ! Fresnel reflection coefficients (nchanprof)
    TYPE(rttov_sunglint)              :: sunglint
    TYPE(rttov_ircld)                 :: ircld
    TYPE(rttov_transmission_scatt_ir) :: transmission_scatt_ir

    TYPE(rttov_radiance_aux)          :: auxrad

    TYPE(rttov_coefs),        POINTER :: coefs            => NULL()
    INTEGER(jpim)                     :: nchanprof
    INTEGER(jpim)                     :: nlevels
    INTEGER(jpim)                     :: nlayers
    TYPE(rttov_options)               :: opts
  END TYPE rttov_traj

  !> @internal RTTOV internal state
  TYPE rttov_traj_dyn
!
! Hold RTTOV trajectory; these variables have counterparts in TL, AD, K,
! but their dimensions are unknown when RTTOV starts running (ncolumns)
!
    LOGICAL(jplm)                     :: from_tladk = .FALSE.  ! True if this direct model traj_dyn came from the TL/AD/K
    INTEGER(jpim)                     :: ncolumns = -1  ! This initialisation used to determine alloc status
    TYPE(rttov_radiance_aux)          :: auxrad_column
    TYPE(rttov_transmission_scatt_ir) :: transmission_scatt_ir_dyn
    TYPE(rttov_transmission_aux)      :: transmission_aux
    TYPE(rttov_profile_dom), POINTER  :: profiles_dom_thermal(:,:) => NULL()
    TYPE(rttov_profile_dom), POINTER  :: profiles_dom_solar(:,:)   => NULL()
    TYPE(rttov_dom_state),   POINTER  :: dom_state_thermal(:)      => NULL()
    TYPE(rttov_dom_state),   POINTER  :: dom_state_solar(:)        => NULL()
    TYPE(rttov_mfasis_refl), POINTER  :: mfasis_refl(:,:)          => NULL()
  END TYPE rttov_traj_dyn

  !> @internal RTTOV internal state (unmodified regressed layer optical depths for v13 pred)
  TYPE rttov_opdp_ref_coef
    ! Dimensions are (nlayers,nchanprof) where nlayers refers to coef layers
    REAL(jprb),     POINTER :: od_mg_ref(:,:)     => NULL()
    REAL(jprb),     POINTER :: od_wv_ref(:,:)     => NULL()
    REAL(jprb),     POINTER :: od_wvcont_ref(:,:) => NULL()
    REAL(jprb),     POINTER :: od_o3_ref(:,:)     => NULL()
    REAL(jprb),     POINTER :: od_co2_ref(:,:)    => NULL()
    REAL(jprb),     POINTER :: od_n2o_ref(:,:)    => NULL()
    REAL(jprb),     POINTER :: od_co_ref(:,:)     => NULL()
    REAL(jprb),     POINTER :: od_ch4_ref(:,:)    => NULL()
    REAL(jprb),     POINTER :: od_so2_ref(:,:)    => NULL()
    REAL(jprb),     POINTER :: od_tot_ref(:,:)    => NULL()
  END TYPE rttov_opdp_ref_coef

  !> @internal RTTOV internal state (optical depths, transmittances)
  TYPE rttov_path_traj_trans
    ! Structure to hold optical depth and transmittance data
    ! within static trajectory.
    REAL(jprb),     POINTER :: tau_level(:,:)       => NULL() ! sat to level transmittance
    REAL(jprb),     POINTER :: tau_surf(:)          => NULL()
    REAL(jprb),     POINTER :: od_level(:,:)        => NULL() ! sat to level optical depth
    REAL(jprb),     POINTER :: od_singlelayer(:,:)  => NULL() ! single layer optical depth
    REAL(jprb),     POINTER :: od_frac(:)           => NULL()

    REAL(jprb),     POINTER   :: opdp_ref_coef(:,:) => NULL() ! layer optical depth before threshold (v789 pred)
    TYPE(rttov_opdp_ref_coef) :: opdp_ref                     ! layer optical depth before threshold (v13 pred)
  END TYPE rttov_path_traj_trans

  !> @internal RTTOV internal state (direct model only)
  TYPE rttov_traj_sta
!
! Hold RTTOV trajectory; these variables do not have counterparts in TL, AD, K
!
    LOGICAL(jplm)          :: do_opdep_calc        ! flag to indicate RTTOV gas optical depth calculation required
    LOGICAL(jplm)          :: do_solar_opdep_calc  ! flag to indicate RTTOV gas optical depth calculation required
                                                   !   for solar channels
    LOGICAL(jplm)          :: do_mfasis            ! flag to indicate MFASIS is being used
    LOGICAL(jplm)          :: do_rayleigh_param    ! flag to indicate Rayleigh parameterisation is active
    LOGICAL(jplm)          :: do_rayleigh_ss       ! flag to indicate Rayleigh single-scattering calculation is active
    LOGICAL(jplm)          :: do_rayleigh_dom      ! flag to indicate Rayleigh DOM calculation is active

    LOGICAL(jplm), POINTER :: thermal(:) => NULL() ! switch for thermal calculations (nchanprof)
    LOGICAL(jplm), POINTER :: solar(:)   => NULL() ! switch for solar calculations (nchanprof)
    LOGICAL(jplm)          :: dothermal            ! flag to indicate thermal calculations required
    LOGICAL(jplm)          :: dosolar              ! flag to indicate solar calculations required

    INTEGER(jpim)          :: dom_nstreams         ! number of streams for DOM
    LOGICAL(jplm)          :: plane_parallel       ! switch for plane-parallel geometry (no atmospheric curvature)

    REAL(jprb),    POINTER :: solar_spec_esd(:) => NULL() ! Solar spectrum adjusted for esd (nchanprof) (mW/m2/cm-1)
    REAL(jprb),    POINTER :: refl_norm(:)      => NULL() ! Normalisation factor for solar surface reflectance

    LOGICAL(jplm), POINTER               :: do_lambertian(:)     => NULL()
    TYPE(rttov_path_traj_trans), POINTER :: thermal_path1        => NULL()
    TYPE(rttov_path_traj_trans), POINTER :: solar_path2          => NULL()
    TYPE(rttov_path_traj_trans), POINTER :: solar_path1          => NULL()
    TYPE(rttov_geometry), POINTER        :: angles(:)            => NULL()
    TYPE(rttov_profile),  POINTER        :: profiles_coef_ref(:) => NULL()
    TYPE(rttov_chanprof), POINTER        :: chanprof_in(:)       => NULL()
    TYPE(rttov_chanprof), POINTER        :: chanprof_pc(:)       => NULL()
    REAL(jprb), POINTER                  :: pc_aer_ref(:,:,:)    => NULL() ! For PC-RTTOV aerosol simulations
    REAL(jprb), POINTER                  :: pc_aer_min(:,:,:)    => NULL() !
    REAL(jprb), POINTER                  :: pc_aer_max(:,:,:)    => NULL() !
  END TYPE rttov_traj_sta

  !> @internal Used for coef testing
  TYPE rttov_lbl_check
    REAL(jprb), POINTER :: atm_layer(:,:)       => NULL()
    REAL(jprb), POINTER :: atm_layer_path2(:,:) => NULL()
    LOGICAL(jplm)       :: plane_geometry
  END TYPE rttov_lbl_check

END MODULE rttov_types
