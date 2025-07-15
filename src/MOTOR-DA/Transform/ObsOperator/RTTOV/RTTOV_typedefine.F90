!!--------------------------------------------------------------------------------------------------
! PROJECT           : rttov_get_utils
! AFFILIATION       : GBA-MWF(SIMI)
! AUTOHR(S)         : Yali Wu
! VERSION           : V 0.0
! HISTORY           :
!   Created by Yali Wu, 2021/8/27
! Information       : Add interfaces to call rttov N.L..
!                   : Only the clear-sky capability is enabled so far.
!!--------------------------------------------------------------------------------------------------
!> @brief
!!
MODULE rttov_typedefine_m

  USE rttov_CONST
  USE parkind1, ONLY: jpim
  USE rttov_TYPES, ONLY: rttov_options, rttov_profile, rttov_coefs, &
                         rttov_radiance, rttov_radiance2, rttov_transmission, rttov_emissivity, &
                         rttov_reflectance, rttov_chanprof, rttov_predictors

  USE kinds_m, ONLY: r_kind, i_kind

  TYPE rttov_in_fixvars_t
    ! Define those mandatory input variables, including 3D profiles and 2D sfc variables
    REAL(r_kind)         ::  seaice, landmask, snowc, sfctype, vegtype, landtype, lai, soiltype   ! surface types
    REAL(r_kind)         ::  psfc, tsfc, q2m, t2m, u10m, v10m, tskin ! 2D surface variables
    REAL(r_kind)         ::  elevation, latitude, longitude, zenangle, azangle, sunzenangle, sunazangle
    REAL(r_kind), ALLOCATABLE, DIMENSION(:)    ::  pres, temp, qvapor, uwind, vwind    ! 3D profiles; POINTER is not working
    REAL(r_kind), ALLOCATABLE, DIMENSION(:)    ::  tsoil, qsoil      ! soil temperature and moisture
    REAL(r_kind), ALLOCATABLE, DIMENSION(:)    ::  rttov_vlevel       ! ???Also set 2D variables to 3D, but with nlev=1???
  END TYPE rttov_in_fixvars_t

  TYPE rttov_in_cldvars_t
    ! Defines those for calculating cloudy radiances
    REAL(r_kind), ALLOCATABLE, DIMENSION(:)    ::  qcloud, qrain, qice, qsnow, qgraup      ! Hydrometeors
    REAL(r_kind), ALLOCATABLE, DIMENSION(:)    ::  re_cloud, re_rain, re_ice, re_snow, re_graup  ! Effective radius
    REAL(r_kind), ALLOCATABLE, DIMENSION(:)    ::  cldfrac      ! Cloud fraction
  END TYPE rttov_in_cldvars_t

  TYPE rttov_in_gasvars_t
    REAL(r_kind), ALLOCATABLE, DIMENSION(:)    ::  h2o, co2, o3, ch4, co, so2            ! gases
    INTEGER(i_kind)       :: gas_ids
    LOGICAL               :: ppmv_gasunit = .FALSE. ! if convert mixing ratio to ppmv
    INTEGER(i_kind)       :: nlevels, nmixed, nwater, nozone, nwvcont, nco2, nn2o, nco, nch4, nso2, npmc
    ! CHARACTER(len=12) :: RTTOV_Absorbers(ngases_max) = (/ 'Mixed_gases ', 'Water_vapour', 'Ozone       ', &
    ! 'WV_Continuum', 'CO2         ', 'N2O         ', &
    ! 'CO          ', 'CH4         ', 'SO2         '/)
    ! INTEGER(i_kind) :: RTTOV_Absorbers_ID(ngases_max) = (/gas_id_mixed, gas_id_watervapour, &
    ! gas_id_ozone, gas_id_wvcont, gas_id_co2, &
    ! gas_id_n2o, gas_id_co, gas_id_ch4, gas_id_so2/)
  END TYPE rttov_in_gasvars_t

  TYPE rttov_instrument_t
    INTEGER(i_kind)    :: chans(1)
    CHARACTER(len=50), ALLOCATABLE, DIMENSION(:)    :: gas_name
    INTEGER(i_kind)    :: nchans = 1
    INTEGER(i_kind)    :: ngases
    LOGICAL               :: rttov_clouds = .FALSE.
    LOGICAL               :: rttov_aerosols = .FALSE.
    LOGICAL               :: use_kmodel = .TRUE.
    LOGICAL               :: extra_qc = .FALSE.
    LOGICAL               :: write_profs = .FALSE.
    LOGICAL               :: write_jac = .FALSE.
    LOGICAL               :: only_over_sea = .FALSE.
  END TYPE rttov_instrument_t

  TYPE rttov_in_userdefs_t
    CHARACTER(len=50)   :: inst_name, platform_name, sensor_name
    ! Note to follow the naming rules of each instrument/platform/sensor
    INTEGER(i_kind)    :: nplats, ninsts, nsensors
    TYPE(rttov_instrument_t)  :: rttov_insts(1)
    LOGICAL               :: debug_dev = .FALSE.
    INTEGER(i_kind)    :: interp_method = 1
  END TYPE rttov_in_userdefs_t

  TYPE rttov_setup_input_t
    TYPE(rttov_in_fixvars_t), ALLOCATABLE   :: rttov_in_fixvars
    TYPE(rttov_in_cldvars_t), ALLOCATABLE   :: rttov_in_cldvars
    TYPE(rttov_in_gasvars_t), ALLOCATABLE   :: rttov_in_gasvars
    TYPE(rttov_in_userdefs_t), ALLOCATABLE  :: rttov_in_userdefs
  END TYPE rttov_setup_input_t

  TYPE rttov_in_calc_t
    LOGICAL, POINTER                   :: calcemis(:) => NULL()  ! in
    LOGICAL, POINTER                   :: calcrefl(:) => NULL()
    TYPE(rttov_reflectance), POINTER   :: reflectance(:) => NULL()
    TYPE(rttov_reflectance), POINTER   :: reflectance_k(:) => NULL()
    TYPE(rttov_reflectance), POINTER   :: reflectance_tl(:) => NULL()
    TYPE(rttov_reflectance), POINTER   :: reflectance_ad(:) => NULL()
    TYPE(rttov_radiance)               :: radiance, radiance_k        ! out
    TYPE(rttov_radiance)               :: radiance_tl, radiance_ad        ! out
    TYPE(rttov_radiance2)              :: radiance2                   ! out
    TYPE(rttov_transmission)           :: transm, transm_k    ! out
    TYPE(rttov_transmission)           :: transm_tl, transm_ad    ! out
    TYPE(rttov_emissivity), POINTER    :: emis(:) => NULL(), emis_k(:) => NULL()  ! inout
    TYPE(rttov_emissivity), POINTER    :: emis_tl(:) => NULL(), emis_ad(:) => NULL()  ! inout
    TYPE(rttov_profile), POINTER       :: profs(:) => NULL(), profs_k(:) => NULL()    ! in
    TYPE(rttov_profile), POINTER       :: profs_tl(:) => NULL(), profs_ad(:) => NULL()
    !TYPE(rttov_profile), ALLOCATABLE  :: profs(:) , profs_k(:)    ! in
    TYPE(rttov_chanprof), POINTER      :: chanprof(:) => NULL()   ! in
  END TYPE rttov_in_calc_t

  TYPE rttov_diag_output
    ! For pressure levels
    INTEGER :: nvars
    REAL(r_kind), ALLOCATABLE, DIMENSION(:) :: pres
    ! For BT and BTclr, unit ((mW/m2/sr/cm-1))
    REAL(r_kind) :: bt, btclr
    ! For transmittance output, unit (0-1)
    REAL(r_kind) :: tau_total, tau_total_cld
    REAL(r_kind), ALLOCATABLE, DIMENSION(:) :: tau_levels, tau_levels_cld
    ! For up/down radiance/emission output, unit (mW/m2/sr/cm-1)
    REAL(r_kind) :: upclr, dnclr
    REAL(r_kind), ALLOCATABLE, DIMENSION(:) :: up, down, surf
    ! For Jacobian output
    REAL(r_kind), ALLOCATABLE, DIMENSION(:) :: t, q, qc, qi, jac_t, jac_q, jac_qc, jac_qi
    REAL(r_kind), ALLOCATABLE, DIMENSION(:) :: cldfrac
    REAL(r_kind), ALLOCATABLE, DIMENSION(:) :: dtaudlnp
    REAL(r_kind), ALLOCATABLE, DIMENSION(:, :) :: outvar
    ! For emissivity output
    REAL(r_kind) :: emis_out
  END TYPE rttov_diag_output

END MODULE rttov_typedefine_m
