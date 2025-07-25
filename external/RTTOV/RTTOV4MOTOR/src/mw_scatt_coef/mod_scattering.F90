module mod_scattering

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

! Module for use with create_tables_spectra and associated subroutines

! Current Code Owner: SAF NWP

! History:
! Version   Date        Comment
! -------   ----        -------
!           04/09/2008   transferred to fcm (Amy Doherty)
!           12/03/2010   restructuring (Alan Geer)
!           03/03/2011   New snow options + cleaning (Alan Geer)
!           15/03/2013   Fully-flexible PSD, shape & density options (Alan Geer)
!           01/04/2014   Flexible single-scattering database inputs (Alan Geer)
!           31/10/2017   Adapted for ARTS-SSDB use (Jana Mendrok)
!           10/01/2018   add option to select liquid water permittivity model (Katrin Lonitz)
!           20/03/2020   SI and more PSD flexibility (Alan Geer)

use parkind1,  only: jprb, jpim, jplm

implicit none


!* constants
real (kind=jprb)  , parameter :: pi          = 3.141592654_jprb

real (kind=jprb)  , parameter :: not_available = 1.0e20_jprb

! Density [kg m^-3]
real (kind=jprb)  , parameter :: density_wat = 1e3_jprb
real (kind=jprb)  , parameter :: density_air = 1.225_jprb
real (kind=jprb)  , parameter :: density_ice = 917.0_jprb
real (kind=jprb)  , parameter :: density_lmt = 50.0_jprb 

complex (kind=jprb), parameter :: perm_air    = (1.0006_jprb, 0.0_jprb)

!* constants for melting routines - not converted to SI as expected to be replaced by DDA
!* models in the near future
integer (kind=jpim), parameter :: max_iter = 1000
real    (kind=jprb), parameter :: eps      = 0.001_jprb
      
real (kind=jprb)  , parameter :: onethird = 0.3333_jprb

real (kind=jprb)  , parameter :: t_0c = 273.15_jprb      ! T in K at 0 degrees Celcius 

real (kind=jprb)  , parameter :: t_ml = 275.0_jprb       ! T at bottom of melting layer in K
real (kind=jprb)  , parameter :: t_fl = 273.15_jprb      ! temperature at FL in k
real (kind=jprb)  , parameter :: rh   = 0.98_jprb        ! relative humidity

real (kind=jprb)  , parameter :: ka = 5.69e-05_jprb      ! heat conductivity in cal/cm/s/c
real (kind=jprb)  , parameter :: dv = 0.211_jprb         ! diffusivity of vapor in air in cm*cm/s
real (kind=jprb)  , parameter :: mw = 18.0_jprb          ! molecular weight of water in g/mol
real (kind=jprb)  , parameter :: le = 797.3_jprb         ! latent heat of evaporation in cal/g
real (kind=jprb)  , parameter :: lm = 79.7_jprb          ! latent heat of melting in cal/g
real (kind=jprb)  , parameter :: r  = 8.314e+07_jprb     ! gas constant in erg/mol/k
real (kind=jprb)  , parameter :: va = 0.1346_jprb        ! kinematic viscosity of air in cm*cm/s

real (kind=jprb)  , parameter :: g       = 981.0_jprb    ! gravity [cm/s^2]
real (kind=jprb)  , parameter :: cd_snow = 1.2_jprb      ! drag coeff. snow []
real (kind=jprb)  , parameter :: cd_grpl = 1.0_jprb      ! drag coeff. graupel []
real (kind=jprb)  , parameter :: cd_hail = 0.6_jprb      ! drag coeff. hail []

!* parameterizations
real (kind=jprb), parameter :: lwc2rr_a = 20.95_jprb ! non-SI units (see n0_t.F90)
real (kind=jprb), parameter :: lwc2rr_b = 1.12_jprb

real    (kind=jprb), parameter :: wc_slope  = 0.01_jprb
integer (kind=jpim), parameter :: wc_offset = 301

integer (kind=jpim), parameter :: lmax    = 500 !Maximum length of filename/path and strings in text files
integer(kind=jpim),  parameter :: nmaxconfig = 100 !Max length of config section

!* Hydro-tables
!*************
!*
integer (kind=jpim), parameter :: n_temp = 70
integer (kind=jpim), parameter :: n_lwc  = 401
integer (kind=jpim), parameter :: n_dia  = 100

!* Hydrometeor types
!********************
!*
integer (kind=jpim), parameter :: n_type = 7
integer (kind=jpim), parameter :: i_rain      = 1
integer (kind=jpim), parameter :: i_snow      = 2
integer (kind=jpim), parameter :: i_graupel   = 3
integer (kind=jpim), parameter :: i_aggregate = 4
integer (kind=jpim), parameter :: i_clw       = 5
integer (kind=jpim), parameter :: i_ciw       = 6
integer (kind=jpim), parameter :: i_totalice  = 7

! PSD options
integer (kind=jpim), parameter :: psd_gamma           = 1
integer (kind=jpim), parameter :: psd_marshall_palmer = 2
integer (kind=jpim), parameter :: psd_field_2005      = 3
integer (kind=jpim), parameter :: psd_field_2007      = 4
integer (kind=jpim), parameter :: psd_modified_gamma  = 5
integer (kind=jpim), parameter :: psd_heymsfield_2013 = 6
integer (kind=jpim), parameter :: psd_mh97            = 7

! Parameters of the modified gamma PSD (see Petty and Huang, 2011, JAS)
integer (kind=jpim), parameter :: modified_gamma_n0     = 1
integer (kind=jpim), parameter :: modified_gamma_mu     = 2
integer (kind=jpim), parameter :: modified_gamma_lambda = 3
integer (kind=jpim), parameter :: modified_gamma_gamma  = 4
integer (kind=jpim), parameter :: modified_gamma_fixed  = -1

! Default Marshall-Palmer PSD
real (kind=jprb)   , parameter :: n0_liquid = 8e6_jprb ! [m^-4]
real (kind=jprb)   , parameter :: n0_frozen = 4e6_jprb ! [m^-4]

! Field parameters
real (kind=jprb)   , parameter :: d_field_min = 1e-4_jprb     ! [m]
logical (kind=jplm), parameter :: lfield_ignore_temp = .false.! Ignore F07 temperature dependence 
real (kind=jprb)   , parameter :: f07_fixed_temp = 253.0_jprb ! [K] fixed temperature to use with above

! Density options
integer (kind=jpim), parameter :: dens_wilson_ballard_1999 = 1
integer (kind=jpim), parameter :: dens_jones_1995 = 2
integer (kind=jpim), parameter :: dens_brown_francis_1995 = 3
integer (kind=jpim), parameter :: dens_fixed = 4

! Liquid water permittivity options
integer (kind=jpim), parameter :: permwat_liebe_89 = 1
integer (kind=jpim), parameter :: permwat_rosenkranz_15 = 2
integer (kind=jpim), parameter :: permwat_TKC_16 = 3

! Pure ice permittivity options
integer (kind=jpim), parameter :: permice_pre_v13 = 1
integer (kind=jpim), parameter :: permice_maetzler_2006 = 2

integer (kind=jpim), parameter :: ctype_max = 11
character (len=ctype_max) :: ctype(n_type) = (/&
& 'rain       ',&
& 'snow       ',&
& 'graupel    ',&
& 'aggregate  ',&
& 'cloud-water',&
& 'cloud-ice  ',&
& 'total_ice  '/)

! Scattering computation types
integer (kind=jpim), parameter :: n_scat   = 3
integer (kind=jpim), parameter :: i_mie    = 1
integer (kind=jpim), parameter :: i_liu    = 2
integer (kind=jpim), parameter :: i_arts   = 3

!* Density, [g/cm3], min/max diameters [m]
real    (kind=jprb) :: density_sphere (n_type) = (/ 1e3,     1e2,    4e2,    9e2,    1e3,     9e2,    5e2    /)
real    (kind=jprb) :: d_min_sphere   (n_type) = (/ 1e-4,    1e-4,   5e-4,   5e-3,   5e-6,    5e-6,   5e-6   /)
real    (kind=jprb) :: d_max_sphere   (n_type) = (/ 1e-2,    2e-2,   3e-3,   2e-2,   1e-4,    1e-4,   2e-2   /)
real    (kind=jprb) :: temp_offset    (n_type) = (/ 233.0,   203.0,  203.0,  203.0,  233.0,   203.0,  203.0  /)
logical (kind=jplm) :: is_frozen      (n_type) = (/ .false., .true., .true., .true., .false., .true., .true. /)

! see ice_density.F90:
real (kind=jprb), dimension (n_type) :: ice_a = (/ 1.00,   0.115,  0.0112, 0.90,   1.0,    0.917,  0.5    /)
real (kind=jprb), dimension (n_type) :: ice_b = (/ 0.0,    0.863,  0.815,  0.0,    0.0,    0.0,    0.0    /)

!* Min / max diameters [m] appropriate to Liu (2008) shapes but with a  
!* 1e-4 m minimum to avoid extrapolating Field et al. 2007 DSD
real (kind=jprb), dimension (0:10) :: d_liu_min = &
 & (/0.00015_jprb, d_field_min, d_field_min, d_field_min, 0.00015_jprb, d_field_min, &
 &  d_field_min, d_field_min, d_field_min, d_field_min, d_field_min/)
real (kind=jprb), dimension (0:10) :: d_liu_max = &
 & (/0.0048,     0.0033,      0.0025,      0.0032,      0.0050,       0.01, &
 &  0.01,        0.01,        0.01,        0.01,        0.01  /)

type modified_gamma_data
  integer (kind=jpim) :: i_mtype
  integer (kind=jpim), allocatable :: i_row(:)
  real    (kind=jprb), allocatable :: params(:,:)
  integer (kind=jpim), allocatable :: i_free(:)
end type modified_gamma_data

type diag_data
  character(len=lmax) :: filename
  integer (kind=jpim) :: itype
  real    (kind=jprb) :: iwc, freq, temp
  real    (kind=jprb) :: integrand(n_dia)
  real    (kind=jprb) :: renorm
end type diag_data

contains

! Returns water content (ice or liquid) in [kg/m3]
function get_lwc(iwc)
  real (kind=jprb) :: get_lwc
  integer (kind=jpim), intent(in) :: iwc
  get_lwc = 1e-3_jprb * 10.0_jprb ** (wc_slope * (iwc - wc_offset))
end function get_lwc

function get_d_dia(dia,n_dia)
  real (kind=jprb) :: get_d_dia
  integer (kind=jpim), intent(in) :: n_dia
  real    (kind=jprb), intent(in) :: dia(n_dia)
  get_d_dia = (dia(n_dia) - dia(1)) / (n_dia - 1)
end function get_d_dia

end module mod_scattering
