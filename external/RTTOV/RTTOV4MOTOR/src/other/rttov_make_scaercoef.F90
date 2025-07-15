! Description:
!> @file
!!   Generate an aerosol optical property file for use in RTTOV visible/IR
!!   scattering simulations.
!
!> @brief
!!   Generate an aerosol optical property file for use in RTTOV visible/IR
!!   scattering simulations.
!!
!! @details
!!   This executable allows users to generate their own "scaercoef" optical
!!   property files.
!!
!!   Usage: rttov_make_scaercoef.exe --config_file     STRING
!!                                   --rtcoef_file     STRING
!!                                   --output_file     STRING
!!                                   --normalise_size           (OPTIONAL)
!!                                   --nthreads        INTEGER  (OPTIONAL)
!!                                   --chou_only                (OPTIONAL)
!!                                   --max_nmom        INTEGER  (OPTIONAL)
!!
!!   where config_file    - the path to a file which defines the aerosol
!!                          properties (see below)
!!         rtcoef_file    - the path to an RTTOV optical depth coefficient file
!!         output_file    - the output aerosol property file name
!!         normalise_size - if this argument is present the size distributions
!!                          are explicitly normalised to 1 particle/cm^3. By
!!                          default this is not done.
!!         nthreads       - the number of threads to use (optional, defaults
!!                          to one, recommended to compile with OpenMP
!!                          enabled and to make use of multiple threads)
!!         chou_only      - if specified the output file will only contain
!!                          optical properties required for IR Chou-scaling
!!                          simulations: much smaller file, but cannot be
!!                          used with the DOM solver or for solar simulations
!!                          (default false i.e. generate optical properties
!!                          for all solvers)
!!         max_nmom       - the maximum number of Legendre coefficients to
!!                          calculate for each phase function; this defines
!!                          the max number of DOM streams that can be used in
!!                          RTTOV with these optical properties, default 128
!!
!!   Optical properties are computed using Mie theory. Optical properties can
!!   be calculated for one or more particle types. Each particle type may be
!!   hydrophobic (no dependence on relative humidity) or hydrophilic. In the
!!   latter case optical properties are calculated for multiple values of
!!   relative humidity (RH): within RTTOV these are interpolated to the RH of
!!   the layer.
!!
!!   Optical properties are computed for all channels in the rtcoef_file. If
!!   any channels are solar-enabled (applies only to "v9 predictor" coef files)
!!   then phase functions will be stored for these channels unless the
!!   chou_only argument is passed. If you only need optical properties for a
!!   subset of channels then use the rttov_conv_coef.exe executable to extract
!!   these channels first. The format of the output aerosol file is the same as
!!   that of the supplied rtcoef_file (ASCII, binary or HDF5).
!!
!!   Note that the number of levels of the coefficient file and the predictor
!!   version do not affect the calculated optical properties (except in the
!!   case of solar-enabled channels in v9 predictor files). Therefore you can
!!   use the same aerosol file with v7, v8 or v9 IASI coefficients on 54L or
!!   101L for non-solar calculations. If you want to include solar radiation
!!   you would need to generate an aerosol file using a v9 predictor IASI file.
!!
!!   The inputs required for each particle type (and for each RH value for
!!   hydrophilic species) are:
!!   - complex refractive index defined on some wavelength grid covering the
!!     spectral range of the channels in the rtcoef_file
!!   - particle size distribution defined on a suitable radius grid,
!!     normalised to 1 particle per unit volume
!!   - information about the particle density (see below)
!!
!!   This information is supplied via the config_file. A sample file is
!!   provided in data/example_rttov_make_scaercoef_config.txt.
!!   This defines the following:
!!   - number of particle types
!!   - for each particle type:
!!     - particle name
!!     - the number of relative humidity values
!!     - the relative humidity values
!!     - the factor used to convert mass mixing ratio to number density
!!       (usually for RH=0%); if <=0 it is computed from the size distribution
!!       for RH=0% and the supplied mass density value
!!     - mass density in g.cm^-3; this value is only used if the unit
!!       conversion factor is <=0
!!     - for each relative humidity value:
!!       - the path to a file containing the refractive indices
!!       - the path to a file containing the size distribution
!!
!!   There is no limit on the number of particle types.
!!
!!   The particle name is a four character field used to identify the particle.
!!
!!   For hydrophobic species the number of RH values should be set to 1 and
!!   the RH value should be 0%.
!!
!!   For hydrophilic species there is no limit to the number of RH values, but
!!   more values imply larger files. RH is specified as integer percent values
!!   starting at 0%, monotonically ascending, and not exceeding 99%. As an
!!   example, the OPAC aerosols use 8 RH values: 0, 50, 70, 80, 90, 95, 98, 99.
!!
!!   The complex refractive index data should be provided for a range of
!!   wavelengths covering the central wavenumbers of the channels in the
!!   supplied rtcoef_file. A warning is printed out if channels lie beyond the
!!   range of the data: optical properties are computed using constant-value
!!   extrapolation of the refractive indices in this case. The wavelength grid
!!   does not need to be evenly spaced. Example refractive index input files
!!   are provided in the data/ directory. The imaginary part of the refractive
!!   index (k) may be positive or negative, but is always treated as if it is
!!   negative i.e. the value -ABS(k) is used in the calculations.
!!
!!   The size distribution is specified as number density values normalised to
!!   one particle per unit volume for a range of particle radius values. The
!!   radius units are microns (um) and the normalised size distribution units
!!   are um^-1. The optical properties stored in the scaercoef files are
!!   computed per [particle.cm^-3]. Inside RTTOV they are then multiplied by
!!   the particle number density in [particles.cm^-3] to compute optical
!!   depths. If the --normalise_size argument is passed to the executable then
!!   the size distributions will be explicitly normalised by the code. The
!!   radius grid does not need to be evenly spaced. Example size distribution
!!   input files are provided in the data/ directory.
!!
!!   NB In all input files, comment lines begin with ! and are optional, and
!!   file paths must be enclosed in quotes "".
!!
!!   The conversion factor is used inside RTTOV to convert input aerosol mass
!!   mixing ratio values to particle number density. It is not used if the
!!   input aerosol units are number density so if you are only going to supply
!!   aerosol concentrations in number density then the value of this conversion
!!   factor does not matter (e.g. set it to 1).
!!
!!   The units of the conversion factor are [g.m^-3]/[particle.cm^-3] and it
!!   represents the mass density per particle density. Given aerosol mass
!!   density rho, the conversion factor is calculated as:
!!
!!          rho * 4pi/3 * integral[r^3.size_dist(r).dr]
!!
!!   where the integral is over the range of radius values for which the size
!!   distribution is specified in microns.
!!
!!   The conversion factor may be calculated and specified explicitly in the
!!   config_file. For example, it may be desirable to limit the calculation to
!!   a maximum particle size: this is done for the OPAC aerosols. Alternatively
!!   if the conversion factor in the config_file is set to a value less than or
!!   equal to zero then a value for the density (rho) of the aerosol material
!!   must be supplied in units of [g.cm^-3]. The conversion factor is then
!!   computed from the size distribution for RH=0% and this density value.
!!
!!   NB The value of confac stored in the scaercoef file is the reciprocal of
!!   the value given by the expression above, or, if greater than zero, the
!!   reciprocal of the value explicitly specified in the config_file.
!!
!!   You can test this executable using the example input data from the
!!   top-level RTTOV directory as follows:
!!
!!   $ bin/rttov_make_scaercoef.exe \
!!   > --config_file data/example_rttov_make_scaercoef_config.txt \
!!   > --rtcoef_file rtcoef_rttov13/rttov7pred54L/rtcoef_msg_3_seviri.dat \
!!   > --output_file scaercoef_msg_3_seviri_example.dat \
!!   > --normalise_size \
!!   > --nthreads 8 \
!!   > --chou_only
!!
!!   You can then compare your output file with the reference one:
!!     data/scaercoef_msg_3_seviri_example.dat
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
!    Copyright 2018, EUMETSAT, All Rights Reserved.
!
PROGRAM rttov_make_scaercoef

#include "throw.h"

  USE rttov_types, ONLY : rttov_coefs, rttov_optp_data

  USE rttov_const, ONLY :   &
    deg2rad,                &
    nphangle_lores,         &
    phangle_lores,          &
    nphangle_hires,         &
    phangle_hires,          &
    aer_id_user

  USE parkind1, ONLY : jprb, jpim, jplm

  USE rttov_coef_io_mod, ONLY : getlun, closelun

  USE rttov_getoptions, ONLY : initoptions, getoption, checkoptions

  USE rttov_unix_env, ONLY: rttov_exit

  USE rttov_scattering_mod, ONLY : &
    gauss_quad,                    &
    spline_interp,                 &
    normalise,                     &
    calc_legendre_coef_gauss,      &
    integrate

  USE rttov_mie_params_mod, ONLY : &
    interpolate_ref_index,         &
    calc_mie_ext_sca,              &
    calc_mie_phasefn,              &
    calc_mass,                     &
    init_new_optp

#ifdef _RTTOV_HDF
  USE hdf5
  USE rttov_hdf_mod
#endif

  IMPLICIT NONE

#include "rttov_read_ascii_coef.interface"
#include "rttov_read_binary_coef.interface"
#include "rttov_init_coef.interface"
#include "rttov_nullify_coefs.interface"
#include "rttov_write_ascii_scaercoef.interface"
#include "rttov_write_binary_scaercoef.interface"
#include "rttov_dealloc_coefs.interface"
#include "rttov_errorreport.interface"
#include "rttov_cmpuc.interface"
#include "rttov_skipcommentline.interface"
#include "rttov_bpr_init.interface"
#include "rttov_bpr_calc.interface"
#include "rttov_bpr_dealloc.interface"

#ifdef _RTTOV_HDF
#include "rttov_hdf_load.interface"
#include "rttov_hdf_save.interface"
#endif

  INTEGER(jpim), PARAMETER :: version = 1            ! Version of file format

  INTEGER(jpim), PARAMETER :: default_max_nmom = 128 ! Default max number of Leg. moments
  INTEGER(jpim), PARAMETER :: min_nmom = 32          ! Minimum number of Leg. moments
  INTEGER(jpim), PARAMETER :: ngauss = 1000          ! Size of Gaussian quadrature

  INTEGER(jpim), PARAMETER :: fid_config = 60        ! Logical unit used to read config file
  INTEGER(jpim), PARAMETER :: fid_data   = 62        ! Logical unit used to read/write input data files

  CHARACTER(LEN=512) :: config_file, rtcoef_file, output_file
  INTEGER(jpim)      :: nthreads, max_nmom
  LOGICAL(jplm)      :: chou_only, normalise_size

  INTEGER(jpim)      :: err, fid_coef
  TYPE(rttov_coefs)  :: coefs
  LOGICAL(jplm)      :: exists
  CHARACTER(LEN=32)  :: file_format, format_in
  CHARACTER(LEN=256) :: fname, f1
  TYPE(rttov_optp_data), POINTER :: optp_data

  INTEGER(jpim)      :: i, iaer, irh, iwav
  INTEGER(jpim)      :: ntypes, nrh, nwvl_data, nrarr
  INTEGER(jpim)      :: ifail
  CHARACTER(LEN=2)   :: rhumstr
  REAL(jprb)         :: error
  REAL(jprb)         :: ecoef, scoef, size_int, mass, rho
  REAL(jprb)         :: refindreal, refindimag
  REAL(jprb)         :: q(ngauss), w(ngauss), phase_interp(ngauss)
  INTEGER(jpim)      :: thisnphangle
  REAL(jprb)         :: thisphangle(nphangle_hires)

  REAL(jprb),    ALLOCATABLE :: rarr(:)             ! Radii in input data file
  REAL(jprb),    ALLOCATABLE :: n(:)                ! Size distribution in input data file
  REAL(jprb),    ALLOCATABLE :: wvl_data(:)         ! Wavelengths in input data file
  COMPLEX(jprb), ALLOCATABLE :: refind_data(:)      ! Refrac. indices in input data file
  REAL(jprb),    ALLOCATABLE :: wvl_chan(:)         ! Instrument channel wavelengths
  COMPLEX(jprb), ALLOCATABLE :: refind_chan(:)      ! Refrac. indices at chan wavelengths
  REAL(jprb),    ALLOCATABLE :: coef_phase(:,:,:,:) ! Phase functions for all coef channels

  ! ---------------------------------------------------------------------------

  TRY

  !----------------------------------------------------------------------------
  ! Some initial set up
  !----------------------------------------------------------------------------

  CALL initoptions("This program generates aerosol optical property (scaercoef) files")

  CALL getoption("--config_file", config_file, mnd=.TRUE._jplm, &
                 use="configuration file defining aerosol properties")

  CALL getoption("--rtcoef_file", rtcoef_file, mnd=.TRUE._jplm, &
                 use="input rtcoef file name")

  CALL getoption("--output_file", output_file, mnd=.TRUE._jplm, &
                 use="output file name")

  normalise_size = .FALSE.
  CALL getoption("--normalise_size", normalise_size, &
                 use="optional, if present size distributions are explicitly normalised to 1 part./cm^3")

  nthreads = -1
  CALL getoption("--nthreads", nthreads, &
                 use="integer, optional, default 1")
  IF (nthreads <= 0) nthreads = 1

  chou_only = .FALSE.
  CALL getoption("--chou_only", chou_only, &
                 use="optional, if present only compute parameters for Chou-scaling")

  max_nmom = -1
  CALL getoption("--max_nmom", max_nmom, &
                 use="number of Leg. moments to compute, optional, default 128")
  IF (max_nmom <= 0) max_nmom = default_max_nmom
  IF (max_nmom <= min_nmom) max_nmom = min_nmom

  CALL checkoptions()

  INQUIRE(file=config_file, exist=exists)
  IF (.NOT. exists) THEN
    PRINT *, 'Cannot find config file: '//TRIM(config_file)
    STOP
  ENDIF
  INQUIRE(file=rtcoef_file, exist=exists)
  IF (.NOT. exists) THEN
    PRINT *, 'Cannot find rtcoef file: '//TRIM(rtcoef_file)
    STOP
  ENDIF

  ! Gaussian quadrature for phase fn Leg. coef integration
  CALL gauss_quad(-1._jprb, 1._jprb, q, w)

  ! We only need bpr for the low-res phase angle array
  CALL rttov_bpr_init(err, phangle_lores(:))
  THROWM(err.NE.0, 'Error initialising bpr calculation tables')

  !----------------------------------------------------------------------------
  ! Read optical depth coefficient file
  !----------------------------------------------------------------------------

  CALL rttov_nullify_coefs(coefs)

  CALL getlun(err, fid_coef, f1, file_format, .FALSE._jplm, 'rtcoef', f=TRIM(rtcoef_file))
  THROW(err.NE.0)

  IF (rttov_cmpuc(file_format, 'unformatted')) THEN
    CALL rttov_read_binary_coef(err, coefs%coef, fid_coef)
    THROWM(err.NE.0, 'Cannot open binary coefficient file '//TRIM(rtcoef_file))
    format_in = 'unformatted'
  ELSE IF (rttov_cmpuc(file_format, 'formatted')) THEN
    CALL rttov_read_ascii_coef(err, coefs%coef, fid_coef)
    THROWM(err.NE.0, 'Cannot open ASCII coefficient file '//TRIM(rtcoef_file))
    format_in = 'formatted'
  ELSE IF (rttov_cmpuc(file_format, 'hdf5')) THEN
#ifndef _RTTOV_HDF
    err = errorstatus_fatal
    THROWM(err.NE.0, 'This program is not compiled with HDF5 capability; use RTTOV_HDF=1 with Makefile.PL')
#else
    CALL open_hdf(.TRUE._jplm, err)
    THROWM(err.NE.0, 'Error opening HDF5 interface')

    CALL rttov_hdf_load(err, f1, "/COEF", coef=coefs%coef)
    THROWM(err.NE.0, 'Cannot open HDF5 coefficient file '//TRIM(rtcoef_file))
    format_in = 'hdf5'

    CALL close_hdf(err)
    THROWM(err.NE.0, 'Error closing HDF5 interface')
#endif
  ELSE
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Unknown coefficient file format '//TRIM(file_format))
  ENDIF

  CALL closelun(err, fid_coef)
  THROW(err.NE.0)

  CALL rttov_init_coef(err, coefs%coef)

  !----------------------------------------------------------------------------
  ! Populate the meta-data in the aerosol params structure
  !----------------------------------------------------------------------------

  ! Open config file and read number of particle types
  OPEN(fid_config, file=config_file, form='formatted', iostat=err)
  THROWM(err.NE.0,'Error opening config file')

  ! Number of aerosol types
  CALL rttov_skipcommentline(fid_config, err)
  THROWM(err.NE.0,'Error reading config file')
  READ(fid_config, *, iostat=err) ntypes
  THROWM(err.NE.0,'Error reading config file, number of aerosol types')

  CALL init_new_optp(coefs%coef, coefs%coef_scatt%optp_aer, ntypes, chou_only)

  ! Version and ID numbers
  coefs%coef_scatt%optp_aer%version = version
  coefs%coef_scatt%optp_aer%id      = aer_id_user

  IF (coefs%coef_scatt%optp_aer%nchan_pha > 0) THEN
    ALLOCATE(coefs%coef_scatt%optp_aer%phangle(nphangle_hires))
    coefs%coef_scatt%optp_aer%phangle(:) = phangle_hires(:)
    coefs%coef_scatt%optp_aer%nphangle = nphangle_hires
  ENDIF

  ALLOCATE(wvl_chan(coefs%coef%fmv_chn), refind_chan(coefs%coef%fmv_chn))
  wvl_chan(:) = 10000._jprb / coefs%coef%ff_cwn(:)

  !----------------------------------------------------------------------------
  ! Process each aerosol type in turn
  !----------------------------------------------------------------------------

  DO iaer = 1, ntypes

    PRINT *,'Processing aerosol number ', iaer

    optp_data => coefs%coef_scatt%optp_aer%data(iaer)
    optp_data%ndeff = 1
    ALLOCATE(optp_data%deff(optp_data%ndeff))
    optp_data%deff = 0.

    !-----------------------------------------------------------
    ! Read information from config file
    !-----------------------------------------------------------

    ! Particle name
    CALL rttov_skipcommentline(fid_config, err)
    THROWM(err.NE.0,'Error reading config file')
    READ(fid_config, *, iostat=err) optp_data%name
    THROWM(err.NE.0,'Error reading config file, aerosol name')

    PRINT *,'Aerosol name ', optp_data%name

    ! Number of RH values
    CALL rttov_skipcommentline(fid_config, err)
    THROWM(err.NE.0,'Error reading config file')
    READ(fid_config, *, iostat=err) nrh
    THROWM(err.NE.0,'Error reading config file, number of RH values')
    IF (nrh < 1) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0,'Number of RH values must be >=1')
    ENDIF
    optp_data%nrelhum = nrh
    ALLOCATE(optp_data%relhum(nrh))

    ! RH values
    CALL rttov_skipcommentline(fid_config, err)
    THROWM(err.NE.0,'Error reading config file')
    READ(fid_config, *, iostat=err) optp_data%relhum
    THROWM(err.NE.0,'Error reading config file, RH values')
    IF (ANY(optp_data%relhum < 0._jprb) .OR. ANY(optp_data%relhum > 99._jprb)) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0,'RH values must lie in the interval [0,99]')
    ENDIF

    ! Unit conversion factor and density
    CALL rttov_skipcommentline(fid_config, err)
    THROWM(err.NE.0,'Error reading config file')
    READ(fid_config, *, iostat=err) optp_data%confac
    THROWM(err.NE.0,'Error reading config file, conversion factor')
    CALL rttov_skipcommentline(fid_config, err)
    THROWM(err.NE.0,'Error reading config file')
    READ(fid_config, *, iostat=err) rho
    THROWM(err.NE.0,'Error reading config file, density')

    ! Allocate aerosol coef structure arrays for optical parameters
    ALLOCATE(optp_data%abs(nrh,1,coefs%coef%fmv_chn))
    ALLOCATE(optp_data%sca(nrh,1,coefs%coef%fmv_chn))
    ALLOCATE(optp_data%bpr(nrh,1,coefs%coef%fmv_chn))
    optp_data%abs(:,:,:) = 0._jprb
    optp_data%sca(:,:,:) = 0._jprb
    optp_data%bpr(:,:,:) = 0._jprb
    ALLOCATE(optp_data%nmom(nrh,coefs%coef%fmv_chn))
    IF (chou_only) THEN
      ALLOCATE(optp_data%legcoef(1,nrh,1,coefs%coef%fmv_chn))
      optp_data%nmom = 0
      optp_data%legcoef = 0
    ELSE
      ALLOCATE(optp_data%legcoef(1:max_nmom+1,nrh,1,coefs%coef%fmv_chn))
      optp_data%nmom(:,:) = 0
      optp_data%legcoef(:,:,:,:) = 0._jprb
      IF (coefs%coef_scatt%optp_aer%nchan_pha > 0) THEN
        ALLOCATE(optp_data%pha(nphangle_hires,nrh,1,coefs%coef_scatt%optp_aer%nchan_pha))
        optp_data%pha(:,:,:,:) = 0._jprb
      ENDIF
    ENDIF

    ! Allocate temporary storage for phase functions for all channels
    ALLOCATE(coef_phase(nphangle_hires,nrh,1,coefs%coef%fmv_chn))
    coef_phase(:,:,:,:) = 0._jprb

    DO irh = 1, nrh

      WRITE(rhumstr, '(i2.2)') INT(optp_data%relhum(irh))
      PRINT *,'RH = '//rhumstr

      !-----------------------------------------------------------
      ! Read refractive index data, interpolate onto channels
      !-----------------------------------------------------------

      CALL rttov_skipcommentline(fid_config, err)
      THROWM(err.NE.0,'Error reading config file')
      READ(fid_config, *, iostat=err) fname
      THROWM(err.NE.0,'Error reading config file, refractive index filename')

      OPEN(fid_data, file=TRIM(ADJUSTL(fname)), form='formatted', iostat=err)
      THROWM(err.NE.0,'Error opening refractive index file')

      CALL rttov_skipcommentline(fid_data, err)
      THROWM(err.NE.0,'Error reading refractive index file')
      READ(fid_data, *, iostat=err) nwvl_data
      THROWM(err.NE.0,'Error reading refractive index file, number of wavelengths')

      ALLOCATE(wvl_data(nwvl_data), refind_data(nwvl_data))

      CALL rttov_skipcommentline(fid_data, err)
      THROWM(err.NE.0,'Error reading refractive index file')
      DO i = 1, nwvl_data
        READ(fid_data, *, iostat=err) wvl_data(i), refindreal, refindimag
        THROWM(err.NE.0,'Error reading refractive index file, data')
        refind_data(i) = CMPLX(refindreal, ABS(refindimag), KIND=jprb)

        IF (i > 1) THEN
          IF (wvl_data(i) <= wvl_data(i-1)) THEN
            err = errorstatus_fatal
            THROWM(err.NE.0,'Error reading ref. index file, wavelengths must be monotonically increasing')
          ENDIF
        ENDIF
      ENDDO

      IF (wvl_data(1) > MINVAL(wvl_chan) .OR. wvl_data(nwvl_data) < MAXVAL(wvl_chan)) THEN
        WARN('WARNING: wavelengths in ref. index file do not cover spectral range of channels')
      ENDIF

      CALL interpolate_ref_index(nwvl_data, wvl_data, refind_data, wvl_chan, refind_chan)

      DEALLOCATE(wvl_data, refind_data)

      CLOSE(fid_data)

      !-----------------------------------------------------------
      ! Read size distribution data
      !-----------------------------------------------------------

      CALL rttov_skipcommentline(fid_config, err)
      THROWM(err.NE.0,'Error reading config file')
      READ(fid_config, *, iostat=err) fname
      THROWM(err.NE.0,'Error reading config file, size distribution filename')

      OPEN(fid_data, file=TRIM(ADJUSTL(fname)), form='formatted', iostat=err)
      THROWM(err.NE.0,'Error opening size distribution file')

      CALL rttov_skipcommentline(fid_data, err)
      THROWM(err.NE.0,'Error reading size distribution file')
      READ(fid_data, *, iostat=err) nrarr
      THROWM(err.NE.0,'Error reading size distribution file, number of data values')

      ALLOCATE(rarr(nrarr), n(nrarr))

      CALL rttov_skipcommentline(fid_data, err)
      THROWM(err.NE.0,'Error reading size distribution file')
      DO i = 1, nrarr
        READ(fid_data, *, iostat=err) rarr(i), n(i)
        THROWM(err.NE.0,'Error reading size distribution file, data')

        IF (i > 1) THEN
          IF (rarr(i) <= rarr(i-1)) THEN
            err = errorstatus_fatal
            THROWM(err.NE.0,'Error reading size dist. file, radii must be monotonically increasing')
          ENDIF
        ENDIF
      ENDDO

      IF (normalise_size) THEN
        CALL integrate(nrarr, rarr, n, size_int, error, ifail)
        n = n / size_int
      ENDIF

      CLOSE(fid_data)

      !-----------------------------------------------------------
      ! If necessary compute unit conversion factor
      !-----------------------------------------------------------

      IF (irh == 1) THEN
        IF (optp_data%confac <= 0._jprb) THEN
          IF (rho <= 0._jprb) THEN
            err = errorstatus_fatal
            THROWM(err.NE.0,'Valid conversion factor or density value must be supplied')
          ENDIF
          CALL calc_mass(nrarr, rarr, n, rho, mass)
          ! mass = density * volume
          !       [g.cm^-3 . um^3] * 10^-6 => [g.m^-3/cm^-3]
          optp_data%confac = 1.E6_jprb / mass
        ELSE
          optp_data%confac = 1._jprb / optp_data%confac
        ENDIF
      ENDIF

      !-----------------------------------------------------------
      ! Compute optical properties
      !-----------------------------------------------------------

      PRINT *, 'Computing optical properties for ', coefs%coef%fmv_chn, ' channels'
      DO iwav = 1, coefs%coef%fmv_chn

        CALL calc_mie_ext_sca(nrarr, rarr, n, wvl_chan(iwav), refind_chan(iwav), ecoef, scoef, nthreads)

        optp_data%abs(irh,1,iwav) = (ecoef - scoef) * 0.001_jprb
        optp_data%sca(irh,1,iwav) = scoef * 0.001_jprb

        !------------------------------------------
        ! Calculate the phase function
        !------------------------------------------

        IF (coefs%coef%ss_val_chn(iwav) == 2) THEN
          thisnphangle = nphangle_hires
          thisphangle(1:thisnphangle) = phangle_hires(:)
        ELSE
          thisnphangle = nphangle_lores
          thisphangle(1:thisnphangle) = phangle_lores(:)
        ENDIF

        CALL calc_mie_phasefn(thisnphangle, thisphangle(1:thisnphangle), nrarr, rarr, &
                              n, wvl_chan(iwav), refind_chan(iwav), scoef, &
                              coef_phase(1:thisnphangle,irh,1,iwav), nthreads)

        IF (.NOT. chou_only) THEN
          !------------------------------------------
          ! Calculate the Legendre coefficients for the phase fn
          !------------------------------------------
          CALL spline_interp(thisnphangle, COS(thisphangle(thisnphangle:1:-1) * deg2rad), &
                             coef_phase(thisnphangle:1:-1,irh,1,iwav), ngauss, q, phase_interp)
          CALL normalise(ngauss, w, phase_interp)
          CALL calc_legendre_coef_gauss(q, w, phase_interp, min_nmom, max_nmom, &
                                        optp_data%nmom(irh,iwav), optp_data%legcoef(:,irh,1,iwav))
        ENDIF

        IF (coefs%coef%ss_val_chn(iwav) < 2) THEN
          !------------------------------------------
          ! Calculate the backscattering parameter bpr
          !------------------------------------------

          CALL rttov_bpr_calc(err, coef_phase(1:thisnphangle,irh,1,iwav), &
                              thisphangle, ecoef, nthreads=nthreads)
          optp_data%bpr(irh,1,iwav) = ecoef
        ENDIF

        IF (.NOT. chou_only) THEN
          !------------------------------------------
          ! Calculate the phase function for output
          !------------------------------------------

          ! Solar channels were already calculated on phangle_hires grid so only recompute for mixed channels
          IF (coefs%coef%ss_val_chn(iwav) == 1) THEN

            CALL calc_mie_phasefn(nphangle_hires, phangle_hires, nrarr, rarr, &
                                  n, wvl_chan(iwav), refind_chan(iwav), scoef, &
                                  coef_phase(:,irh,1,iwav), nthreads)

          ENDIF
        ENDIF

      ENDDO !iwav

      DEALLOCATE(rarr, n)

      !-----------------------------------------------------------
      ! Store phase function data in coefs structure
      !-----------------------------------------------------------

      ! coef_phase(:,:,:) contains phase fns in absolute channel indexes.
      ! Copy just the phase fns for solar-affected channels into the coef pha(:,:,:) array
      IF (coefs%coef_scatt%optp_aer%nchan_pha > 0) THEN
        DO i = 1, coefs%coef_scatt%optp_aer%nchan_pha
          optp_data%pha(:,irh,1,i) = coef_phase(:,irh,1,coefs%coef_scatt%optp_aer%chan_pha(i))
        ENDDO
      ENDIF

    ENDDO ! relhum

    !-----------------------------------------------------------
    ! Tidy up
    !-----------------------------------------------------------

    DEALLOCATE(coef_phase)

  ENDDO ! particle types

  IF (.NOT. chou_only) coefs%coef_scatt%optp_aer%maxnmom = max_nmom

  CLOSE(fid_config)

  !-------------------------------------------------------------------------------
  ! Write out aerosol coefficients
  !-------------------------------------------------------------------------------

  CALL getlun(err, fid_coef, f1, file_format, .true._jplm, 'scaercoef', f=output_file, form=format_in)
  THROW(err.NE.0)

  IF (TRIM(file_format) == 'formatted') THEN
    CALL rttov_write_ascii_scaercoef(err, coefs%coef, coefs%coef_scatt, fid_coef, verbose=.TRUE._jplm)
    THROWM(err.NE.0, 'Error writing ASCII aerosol coefficients')
  ELSE IF (TRIM(file_format) == 'unformatted') THEN
    CALL rttov_write_binary_scaercoef(err, coefs%coef, coefs%coef_scatt, fid_coef, verbose=.TRUE._jplm)
    THROWM(err.NE.0, 'Error writing binary aerosol coefficients')
  ELSE IF (TRIM(file_format) == 'hdf5') THEN
#ifndef _RTTOV_HDF
    err = errorstatus_fatal
    THROWM(err.NE.0, 'This program is not compiled with HDF5 capability; use RTTOV_HDF=1 with Makefile.PL')
#else
    CALL open_hdf(.TRUE._jplm, err)
    THROWM(err.NE.0, 'Error opening HDF5 interface')

    CALL rttov_hdf_save(err, f1, "/SCAER", create=.TRUE._jplm, &
                        scaercoef=coefs%coef_scatt, compress=.TRUE._jplm)
    THROWM(err.NE.0, 'Error writing HDF5 aerosol coefficients')

    CALL close_hdf(err)
    THROWM(err.NE.0, 'Error closing HDF5 interface')
#endif
  ELSE
    THROWM(err.NE.0, 'Unrecognised format for output')
  ENDIF

  CALL closelun(err, fid_coef)
  THROW(err.NE.0)

  !-------------------------------------------------------------------------------
  ! Clean up
  !-------------------------------------------------------------------------------

  DEALLOCATE(wvl_chan, refind_chan)

  CALL rttov_bpr_dealloc(err)
  THROWM(err.NE.0, 'Error deallocating bpr calculation tables')

  CALL rttov_dealloc_coefs(err, coefs)
  THROWM(err.NE.0, 'Error deallocating coefficients')

  PCATCH

END PROGRAM rttov_make_scaercoef
