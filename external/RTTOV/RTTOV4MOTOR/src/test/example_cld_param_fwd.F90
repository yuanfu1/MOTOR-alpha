! Description:
!> @file
!!   Example calling direct model visible/IR cloud simulation using
!!   explicitly input particle optical properties.
!
!> @brief
!!   Example calling direct model visible/IR cloud simulation using
!!   explicitly input particle optical properties.
!!
!! @details
!!   This example is most easily run using the run_example_cld_param_fwd.sh
!!   script (see the user guide).
!!
!!   This program requires the following files:
!!     the file containing input profiles (e.g. prof.dat)
!!     the file containing cloud optical parameter profiles for each channel
!!     the RTTOV optical depth coefficient file
!!
!!   The output is written to a file called example_cld_param_fwd_output.dat
!!
!!   This program is very similar to example_fwd.F90 with the addition
!!   of the cloud scattering.
!!
!!   You may wish to base your own code on this example in which case you
!!   should edit in particular the sections of code marked as follows:
!!       !================================
!!       !======Read =====start===========
!!            code to be modified
!!       !======Read ===== end ===========
!!       !================================
!!
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
!    Copyright 2017, EUMETSAT, All Rights Reserved.
!
PROGRAM example_cld_param_fwd

  ! rttov_const contains useful RTTOV constants
  USE rttov_const, ONLY :     &
         errorstatus_success, &
         errorstatus_fatal,   &
         platform_name,       &
         inst_name

  ! rttov_types contains definitions of all RTTOV data types
  USE rttov_types, ONLY :     &
         rttov_options,       &
         rttov_coefs,         &
         rttov_profile,       &
         rttov_transmission,  &
         rttov_radiance,      &
         rttov_chanprof,      &
         rttov_emissivity,    &
         rttov_reflectance,   &
         rttov_opt_param

  ! jpim, jprb and jplm are the RTTOV integer, real and logical KINDs
  USE parkind1, ONLY : jpim, jprb, jplm

  USE rttov_unix_env, ONLY : rttov_exit

  IMPLICIT NONE

#include "rttov_direct.interface"
#include "rttov_parallel_direct.interface"
#include "rttov_read_coefs.interface"
#include "rttov_dealloc_coefs.interface"
#include "rttov_alloc_direct.interface"
#include "rttov_init_emis_refl.interface"
#include "rttov_init_opt_param.interface"
#include "rttov_bpr_init.interface"
#include "rttov_bpr_calc.interface"
#include "rttov_bpr_dealloc.interface"
#include "rttov_legcoef_calc.interface"
#include "rttov_user_options_checkinput.interface"
#include "rttov_print_opts.interface"
#include "rttov_print_profile.interface"
#include "rttov_skipcommentline.interface"

  !--------------------------
  !
  INTEGER(KIND=jpim), PARAMETER :: iup   = 20   ! unit for input profile file
  INTEGER(KIND=jpim), PARAMETER :: ioout = 21   ! unit for output

  ! RTTOV variables/structures
  !====================
  TYPE(rttov_options)              :: opts                     ! Options structure
  TYPE(rttov_coefs)                :: coefs                    ! Coefficients structure
  TYPE(rttov_chanprof),    POINTER :: chanprof(:)    => NULL() ! Input channel/profile list
  LOGICAL(KIND=jplm),      POINTER :: calcemis(:)    => NULL() ! Flag to indicate calculation of emissivity within RTTOV
  TYPE(rttov_emissivity),  POINTER :: emissivity(:)  => NULL() ! Input/output surface emissivity
  LOGICAL(KIND=jplm),      POINTER :: calcrefl(:)    => NULL() ! Flag to indicate calculation of BRDF within RTTOV
  TYPE(rttov_reflectance), POINTER :: reflectance(:) => NULL() ! Input/output surface BRDF
  TYPE(rttov_profile),     POINTER :: profiles(:)    => NULL() ! Input profiles
  TYPE(rttov_opt_param)            :: cld_opt_param            ! Input cloud optical parameters
  TYPE(rttov_transmission)         :: transmission             ! Output transmittances
  TYPE(rttov_radiance)             :: radiance                 ! Output radiances

  INTEGER(KIND=jpim)               :: errorstatus              ! Return error status of RTTOV subroutine calls

  INTEGER(KIND=jpim) :: alloc_status
  CHARACTER(LEN=21)  :: NameOfRoutine = 'example_cld_param_fwd'

  ! variables for input
  !====================
  CHARACTER(LEN=256) :: coef_filename
  CHARACTER(LEN=256) :: prof_filename
  CHARACTER(LEN=256) :: file_cld_opt_param
  INTEGER(KIND=jpim) :: nthreads
  INTEGER(KIND=jpim) :: dosolar
  INTEGER(KIND=jpim) :: nlevels
  INTEGER(KIND=jpim) :: nprof
  INTEGER(KIND=jpim) :: nchannels
  INTEGER(KIND=jpim) :: nchanprof
  INTEGER(KIND=jpim), ALLOCATABLE :: channel_list(:)
  INTEGER(KIND=jpim) :: nphangle
  ! loop variables
  INTEGER(KIND=jpim) :: i, j, jch
  INTEGER(KIND=jpim) :: nch
  INTEGER(KIND=jpim) :: iprof, joff
  INTEGER            :: ios

  !- End of header --------------------------------------------------------

  ! The usual steps to take when running RTTOV are as follows:
  !   1. Specify required RTTOV options
  !   2. Read coefficients
  !   3. Allocate RTTOV input and output structures
  !   4. Set up the chanprof array with the channels/profiles to simulate
  !   5. Read input profile(s)
  !   6. Set up surface emissivity and/or reflectance
  !   7. Call rttov_direct and store results
  !   8. Deallocate all structures and arrays

  ! In this example we supply a profile of the optical parameters for each
  ! channel. Therefore a cloud coefficient file is not required.

  ! The optical parameters required in the input file are:
  !   - the list of angles for which the phase function is defined
  ! and for each input profile:
  !   - absorption coefficients for all channels for each layer
  !   - scattering coefficients for all channels for each layer
  !   - the phase function over all angles for each channel for each layer

  ! The phase functions are always required. They are used:
  ! - to calculate the "b parameters" required by the Chou-scaling parameterisation
  ! - to calculate the Legendre coefficients required by the DOM solver
  ! - explicitly for solar scattering calculations

  ! The calculations of the Legendre decomposition and also, in particular,
  ! the "b parameters" are relatively slow so you may want to consider calculating
  ! these inputs off-line and storing them in a file if possible.

  ! The additional steps required when running with optical parameters are:
  !   1. Allocate the optical parameter structure
  !   2. Read the optical parameters from the input file
  !   3. Calculate the b parameters from the phase functions (for Chou-scaling only)
  !   4. Calculate Legendre coefficients from the phase functions (for DOM only)
  !   5. Call rttov_init_opt_params to pre-calculate some phase angle data
  !      (for solar calculations only)

  ! If nthreads is greater than 1 the parallel RTTOV interface is used.
  ! To take advantage of multi-threaded execution you must have compiled
  ! RTTOV with openmp enabled. See the user guide and the compiler flags.

  errorstatus = 0_jpim

  !=====================================================
  !========== Interactive inputs == start ==============

  WRITE(0,*) 'enter path of coefficient file'
  READ(*,*) coef_filename
  WRITE(0,*) 'enter path of file containing profile data'
  READ(*,*) prof_filename
  WRITE(0,*) 'enter path of cloud optical parameter file'
  READ(*,*) file_cld_opt_param
  WRITE(0,*) 'enter number of angles for which phase function is defined'
  READ(*,*) nphangle
  WRITE(0,*) 'enter number of profiles'
  READ(*,*) nprof
  WRITE(0,*) 'enter number of profile levels'
  READ(*,*) nlevels
  WRITE(0,*) 'turn on solar simulations? (0=no, 1=yes)'
  READ(*,*) dosolar
  WRITE(0,*) 'enter number of channels to simulate per profile'
  READ(*,*) nchannels
  ALLOCATE(channel_list(nchannels))
  WRITE(0,*) 'enter space-separated channel list'
  READ(*,*,iostat=ios) channel_list(:)
  WRITE(0,*) 'enter number of threads to use'
  READ(*,*) nthreads
  IF (nthreads < 1) nthreads = 1


  ! --------------------------------------------------------------------------
  ! 1. Initialise RTTOV options structure
  ! --------------------------------------------------------------------------

  IF (dosolar == 1) THEN
    opts % rt_ir % addsolar = .TRUE.           ! Include solar radiation
  ELSE
    opts % rt_ir % addsolar = .FALSE.          ! Do not include solar radiation
  ENDIF
  opts % interpolation % addinterp   = .TRUE.  ! Allow interpolation of input profile
  opts % interpolation % interp_mode = 1       ! Set interpolation method
  opts % rt_all % addrefrac          = .TRUE.  ! Include refraction in path calc
  opts % rt_ir % addaerosl           = .FALSE. ! Don't include aerosol effects

  opts % rt_ir % addclouds           = .TRUE.  ! Include cloud effects
  opts % rt_ir % user_cld_opt_param  = .TRUE.  ! Supply optical parameters explictly

  opts % rt_ir % ir_scatt_model      = 2       ! Scattering model for emission source term:
                                               !   1 => DOM; 2 => Chou-scaling
  opts % rt_ir % vis_scatt_model     = 1       ! Scattering model for solar source term:
                                               !   1 => DOM; 2 => single-scattering
  opts % rt_ir % dom_nstreams        = 8       ! Number of streams for Discrete Ordinates (DOM)

  opts % rt_all % ozone_data         = .FALSE. ! Set the relevant flag to .TRUE.
  opts % rt_all % co2_data           = .FALSE. !   when supplying a profile of the
  opts % rt_all % n2o_data           = .FALSE. !   given trace gas (ensure the
  opts % rt_all % ch4_data           = .FALSE. !   coef file supports the gas)
  opts % rt_all % co_data            = .FALSE. !
  opts % rt_all % so2_data           = .FALSE. !
  opts % rt_mw % clw_data            = .FALSE. !

  opts % config % verbose            = .TRUE.  ! Enable printing of warnings

  !========== Interactive inputs == end ==============
  !===================================================


  ! --------------------------------------------------------------------------
  ! 2. Read coefficients
  ! --------------------------------------------------------------------------
  CALL rttov_read_coefs(errorstatus, coefs, opts, file_coef=coef_filename)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'fatal error reading coefficients'
    CALL rttov_exit(errorstatus)
  ENDIF

  ! Ensure input number of channels is not higher than number stored in coefficient file
  IF (nchannels > coefs % coef % fmv_chn) THEN
    nchannels = coefs % coef % fmv_chn
  ENDIF

  ! Ensure the options and coefficients are consistent
  CALL rttov_user_options_checkinput(errorstatus, opts, coefs)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'error in rttov options'
    CALL rttov_exit(errorstatus)
  ENDIF


  ! --------------------------------------------------------------------------
  ! 3. Allocate RTTOV input and output structures
  ! --------------------------------------------------------------------------

  ! Determine the total number of radiances to simulate (nchanprof).
  ! In this example we simulate all specified channels for each profile, but
  ! in general one can simulate a different number of channels for each profile.

  nchanprof = nchannels * nprof

  ! Allocate structures for rttov_direct
  CALL rttov_alloc_direct( &
        errorstatus,                         &
        1_jpim,                              &  ! 1 => allocate
        nprof,                               &
        nchanprof,                           &
        nlevels,                             &
        chanprof,                            &
        opts,                                &
        profiles,                            &
        coefs,                               &
        transmission,                        &
        radiance,                            &
        calcemis=calcemis,                   &
        emissivity=emissivity,               &
        calcrefl=calcrefl,                   &
        reflectance=reflectance,             &
        cld_maxnmom=opts%rt_ir%dom_nstreams, &
        cld_nphangle=nphangle,               &
        cld_opt_param=cld_opt_param,         &
        init=.TRUE._jplm)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'allocation error for rttov_direct structures'
    CALL rttov_exit(errorstatus)
  ENDIF


  ! --------------------------------------------------------------------------
  ! 4. Build the list of profile/channel indices in chanprof
  ! --------------------------------------------------------------------------

  nch = 0_jpim
  DO j = 1, nprof
    DO jch = 1, nchannels
      nch = nch + 1_jpim
      chanprof(nch)%prof = j
      chanprof(nch)%chan = channel_list(jch)
    ENDDO
  ENDDO


  ! --------------------------------------------------------------------------
  ! 5. Read profile data
  ! --------------------------------------------------------------------------

  !===============================================
  !========== Read profiles == start =============

  OPEN(iup, file=TRIM(prof_filename), status='old', iostat=ios)
  IF (ios /= 0) THEN
    WRITE(*,*) 'error opening profile file ios= ', ios
    CALL rttov_exit(errorstatus_fatal)
  ENDIF
  CALL rttov_skipcommentline(iup, errorstatus)

  ! Read gas units for profiles
  READ(iup,*) profiles(1) % gas_units
  profiles(:) % gas_units = profiles(1) % gas_units
  CALL rttov_skipcommentline(iup, errorstatus)

  ! Loop over all profiles and read data for each one
  DO iprof = 1, nprof

    ! Read pressure (hPa), temp (K), WV, O3 (gas units ppmv or kg/kg - as read above)
    READ(iup,*) profiles(iprof) % p(:)
    CALL rttov_skipcommentline(iup, errorstatus)
    READ(iup,*) profiles(iprof) % t(:)
    CALL rttov_skipcommentline(iup, errorstatus)
    READ(iup,*) profiles(iprof) % q(:)
    CALL rttov_skipcommentline(iup, errorstatus)
    ! Ozone profile is commented out in input profile data
!     READ(iup,*) profiles(iprof) % o3(:)
!     CALL rttov_skipcommentline(iup, errorstatus)

    ! 2 meter air variables
    READ(iup,*) profiles(iprof) % s2m % t, &
                profiles(iprof) % s2m % q, &
                profiles(iprof) % s2m % p, &
                profiles(iprof) % s2m % u, &
                profiles(iprof) % s2m % v, &
                profiles(iprof) % s2m % wfetc
    CALL rttov_skipcommentline(iup, errorstatus)

    ! Skin variables
    READ(iup,*) profiles(iprof) % skin % t
    CALL rttov_skipcommentline(iup, errorstatus)

    ! Surface type and water type
    READ(iup,*) profiles(iprof) % skin % surftype, &
                profiles(iprof) % skin % watertype
    CALL rttov_skipcommentline(iup, errorstatus)

    ! Elevation, latitude and longitude
    READ(iup,*) profiles(iprof) % elevation, &
                profiles(iprof) % latitude,  &
                profiles(iprof) % longitude
    CALL rttov_skipcommentline(iup, errorstatus)

    ! Satellite and solar angles
    READ(iup,*) profiles(iprof) % zenangle,    &
                profiles(iprof) % azangle,     &
                profiles(iprof) % sunzenangle, &
                profiles(iprof) % sunazangle
    CALL rttov_skipcommentline(iup, errorstatus)

    ! Cloud variables for simple cloud scheme, set cfraction to 0. to turn this off (VIS/IR only)
    ! (Ignored for full cloud scattering simulations)
    READ(iup,*) profiles(iprof) % ctp, &
                profiles(iprof) % cfraction
    CALL rttov_skipcommentline(iup, errorstatus)

  ENDDO
  CLOSE(iup)

  ! --------------------------------------------------------------------------
  ! Read the cloud optical parameter data for a single atmospheric profile
  ! --------------------------------------------------------------------------

  ! The parameters are supplied per channel so this file is specific to the instrument
  OPEN(iup, file=file_cld_opt_param, status='old', iostat=ios)
  IF (ios /= 0) THEN
    WRITE(*,*) 'error opening profile file ios= ', ios
    CALL rttov_exit(errorstatus_fatal)
  ENDIF
  CALL rttov_skipcommentline(iup, errorstatus)

  ! Read in the angles over which the phase function is defined
  READ(iup,*) cld_opt_param % phangle(:)
  CALL rttov_skipcommentline(iup, errorstatus)

  DO iprof = 1, nprof
    joff = nchannels * (iprof-1) + 1

    ! Read cloud fraction for each layer
    DO i = 1, nlevels-1_jpim
      READ(iup,*) profiles(iprof) % cfrac(i)
    ENDDO
    CALL rttov_skipcommentline(iup, errorstatus)

    ! Read absorption and scattering coefficients for each layer for each channel
    DO i = 1, nlevels-1_jpim
      READ(iup,*) cld_opt_param % abs(i,joff:joff+nchannels-1)
    ENDDO
    CALL rttov_skipcommentline(iup, errorstatus)

    DO i = 1, nlevels-1_jpim
      READ(iup,*) cld_opt_param % sca(i,joff:joff+nchannels-1)
    ENDDO
    CALL rttov_skipcommentline(iup, errorstatus)

    ! Read in phase function: each full phase function read for each layer for each channel
    DO i = 1, nlevels-1_jpim
      DO j = 1, nchannels
        READ(iup,*) cld_opt_param % pha(:,i,joff+j-1)
      ENDDO
    ENDDO
    CALL rttov_skipcommentline(iup, errorstatus)
  ENDDO
  CLOSE(iup)


  IF (opts % rt_ir % ir_scatt_model == 2) THEN
    ! For Chou-scaling the "b parameters" must be calculated from the phase functions:
    ! this could be done off-line and the b parameters could be stored to avoid repeating
    ! this calculation for the same phase functions as it is relatively expensive.
    ! NB if the phase function is the same for multiple layers it does not need
    ! re-calculating for every layer.

    WRITE(*,*) 'calculating b parameters...'

    ! Initialise the data for the bpr calculation
    CALL rttov_bpr_init(errorstatus, cld_opt_param % phangle(:))
    IF (errorstatus /= errorstatus_success) THEN
      WRITE(*,*) 'error calling rttov_bpr_init'
      CALL rttov_exit(errorstatus)
    ENDIF

    ! Calculate the b parameterfor each scattering layer for each channel
    DO iprof = 1, nprof
      DO j = 1, nchannels
        joff = nchannels * (iprof-1)
        DO i = 1, nlevels-1_jpim

          ! Only do calculation for layers containing cloud
          IF (ANY(cld_opt_param % pha(:,i,j+joff) > 0._jprb)) THEN

            ! If this phase function is same as previous layer copy the bpr
            IF (i > 1 .AND. ALL(cld_opt_param % pha(:,i,j+joff) == cld_opt_param % pha(:,MAX(i-1,1),j+joff))) THEN
              cld_opt_param % bpr(i,j+joff) = cld_opt_param % bpr(i-1,j+joff)
            ELSE
              CALL rttov_bpr_calc(errorstatus, &
                                  cld_opt_param % pha(:,i,j+joff), &! input phase function
                                  cld_opt_param % phangle(:),      &! input phase angles
                                  cld_opt_param % bpr(i,j+joff),   &! output b parameter
                                  nthreads)
              IF (errorstatus /= errorstatus_success) THEN
                WRITE(*,*) 'error calculating b parameters'
                CALL rttov_exit(errorstatus)
              ENDIF
            ENDIF

          ELSE
            cld_opt_param % bpr(i,j+joff) = 0._jprb
          ENDIF

        ENDDO
      ENDDO
    ENDDO

    ! Deallocate the data for the bpr calculation
    CALL rttov_bpr_dealloc(errorstatus)
    IF (errorstatus /= errorstatus_success) THEN
      WRITE(*,*) 'error calling rttov_bpr_dealloc'
    ENDIF

    WRITE(*,*) 'done calculating b parameters'
  ENDIF

  IF ( opts % rt_ir % ir_scatt_model  == 1 .OR. &
      (opts % rt_ir % vis_scatt_model == 1 .AND. opts % rt_ir % addsolar)) THEN

    ! For DOM simulations the first dom_nstreams Legendre coefficients must be
    ! calculated for each phase function.

    WRITE(*,*) 'calculating Legendre coefficients...'

    ! Calculate the Legendre coefficients for each scattering layer for each channel
    DO iprof = 1, nprof
      DO j = 1, nchannels
        joff = nchannels * (iprof-1)
        DO i = 1, nlevels-1_jpim

          ! Only do calculation for layers containing cloud
          IF (ANY(cld_opt_param % pha(:,i,j+joff) > 0._jprb)) THEN

            ! If this phase function is same as previous layer copy the coefficients
            IF (i > 1 .AND. ALL(cld_opt_param % pha(:,i,j+joff) == cld_opt_param % pha(:,MAX(i-1,1),j+joff))) THEN
              cld_opt_param % legcoef(:,i,j+joff) = cld_opt_param % legcoef(:,i-1,j+joff)
            ELSE
              CALL rttov_legcoef_calc(errorstatus, &
                                  cld_opt_param % pha(:,i,j+joff),   &! input phase function
                                  cld_opt_param % phangle(:),        &! input phase angles
                                  opts % rt_ir % dom_nstreams,       &! number of Legendre coefficients to calculate
                                  cld_opt_param % legcoef(:,i,j+joff))! output Legendre coefficients
              IF (errorstatus /= errorstatus_success) THEN
                WRITE(*,*) 'error calculating Legendre coefficients'
                CALL rttov_exit(errorstatus)
              ENDIF
            ENDIF

          ELSE
            cld_opt_param % legcoef(:,i,j+joff) = 0._jprb
          ENDIF

        ENDDO
      ENDDO
    ENDDO

    WRITE(*,*) 'done calculating Legendre coefficients'
  ENDIF

  ! If doing solar calculations pre-calculate some phase angle data for scattering calculations
  ! This needs to be re-run if the phase *angles* change (but not if the phase *functions* change)
  IF (opts % rt_ir % addsolar) THEN
    CALL rttov_init_opt_param(errorstatus, opts, cld_opt_param)
    IF (errorstatus /= errorstatus_success) THEN
      WRITE(*,*) 'error initialising cld_opt_param'
      CALL rttov_exit(errorstatus)
    ENDIF
  ENDIF

  !========== Read profiles == end =============
  !=============================================


  ! --------------------------------------------------------------------------
  ! 6. Specify surface emissivity and reflectance
  ! --------------------------------------------------------------------------

  ! In this example we have no values for input emissivities or reflectances
  ! so we initialise all inputs to zero
  CALL rttov_init_emis_refl(emissivity, reflectance)

  ! Calculate emissivity within RTTOV where the input emissivity value is
  ! zero or less (all channels in this case)
  calcemis(:) = (emissivity(:) % emis_in <= 0._jprb)

  ! Calculate reflectances within RTTOV where the input BRDF value is zero or
  ! less (all channels in this case)
  calcrefl(:) = (reflectance(:) % refl_in <= 0._jprb)


  ! --------------------------------------------------------------------------
  ! 7. Call RTTOV forward model
  ! --------------------------------------------------------------------------
  IF (nthreads <= 1) THEN
    CALL rttov_direct(                   &
            errorstatus,                 &! out   error flag
            chanprof,                    &! in    channel and profile index structure
            opts,                        &! in    options structure
            profiles,                    &! in    profile array
            coefs,                       &! in    coefficients structure
            transmission,                &! inout computed transmittances
            radiance,                    &! inout computed radiances
            calcemis      = calcemis,    &! in    flag for internal emissivity calcs
            emissivity    = emissivity,  &! inout input/output emissivities per channel
            calcrefl      = calcrefl,    &! in    flag for internal BRDF calcs
            reflectance   = reflectance, &! inout input/output BRDFs per channel
            cld_opt_param = cld_opt_param)! in    cloud optical parameters
  ELSE
    CALL rttov_parallel_direct(           &
            errorstatus,                  &! out   error flag
            chanprof,                     &! in    channel and profile index structure
            opts,                         &! in    options structure
            profiles,                     &! in    profile array
            coefs,                        &! in    coefficients structure
            transmission,                 &! inout computed transmittances
            radiance,                     &! inout computed radiances
            calcemis      = calcemis,     &! in    flag for internal emissivity calcs
            emissivity    = emissivity,   &! inout input/output emissivities per channel
            calcrefl      = calcrefl,     &! in    flag for internal BRDF calcs
            reflectance   = reflectance,  &! inout input/output BRDFs per channel
            cld_opt_param = cld_opt_param,&! in    cloud optical parameters
            nthreads      = nthreads)      ! in    number of threads to use
  ENDIF

  IF (errorstatus /= errorstatus_success) THEN
    WRITE (*,*) 'rttov_direct error'
    CALL rttov_exit(errorstatus)
  ENDIF


  !=====================================================
  !============== Output results == start ==============

  ! Open output file where results are written
  OPEN(ioout, file='output_'//NameOfRoutine//'.dat', status='unknown', form='formatted', iostat=ios)
  IF (ios /= 0) THEN
    WRITE(*,*) 'error opening the output file ios= ', ios
    CALL rttov_exit(errorstatus_fatal)
  ENDIF

  WRITE(ioout,*)' -----------------'
  WRITE(ioout,*)' Instrument ', inst_name(coefs % coef % id_inst)
  WRITE(ioout,*)' -----------------'
  WRITE(ioout,*)' '
  CALL rttov_print_opts(opts, lu=ioout)

  DO iprof = 1, nprof

    joff = (iprof-1_jpim) * nchannels

    WRITE(ioout,*)' '
    WRITE(ioout,*)' Profile ', iprof

    CALL rttov_print_profile(profiles(iprof), lu=ioout)

    WRITE(ioout,777)'CHANNELS PROCESSED FOR SAT ', platform_name(coefs % coef % id_platform), coefs % coef % id_sat
    WRITE(ioout,111) (chanprof(j) % chan, j = 1+joff, nchannels+joff)
    WRITE(ioout,*)' '
    WRITE(ioout,*)'CALCULATED BRIGHTNESS TEMPERATURES (K):'
    WRITE(ioout,222) (radiance % bt(j), j = 1+joff, nchannels+joff)
    IF (opts % rt_ir % addsolar) THEN
      WRITE(ioout,*)' '
      WRITE(ioout,*)'CALCULATED SATELLITE REFLECTANCES (BRF):'
      WRITE(ioout,444) (radiance % refl(j), j = 1+joff, nchannels+joff)
    ENDIF
    WRITE(ioout,*)' '
    WRITE(ioout,*)'CALCULATED RADIANCES (mW/m2/sr/cm-1):'
    WRITE(ioout,222) (radiance % total(j), j = 1+joff, nchannels+joff)
    WRITE(ioout,*)' '
    WRITE(ioout,*)'CALCULATED SURFACE EMISSIVITIES:'
    WRITE(ioout,444) (emissivity(j) % emis_out, j = 1+joff, nchannels+joff)
    IF (opts % rt_ir % addsolar) THEN
      WRITE(ioout,*)' '
      WRITE(ioout,*)'CALCULATED SURFACE BRDF:'
      WRITE(ioout,444) (reflectance(j) % refl_out, j = 1+joff, nchannels+joff)
    ENDIF
  ENDDO

  ! Close output file
  CLOSE(ioout, iostat=ios)
  IF (ios /= 0) THEN
    WRITE(*,*) 'error closing the output file ios= ', ios
    CALL rttov_exit(errorstatus_fatal)
  ENDIF

  !============== Output results == end ==============
  !=====================================================


  ! --------------------------------------------------------------------------
  ! 8. Deallocate all RTTOV arrays and structures
  ! --------------------------------------------------------------------------
  DEALLOCATE (channel_list, stat=alloc_status)
  IF (alloc_status /= 0) THEN
    WRITE(*,*) 'mem dellocation error'
  ENDIF

  CALL rttov_alloc_direct( &
        errorstatus,                         &
        0_jpim,                              &  ! 0 => deallocate
        nprof,                               &
        nchanprof,                           &
        nlevels,                             &
        chanprof,                            &
        opts,                                &
        profiles,                            &
        coefs,                               &
        transmission,                        &
        radiance,                            &
        calcemis=calcemis,                   &
        emissivity=emissivity,               &
        calcrefl=calcrefl,                   &
        reflectance=reflectance,             &
        cld_maxnmom=opts%rt_ir%dom_nstreams, &
        cld_nphangle=nphangle,               &
        cld_opt_param=cld_opt_param)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'deallocation error for rttov_direct structures'
    CALL rttov_exit(errorstatus)
  ENDIF

  CALL rttov_dealloc_coefs(errorstatus, coefs)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'coefs deallocation error'
  ENDIF


! Format definitions for output
111  FORMAT(1X,10I8)
222  FORMAT(1X,10F8.2)
444  FORMAT(1X,10F8.3)
777  FORMAT(/,A,A9,I3)

END PROGRAM example_cld_param_fwd
