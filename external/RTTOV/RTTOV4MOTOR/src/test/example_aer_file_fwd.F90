! Description:
!> @file
!!   Example calling direct model visible/IR aerosol simulation using
!!   pre-defined particle optical properties.
!
!> @brief
!!   Example calling direct model visible/IR aerosol simulation using
!!   pre-defined particle optical properties.
!!
!! @details
!!   This example is most easily run using the run_example_aer_file_fwd.sh
!!   script (see the user guide).
!!
!!   This program requires the following files:
!!     the file containing input profiles (e.g. prof.dat)
!!     the file containing input aerosol data
!!     the RTTOV optical depth coefficient file
!!     the RTTOV aerosol coefficient file
!!
!!   The output is written to a file called example_aer_file_fwd_output.dat
!!
!!   This program is very similar to example_fwd.F90 with the addition
!!   of the aerosol scattering.
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
PROGRAM example_aer_file_fwd

  ! rttov_const contains useful RTTOV constants
  USE rttov_const, ONLY :     &
         errorstatus_success, &
         errorstatus_fatal,   &
         platform_name,       &
         inst_name,           &
         naer_opac

  ! rttov_types contains definitions of all RTTOV data types
  USE rttov_types, ONLY :     &
         rttov_options,       &
         rttov_coefs,         &
         rttov_profile,       &
         rttov_transmission,  &
         rttov_radiance,      &
         rttov_chanprof,      &
         rttov_emissivity,    &
         rttov_reflectance

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
#include "rttov_user_options_checkinput.interface"
#include "rttov_print_opts.interface"
#include "rttov_print_profile.interface"
#include "rttov_skipcommentline.interface"
#include "rttov_aer_clim_prof.interface"

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
  TYPE(rttov_transmission)         :: transmission             ! Output transmittances
  TYPE(rttov_radiance)             :: radiance                 ! Output radiances

  INTEGER(KIND=jpim)               :: errorstatus              ! Return error status of RTTOV subroutine calls

  INTEGER(KIND=jpim) :: alloc_status
  CHARACTER(LEN=20)  :: NameOfRoutine = 'example_aer_file_fwd'

  ! variables for input
  !====================
  CHARACTER(LEN=256) :: coef_filename
  CHARACTER(LEN=256) :: aer_coef_filename
  CHARACTER(LEN=256) :: prof_filename
  CHARACTER(LEN=256) :: aer_prof_filename
  INTEGER(KIND=jpim) :: nthreads
  INTEGER(KIND=jpim) :: dosolar
  INTEGER(KIND=jpim) :: nlevels
  INTEGER(KIND=jpim) :: nprof
  INTEGER(KIND=jpim) :: nchannels
  INTEGER(KIND=jpim) :: nchanprof
  INTEGER(KIND=jpim), ALLOCATABLE :: channel_list(:)
  INTEGER(KIND=jpim) :: aer_clim_index
  REAL(KIND=jprb), ALLOCATABLE :: aerprof(:,:,:)
  ! loop variables
  INTEGER(KIND=jpim) :: j, jch, lay
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

  ! If nthreads is greater than 1 the parallel RTTOV interface is used.
  ! To take advantage of multi-threaded execution you must have compiled
  ! RTTOV with openmp enabled. See the user guide and the compiler flags.

  errorstatus = 0_jpim

  !=====================================================
  !========== Interactive inputs == start ==============

  WRITE(0,*) 'enter path of coefficient file'
  READ(*,*) coef_filename
  WRITE(0,*) 'enter path of aerosol coefficient file'
  READ(*,*) aer_coef_filename
  WRITE(0,*) 'enter path of file containing (clear-sky) profile data'
  READ(*,*) prof_filename
  WRITE(0,*) 'enter path of file containing aerosol profile data'
  READ(*,*) aer_prof_filename
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
  opts % rt_ir % addclouds           = .FALSE. ! Don't include cloud effects

  opts % rt_ir % addaerosl           = .TRUE.  ! Include aerosol effects

  opts % rt_ir % ir_scatt_model      = 2       ! Scattering model for emission source term:
                                               !   1 => DOM; 2 => Chou-scaling
  opts % rt_ir % vis_scatt_model     = 1       ! Scattering model for solar source term:
                                               !   1 => DOM; 2 => single-scattering; 3 => MFASIS
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

  ! Read optical depth and aerosol coefficient files together
  CALL rttov_read_coefs(errorstatus, coefs, opts, &
                        file_coef=coef_filename, file_scaer=aer_coef_filename)
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
        errorstatus,             &
        1_jpim,                  &  ! 1 => allocate
        nprof,                   &
        nchanprof,               &
        nlevels,                 &
        chanprof,                &
        opts,                    &
        profiles,                &
        coefs,                   &
        transmission,            &
        radiance,                &
        calcemis=calcemis,       &
        emissivity=emissivity,   &
        calcrefl=calcrefl,       &
        reflectance=reflectance, &
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
    ! (Ignored for aerosol simulations)
    READ(iup,*) profiles(iprof) % ctp, &
                profiles(iprof) % cfraction
    CALL rttov_skipcommentline(iup, errorstatus)

  ENDDO
  CLOSE(iup)


  ! -----------------------------------
  ! Specify the aerosol input profiles
  ! -----------------------------------

  ! The aerosol coefficient files contain optical parameters for 13 pre-defined aerosol types.
  ! Aerosol simulations are performed by specifying number density profiles for one or more
  ! of these particle types in the profile%aerosols(:,:) array.

  ! RTTOV has a subroutine rttov_aer_clim_prof which can be used to generate climatological
  ! aerosol profiles (see user guide for more details).

  ! The aerosol profile input file can either provide an index indicating the climatological
  ! aerosol profile type (integer from 1-10) or otherwise it can provide number density profiles.

  OPEN(iup, file=TRIM(aer_prof_filename), status='old', iostat=ios)
  IF (ios /= 0) THEN
    WRITE(*,*) 'error opening profile file ios= ', ios
    CALL rttov_exit(errorstatus_fatal)
  ENDIF
  CALL rttov_skipcommentline(iup, errorstatus)

  ! Read flag to indicate aerosol units (T => kg/kg; F => number density cm^-3)
  READ(iup,*) profiles(1)%mmr_cldaer
  profiles(:)%mmr_cldaer = profiles(1)%mmr_cldaer
  CALL rttov_skipcommentline(iup, errorstatus)

  ALLOCATE(aerprof(nlevels-1,10,naer_opac))

  ! Loop over all profiles and read data for each one
  DO iprof = 1, nprof

    ! Read aer_clim_index: if the value is 1-10, use rttov_aer_clim_prof.F90;
    ! otherwise, read aerosol profiles from the input file.

    READ(iup,*) aer_clim_index
    CALL rttov_skipcommentline(iup, errorstatus)

    IF (aer_clim_index >= 1 .AND. aer_clim_index <= 10) THEN

      ! Use climatological aerosol profile

      CALL rttov_aer_clim_prof(       &
        profiles(iprof) % p,          &
        profiles(iprof) % t,          &
        profiles(iprof) % q,          &
        profiles(iprof) % gas_units,  &
        profiles(iprof) % mmr_cldaer, &
        nlevels,                      & ! Profile level corresponding to surface
        profiles(iprof) % latitude,   &
        profiles(iprof) % elevation,  &
        1.0_jprb,                     & ! Scale factor for output profiles
        aerprof)                        ! Output profiles

      DO lay = 1, nlevels-1
        profiles(iprof) % aerosols(:,lay) = aerprof(lay,aer_clim_index,:)
      ENDDO

    ELSE

      ! Read aerosol number density profiles
      READ(iup,*) profiles(iprof) % aerosols(:,:)
      CALL rttov_skipcommentline(iup, errorstatus)

    ENDIF

  ENDDO
  CLOSE(iup)

  DEALLOCATE(aerprof)

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
    CALL rttov_direct(                &
            errorstatus,              &! out   error flag
            chanprof,                 &! in    channel and profile index structure
            opts,                     &! in    options structure
            profiles,                 &! in    profile array
            coefs,                    &! in    coefficients structure
            transmission,             &! inout computed transmittances
            radiance,                 &! inout computed radiances
            calcemis    = calcemis,   &! in    flag for internal emissivity calcs
            emissivity  = emissivity, &! inout input/output emissivities per channel
            calcrefl    = calcrefl,   &! in    flag for internal BRDF calcs
            reflectance = reflectance) ! inout input/output BRDFs per channel
  ELSE
    CALL rttov_parallel_direct(       &
            errorstatus,              &! out   error flag
            chanprof,                 &! in    channel and profile index structure
            opts,                     &! in    options structure
            profiles,                 &! in    profile array
            coefs,                    &! in    coefficients structure
            transmission,             &! inout computed transmittances
            radiance,                 &! inout computed radiances
            calcemis    = calcemis,   &! in    flag for internal emissivity calcs
            emissivity  = emissivity, &! inout input/output emissivities per channel
            calcrefl    = calcrefl,   &! in    flag for internal BRDF calcs
            reflectance = reflectance,&! inout input/output BRDFs per channel
            nthreads    = nthreads)    ! in    number of threads to use
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

  ! Deallocate structures for rttov_direct
  CALL rttov_alloc_direct( &
        errorstatus,             &
        0_jpim,                  &  ! 0 => deallocate
        nprof,                   &
        nchanprof,               &
        nlevels,                 &
        chanprof,                &
        opts,                    &
        profiles,                &
        coefs,                   &
        transmission,            &
        radiance,                &
        calcemis=calcemis,       &
        emissivity=emissivity,   &
        calcrefl=calcrefl,       &
        reflectance=reflectance)
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

END PROGRAM example_aer_file_fwd
