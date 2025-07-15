! Description:
!> @file
!!   Example calling PC-RTTOV direct model clear-sky simulation.
!
!> @brief
!!   Example calling PC-RTTOV direct model clear-sky simulation.
!!
!! @details
!!   This example is most easily run using the run_example_pc_fwd.sh
!!   script (see the user guide).
!!
!!   This program requires the following files:
!!     the file containing input profiles (e.g. prof.dat)
!!     the file containing the list of channels for which
!!       to reconstruct radiances (e.g. radrec.dat) - optional
!!     the RTTOV "rtcoef" and "pccoef" coefficient files
!!
!!   The output is written to a file called example_pc_fwd_output.dat
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
PROGRAM example_pc_fwd

  ! rttov_const contains useful RTTOV constants
  USE rttov_const, ONLY :     &
         errorstatus_success, &
         errorstatus_fatal,   &
         inst_name,           &
         surftype_sea

  ! rttov_types contains definitions of all RTTOV data types
  USE rttov_types, ONLY :     &
         rttov_options,       &
         rttov_coefs,         &
         rttov_profile,       &
         rttov_transmission,  &
         rttov_radiance,      &
         rttov_chanprof,      &
         rttov_emissivity,    &
         rttov_pccomp

  ! jpim, jprb and jplm are the RTTOV integer, real and logical KINDs
  USE parkind1, ONLY : jpim, jplm !, jprb (jprb currently unused here)

  USE rttov_unix_env, ONLY : rttov_exit

  IMPLICIT NONE

#include "rttov_direct.interface"
#include "rttov_parallel_direct.interface"
#include "rttov_read_coefs.interface"
#include "rttov_dealloc_coefs.interface"
#include "rttov_alloc_direct.interface"
#include "rttov_init_emis_refl.interface"
#include "rttov_user_options_checkinput.interface"
#include "rttov_get_pc_predictindex.interface"
#include "rttov_print_opts.interface"
#include "rttov_print_profile.interface"
#include "rttov_skipcommentline.interface"

  !--------------------------
  !
  INTEGER(KIND=jpim), PARAMETER :: iup   = 20   ! unit for input profile file
  INTEGER(KIND=jpim), PARAMETER :: ioout = 21   ! unit for output

  ! RTTOV variables/structures
  !====================
  TYPE(rttov_options)             :: opts                      ! Options structure
  TYPE(rttov_coefs)               :: coefs                     ! Coefficients structure
  TYPE(rttov_chanprof),   POINTER :: chanprof(:)    => NULL()  ! Input channel/profile list
  LOGICAL(KIND=jplm),     POINTER :: calcemis(:)    => NULL()  ! Flag to indicate calculation of emissivity within RTTOV
  TYPE(rttov_emissivity), POINTER :: emissivity(:)  => NULL()  ! Input/output surface emissivity
  TYPE(rttov_profile),    POINTER :: profiles(:)    => NULL()  ! Input profiles
  TYPE(rttov_transmission)        :: transmission              ! Output transmittances
  TYPE(rttov_radiance)            :: radiance                  ! Output radiances
  TYPE(rttov_pccomp)              :: pccomp                    ! Output PC structure
  INTEGER(KIND=jpim),     POINTER :: channels_rec(:) => NULL() ! Reconstructed radiance channel list

  INTEGER(KIND=jpim)              :: errorstatus               ! Return error status of RTTOV subroutine calls

  INTEGER(KIND=jpim) :: alloc_status
  CHARACTER(LEN=14)  :: NameOfRoutine = 'example_pc_fwd'

  ! variables for input
  !====================
  CHARACTER(LEN=256) :: coef_filename, pccoef_filename
  CHARACTER(LEN=256) :: prof_filename, radrec_filename
  INTEGER(KIND=jpim) :: nthreads
  INTEGER(KIND=jpim) :: nlevels
  INTEGER(KIND=jpim) :: nprof
  INTEGER(KIND=jpim) :: ipcbnd, ipcreg, npcscores
  INTEGER(KIND=jpim), POINTER :: predictindex(:)
  INTEGER(KIND=jpim) :: nchannels
  INTEGER(KIND=jpim) :: nchannels_rec
  INTEGER(KIND=jpim) :: nchanprof
  LOGICAL(KIND=jplm) :: exists
  ! loop variables
  INTEGER(KIND=jpim) :: j, lo, hi
  INTEGER(KIND=jpim) :: iprof
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

  ! In addition the following considerations apply to PC-RTTOV:
  !   - channels for PC predictors are stored in chanprof and are fixed
  !     according to the instrument and number of predictors
  !   - number of predictors is a function of the regression set (ipcbnd, ipcreg)
  !   - RTTOV optical depth coefficient file should be compatible with PCs
  !     (requires v9 predictors)
  !   - some PC coefficients can be used over all surface types - filenames include "landsea";
  !     otherwise they can only be used for profiles over sea
  !   - sea surface emissivity must be calculated within RTTOV (calcemis TRUE), land surface
  !     emissivity should be provided by the UW IR emissivity atlas
  !   - some PC coefficients optionally allow the NLTE bias correction - filenames include "nlte";
  !     otherwise they cannot be used with the NLTE correction
  !   - some PC coefficients allow all variable trace gases (except SO2) - filenames include "trace";
  !     otherwise the only variable trace gases allowed are water vapour and ozone
  !   - some PC coefficients allow aerosol simulations - filenames include "aer";
  !     otherwise the simulations must be clear-sky
  ! See user guide for more details on these points.

  ! If nthreads is greater than 1 the parallel RTTOV interface is used.
  ! To take advantage of multi-threaded execution you must have compiled
  ! RTTOV with openmp enabled. See the user guide and the compiler flags.

  errorstatus = 0

  !=====================================================
  !========== Interactive inputs == start ==============

  WRITE(0,*) 'enter path of rtcoef coefficient file'
  READ(*,*) coef_filename
  WRITE(0,*) 'enter path of PC coefficient file'
  READ(*,*) pccoef_filename
  WRITE(0,*) 'enter path of file containing profile data'
  READ(*,*) prof_filename
  WRITE(0,*) 'enter path of file containing reconstructed channel list'
  READ(*,*) radrec_filename
  WRITE(0,*) 'enter number of profiles'
  READ(*,*) nprof
  WRITE(0,*) 'enter number of profile levels'
  READ(*,*) nlevels
  WRITE(0,*) 'enter ipcbnd (AIRS: 1, IASI: 1-3)'
  READ(*,*) ipcbnd
  WRITE(0,*) 'enter ipcreg (1-4)'
  READ(*,*) ipcreg
  WRITE(0,*) 'enter number of pcscores (1-400)'
  READ(*,*) npcscores
  WRITE(0,*) 'enter number of threads to use'
  READ(*,*) nthreads


  ! --------------------------------------------------------------------------
  ! 1. Initialise RTTOV options structure
  ! --------------------------------------------------------------------------

  opts % rt_ir % pc % addpc          = .TRUE.  ! Use PC-RTTOV
  opts % rt_ir % pc % ipcbnd         = ipcbnd
  opts % rt_ir % pc % ipcreg         = ipcreg
  opts % rt_ir % pc % npcscores      = npcscores

  ! In this example we reconstruct radiances if there is an input file
  ! containing a channel list
  INQUIRE(file=radrec_filename, exist=exists)
  opts % rt_ir % pc % addradrec      = exists

  opts % interpolation % addinterp   = .TRUE.  ! Allow interpolation of input profile
  opts % interpolation % interp_mode = 1       ! Set interpolation method
  opts % rt_all % addrefrac          = .TRUE.  ! Include refraction in path calc (always for PC)
  opts % rt_ir % addclouds           = .FALSE. ! Don't include cloud effects
  opts % rt_ir % addaerosl           = .FALSE. ! Don't include aerosol effects
  opts % rt_ir % addsolar            = .FALSE. ! Do not include solar radiation

  opts % rt_all % ozone_data         = .TRUE.  ! We have an ozone profile
  opts % rt_all % co2_data           = .FALSE. ! We do not have profiles
  opts % rt_all % n2o_data           = .FALSE. !   for any other constituents
  opts % rt_all % ch4_data           = .FALSE. !
  opts % rt_all % co_data            = .FALSE. !
  opts % rt_all % so2_data           = .FALSE. !
  opts % rt_mw % clw_data            = .FALSE. !

  opts % config % verbose            = .TRUE.  ! Enable printing of warnings

  !========== Interactive inputs == end ==============
  !===================================================


  ! --------------------------------------------------------------------------
  ! 2. Read coefficients
  ! --------------------------------------------------------------------------
  CALL rttov_read_coefs(errorstatus, coefs, opts, &
                        file_coef=coef_filename, file_pccoef=pccoef_filename)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'fatal error reading coefficients'
    CALL rttov_exit(errorstatus)
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

  ! Determine the size of the chanprof array (the number of PC-RTTOV predictor channels).
  ! This is determined by ipcbnd and ipcreg and is the same for all simulated profiles.
  ! The channel list in chanprof is the set of predictor channels for the PC scheme.
  ! The channel list may be obtained using the rttov_get_pc_predictindex subroutine:

  NULLIFY(predictindex)
  CALL rttov_get_pc_predictindex(errorstatus, opts, predictindex, file_pccoef=pccoef_filename)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'rttov_get_pc_predictindex fatal error'
    CALL rttov_exit(errorstatus)
  ENDIF
  nchannels = SIZE(predictindex)
  nchanprof = nchannels * nprof  ! Size of chanprof array is total number of channels over all profiles


  ! Determine the number of reconstructed radiances per profile (nchannels_rec)
  nchannels_rec = 0
  IF (opts % rt_ir % pc % addradrec) THEN

    ! Read list of channels for which to reconstruct radiances
    OPEN(iup, file=TRIM(radrec_filename), status='old', iostat=ios)
    IF (ios /= 0) THEN
      WRITE(*,*) 'error opening radrec file ios= ', ios
      CALL rttov_exit(errorstatus_fatal)
    ENDIF
    CALL rttov_skipcommentline(iup, errorstatus)

    READ(iup,*) nchannels_rec
    CALL rttov_skipcommentline(iup, errorstatus)


    IF (nchannels_rec < 0) THEN
      ! If the number of channels is negative, don't reconstruct radiances at all

      opts % rt_ir % pc % addradrec = .FALSE.

    ELSE IF (nchannels_rec == 0) THEN
      ! If the number of channels is set to 0 then reconstruct all instrument channels

      nchannels_rec = coefs % coef % fmv_chn
      ALLOCATE(channels_rec(nchannels_rec))
      channels_rec(:) = (/ (j, j = 1, nchannels_rec) /)

    ELSE
      ! Otherwise read the channel list from the file

      ALLOCATE(channels_rec(nchannels_rec))
      READ(iup,*) channels_rec(:)

    ENDIF

    CLOSE(iup)
  ENDIF

  ! Ensure we don't have unassociated pointers below when addradrec is FALSE
  IF (nchannels_rec <= 0) ALLOCATE(channels_rec(0))


  ! Allocate structures for rttov_direct
  CALL rttov_alloc_direct(                   &
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
        npcscores=npcscores * nprof,         &
        nchannels_rec=nchannels_rec * nprof, &
        pccomp=pccomp,                       &
        init=.TRUE._jplm)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'allocation error for rttov_direct structures'
    CALL rttov_exit(errorstatus)
  ENDIF


  ! --------------------------------------------------------------------------
  ! 4. Build the list of profile/channel indices in chanprof
  ! --------------------------------------------------------------------------

  ! Populate chanprof using the channel list obtained above in predictindex(:)
  DO j = 1, nprof
    lo = (j - 1) * nchannels + 1
    hi = lo + nchannels - 1
    chanprof(lo:hi)%prof = j
    chanprof(lo:hi)%chan = predictindex(:)
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
    READ(iup,*) profiles(iprof) % o3(:)
    CALL rttov_skipcommentline(iup, errorstatus)

    ! 2 meter air variables
    READ(iup,*) profiles(iprof) % s2m % t, &
                profiles(iprof) % s2m % q, &
                profiles(iprof) % s2m % p, &
                profiles(iprof) % s2m % u, &
                profiles(iprof) % s2m % v
    CALL rttov_skipcommentline(iup, errorstatus)

    ! Skin temperature
    READ(iup,*) profiles(iprof) % skin % t
    CALL rttov_skipcommentline(iup, errorstatus)

    ! Elevation and latitude
    READ(iup,*) profiles(iprof) % elevation, &
                profiles(iprof) % latitude
    CALL rttov_skipcommentline(iup, errorstatus)

    ! Satellite zenith angle
    READ(iup,*) profiles(iprof) % zenangle
    CALL rttov_skipcommentline(iup, errorstatus)

    ! Set surface type (for land surfaces the UW IR emissivity
    ! atlas is recommended - see example_atlas_fwd.F90)
    profiles(iprof) % skin % surftype = surftype_sea

  ENDDO
  CLOSE(iup)

  !========== Read profiles == end =============
  !=============================================


  ! --------------------------------------------------------------------------
  ! 6. Specify surface emissivity
  ! --------------------------------------------------------------------------

  ! For PC-RTTOV calcemis must be true for sea surfaces
  CALL rttov_init_emis_refl(emissivity)
  calcemis(:) = .TRUE.


  ! --------------------------------------------------------------------------
  ! 7. Call RTTOV forward model
  ! --------------------------------------------------------------------------
  IF (nthreads <= 1) THEN
    CALL rttov_direct(                  &
            errorstatus,                &! out   error flag
            chanprof,                   &! in    channel and profile index structure
            opts,                       &! in    options structure
            profiles,                   &! in    profile array
            coefs,                      &! in    coefficients structure
            transmission,               &! inout computed transmittances for predictors
            radiance,                   &! inout computed radiances for predictors
            calcemis     = calcemis,    &! in    flag for internal emissivity calcs
            emissivity   = emissivity,  &! inout input/output emissivities per channel
            pccomp       = pccomp,      &! inout computed PC scores
            channels_rec = channels_rec) ! in    reconstructed channel list
  ELSE
    CALL rttov_parallel_direct(         &
            errorstatus,                &! out   error flag
            chanprof,                   &! in    channel and profile index structure
            opts,                       &! in    options structure
            profiles,                   &! in    profile array
            coefs,                      &! in    coefficients structure
            transmission,               &! inout computed transmittances for predictors
            radiance,                   &! inout computed radiances for predictors
            calcemis     = calcemis,    &! in    flag for internal emissivity calcs
            emissivity   = emissivity,  &! inout input/output emissivities per channel
            pccomp       = pccomp,      &! inout computed PC scores
            channels_rec = channels_rec,&! in    reconstructed channel list
            nthreads     = nthreads)     ! in    number of threads to use
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

    WRITE(ioout,*)' '
    WRITE(ioout,*)' Profile ', iprof

    CALL rttov_print_profile(profiles(iprof), lu=ioout)

    WRITE(ioout,*)' '
    WRITE(ioout,*)' PC Scores '
    lo = 1 + (iprof-1) * npcscores
    hi = iprof * npcscores

    WRITE(ioout, '(A," = (")' ) 'PCCOMP%TOTAL_PCSCORES'
    WRITE(ioout, '(10(G15.6E3))') pccomp%total_pcscores(lo:hi)
    WRITE(ioout,*) ')'

    IF ( opts % rt_ir % pc % addradrec ) THEN
      lo = 1 + (iprof-1) * nchannels_rec
      hi = iprof * nchannels_rec
      WRITE(ioout,*)
      WRITE(ioout,*) 'Radiances reconstructed using principal components'
      WRITE(ioout, '(A," = (")' ) 'PCCOMP%TOTAL_PCCOMP'
      WRITE(ioout, '(10(G15.6E3))') pccomp % total_pccomp (lo:hi)
      WRITE(ioout,*) ')'

      WRITE(ioout,*)' '
      WRITE(ioout,*) 'Brightness temp equivalent to radiances reconstructed using PCs'
      WRITE(ioout, '(A," = (")' ) 'PCCOMP%BT_PCCOMP'
      WRITE(ioout, '(10(G15.6E3))') pccomp % bt_pccomp (lo:hi)
      WRITE(ioout,*) ')'
      WRITE(ioout,*)
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
  IF (ASSOCIATED(predictindex)) DEALLOCATE (predictindex, stat=alloc_status)
  IF (alloc_status /= 0) THEN
    WRITE(*,*) 'mem dellocation error'
  ENDIF

  IF (ASSOCIATED(channels_rec)) DEALLOCATE (channels_rec, stat=alloc_status)
  IF (alloc_status /= 0) THEN
    WRITE(*,*) 'mem dellocation error'
  ENDIF

  ! Deallocate structures for rttov_direct
  CALL rttov_alloc_direct(                   &
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
        npcscores=npcscores * nprof,         &
        nchannels_rec=nchannels_rec * nprof, &
        pccomp=pccomp)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'deallocation error for rttov_direct structures'
    CALL rttov_exit(errorstatus)
  ENDIF

  CALL rttov_dealloc_coefs(errorstatus, coefs)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'coefs deallocation error'
  ENDIF

END PROGRAM example_pc_fwd
