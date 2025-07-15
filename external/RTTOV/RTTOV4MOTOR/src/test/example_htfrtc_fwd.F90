! Description:
!> @file
!!   Example calling HTFRTC RT model through RTTOV for a direct model
!!   clear-sky simulation.
!
!> @brief
!!   Example calling HTFRTC RT model through RTTOV for a direct model
!!   clear-sky simulation.
!!
!! @details
!!   This example is most easily run using the run_example_htfrtc_fwd.sh
!!   script.
!!
!!   This program requires the following files:
!!     the file containing input profiles (e.g. prof.dat)
!!     the file containing the list of channels for which
!!       to reconstruct radiances (e.g. radrec.dat) - optional
!!     the HTFRTC static input file (htfrtc_coef_static.nc)
!!     the HTFRTC input file for the sensor to be simulated
!!
!!   The output is written to a file called example_htfrtc_fwd_output.dat
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
!    Copyright 2018, EUMETSAT, All Rights Reserved.
!
PROGRAM example_htfrtc_fwd

  ! rttov_const contains useful RTTOV constants
  USE rttov_const, ONLY :     &
         errorstatus_success, &
         errorstatus_fatal

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

   ! The rttov_emis_atlas_data type must be imported separately
   USE mod_rttov_emis_atlas, ONLY : &
         rttov_emis_atlas_data, &
         atlas_type_ir

  ! jpim, jprb and jplm are the RTTOV integer, real and logical KINDs
  USE parkind1, ONLY : jpim, jplm, jprb

  USE rttov_unix_env, ONLY : rttov_exit

  IMPLICIT NONE

#include "rttov_direct.interface"
#include "rttov_parallel_direct.interface"
#include "rttov_read_coefs_htfrtc.interface"
#include "rttov_dealloc_coefs.interface"
#include "rttov_alloc_direct.interface"
#include "rttov_print_opts.interface"
#include "rttov_print_profile.interface"
#include "rttov_skipcommentline.interface"

! Use emissivity atlas
#include "rttov_setup_emis_atlas.interface"
#include "rttov_get_emis.interface"
#include "rttov_deallocate_emis_atlas.interface"

  !--------------------------
  !
  INTEGER(KIND=jpim), PARAMETER :: iup   = 20   ! unit for input profile file
  INTEGER(KIND=jpim), PARAMETER :: ioout = 21   ! unit for output

  ! RTTOV variables/structures
  !====================
  TYPE(rttov_options)             :: opts                      ! Options structure
  TYPE(rttov_coefs)               :: coefs                     ! Coefficients structure
  TYPE(rttov_chanprof),   POINTER :: chanprof(:)    => NULL()  ! Required argument, unused by HTFRTC
  LOGICAL(KIND=jplm),     POINTER :: calcemis(:)    => NULL()  ! Optional, flag to indicate calculation of emissivity within HTFRTC
  TYPE(rttov_emissivity), POINTER :: emissivity(:)  => NULL()  ! Optional, input/output surface emissivity
  TYPE(rttov_profile),    POINTER :: profiles(:)    => NULL()  ! Input profiles
  TYPE(rttov_transmission)        :: transmission              ! Required argument, unused by HTFRTC
  TYPE(rttov_radiance)            :: radiance                  ! Required argument, unused by HTFRTC
  TYPE(rttov_pccomp)              :: pccomp                    ! Output PC structure
  INTEGER(KIND=jpim),     POINTER :: channels_rec(:) => NULL() ! Reconstructed radiance channel list

  TYPE(rttov_emis_atlas_data)     :: emis_atlas                ! Data structure for emissivity atlas

  INTEGER(KIND=jpim)              :: errorstatus               ! Return error status of RTTOV subroutine calls

  CHARACTER(LEN=18)  :: NameOfRoutine = 'example_htfrtc_fwd'

  ! variables for input
  !====================
  CHARACTER(LEN=256) :: fname_coef, fname_sensor
  CHARACTER(LEN=256) :: prof_filename, radrec_filename
  INTEGER(KIND=jpim) :: nthreads
  INTEGER(KIND=jpim) :: imonth
  INTEGER(KIND=jpim) :: nlevels
  INTEGER(KIND=jpim) :: nprof
  INTEGER(KIND=jpim) :: npcscores
  INTEGER(KIND=jpim) :: nchannels_rec
  INTEGER(KIND=jpim) :: nchanprof
  LOGICAL(KIND=jplm) :: exists
  ! loop variables
  INTEGER(KIND=jpim) :: lo, hi
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

  ! However, when running HTFRTC step 3 is simplified and step 4 is omitted as
  ! HTFRTC has its own training data and does not require predictors channels 
  ! (and as such has no need for RTTOV to allocate memory for them).
  ! Consequently, it is only necessary to allocate the profile structure
  ! and pccomp, where the output will be stored. Step 6 is optional:
  ! by default (if the calcemis and emissivity arguments are omitted) HTFRTC
  ! calculates sea surface emissivities using IREMIS or the PC-RTTOV emissivity
  ! model and uses 0.98 and 0.99 for land and sea-ice respectively. Alternatively
  ! it is possible to pass emissivities in as for RTTOV. The calcemis and
  ! emissivity arguments must be of size nchanprof = n_f * nprofiles where n_f
  ! is the number of predictor frequencies used by HTFRTC, available in
  ! coefs%coef_htfrtc%n_f (note this is analogous to PC-RTTOV). The associated
  ! wavenumbers are given in coefs%coef_htfrtc%freq(:).
  ! To pass emissivities in for a given profile you must set the
  ! elements of calcemis corresponding to that profile to false for all 
  ! frequencies. Note that passing the emissivity array provides access to the
  ! emissivities used/computed by HTFRTC for the centroid (predictor)
  ! frequencies in the emis_out member, similar to PC-RTTOV.

  ! If nthreads is greater than 1 the parallel RTTOV interface is used.
  ! To take advantage of multi-threaded execution you must have compiled
  ! RTTOV with openmp enabled. See the user guide and the compiler flags.

  errorstatus = 0_jpim

  !=====================================================
  !========== Interactive inputs == start ==============

  WRITE(0,*) 'enter path of HTFRTC static coefficient file'
  READ(*,*) fname_coef
  WRITE(0,*) 'enter path of HTFRTC sensor coefficient file'
  READ(*,*) fname_sensor
  WRITE(0,*) 'enter path of file containing profile data'
  READ(*,*) prof_filename
  WRITE(0,*) 'enter path of file containing reconstructed channel list'
  READ(*,*) radrec_filename
  WRITE(0,*) 'enter number of profiles'
  READ(*,*) nprof
  WRITE(0,*) 'enter number of profile levels'
  READ(*,*) nlevels
  WRITE(0,*) 'enter profile month 1-12 (for emis atlas) or 0 for no atlas'
  READ(*,*) imonth
  WRITE(0,*) 'enter number of pcscores (1-300)'
  READ(*,*) npcscores
  WRITE(0,*) 'enter number of threads to use'
  READ(*,*) nthreads

  ! --------------------------------------------------------------------------
  ! 1. Initialise HTFRTC specific RTTOV options
  ! --------------------------------------------------------------------------

  opts % htfrtc_opts % htfrtc = .TRUE.        ! Select HTFRTC
  opts % htfrtc_opts % n_pc_in = npcscores

  opts % htfrtc_opts % simple_cloud = .FALSE. ! Set to true to enable simple cloud scheme
  opts % htfrtc_opts % overcast     = .FALSE. ! Set to true to enable overcast radiance calculations

  opts % rt_all % ozone_data        = .TRUE.  ! Set the relevant flag to .TRUE.
  opts % rt_all % co2_data          = .FALSE. !   when supplying a profile of the
  opts % rt_all % n2o_data          = .FALSE. !   given trace gas
  opts % rt_all % ch4_data          = .FALSE. !
  opts % rt_all % co_data           = .FALSE. !
  opts % rt_all % so2_data          = .FALSE. !

  ! In this example we reconstruct radiances if there is an input file
  ! containing a channel list
  INQUIRE(file=radrec_filename, exist=exists)
  opts % htfrtc_opts % reconstruct = exists

  !========== Interactive inputs == end ==============
  !===================================================

  ! --------------------------------------------------------------------------
  ! Determine the number of reconstructed radiances per profile (nchannels_rec)
  ! --------------------------------------------------------------------------

  nchannels_rec = 0
  IF (opts % htfrtc_opts % reconstruct) THEN

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

      opts % htfrtc_opts % reconstruct = .FALSE.

    ELSE IF (nchannels_rec == 0) THEN
      ! If the number of channels is set to 0 then reconstruct all instrument channels
      ! in HTFRTC file.

      ! To do this we do not supply the channels_rec argument to rttov_read_coefs_htfrtc,
      ! but we do need to obtain the number of sensor channels in order to allocate the
      ! RTTOV structures below: this is done below after the coefficients have been read in.

    ELSE IF (nchannels_rec > 0) THEN
      ! Otherwise read the channel list from the file

      ALLOCATE(channels_rec(nchannels_rec))
      READ(iup,*) channels_rec(:)

    ENDIF

    CLOSE(iup)
  ENDIF

  ! Ensure we don't have unassociated pointers below when reconstruct is FALSE
  IF (nchannels_rec <= 0) ALLOCATE(channels_rec(0))


  ! --------------------------------------------------------------------------
  ! 2. Read coefficients
  ! --------------------------------------------------------------------------
  IF (opts % htfrtc_opts % reconstruct .AND. nchannels_rec > 0) THEN
    CALL rttov_read_coefs_htfrtc(errorstatus, coefs, fname_coef, fname_sensor, channels_rec)
  ELSE
    CALL rttov_read_coefs_htfrtc(errorstatus, coefs, fname_coef, fname_sensor)
  ENDIF
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'fatal error reading coefficients'
    CALL rttov_exit(errorstatus)
  ENDIF

  ! Obtain the number of sensor channels if we are reconstructing all radiances
  IF (opts % htfrtc_opts % reconstruct .AND. nchannels_rec == 0) THEN
    nchannels_rec = coefs % coef_htfrtc % n_ch
  ENDIF


  ! --------------------------------------------------------------------------
  ! 3. Allocate RTTOV input and output structures
  ! --------------------------------------------------------------------------

  ! NB If passing calcemis/emissivity arguments to RTTOV, then allocate chanprof
  !    to be the same size to satisfy the interface (it does not need to be 
  !    filled with any values). For HTFRTC these arrays are n_f * nprofiles
  !    in size as described above. If you wish to allow HTFRTC to provide all
  !    emissivities then you can omit the calcemis and emissivity arguments
  !    (and there is no need to allocate them).
  nchanprof = coefs % coef_htfrtc % n_f * nprof
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

  IF (imonth >= 1 .AND. imonth <= 12) THEN
    CALL rttov_setup_emis_atlas(          &
                errorstatus,              &
                opts,                     &
                imonth,                   &
                atlas_type_ir,            &
                emis_atlas,               &
                path = '../../emis_data', & ! The default path to atlas data
                coefs = coefs) ! Always supply this for HTFRTC as the predictor
                               ! frequencies are the same for all sensors: this
                               ! makes the atlas much faster to access
    IF (errorstatus /= errorstatus_success) THEN
      WRITE(*,*) 'error initialising emissivity atlas'
      CALL rttov_exit(errorstatus)
    ENDIF
  ENDIF

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

    ! Read pressure (hPa), temp (K), WV, O3 (gas units as read above)
    READ(iup,*) profiles(iprof) % p(:)
    CALL rttov_skipcommentline(iup, errorstatus)
    READ(iup,*) profiles(iprof) % t(:)
    CALL rttov_skipcommentline(iup, errorstatus)
    READ(iup,*) profiles(iprof) % q(:)
    CALL rttov_skipcommentline(iup, errorstatus)
    READ(iup,*) profiles(iprof) % o3(:)
    CALL rttov_skipcommentline(iup, errorstatus)

    ! Near surface variables
    READ(iup,*) profiles(iprof) % s2m % t, &
                profiles(iprof) % s2m % q, &
                profiles(iprof) % s2m % p, &
                profiles(iprof) % s2m % u, &
                profiles(iprof) % s2m % v
    CALL rttov_skipcommentline(iup, errorstatus)

    ! Skin temperature
    READ(iup,*) profiles(iprof) % skin % t
    CALL rttov_skipcommentline(iup, errorstatus)

    ! Surface type
    READ(iup,*) profiles(iprof) % skin % surftype
    CALL rttov_skipcommentline(iup, errorstatus)

    ! Elevation
    READ(iup,*) profiles(iprof) % elevation, &
                profiles(iprof) % latitude,  &
                profiles(iprof) % longitude
    CALL rttov_skipcommentline(iup, errorstatus)

    ! Satellite zenith angle
    READ(iup,*) profiles(iprof) % zenangle
    CALL rttov_skipcommentline(iup, errorstatus)

    ! Cloud variables for simple cloud scheme (requires HTFRTC simple_cloud flag to be true)
    READ(iup,*) profiles(iprof) % ctp, &
                profiles(iprof) % cfraction
    CALL rttov_skipcommentline(iup, errorstatus)

  ENDDO
  CLOSE(iup)

  !========== Read profiles == end =============
  !=============================================

  ! --------------------------------------------------------------------------
  ! 6. Specify surface emissivity
  ! --------------------------------------------------------------------------

  IF (imonth >= 1 .AND. imonth <= 12) THEN
    ! Use emissivity atlas: emissivities are returned for the centroid
    ! (predictor) wavenumbers as required by HTFRTC
    CALL rttov_get_emis(             &
              errorstatus,           &
              opts,                  &
              chanprof,              &
              profiles,              &
              coefs,                 &
              emis_atlas,            &
              emissivity(:) % emis_in)
    IF (errorstatus /= errorstatus_success) THEN
      WRITE(*,*) 'error reading emissivity atlas'
      CALL rttov_exit(errorstatus)
    ENDIF

    ! Calculate emissivity within RTTOV where the atlas emissivity value is
    ! zero or less
    calcemis(:) = (emissivity(:) % emis_in <= 0._jprb)
  ELSE
    emissivity(:) % emis_in = 0._jprb
    calcemis(:) = .TRUE.
  ENDIF

  ! --------------------------------------------------------------------------
  ! 7. Call RTTOV forward model
  ! --------------------------------------------------------------------------
  IF (nthreads <= 1) THEN
    CALL rttov_direct(               &
           errorstatus,              &! out   error flag
           chanprof,                 &! in    not used by HTFRTC
           opts,                     &! in    options structure
           profiles,                 &! in    profile array
           coefs,                    &! in    not used by HTFRTC
           transmission,             &! inout not used by HTFRTC
           radiance,                 &! inout not used by HTFRTC
           calcemis   = calcemis,    &! in    flag for internal emissivity calcs, optional
           emissivity = emissivity,  &! inout input/output emissivities per channel, optional
           pccomp     = pccomp)       ! inout computed PC scores
  ELSE
    CALL rttov_parallel_direct(      &
           errorstatus,              &! out   error flag
           chanprof,                 &! in    not used by HTFRTC
           opts,                     &! in    options structure
           profiles,                 &! in    profile array
           coefs,                    &! in    not used by HTFRTC
           transmission,             &! inout not used by HTFRTC
           radiance,                 &! inout not used by HTFRTC
           calcemis   = calcemis,    &! in    flag for internal emissivity calcs, optional
           emissivity = emissivity,  &! inout input/output emissivities per channel, optional
           pccomp     = pccomp,      &! inout computed PC scores
           nthreads   = nthreads)     ! in    number of threads to use
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
  WRITE(ioout,*)' Instrument '//TRIM(fname_sensor)
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
    WRITE(ioout, '(10(G15.6E3))') pccomp % total_pcscores(lo:hi)
    WRITE(ioout,*) ')'

    IF ( opts % htfrtc_opts % reconstruct ) THEN
      lo = 1 + (iprof-1) * nchannels_rec
      hi = iprof * nchannels_rec
      WRITE(ioout,*)
      WRITE(ioout,*) 'Radiances reconstructed using principal components'
      WRITE(ioout, '(A," = (")' ) 'PCCOMP%TOTAL_PCCOMP'
      WRITE(ioout, '(10(G15.6E3))') pccomp % total_pccomp(lo:hi)
      WRITE(ioout,*) ')'

      WRITE(ioout,*)' '
      WRITE(ioout,*) 'Brightness temp equivalent to radiances reconstructed using PCs'
      WRITE(ioout, '(A," = (")' ) 'PCCOMP%BT_PCCOMP'
      WRITE(ioout, '(10(G15.6E3))') pccomp % bt_pccomp(lo:hi)
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
  IF (imonth >= 1 .AND. imonth <= 12) THEN
    ! Deallocate emissivity atlas
    CALL rttov_deallocate_emis_atlas(emis_atlas)
  ENDIF

  CALL rttov_alloc_direct(                   &
        errorstatus,                         &
        0_jpim,                              &  ! 1 => allocate
        nprof,                               &
        nchanprof,                           &
        nlevels,                             &
        chanprof,                            &
        opts,                                &
        profiles,                            &
        coefs,                               &
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

  DEALLOCATE(channels_rec, stat=errorstatus)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'deallocation error'
    CALL rttov_exit(errorstatus)
  ENDIF

  CALL rttov_dealloc_coefs(errorstatus, coefs)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'deallocation error for coefs structure'
    CALL rttov_exit(errorstatus)
  ENDIF

END PROGRAM example_htfrtc_fwd
