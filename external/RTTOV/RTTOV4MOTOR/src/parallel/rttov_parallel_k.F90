! Description:
!> @file
!!   Runs the RTTOV Jacobian model with multiple threads under OpenMP.
!
!> @brief
!!   Runs the RTTOV Jacobian model with multiple threads under OpenMP.
!!
!! @details
!!   To take advantage of this RTTOV must be compiled with OpenMP
!!   compiler flags.
!!
!!   The arguments are identical to those for rttov_k.F90 though
!!   the traj and traj_k arguments cannot be used with the parallel interface.
!!
!!   The only additional optional final argument is nthreads which is
!!   used to specify the number of threads.
!!
!!   The additional "strategy" and "debug" arguments are for debugging
!!   purposes and are not intended for general use.
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
SUBROUTINE rttov_parallel_k(      &
          errorstatus      &
        , chanprof         &
        , opts             &
        , profiles         &
        , profiles_k       &
        , coefs            &
        , transmission     &
        , transmission_k   &
        , radiance         &
        , radiance_k       &
        , radiance2        &
        , calcemis         &
        , emissivity       &
        , emissivity_k     &
        , calcrefl         &
        , reflectance      &
        , reflectance_k    &
        , aer_opt_param    &
        , aer_opt_param_k  &
        , cld_opt_param    &
        , cld_opt_param_k  &
        , traj             &
        , traj_k           &
        , pccomp           &
        , pccomp_k         &
        , profiles_k_pc    &
        , profiles_k_rec   &
        , channels_rec     &
        , nthreads         &
        , strategy         &
        , debug)

!INTF_OFF
  USE rttov_const, ONLY: &
         errorstatus_success, &
         errorstatus_fatal,   &
         sensor_id_po
!INTF_ON
  USE rttov_types, ONLY:        &
          rttov_options,        &
          rttov_chanprof,       &
          rttov_emissivity,     &
          rttov_reflectance,    &
          rttov_opt_param,      &
          rttov_pccomp,         &
          rttov_coefs,          &
          rttov_profile,        &
          rttov_transmission,   &
          rttov_radiance,       &
          rttov_radiance2,      &
          rttov_traj

  USE parkind1, ONLY:jpim, jplm

  IMPLICIT NONE

  INTEGER(KIND=jpim),              INTENT(OUT)   :: errorstatus
  TYPE(rttov_chanprof),            INTENT(IN)    :: chanprof(:)
  TYPE(rttov_options),             INTENT(IN)    :: opts
  TYPE(rttov_profile),             INTENT(IN)    :: profiles(:)
  TYPE(rttov_coefs),               INTENT(IN)    :: coefs
  TYPE(rttov_transmission),        INTENT(INOUT) :: transmission
  TYPE(rttov_radiance),            INTENT(INOUT) :: radiance
  TYPE(rttov_radiance2), OPTIONAL, INTENT(INOUT) :: radiance2
  TYPE(rttov_profile),             INTENT(INOUT) :: profiles_k(:)
  TYPE(rttov_transmission),        INTENT(INOUT) :: transmission_k
  TYPE(rttov_radiance),            INTENT(INOUT) :: radiance_k

  LOGICAL(KIND=jplm),      OPTIONAL, TARGET, INTENT(IN)    :: calcemis(:)
  TYPE(rttov_emissivity),  OPTIONAL, TARGET, INTENT(INOUT) :: emissivity(:)
  TYPE(rttov_emissivity),  OPTIONAL, TARGET, INTENT(INOUT) :: emissivity_k(:)

  LOGICAL(KIND=jplm),      OPTIONAL, TARGET, INTENT(IN)    :: calcrefl(:)
  TYPE(rttov_reflectance), OPTIONAL, TARGET, INTENT(INOUT) :: reflectance(:)
  TYPE(rttov_reflectance), OPTIONAL, TARGET, INTENT(INOUT) :: reflectance_k(:)

  TYPE(rttov_opt_param), OPTIONAL, INTENT(IN)  :: aer_opt_param
  TYPE(rttov_opt_param), OPTIONAL, INTENT(IN)  :: aer_opt_param_k

  TYPE(rttov_opt_param), OPTIONAL, INTENT(IN)  :: cld_opt_param
  TYPE(rttov_opt_param), OPTIONAL, INTENT(IN)  :: cld_opt_param_k

  TYPE(rttov_traj),    OPTIONAL, INTENT(INOUT) :: traj
  TYPE(rttov_traj),    OPTIONAL, INTENT(INOUT) :: traj_k
  TYPE(rttov_pccomp),  OPTIONAL, INTENT(INOUT) :: pccomp
  TYPE(rttov_pccomp),  OPTIONAL, INTENT(INOUT) :: pccomp_k
  TYPE(rttov_profile), OPTIONAL, INTENT(INOUT) :: profiles_k_pc(:)
  TYPE(rttov_profile), OPTIONAL, INTENT(INOUT) :: profiles_k_rec(:)
  INTEGER(KIND=jpim),  OPTIONAL, INTENT(IN)    :: channels_rec(:)


  INTEGER(KIND=jpim),  OPTIONAL, INTENT(IN)    :: nthreads
  ! Strategy: 0 = no strategy (RTTOV computes band limits), 1 = profile-wise, 2 = channel-wise
  INTEGER(KIND=jpim),  OPTIONAL, INTENT(IN)    :: strategy
  LOGICAL(KIND=jplm),  OPTIONAL, INTENT(IN)    :: debug
!INTF_END

#include "rttov_k.interface"
#include "rttov_errorreport.interface"

  CHARACTER(LEN=*), PARAMETER :: NameOfRoutine = &
    "rttov_parallel_k"

  INTEGER(KIND=jpim) :: iband, nbands
  INTEGER(KIND=jpim) :: mthreads
  INTEGER(KIND=jpim) :: nsplits
  INTEGER(KIND=jpim) :: ichannelk, iprofilek
  INTEGER(KIND=jpim) :: nchannels
  INTEGER(KIND=jpim) :: nprofiles
  INTEGER(KIND=jpim) :: strategy1
  LOGICAL(KIND=jplm) :: debug1
  LOGICAL(KIND=jplm) :: do_pc, do_pcrttov_radrec, do_htfrtc_radrec, do_radrec
  INTEGER            :: mstat
  INTEGER            :: npcscores, nchan_rec
  INTEGER(KIND=jpim),       ALLOCATABLE :: errorstatus_x(:)                 ! errorstatus/bands
  INTEGER(KIND=jpim),       ALLOCATABLE :: ichannel1(:), ichannel2(:)       ! channel limits
  INTEGER(KIND=jpim),       ALLOCATABLE :: iprofile1(:), iprofile2(:)       ! profile limits
  INTEGER(KIND=jpim),       ALLOCATABLE :: ipc1(:), ipc2(:)                 ! PCscores limits
  INTEGER(KIND=jpim),       ALLOCATABLE :: ichannelrec1(:), ichannelrec2(:) ! PC reconstructed channel limits
  TYPE(rttov_chanprof),     ALLOCATABLE :: chanprof_x(:,:)                  ! chanprof/bands

  TYPE(rttov_opt_param),    ALLOCATABLE :: aer_opt_param_x(:)               ! user aerosol opt params
  TYPE(rttov_opt_param),    ALLOCATABLE :: aer_opt_param_k_x(:)

  TYPE(rttov_opt_param),    ALLOCATABLE :: cld_opt_param_x(:)               ! user cloud opt params
  TYPE(rttov_opt_param),    ALLOCATABLE :: cld_opt_param_k_x(:)

  LOGICAL(KIND=jplm),       POINTER :: calcemis_x(:)             ! switch for emis/bands (pointer assoc)
  TYPE(rttov_emissivity),   POINTER :: emissivity_x(:)           ! surface emis/bands (pointer assoc)
  TYPE(rttov_emissivity),   POINTER :: emissivity_k_x(:)

  LOGICAL(KIND=jplm),       POINTER :: calcrefl_x(:)             ! switch for refl/bands (pointer assoc)
  TYPE(rttov_reflectance),  POINTER :: reflectance_x(:)          ! surface refl/bands (pointer assoc)
  TYPE(rttov_reflectance),  POINTER :: reflectance_k_x(:)

  TYPE(rttov_transmission), ALLOCATABLE :: transmission_x(:)         ! transmission/bands (pointer assoc)
  TYPE(rttov_radiance),     ALLOCATABLE :: radiance_x(:)             ! radiance/bands (pointer assoc)
  TYPE(rttov_pccomp),       ALLOCATABLE :: pccomp_x(:)               ! pccomp/bands (pointer assoc)

  LOGICAL(KIND=jplm)                    :: calc_rad2
  TYPE(rttov_radiance2),    ALLOCATABLE :: radiance2_x(:)
  TYPE(rttov_transmission), ALLOCATABLE :: transmission_k_x(:)
  TYPE(rttov_radiance),     ALLOCATABLE :: radiance_k_x(:)
  TYPE(rttov_pccomp),       ALLOCATABLE :: pccomp_k_x(:)


! Activated if openmp
!
!$  INTEGER, External :: omp_get_thread_num
!$  INTEGER, External :: omp_get_max_threads
!

  do_pc = opts%rt_ir%pc%addpc .OR. opts%htfrtc_opts%htfrtc
  do_pcrttov_radrec = opts%rt_ir%pc%addpc .AND. opts%rt_ir%pc%addradrec
  do_htfrtc_radrec = opts%htfrtc_opts%htfrtc .AND. opts%htfrtc_opts%reconstruct
  do_radrec = do_pcrttov_radrec .OR. do_htfrtc_radrec

  nprofiles = SIZE(profiles)
  IF (opts%htfrtc_opts%htfrtc) THEN
    nchannels = coefs%coef_htfrtc%n_f * nprofiles
  ELSE
    nchannels = SIZE(chanprof)
  ENDIF

  debug1 = .FALSE.
  IF (PRESENT(debug)) debug1 = debug

  IF (PRESENT(traj) .AND. opts%config%verbose) THEN
    CALL rttov_errorreport(errorstatus_success, "traj argument will not be used !!", NameOfRoutine)
  ENDIF
  IF (PRESENT(traj_k) .AND. opts%config%verbose) THEN
    CALL rttov_errorreport(errorstatus_success, "traj_k argument will not be used !!", NameOfRoutine)
  ENDIF

  calc_rad2 = PRESENT(radiance2)

  strategy1 = 0
  IF (PRESENT(strategy)) THEN
    strategy1 = strategy
  ENDIF
  ! Must distribute work profile-wise among threads for PC-RTTOV, HTFRTC and polarised sensors
  IF (do_pc .OR. coefs%coef%id_sensor == sensor_id_po) strategy1 = 1

  IF (PRESENT(nthreads)) THEN
    mthreads = nthreads
  ELSE
! Default value for mthreads
    mthreads = 1_jpim
! Activated if openmp
!
!$  mthreads = omp_get_max_threads()
!
  ENDIF

! Check we do not have more bands than channels

  nsplits = 1_jpim
  nbands = mthreads * nsplits

  IF (strategy1 == 0) THEN
    IF (nbands > nchannels) nbands = nchannels
  ELSE IF (strategy1 == 1) THEN
    nbands = nprofiles
  ELSE IF (strategy1 == 2) THEN
    nbands = nchannels
  ELSE
    errorstatus = errorstatus_fatal
    CALL rttov_errorreport(errorstatus_fatal, "Unknown strategy", NameOfRoutine)
    RETURN
  ENDIF


  IF (debug1) THEN
    WRITE (*, '(A20," = ",I8)') "nprofiles", nprofiles
    WRITE (*, '(A20," = ",I8)') "nchannels", nchannels
    WRITE (*, '(A20," = ",I5)') "nbands", nbands
  ENDIF


  ALLOCATE(ichannel1(nbands), stat = mstat)
  IF (mstat .NE. 0) GOTO 100
  ALLOCATE(ichannel2(nbands), stat = mstat)
  IF (mstat .NE. 0) GOTO 100

  IF (do_pc) THEN
    ALLOCATE(ipc1(nbands), stat = mstat)
    IF (mstat .NE. 0) GOTO 100
    ALLOCATE(ipc2(nbands), stat = mstat)
    IF (mstat .NE. 0) GOTO 100
    IF (do_radrec) THEN
      ALLOCATE(ichannelrec1(nbands), stat = mstat)
      IF (mstat .NE. 0) GOTO 100
      ALLOCATE(ichannelrec2(nbands), stat = mstat)
      IF (mstat .NE. 0) GOTO 100
    ENDIF
  ENDIF

  ALLOCATE(iprofile1(nbands), stat = mstat)
  IF (mstat .NE. 0) GOTO 100
  ALLOCATE(iprofile2(nbands), stat = mstat)
  IF (mstat .NE. 0) GOTO 100

  ALLOCATE(transmission_x(nbands), stat = mstat)
  IF (mstat .NE. 0) GOTO 100
  ALLOCATE(radiance_x(nbands), stat = mstat)
  IF (mstat .NE. 0) GOTO 100
  IF (do_pc) THEN
    ALLOCATE(pccomp_x(nbands), stat = mstat)
    IF (mstat .NE. 0) GOTO 100
  ENDIF

  ALLOCATE(aer_opt_param_x(nbands), stat = mstat)
  IF (mstat .NE. 0) GOTO 100
  ALLOCATE(aer_opt_param_k_x(nbands), stat = mstat)
  IF (mstat .NE. 0) GOTO 100

  ALLOCATE(cld_opt_param_x(nbands), stat = mstat)
  IF (mstat .NE. 0) GOTO 100
  ALLOCATE(cld_opt_param_k_x(nbands), stat = mstat)
  IF (mstat .NE. 0) GOTO 100

  IF (calc_rad2) THEN
    ALLOCATE(radiance2_x(nbands), stat = mstat)
    IF (mstat .NE. 0) GOTO 100
  ENDIF
  ALLOCATE(transmission_k_x(nbands), stat = mstat)
  IF (mstat .NE. 0) GOTO 100
  ALLOCATE(radiance_k_x(nbands), stat = mstat)
  IF (mstat .NE. 0) GOTO 100
  IF (do_pc) THEN
    ALLOCATE(pccomp_k_x(nbands), stat = mstat)
    IF (mstat .NE. 0) GOTO 100
  ENDIF

  IF (strategy1 == 1) THEN
    IF (do_pc) THEN
      ! Always the same number of channels per profile for PCs (nbands == nprofiles)
      ALLOCATE(chanprof_x(nchannels / nbands, nbands), stat = mstat)
    ELSE IF (coefs%coef%id_sensor == sensor_id_po) THEN
      ! Account for extreme variations in channels per profile
      ALLOCATE(chanprof_x(coefs%coef%fmv_chn, nbands), stat = mstat)
    ELSE
      ! In the general case ensure chanprof_x is large enough if one profile
      ! is run for every channel and another is run for only one channel.
      ! The only way to do this without explicitly finding the maximum number of
      ! channels per profile is to size chanprof_x by chanprof. This shouldn't be a
      ! problem since the prof-by-prof strategy is intended for debugging rather than
      ! for general use.
      ALLOCATE(chanprof_x(nchannels, nbands), stat = mstat)
    ENDIF
  ELSE IF (strategy1 == 2) THEN
    ! One channel per band (i.e. nbands == nchannels)
    ALLOCATE(chanprof_x(1, nbands), stat = mstat)
  ELSE
    ALLOCATE(chanprof_x(nchannels / nbands + nbands, nbands), stat = mstat)
  ENDIF
  IF (mstat .NE. 0) GOTO 100
  ALLOCATE(errorstatus_x(nbands), stat = mstat)
  IF (mstat .NE. 0) GOTO 100

  ichannel1(:) = 0_jpim
  ichannel2(:) = 0_jpim

  IF (do_pc) THEN
    ipc1(:) = 0_jpim
    ipc2(:) = 0_jpim
    IF (do_radrec) THEN
      ichannelrec1(:) = 0_jpim
      ichannelrec2(:) = 0_jpim
    ENDIF
  ENDIF

  iprofile1(:) = 0_jpim
  iprofile2(:) = 0_jpim

  !
  ! Build the array of lower/upper bounds of bands
  !

  IF (strategy1 == 1) THEN
    IF (do_pc) THEN
      ! Always test for HTFRTC before PC-RTTOV as HTFRTC takes precedence
      IF (opts%htfrtc_opts%htfrtc) THEN
        npcscores = MIN(MAX(opts%htfrtc_opts%n_pc_in, 1), coefs%coef_htfrtc%n_pc)
      ELSE
        npcscores = MIN(MAX(opts%rt_ir%pc%npcscores, 1), coefs%coef_pccomp%fmv_pc_mnum)
      ENDIF
      IF (do_htfrtc_radrec) THEN
        nchan_rec = coefs%coef_htfrtc%n_ch
      ELSEIF (do_pcrttov_radrec) THEN
        nchan_rec = SIZE(channels_rec)
      ENDIF
    ENDIF
    IF (do_pc) THEN
      ! With one profile per band, and the same number of channels per profile,
      ! splitting the channels, pcscores and reconstructed channels among bands
      ! is quite simple
      DO iband = 1, nbands
        ichannel1(iband) = (iband-1_jpim) * nchannels / nprofiles + 1_jpim
        ichannel2(iband) = ichannel1(iband) + nchannels / nprofiles - 1_jpim
        ipc1(iband) = (iband-1_jpim) * npcscores + 1_jpim
        ipc2(iband) = ipc1(iband) + npcscores - 1_jpim
        IF (do_radrec) THEN
          ichannelrec1(iband) = (iband-1_jpim) * nchan_rec + 1_jpim
          ichannelrec2(iband) = ichannelrec1(iband) + nchan_rec - 1_jpim
        ENDIF
      ENDDO
    ELSE
      ichannelk = 1_jpim
      DO iband = 1, nbands
        ichannel1(iband) = ichannelk
        iprofilek = chanprof(ichannelk)%prof
        DO
          IF (ichannelk > nchannels) EXIT
          IF (iprofilek > nprofiles) EXIT
          IF (iprofilek .NE. chanprof(ichannelk)%prof) EXIT
          ichannelk = ichannelk + 1_jpim
        ENDDO
        ichannel2(iband) = ichannelk - 1_jpim
      ENDDO
    ENDIF
  ELSE
    !
    ! For IR & MW sensors, we try to make bands
    ! with the same number of channels
    !
    ichannelk = 1_jpim
    DO iband = 1, nbands
      ichannel1(iband) = ichannelk
      ichannelk = ichannelk + nchannels / nbands - 1_jpim
      IF (iband == nbands) ichannelk = nchannels
      ichannel2(iband) = ichannelk
      ichannelk = ichannelk + 1_jpim
    ENDDO
  ENDIF

  chanprof_x(:,:)%chan = 0_jpim
  chanprof_x(:,:)%prof = 0_jpim

  IF (opts%htfrtc_opts%htfrtc) THEN
    DO iband = 1, nbands
      iprofile1(iband) = iband
      iprofile2(iband) = iband
    ENDDO
  ELSE
    DO iband = 1, nbands
      iprofile1(iband) = chanprof(ichannel1(iband))%prof
      iprofile2(iband) = chanprof(ichannel2(iband))%prof
      chanprof_x(1:ichannel2(iband) - ichannel1(iband) + 1, iband)%chan = &
              chanprof(ichannel1(iband):ichannel2(iband))%chan
      DO ichannelk = ichannel1(iband), ichannel2(iband)
        chanprof_x(ichannelk - ichannel1(iband) + 1, iband)%prof &
            = chanprof(ichannelk)%prof - iprofile1(iband) + 1
      ENDDO
    ENDDO
  ENDIF




  IF (debug1) THEN
    DO iband = 1, nbands
      WRITE (*, *)
      WRITE (*, '(A20," = ",I5)') "band", iband
      WRITE (*, '(A20," = ",100(I5))') "channels", chanprof(ichannel1(iband):ichannel2(iband))%chan
      WRITE (*, '(A20," = ",100(I5))') "lprofiles", chanprof(ichannel1(iband):ichannel2(iband))%prof
      WRITE (*, '(A20," = ",100(I5))') "chanprof_x%chan", chanprof_x(:,iband)%chan
      WRITE (*, '(A20," = ",100(I5))') "chanprof_x%prof", chanprof_x(:,iband)%prof
      WRITE (*, *)
    ENDDO
  ENDIF


  IF (PRESENT(calcemis)) THEN
    calcemis_x => calcemis
  ELSE
    ALLOCATE(calcemis_x(nchannels), stat = mstat)
    IF (mstat .NE. 0) GOTO 100
    calcemis_x(:) = .TRUE.
  ENDIF

  IF (PRESENT(emissivity)) THEN
    emissivity_x => emissivity
  ELSE
    ALLOCATE(emissivity_x(nchannels), stat = mstat)
    IF (mstat .NE. 0) GOTO 100
  ENDIF

  IF (PRESENT(emissivity_k)) THEN
    emissivity_k_x => emissivity_k
  ELSE
    ALLOCATE(emissivity_k_x(nchannels), stat = mstat)
    IF (mstat .NE. 0) GOTO 100
  ENDIF

  IF (PRESENT(calcrefl)) THEN
    calcrefl_x => calcrefl
  ELSE
    ALLOCATE(calcrefl_x(nchannels), stat = mstat)
    IF (mstat .NE. 0) GOTO 100
    calcrefl_x(:) = .TRUE.
  ENDIF

  IF (PRESENT(reflectance)) THEN
    reflectance_x => reflectance
  ELSE
    ALLOCATE(reflectance_x(nchannels), stat = mstat)
    IF (mstat .NE. 0) GOTO 100
  ENDIF

  IF (PRESENT(reflectance_k)) THEN
    reflectance_k_x => reflectance_k
  ELSE
    ALLOCATE(reflectance_k_x(nchannels), stat = mstat)
    IF (mstat .NE. 0) GOTO 100
  ENDIF


  !
  ! we make the correct pointers associations between temporary transmissions
  ! and radiances and arrays passed as input
  !
  DO iband = 1, nbands
    IF (.NOT. opts%htfrtc_opts%htfrtc) THEN
      CALL AssociateTransmission(ichannel1(iband), ichannel2(iband), &
         transmission_x(iband), transmission)
      CALL AssociateRadiance(ichannel1(iband), ichannel2(iband), &
         radiance_x(iband), radiance)
    ENDIF

    IF (do_pc) THEN
      IF (do_radrec) THEN
        CALL AssociatePCcomp(ipc1(iband), ipc2(iband), &
          pccomp_x(iband), pccomp, ichannelrec1(iband), ichannelrec2(iband))
      ELSE
        CALL AssociatePCcomp(ipc1(iband), ipc2(iband), &
          pccomp_x(iband), pccomp)
      ENDIF
    ENDIF

    IF (PRESENT(aer_opt_param) .AND. opts%rt_ir%user_aer_opt_param) THEN
      CALL AssociateOptParam(ichannel1(iband), ichannel2(iband), &
        aer_opt_param_x(iband), aer_opt_param)
    ENDIF
    IF (PRESENT(aer_opt_param_k) .AND. opts%rt_ir%user_aer_opt_param) THEN
      CALL AssociateOptParam(ichannel1(iband), ichannel2(iband), &
        aer_opt_param_k_x(iband), aer_opt_param_k)
    ELSE
      NULLIFY(aer_opt_param_k_x(iband)%abs,     &
              aer_opt_param_k_x(iband)%sca,     &
              aer_opt_param_k_x(iband)%bpr,     &
              aer_opt_param_k_x(iband)%legcoef, &
              aer_opt_param_k_x(iband)%pha,     &
              aer_opt_param_k_x(iband)%phangle)
    ENDIF

    IF (PRESENT(cld_opt_param) .AND. opts%rt_ir%user_cld_opt_param) THEN
      CALL AssociateOptParam(ichannel1(iband), ichannel2(iband), &
        cld_opt_param_x(iband), cld_opt_param)
    ENDIF
    IF (PRESENT(cld_opt_param_k) .AND. opts%rt_ir%user_cld_opt_param) THEN
      CALL AssociateOptParam(ichannel1(iband), ichannel2(iband), &
        cld_opt_param_k_x(iband), cld_opt_param_k)
    ELSE
      NULLIFY(cld_opt_param_k_x(iband)%abs,     &
              cld_opt_param_k_x(iband)%sca,     &
              cld_opt_param_k_x(iband)%bpr,     &
              cld_opt_param_k_x(iband)%legcoef, &
              cld_opt_param_k_x(iband)%pha,     &
              cld_opt_param_k_x(iband)%phangle)
    ENDIF

    IF (calc_rad2) THEN
      CALL AssociateRadiance2(ichannel1(iband), ichannel2(iband), &
         radiance2_x(iband), radiance2)
    ENDIF
    IF (.NOT. opts%htfrtc_opts%htfrtc) THEN
      CALL AssociateTransmission(ichannel1(iband), ichannel2(iband), &
         transmission_k_x(iband), transmission_k)
      CALL AssociateRadiance(ichannel1(iband), ichannel2(iband), &
         radiance_k_x(iband), radiance_k)
    ENDIF

    IF (do_pc) THEN
      IF (do_radrec) THEN
        CALL AssociatePCcomp(ipc1(iband), ipc2(iband), &
          pccomp_k_x(iband), pccomp_k, ichannelrec1(iband), ichannelrec2(iband))
      ELSE
        CALL AssociatePCcomp(ipc1(iband), ipc2(iband), &
          pccomp_k_x(iband), pccomp_k)
      ENDIF
    ENDIF

  ENDDO


!$OMP PARALLEL DO PRIVATE(iband) NUM_THREADS(mthreads) SCHEDULE(DYNAMIC)
DO iband = 1, nbands

  IF (debug1) THEN
!    PRINT *, "================================================"
!    PRINT *, "THREAD = ", OMP_GET_THREAD_NUM(), " BAND = ", iband, &
!    " ichannel1 = ", ichannel1(iband), " ichannel2 = ", ichannel2(iband)

    WRITE (*, *)
    WRITE (*, '(A)') "=========================== rttov_k ============================"
    WRITE (*, '(A20)') "errorstatus"
    WRITE (*, '(A20," = ",I5)') "nprofiles", iprofile2(iband) - iprofile1(iband) + 1
    WRITE (*, '(A20," = ",I5)') "nchannels", ichannel2(iband) - ichannel1(iband) + 1
    WRITE (*, '(A20," = ",100(I5))') "channels", chanprof_x(:,iband)%chan
    WRITE (*, '(A20," = ",100(I5))') "lprofiles", chanprof_x(:,iband)%prof
    WRITE (*, '(A20," = ",A20," (",I5,":",I5,") ")') "profiles",          &
     "profiles", iprofile1(iband), iprofile2(iband)
    WRITE (*, '(A20," = ",A20," (",I5,":",I5,") ")') "profiles_k",        &
     "profiles_k", ichannel1(iband), ichannel2(iband)
    WRITE (*, '(A20)') "coefs"

#ifndef RTTOV_NAG53
    IF (PRESENT(calcemis)) THEN
      WRITE (*, '(A20," = ",A20," (",I5,":",I5,") ")') "calcemis",        &
       "calcemis", ichannel1(iband), ichannel2(iband)
    ENDIF
    IF (PRESENT(emissivity)) THEN
      WRITE (*, '(A20," = ",A20," (",I5,":",I5,") ")') "emis%emis_in",    &
       "emis%emis_in", ichannel1(iband), ichannel2(iband)
      WRITE (*, '(A20," = ",A20," (",I5,":",I5,") ")') "emis%emis_out",   &
       "emis%emis_out", ichannel1(iband), ichannel2(iband)
    ENDIF
    IF (PRESENT(emissivity_k)) THEN
      WRITE (*, '(A20," = ",A20," (",I5,":",I5,") ")') "emis_k%emis_in", &
       "emis_k%emis_in", ichannel1(iband), ichannel2(iband)
      WRITE (*, '(A20," = ",A20," (",I5,":",I5,") ")') "emis_k%emis_out", &
       "emis_k%emis_out", ichannel1(iband), ichannel2(iband)
    ENDIF
    IF (PRESENT(calcrefl)) THEN
      WRITE (*, '(A20," = ",A20," (",I5,":",I5,") ")') "calcrefl",        &
       "calcrefl", ichannel1(iband), ichannel2(iband)
    ENDIF
    IF (PRESENT(reflectance)) THEN
      WRITE (*, '(A20," = ",A20," (",I5,":",I5,") ")') "refl%refl_in",    &
       "refl%refl_in", ichannel1(iband), ichannel2(iband)
      WRITE (*, '(A20," = ",A20," (",I5,":",I5,") ")') "refl%refl_out",   &
       "refl%refl_out", ichannel1(iband), ichannel2(iband)
    ENDIF
    IF (PRESENT(reflectance_k)) THEN
      WRITE (*, '(A20," = ",A20," (",I5,":",I5,") ")') "refl_k%refl_in", &
       "refl_k%refl_in", ichannel1(iband), ichannel2(iband)
      WRITE (*, '(A20," = ",A20," (",I5,":",I5,") ")') "refl_k%refl_out", &
       "refl_k%refl_out", ichannel1(iband), ichannel2(iband)
    ENDIF
#endif
    WRITE (*, *)
  ENDIF


  IF (do_htfrtc_radrec) THEN
    CALL &
        rttov_k                                                                     &
        (errorstatus_x(iband)                                                       &
        , chanprof_x(1:ichannel2(iband) - ichannel1(iband) + 1, iband)              &
        , opts                                                                      &
        , profiles(iprofile1(iband):iprofile2(iband))                               &
        , profiles_k(ichannel1(iband):ichannel2(iband))                             &
        , coefs                                                                     &
        , transmission_x(iband)                                                     &
        , transmission_k_x(iband)                                                   &
        , radiance_x(iband)                                                         &
        , radiance_k_x(iband)                                                       &
        , calcemis        = calcemis_x(ichannel1(iband):ichannel2(iband))           &
        , emissivity      = emissivity_x(ichannel1(iband):ichannel2(iband))         &
        , emissivity_k    = emissivity_k_x(ichannel1(iband):ichannel2(iband))       &
        , calcrefl        = calcrefl_x(ichannel1(iband):ichannel2(iband))           &
        , reflectance     = reflectance_x(ichannel1(iband):ichannel2(iband))        &
        , reflectance_k   = reflectance_k_x(ichannel1(iband):ichannel2(iband))      &
        , aer_opt_param   = aer_opt_param_x(iband)                                  &
        , aer_opt_param_k = aer_opt_param_k_x(iband)                                &
        , cld_opt_param   = cld_opt_param_x(iband)                                  &
        , cld_opt_param_k = cld_opt_param_k_x(iband)                                &
        , pccomp          = pccomp_x(iband)                                         &
        , pccomp_k        = pccomp_k_x(iband)                                       &
        , profiles_k_pc   = profiles_k_pc(ipc1(iband):ipc2(iband))                  &
        , profiles_k_rec  = profiles_k_rec(ichannelrec1(iband):ichannelrec2(iband)) &
        , channels_rec    = channels_rec                                            &
      )
  ELSEIF (do_pcrttov_radrec) THEN
    CALL &
        rttov_k                                                                     &
        (errorstatus_x(iband)                                                       &
        , chanprof_x(1:ichannel2(iband) - ichannel1(iband) + 1, iband)              &
        , opts                                                                      &
        , profiles(iprofile1(iband):iprofile2(iband))                               &
        , profiles_k(ichannel1(iband):ichannel2(iband))                             &
        , coefs                                                                     &
        , transmission_x(iband)                                                     &
        , transmission_k_x(iband)                                                   &
        , radiance_x(iband)                                                         &
        , radiance_k_x(iband)                                                       &
        , calcemis        = calcemis_x(ichannel1(iband):ichannel2(iband))           &
        , emissivity      = emissivity_x(ichannel1(iband):ichannel2(iband))         &
        , emissivity_k    = emissivity_k_x(ichannel1(iband):ichannel2(iband))       &
        , calcrefl        = calcrefl_x(ichannel1(iband):ichannel2(iband))           &
        , reflectance     = reflectance_x(ichannel1(iband):ichannel2(iband))        &
        , reflectance_k   = reflectance_k_x(ichannel1(iband):ichannel2(iband))      &
        , aer_opt_param   = aer_opt_param_x(iband)                                  &
        , aer_opt_param_k = aer_opt_param_k_x(iband)                                &
        , cld_opt_param   = cld_opt_param_x(iband)                                  &
        , cld_opt_param_k = cld_opt_param_k_x(iband)                                &
        , pccomp          = pccomp_x(iband)                                         &
        , pccomp_k        = pccomp_k_x(iband)                                       &
        , profiles_k_rec  = profiles_k_rec(ichannelrec1(iband):ichannelrec2(iband)) &
        , channels_rec    = channels_rec                                            &
      )
  ELSEIF (do_pc) THEN
    CALL &
        rttov_k                                                                  &
        (errorstatus_x(iband)                                                    &
        , chanprof_x(1:ichannel2(iband) - ichannel1(iband) + 1, iband)           &
        , opts                                                                   &
        , profiles(iprofile1(iband):iprofile2(iband))                            &
        , profiles_k(ichannel1(iband):ichannel2(iband))                          &
        , coefs                                                                  &
        , transmission_x(iband)                                                  &
        , transmission_k_x(iband)                                                &
        , radiance_x(iband)                                                      &
        , radiance_k_x(iband)                                                    &
        , calcemis         = calcemis_x(ichannel1(iband):ichannel2(iband))       &
        , emissivity       = emissivity_x(ichannel1(iband):ichannel2(iband))     &
        , emissivity_k     = emissivity_k_x(ichannel1(iband):ichannel2(iband))   &
        , calcrefl         = calcrefl_x(ichannel1(iband):ichannel2(iband))       &
        , reflectance      = reflectance_x(ichannel1(iband):ichannel2(iband))    &
        , reflectance_k    = reflectance_k_x(ichannel1(iband):ichannel2(iband))  &
        , aer_opt_param    = aer_opt_param_x(iband)                              &
        , aer_opt_param_k  = aer_opt_param_k_x(iband)                            &
        , cld_opt_param    = cld_opt_param_x(iband)                              &
        , cld_opt_param_k  = cld_opt_param_k_x(iband)                            &
        , pccomp           = pccomp_x(iband)                                     &
        , pccomp_k         = pccomp_k_x(iband)                                   &
        , profiles_k_pc    = profiles_k_pc(ipc1(iband):ipc2(iband))              &
        , channels_rec     = channels_rec                                        &
      )
  ELSE
    IF (calc_rad2) THEN
      CALL &
          rttov_k                                                                  &
          (errorstatus_x(iband)                                                    &
          , chanprof_x(1:ichannel2(iband) - ichannel1(iband) + 1, iband)           &
          , opts                                                                   &
          , profiles(iprofile1(iband):iprofile2(iband))                            &
          , profiles_k(ichannel1(iband):ichannel2(iband))                          &
          , coefs                                                                  &
          , transmission_x(iband)                                                  &
          , transmission_k_x(iband)                                                &
          , radiance_x(iband)                                                      &
          , radiance_k_x(iband)                                                    &
          , radiance2_x(iband)                                                     &
          , calcemis         = calcemis_x(ichannel1(iband):ichannel2(iband))       &
          , emissivity       = emissivity_x(ichannel1(iband):ichannel2(iband))     &
          , emissivity_k     = emissivity_k_x(ichannel1(iband):ichannel2(iband))   &
          , calcrefl         = calcrefl_x(ichannel1(iband):ichannel2(iband))       &
          , reflectance      = reflectance_x(ichannel1(iband):ichannel2(iband))    &
          , reflectance_k    = reflectance_k_x(ichannel1(iband):ichannel2(iband))  &
          , aer_opt_param    = aer_opt_param_x(iband)                              &
          , aer_opt_param_k  = aer_opt_param_k_x(iband)                            &
          , cld_opt_param    = cld_opt_param_x(iband)                              &
          , cld_opt_param_k  = cld_opt_param_k_x(iband)                            &
        )
    ELSE
      CALL &
          rttov_k                                                                  &
          (errorstatus_x(iband)                                                    &
          , chanprof_x(1:ichannel2(iband) - ichannel1(iband) + 1, iband)           &
          , opts                                                                   &
          , profiles(iprofile1(iband):iprofile2(iband))                            &
          , profiles_k(ichannel1(iband):ichannel2(iband))                          &
          , coefs                                                                  &
          , transmission_x(iband)                                                  &
          , transmission_k_x(iband)                                                &
          , radiance_x(iband)                                                      &
          , radiance_k_x(iband)                                                    &
          , calcemis         = calcemis_x(ichannel1(iband):ichannel2(iband))       &
          , emissivity       = emissivity_x(ichannel1(iband):ichannel2(iband))     &
          , emissivity_k     = emissivity_k_x(ichannel1(iband):ichannel2(iband))   &
          , calcrefl         = calcrefl_x(ichannel1(iband):ichannel2(iband))       &
          , reflectance      = reflectance_x(ichannel1(iband):ichannel2(iband))    &
          , reflectance_k    = reflectance_k_x(ichannel1(iband):ichannel2(iband))  &
          , aer_opt_param    = aer_opt_param_x(iband)                              &
          , aer_opt_param_k  = aer_opt_param_k_x(iband)                            &
          , cld_opt_param    = cld_opt_param_x(iband)                              &
          , cld_opt_param_k  = cld_opt_param_k_x(iband)                            &
        )
    ENDIF
  ENDIF

ENDDO
!$OMP END PARALLEL DO


  !
  ! we have to construct the global errorstatus array
  ! from what we got from every thread
  !
  errorstatus = errorstatus_success
  IF (ANY(errorstatus_x(:) /= errorstatus_success)) errorstatus = errorstatus_fatal

!  DO iband = 1, nbands
!    DO iprofilek = iprofile1(iband), iprofile2(iband)
!      PRINT *, iprofilek, " <= ", iband, iprofilek - iprofile1(iband) + 1, &
!        errorstatus_x(iprofilek - iprofile1(iband) + 1, iband)
!    ENDDO
!  ENDDO

  ! Set the plane_parallel flag in the radiance_k structure
  radiance_k%plane_parallel = ANY(radiance_k_x(:)%plane_parallel)

  DEALLOCATE(chanprof_x)
  DEALLOCATE(errorstatus_x)

  IF (PRESENT(calcemis)) THEN
    NULLIFY(calcemis_x)
  ELSE
    DEALLOCATE(calcemis_x)
  ENDIF

  IF (PRESENT(emissivity)) THEN
    NULLIFY(emissivity_x)
  ELSE
    DEALLOCATE(emissivity_x)
  ENDIF

  IF (PRESENT(emissivity_k)) THEN
    NULLIFY(emissivity_k_x)
  ELSE
    DEALLOCATE(emissivity_k_x)
  ENDIF

  IF (PRESENT(calcrefl)) THEN
    NULLIFY(calcrefl_x)
  ELSE
    DEALLOCATE(calcrefl_x)
  ENDIF

  IF (PRESENT(reflectance)) THEN
    NULLIFY(reflectance_x)
  ELSE
    DEALLOCATE(reflectance_x)
  ENDIF

  IF (PRESENT(reflectance_k)) THEN
    NULLIFY(reflectance_k_x)
  ELSE
    DEALLOCATE(reflectance_k_x)
  ENDIF

  DEALLOCATE(aer_opt_param_x)
  DEALLOCATE(aer_opt_param_k_x)

  DEALLOCATE(cld_opt_param_x)
  DEALLOCATE(cld_opt_param_k_x)

  DEALLOCATE(transmission_x)
  DEALLOCATE(radiance_x)

  IF (do_pc) DEALLOCATE(pccomp_x)

  IF (calc_rad2) DEALLOCATE(radiance2_x)
  DEALLOCATE(transmission_k_x)
  DEALLOCATE(radiance_k_x)
  IF (do_pc) DEALLOCATE(pccomp_k_x)


  DEALLOCATE(iprofile2)
  DEALLOCATE(iprofile1)

  IF (do_pc) THEN
    DEALLOCATE(ipc2)
    DEALLOCATE(ipc1)
    IF (do_radrec) THEN
      DEALLOCATE(ichannelrec2)
      DEALLOCATE(ichannelrec1)
    ENDIF
  ENDIF

  DEALLOCATE(ichannel2)
  DEALLOCATE(ichannel1)

  RETURN

100 CONTINUE
    CALL rttov_errorreport(errorstatus_fatal, "Memory allocation failed", NameOfRoutine)
    STOP

CONTAINS

  ! Ancillary routines

  !> Set up one radiance structure to point to a subset of another
  !! @param[in]     ichannel1       first channel in range of subset
  !! @param[in]     ichannel2       last channel in range of subset
  !! @param[out]    radiance_x      radiance structure pointing to subset of user structure
  !! @param[in]     radiance        user radiance structure
  SUBROUTINE AssociateRadiance(ichannel1, ichannel2, radiance_x, radiance)
    INTEGER(KIND=jpim),   INTENT(IN)  :: ichannel1, ichannel2
    TYPE(rttov_radiance), INTENT(IN)  :: radiance
    TYPE(rttov_radiance), INTENT(OUT) :: radiance_x

    radiance_x%clear &
      => radiance%clear(ichannel1:ichannel2)

    radiance_x%cloudy &
      => radiance%cloudy(ichannel1:ichannel2)

    radiance_x%total &
      => radiance%total(ichannel1:ichannel2)

    radiance_x%bt &
      => radiance%bt(ichannel1:ichannel2)

    radiance_x%bt_clear &
      => radiance%bt_clear(ichannel1:ichannel2)

    radiance_x%refl &
      => radiance%refl(ichannel1:ichannel2)

    radiance_x%refl_clear &
      => radiance%refl_clear(ichannel1:ichannel2)

    radiance_x%overcast &
      => radiance%overcast(:,ichannel1:ichannel2)

    radiance_x%quality &
      => radiance%quality(ichannel1:ichannel2)

    radiance_x%geometric_height &
      => radiance%geometric_height(:,ichannel1:ichannel2)

  END SUBROUTINE

  !> Set up one secondary radiance structure to point to a subset of another
  !! @param[in]     ichannel1       first channel in range of subset
  !! @param[in]     ichannel2       last channel in range of subset
  !! @param[out]    radiance_x      secondary radiance structure pointing to subset of user structure
  !! @param[in]     radiance        user secondary radiance structure
  SUBROUTINE AssociateRadiance2(ichannel1, ichannel2, radiance_x, radiance)
    INTEGER(KIND=jpim),   INTENT(IN)   :: ichannel1, ichannel2
    TYPE(rttov_radiance2), INTENT(IN)  :: radiance
    TYPE(rttov_radiance2), INTENT(OUT) :: radiance_x

    radiance_x%upclear &
      => radiance%upclear(ichannel1:ichannel2)

    radiance_x%dnclear &
      => radiance%dnclear(ichannel1:ichannel2)

    radiance_x%refldnclear &
      => radiance%refldnclear(ichannel1:ichannel2)

    radiance_x%up &
      => radiance%up(:,ichannel1:ichannel2)

    radiance_x%down &
      => radiance%down(:,ichannel1:ichannel2)

    radiance_x%surf &
      => radiance%surf(:,ichannel1:ichannel2)

  END SUBROUTINE

  !> Set up one transmittance structure to point to a subset of another
  !! @param[in]     ichannel1       first channel in range of subset
  !! @param[in]     ichannel2       last channel in range of subset
  !! @param[out]    transmission_x  transmittance structure pointing to subset of user structure
  !! @param[in]     transmission    user transmittance structure
  SUBROUTINE AssociateTransmission(ichannel1, ichannel2, transmission_x, transmission)
    INTEGER(KIND=jpim),       INTENT(IN)  :: ichannel1, ichannel2
    TYPE(rttov_transmission), INTENT(IN)  :: transmission
    TYPE(rttov_transmission), INTENT(OUT) :: transmission_x

    transmission_x%tau_levels &
      => transmission%tau_levels(:,ichannel1:ichannel2)

    transmission_x%tau_total &
      => transmission%tau_total(ichannel1:ichannel2)

    transmission_x%tausun_levels_path1 &
      => transmission%tausun_levels_path1(:,ichannel1:ichannel2)

    transmission_x%tausun_total_path1 &
      => transmission%tausun_total_path1(ichannel1:ichannel2)

    transmission_x%tausun_levels_path2 &
      => transmission%tausun_levels_path2(:,ichannel1:ichannel2)

    transmission_x%tausun_total_path2 &
      => transmission%tausun_total_path2(ichannel1:ichannel2)

    transmission_x%tau_levels_cld &
      => transmission%tau_levels_cld(:,ichannel1:ichannel2)

    transmission_x%tau_total_cld &
      => transmission%tau_total_cld(ichannel1:ichannel2)

  END SUBROUTINE

  !> Set up one PC score/reconstructed radiance structure to point to a subset of another
  !! @param[in]     ipc1          first PC score in range of subset
  !! @param[in]     ipc2          last PC score in range of subset
  !! @param[out]    pccomp_x      PC score/reconstructed radiance structure pointing to subset of user structure
  !! @param[in]     pccomp        user PC score/reconstructed radiance structure
  !! @param[in]     ichannelrec1  first rec. rad. channel in range of subset, optional
  !! @param[in]     ichannelrec2  last rec. rad. channel in range of subset, optional
  SUBROUTINE AssociatePCcomp(ipc1, ipc2, pccomp_x, pccomp, ichannelrec1, ichannelrec2)
    INTEGER(KIND=jpim), INTENT(IN)           :: ipc1, ipc2
    TYPE(rttov_pccomp), INTENT(IN)           :: pccomp
    TYPE(rttov_pccomp), INTENT(OUT)          :: pccomp_x
    INTEGER(KIND=jpim), INTENT(IN), OPTIONAL :: ichannelrec1, ichannelrec2

    pccomp_x%total_pcscores &
      => pccomp%total_pcscores(ipc1:ipc2)

    IF (ASSOCIATED(pccomp%clear_pcscores)) &
      pccomp_x%clear_pcscores &
        => pccomp%clear_pcscores(ipc1:ipc2)

    IF (ASSOCIATED(pccomp%overcast_pcscores)) &
      pccomp_x%overcast_pcscores &
        => pccomp%overcast_pcscores(:,ipc1:ipc2)

    IF (ASSOCIATED(pccomp%cloudy_pcscores)) &
      pccomp_x%cloudy_pcscores &
        => pccomp%cloudy_pcscores(ipc1:ipc2)

    IF (PRESENT(ichannelrec1)) THEN
      pccomp_x%bt_pccomp &
        => pccomp%bt_pccomp(ichannelrec1:ichannelrec2)

      pccomp_x%total_pccomp &
        => pccomp%total_pccomp(ichannelrec1:ichannelrec2)

      IF (ASSOCIATED(pccomp%bt_clear_pccomp)) &
        pccomp_x%bt_clear_pccomp &
          => pccomp%bt_clear_pccomp(ichannelrec1:ichannelrec2)

      IF (ASSOCIATED(pccomp%clear_pccomp)) &
        pccomp_x%clear_pccomp &
          => pccomp%clear_pccomp(ichannelrec1:ichannelrec2)

      IF (ASSOCIATED(pccomp%overcast_pccomp)) &
        pccomp_x%overcast_pccomp &
          => pccomp%overcast_pccomp(:,ichannelrec1:ichannelrec2)

      IF (ASSOCIATED(pccomp%cloudy_pccomp)) &
        pccomp_x%cloudy_pccomp &
          => pccomp%cloudy_pccomp(ichannelrec1:ichannelrec2)
    ENDIF

  END SUBROUTINE

  !> Set up one cloud/aerosol optical property structure to point to a subset of another
  !! @param[in]     ichannel1     first channel in range of subset
  !! @param[in]     ichannel2     last channel in range of subset
  !! @param[out]    opt_param_x   cloud/aerosol optical properties pointing to subset of user input
  !! @param[in]     opt_param     cloud/aerosol optical properties input by user
  SUBROUTINE AssociateOptParam(ichannel1, ichannel2, opt_param_x, opt_param)
    INTEGER(KIND=jpim),    INTENT(IN)  :: ichannel1, ichannel2
    TYPE(rttov_opt_param), INTENT(IN)  :: opt_param
    TYPE(rttov_opt_param), INTENT(OUT) :: opt_param_x

    opt_param_x%abs &
      => opt_param%abs(:,ichannel1:ichannel2)

    opt_param_x%sca &
      => opt_param%sca(:,ichannel1:ichannel2)

    opt_param_x%bpr &
      => opt_param%bpr(:,ichannel1:ichannel2)

    opt_param_x%nmom = opt_param%nmom

    opt_param_x%legcoef &
      => opt_param%legcoef(:,:,ichannel1:ichannel2)

    opt_param_x%pha &
      => opt_param%pha(:,:,ichannel1:ichannel2)

    opt_param_x%phangle &
      => opt_param%phangle(:)

    opt_param_x%phasefn_int%zminphadiff = opt_param%phasefn_int%zminphadiff

    IF (ASSOCIATED(opt_param%phasefn_int%cosphangle)) &
      opt_param_x%phasefn_int%cosphangle &
        => opt_param%phasefn_int%cosphangle(:)

    IF (ASSOCIATED(opt_param%phasefn_int%iphangle)) &
      opt_param_x%phasefn_int%iphangle &
        => opt_param%phasefn_int%iphangle(:)

  END SUBROUTINE

END SUBROUTINE
