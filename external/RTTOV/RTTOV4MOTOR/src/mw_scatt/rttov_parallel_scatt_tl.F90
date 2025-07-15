! Description:
!> @file
!!   Runs the RTTOV-SCATT tangent linear model with multiple threads under OpenMP.
!
!> @brief
!!   Runs the RTTOV-SCATT tangent linear model with multiple threads under OpenMP.
!!
!! @details
!!   To take advantage of this RTTOV must be compiled with OpenMP
!!   compiler flags.
!!
!!   The arguments are identical to those for rttov_scatt_tl.F90.
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
!    Copyright 2017, EUMETSAT, All Rights Reserved.
!
SUBROUTINE rttov_parallel_scatt_tl( &
          errorstatus          &
        , opts_scatt           &
        , nlevels              &
        , chanprof             &
        , frequencies          &
        , profiles             &
        , cld_profiles         &
        , coef_rttov           &
        , coef_scatt           &
        , calcemis             &
        , emissivity           &
        , profiles_tl          &
        , cld_profiles_tl      &
        , emissivity_tl        &
        , radiance             &
        , radiance_tl          &
        , reflectivity         &
        , reflectivity_tl      &
        , nthreads             &
        , strategy             &
        , debug)

!INTF_OFF
  USE rttov_const, ONLY: &
         errorstatus_success, &
         errorstatus_fatal,   &
         sensor_id_po
!INTF_ON
  USE rttov_types, ONLY:        &
          rttov_options_scatt,  &
          rttov_chanprof,       &
          rttov_emissivity,     &
          rttov_coefs,          &
          rttov_scatt_coef,     &
          rttov_profile,        &
          rttov_profile_cloud,  &
          rttov_radiance,       &
          rttov_reflectivity

  USE parkind1, ONLY: jpim, jplm

  IMPLICIT NONE

  INTEGER(KIND=jpim),                 INTENT(OUT)   :: errorstatus
  TYPE(rttov_options_scatt),          INTENT(IN)    :: opts_scatt
  INTEGER(KIND=jpim),                 INTENT(IN)    :: nlevels
  TYPE(rttov_chanprof),               INTENT(IN)    :: chanprof(:)
  INTEGER(KIND=jpim),                 INTENT(IN)    :: frequencies(SIZE(chanprof))
  TYPE(rttov_profile),                INTENT(IN)    :: profiles(:)
  TYPE(rttov_profile_cloud),          INTENT(IN)    :: cld_profiles(:)
  TYPE(rttov_coefs),                  INTENT(IN)    :: coef_rttov
  TYPE(rttov_scatt_coef),             INTENT(IN)    :: coef_scatt
  LOGICAL(KIND=jplm),                 INTENT(IN)    :: calcemis(:)
  TYPE(rttov_emissivity),             INTENT(INOUT) :: emissivity(:)

  TYPE(rttov_emissivity),             INTENT(INOUT) :: emissivity_tl(:)
  TYPE(rttov_profile),                INTENT(INOUT) :: profiles_tl(:)
  TYPE(rttov_profile_cloud),          INTENT(INOUT) :: cld_profiles_tl(:)

  TYPE(rttov_radiance),               INTENT(INOUT) :: radiance
  TYPE(rttov_radiance),               INTENT(INOUT) :: radiance_tl
  TYPE(rttov_reflectivity), OPTIONAL, INTENT(INOUT) :: reflectivity ! Optional for radar
  TYPE(rttov_reflectivity), OPTIONAL, INTENT(INOUT) :: reflectivity_tl

  INTEGER(KIND=jpim),       OPTIONAL, INTENT(IN)    :: nthreads
  ! Strategy: 0 = no strategy (RTTOV computes band limits), 1 = profile-wise, 2 = channel-wise
  INTEGER(KIND=jpim),       OPTIONAL, INTENT(IN)    :: strategy
  LOGICAL(KIND=jplm),       OPTIONAL, INTENT(IN)    :: debug
!INTF_END

#include "rttov_scatt_tl.interface"
#include "rttov_errorreport.interface"

  CHARACTER(LEN=*), PARAMETER :: NameOfRoutine = &
    "rttov_parallel_scatt_tl"

  INTEGER(KIND=jpim) :: iband, nbands
  INTEGER(KIND=jpim) :: mthreads
  INTEGER(KIND=jpim) :: nsplits
  INTEGER(KIND=jpim) :: ichannelk, iprofilek
  INTEGER(KIND=jpim) :: nchannels
  INTEGER(KIND=jpim) :: nprofiles
  INTEGER(KIND=jpim) :: strategy1
  LOGICAL(KIND=jplm) :: debug1
  INTEGER            :: mstat
  INTEGER(KIND=jpim),        ALLOCATABLE :: errorstatus_x(:)                 ! errorstatus/bands
  INTEGER(KIND=jpim),        ALLOCATABLE :: ichannel1(:), ichannel2(:)       ! channel limits
  INTEGER(KIND=jpim),        ALLOCATABLE :: iprofile1(:), iprofile2(:)       ! profile limits
  TYPE(rttov_chanprof),      ALLOCATABLE :: chanprof_x(:,:)                  ! chanprof/bands
  TYPE(rttov_radiance),      ALLOCATABLE :: radiance_x(:)                    ! radiance/bands (pointer assoc)
  TYPE(rttov_reflectivity),  ALLOCATABLE :: reflectivity_x(:)
  TYPE(rttov_radiance),      ALLOCATABLE :: radiance_tl_x(:)
  TYPE(rttov_reflectivity),  ALLOCATABLE :: reflectivity_tl_x(:)

! Activated if openmp
!
!$  INTEGER, External :: omp_get_thread_num
!$  INTEGER, External :: omp_get_max_threads
!

  nchannels = SIZE(chanprof)
  nprofiles = SIZE(profiles)

  debug1 = .FALSE.
  IF (PRESENT(debug)) debug1 = debug

  strategy1 = 0
  IF (PRESENT(strategy)) THEN
    strategy1 = strategy
  ENDIF

  ! Must distribute work profile-wise among threads for polarised sensors
  IF (coef_rttov%coef%id_sensor == sensor_id_po) strategy1 = 1

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

  ALLOCATE(iprofile1(nbands), stat = mstat)
  IF (mstat .NE. 0) GOTO 100
  ALLOCATE(iprofile2(nbands), stat = mstat)
  IF (mstat .NE. 0) GOTO 100

  ALLOCATE(radiance_x(nbands), stat = mstat)
  IF (mstat .NE. 0) GOTO 100
  IF (PRESENT(reflectivity)) then
    ALLOCATE(reflectivity_x(nbands), stat = mstat)
    IF (mstat .NE. 0) GOTO 100
  ENDIF

  ALLOCATE(radiance_tl_x(nbands), stat = mstat)
  IF (mstat .NE. 0) GOTO 100
  IF (PRESENT(reflectivity_tl)) then
    ALLOCATE(reflectivity_tl_x(nbands), stat = mstat)
    IF (mstat .NE. 0) GOTO 100
  ENDIF

  IF (strategy1 == 1) THEN
    ALLOCATE(chanprof_x(coef_rttov%coef%fmv_chn, nbands), stat = mstat)
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

  iprofile1(:) = 0_jpim
  iprofile2(:) = 0_jpim

  !
  ! Build the array of lower/upper bounds of bands
  !

  IF (strategy1 == 1) THEN
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
  ELSE
    !
    ! Try to make bands with the same number of channels
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


  !
  ! we make the correct pointers associations between temporary
  ! radiances and arrays passed as input
  !
  DO iband = 1, nbands

    CALL AssociateRadiance(ichannel1(iband), ichannel2(iband), &
       radiance_x(iband), radiance)

    CALL AssociateRadiance(ichannel1(iband), ichannel2(iband), &
       radiance_tl_x(iband), radiance_tl)

    IF (PRESENT(reflectivity)) THEN
      CALL AssociateReflectivity(ichannel1(iband), ichannel2(iband), &
         reflectivity_x(iband), reflectivity)
    ENDIF
    IF (PRESENT(reflectivity_tl)) THEN
      CALL AssociateReflectivity(ichannel1(iband), ichannel2(iband), &
         reflectivity_tl_x(iband), reflectivity_tl)
    ENDIF

  ENDDO


!$OMP PARALLEL DO PRIVATE(iband) NUM_THREADS(mthreads) SCHEDULE(DYNAMIC)
  DO iband = 1, nbands

    IF (debug1) THEN
!    PRINT *, "================================================"
!    PRINT *, "THREAD = ", OMP_GET_THREAD_NUM(), " BAND = ", iband, &
!    " ichannel1 = ", ichannel1(iband), " ichannel2 = ", ichannel2(iband)

      WRITE (*, *)
      WRITE (*, '(A)') "=========================== rttov_tl ==========================="
      WRITE (*, '(A20)') "errorstatus"
      WRITE (*, '(A20," = ",I5)') "nprofiles", iprofile2(iband) - iprofile1(iband) + 1
      WRITE (*, '(A20," = ",I5)') "nchannels", ichannel2(iband) - ichannel1(iband) + 1
      WRITE (*, '(A20," = ",100(I5))') "channels", chanprof_x(:,iband)%chan
      WRITE (*, '(A20," = ",100(I5))') "lprofiles", chanprof_x(:,iband)%prof
      WRITE (*, '(A20," = ",A20," (",I5,":",I5,") ")') "profiles",          &
       "profiles", iprofile1(iband), iprofile2(iband)
      WRITE (*, '(A20," = ",A20," (",I5,":",I5,") ")') "profiles_tl",       &
       "profiles_tl", iprofile1(iband), iprofile2(iband)
      WRITE (*, '(A20)') "coefs"

#ifndef RTTOV_NAG53
      WRITE (*, '(A20," = ",A20," (",I5,":",I5,") ")') "calcemis",        &
       "calcemis", ichannel1(iband), ichannel2(iband)
      WRITE (*, '(A20," = ",A20," (",I5,":",I5,") ")') "emis%emis_in",    &
       "emis%emis_in", ichannel1(iband), ichannel2(iband)
      WRITE (*, '(A20," = ",A20," (",I5,":",I5,") ")') "emis%emis_out",   &
       "emis%emis_out", ichannel1(iband), ichannel2(iband)
      WRITE (*, '(A20," = ",A20," (",I5,":",I5,") ")') "emis_tl%emis_in", &
       "emis_tl%emis_in", ichannel1(iband), ichannel2(iband)
      WRITE (*, '(A20," = ",A20," (",I5,":",I5,") ")') "emis_tl%emis_out", &
       "emis_tl%emis_out", ichannel1(iband), ichannel2(iband)
#endif
      WRITE (*, *)
    ENDIF

    IF (PRESENT(reflectivity) .AND. PRESENT(reflectivity_tl)) THEN
      CALL rttov_scatt_tl                                                      &
        (errorstatus_x(iband)                                                  &
        , opts_scatt                                                           &
        , nlevels                                                              &
        , chanprof_x(1:ichannel2(iband) - ichannel1(iband) + 1, iband)         &
        , frequencies(ichannel1(iband):ichannel2(iband))                       &
        , profiles(iprofile1(iband):iprofile2(iband))                          &
        , cld_profiles(iprofile1(iband):iprofile2(iband))                      &
        , coef_rttov                                                           &
        , coef_scatt                                                           &
        , calcemis(ichannel1(iband):ichannel2(iband))                          &
        , emissivity(ichannel1(iband):ichannel2(iband))                        &
        , profiles_tl(iprofile1(iband):iprofile2(iband))                       &
        , cld_profiles_tl(iprofile1(iband):iprofile2(iband))                   &
        , emissivity_tl(ichannel1(iband):ichannel2(iband))                     &
        , radiance_x(iband)                                                    &
        , radiance_tl_x(iband) &
        , reflectivity = reflectivity_x(iband) &
        , reflectivity_tl = reflectivity_tl_x(iband))
    ELSE
      CALL rttov_scatt_tl                                                      &
        (errorstatus_x(iband)                                                  &
        , opts_scatt                                                           &
        , nlevels                                                              &
        , chanprof_x(1:ichannel2(iband) - ichannel1(iband) + 1, iband)         &
        , frequencies(ichannel1(iband):ichannel2(iband))                       &
        , profiles(iprofile1(iband):iprofile2(iband))                          &
        , cld_profiles(iprofile1(iband):iprofile2(iband))                      &
        , coef_rttov                                                           &
        , coef_scatt                                                           &
        , calcemis(ichannel1(iband):ichannel2(iband))                          &
        , emissivity(ichannel1(iband):ichannel2(iband))                        &
        , profiles_tl(iprofile1(iband):iprofile2(iband))                       &
        , cld_profiles_tl(iprofile1(iband):iprofile2(iband))                   &
        , emissivity_tl(ichannel1(iband):ichannel2(iband))                     &
        , radiance_x(iband)                                                    &
        , radiance_tl_x(iband))
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

  ! Set the plane_parallel flag in the radiance_tl structure
  radiance_tl%plane_parallel = ANY(radiance_tl_x(:)%plane_parallel)


  DEALLOCATE(chanprof_x)
  DEALLOCATE(errorstatus_x)

  DEALLOCATE(radiance_x)
  IF (PRESENT(reflectivity)) THEN
    DEALLOCATE(reflectivity_x)
  ENDIF

  DEALLOCATE(radiance_tl_x)
  IF (PRESENT(reflectivity_tl)) THEN
    DEALLOCATE(reflectivity_tl_x)
  ENDIF


  DEALLOCATE(iprofile2)
  DEALLOCATE(iprofile1)

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
    INTEGER(KIND=jpim),       INTENT(IN)  :: ichannel1, ichannel2
    TYPE(rttov_radiance),     INTENT(IN)  :: radiance
    TYPE(rttov_radiance),     INTENT(OUT) :: radiance_x

    radiance_x % clear &
      => radiance % clear(ichannel1:ichannel2)

    radiance_x % cloudy &
      => radiance % cloudy(ichannel1:ichannel2)

    radiance_x % total &
      => radiance % total(ichannel1:ichannel2)

    radiance_x % bt &
      => radiance % bt(ichannel1:ichannel2)

    radiance_x % bt_clear &
      => radiance % bt_clear(ichannel1:ichannel2)

    radiance_x % refl &
      => radiance % refl(ichannel1:ichannel2)

    radiance_x % refl_clear &
      => radiance % refl_clear(ichannel1:ichannel2)

    radiance_x % overcast &
      => radiance % overcast(:,ichannel1:ichannel2)

    radiance_x % quality &
      => radiance % quality(ichannel1:ichannel2)

    radiance_x % geometric_height &
      => radiance % geometric_height(:,ichannel1:ichannel2)

  END SUBROUTINE

  !> Set up one reflectivity structure to point to a subset of another
  !! @param[in]     ichannel1       first channel in range of subset
  !! @param[in]     ichannel2       last channel in range of subset
  !! @param[out]    reflectivity_x  reflectivity structure pointing to subset of user structure
  !! @param[in]     reflectivity    user reflectivity structure
  SUBROUTINE AssociateReflectivity(ichannel1, ichannel2, reflectivity_x, reflectivity)
    INTEGER(KIND=jpim),       INTENT(IN)  :: ichannel1, ichannel2
    TYPE(rttov_reflectivity), INTENT(IN)  :: reflectivity
    TYPE(rttov_reflectivity), INTENT(OUT) :: reflectivity_x

    reflectivity_x % zef &
      => reflectivity % zef(:,ichannel1:ichannel2)

    reflectivity_x % azef &
      => reflectivity % azef(:,ichannel1:ichannel2)

  END SUBROUTINE


END SUBROUTINE
