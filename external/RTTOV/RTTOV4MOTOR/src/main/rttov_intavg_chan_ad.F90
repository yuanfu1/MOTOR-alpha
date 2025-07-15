!
SUBROUTINE rttov_intavg_chan_ad( &
            & opts,        &
            & thermal,     &
            & solar,       &
            & kni,         &
            & kno,         &
            & chanprof,    &
            & ProfIn,      &
            & ProfOut,     &
            & ProfOut_ad,  &
            & OpdepIn,     &
            & OpdepIn_ad,  &
            & OpdepOut_ad)
! Copyright:
!    This software was developed within the context of
!    the EUMETSAT Satellite Application Facility on
!    Numerical Weather Prediction (NWP SAF), under the
!    Cooperation Agreement dated 7 December 2016, between
!    EUMETSAT and the Met Office, UK, by one or more partners
!    within the NWP SAF. The partners in the NWP SAF are
!    the Met Office, ECMWF, DWD and MeteoFrance.
!
!    Copyright 2010, EUMETSAT, All Rights Reserved.
!
! Current Code Owner: SAF NWP
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".
!
!-----------------------------------------------------------------------
! usage of arguments
!              set of optical depths
!              ---------------------------
! nchannels,     number of channels this profile
! solar         solar switch by channel
! kni           number of source levs
! kno           number of target levs
! nprofiles     number of profiles
! lprofiles     channels -> profiles
! ProfIn        for target p-levs & addsolar
! ProfOut       for source p-levs
! ProfOut_ad    for source p-levs
! OpdepIn       direct for source
! OpdepIn_ad    target for interpolator
! OpdepOut_ad   source for interpolator
!
!--------------------------------------------------------------------------
! application within RTTOV (for NWP SAF)
!   - profiles on USER levels -> prediction of gas/atm opdeps on COEF levels
!   - gas/atm opdeps on COEF levels -> transmittances on USER levels for RTE
!---------------------------------------------------------------------------
!     History:
!
!     Version   Date      Comment
!     -------   ----      -------
!
!     1         12/2006   main code intavg restructured Peter Rayer
!
!     2         10/2007   Optimisation Deborah Salmond and Pascal Brunel
!
!     3         11/2007   Extrapolate optical depths near surface so as to
!                         replicate behaviour of rttov_transmit.F90 when
!                         RTTOV levels are used throughout. - Alan Geer
!                         Niels Bormann
!
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!
!--------------------------------------------------------------------
!
! imported type definitions:
  USE rttov_types, ONLY :  &
       & rttov_chanprof, &
       & rttov_options,  &
       & rttov_profile,  &
       & rttov_opdp_path
  USE parkind1, ONLY : jpim, jplm
!INTF_OFF
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jprb
  USE rttov_const, ONLY :    &
    realtol,                 &
    interp_rochon,           &
    interp_loglinear,        &
    interp_rochon_loglinear, &
    interp_rochon_wfn,       &
    interp_rochon_loglinear_wfn
!INTF_ON
  IMPLICIT NONE
! --- Subroutine arguments
  TYPE(rttov_chanprof),  INTENT(IN)    :: chanprof(:)
  TYPE(rttov_options),   INTENT(IN)    :: opts
  LOGICAL(KIND=jplm),    INTENT(IN)    :: thermal(SIZE(chanprof)) ! Thermal calculation flag
  LOGICAL(KIND=jplm),    INTENT(IN)    :: solar(SIZE(chanprof))   ! Solar calculation flag
  INTEGER(KIND=jpim),    INTENT(IN)    :: kni, kno                ! number of levels
  TYPE(rttov_profile),   INTENT(IN)    :: ProfIn    (:)           ! atmospheric profiles
  TYPE(rttov_profile),   INTENT(IN)    :: ProfOut   (SIZE(profin))! atmospheric profiles
  TYPE(rttov_profile),   INTENT(INOUT) :: ProfOut_ad(SIZE(profin))! atmospheric profiles
  TYPE(rttov_opdp_path), INTENT(IN)    :: OpdepIn                 ! optical depths
  TYPE(rttov_opdp_path), INTENT(INOUT) :: OpdepIn_ad
  TYPE(rttov_opdp_path), INTENT(INOUT) :: OpdepOut_ad
!INTF_END
#include "rttov_layeravg.interface"
#include "rttov_layeravg_ad.interface"
! --- local scalars and arrays
!
  INTEGER(KIND=jpim) :: ji, jo, kn, istart(kno, SIZE(ProfIn)), iend(kno, SIZE(ProfIn))
  INTEGER(KIND=jpim) :: prof
  INTEGER(KIND=jpim) :: nprofiles
  INTEGER(KIND=jpim) :: nchannels
  LOGICAL(KIND=jplm) :: llevels_different
  LOGICAL(KIND=jplm) :: top_flag(SIZE(ProfIn))
  REAL(KIND=jprb)    :: pvlev(kni)
  REAL(KIND=jprb)    :: zlnpi(kni), zlnpi2(kni-1)
  REAL(KIND=jprb)    :: zlnpo(kno), zlnpo2(kno-1)
  REAL   (KIND=jprb), ALLOCATABLE :: zpz(:,:,:)
  REAL(KIND=jprb)    :: zdp(kno)
  REAL(KIND=jprb)    :: ppo(kno), ppo_ad(kno,SIZE(ProfOut))
  REAL(KIND=jprb)    :: zdp_ad(kno), zlnpo_ad(kno), zlnpi_ad(kni)
  REAL(KIND=jprb)    :: zlnpo2_ad(kno-1), zlnpi2_ad(kni-1)
  REAL   (KIND=jprb), ALLOCATABLE :: zpz_ad(:,:,:)
  REAL(KIND=jprb)    :: doddpi(kni-1), doddpo(kno-1)
  REAL(KIND=jprb)    :: doddpi_ad(kni-1), doddpo_ad(kno-1)
  INTEGER(KIND=jpim) :: interp_mode       ! Rochon or log-linear
  LOGICAL(KIND=jplm) :: interp_weightfn   ! Interpolate optical depths or weighting functions
  REAL(KIND=jprb)    :: ZHOOK_HANDLE
!----------------------------------------------------------------------
!--------------------------------------------------------------------------
! initialization
!--------------------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_INTAVG_CHAN_AD', 0_jpim, ZHOOK_HANDLE)
  llevels_different = .FALSE.
  top_flag(:)       = .FALSE.
  nprofiles         = SIZE(ProfIn)
  nchannels         = SIZE(chanprof)
  IF (opts%interpolation%lgradp) THEN
    llevels_different = .TRUE.
  ELSE
    DO prof = 1, nprofiles - 1
      DO jo = 1, kno
        IF (ABS(ProfOut(prof)%p(jo) - ProfOut(prof + 1)%p(jo)) > realtol) THEN
          llevels_different = .TRUE.
          EXIT
        ENDIF
      ENDDO
    ENDDO
  ENDIF

  ! Set interpolation mode: default to Rochon for invalid mode
  ! For Rochon/loglinear, we use loglinear for optical depth interpolation

  interp_mode = interp_rochon
  interp_weightfn = .FALSE.
  IF (opts%interpolation%interp_mode == interp_loglinear .OR. &
      opts%interpolation%interp_mode == interp_rochon_loglinear .OR. &
      opts%interpolation%interp_mode == interp_rochon_loglinear_wfn) THEN
    interp_mode = interp_loglinear
  ENDIF

  IF (opts%interpolation%interp_mode == interp_rochon_wfn .OR. &
      opts%interpolation%interp_mode == interp_rochon_loglinear_wfn) THEN
    ALLOCATE(zpz(kni-1, kno-1, nprofiles))
    IF (opts%interpolation%lgradp) ALLOCATE(zpz_ad(kni-1, kno-1, nprofiles))
    interp_weightfn = .TRUE.
    ! If spacetop is false and the output profile extends above the input
    ! profile, set a flag to move the top input pressure to the top output
    ! pressure. This is a convenient trick to make the interpolation work nicely.
    IF (.NOT. opts%interpolation%spacetop) THEN
      DO prof = 1, nprofiles
        top_flag(prof) = ProfIn(prof)%p(1) > ProfOut(prof)%p(1)
      ENDDO
    ENDIF
  ELSE
    ALLOCATE(zpz(kni, kno, nprofiles))
    IF (opts%interpolation%lgradp) ALLOCATE(zpz_ad(kni, kno, nprofiles))
  ENDIF

  DO prof = 1, nprofiles
    IF (llevels_different .OR. prof == 1) THEN
! source levels
      pvlev(1:kni) = ProfIn(prof)%p(1:kni)
      IF (interp_weightfn) THEN
        IF (top_flag(prof)) pvlev(1) = ProfOut(prof)%p(1)
      ENDIF
      zlnpi(1:kni) = LOG(pvlev(1:kni))
! target levels
      ppo(1:kno)   = ProfOut(prof)%p(1:kno)
      zlnpo(1:kno) = LOG(ppo(1:kno))
! interpolation
!--------------------
      IF (interp_weightfn) THEN
        zlnpi2(1:kni-1) = 0.5_jprb * (zlnpi(1:kni-1) + zlnpi(2:kni))
        zlnpo2(1:kno-1) = 0.5_jprb * (zlnpo(1:kno-1) + zlnpo(2:kno))
        CALL rttov_layeravg(           &
              & zlnpo2,                &
              & zlnpi2,                &
              & kno-1_jpim,            &
              & kni-1_jpim,            &
              & zpz(:, :, prof),       &
              & istart(1:kno-1, prof), &
              & iend(1:kno-1, prof),   &
              & interp_mode)
      ELSE
        CALL rttov_layeravg( &
              & zlnpo,           &
              & zlnpi,           &
              & kno,             &
              & kni,             &
              & zpz(:, :, prof), &
              & istart(:, prof), &
              & iend(:, prof),   &
              & interp_mode)
      ENDIF
! The layer averaging tends towards a constant value when target levels are
! below the lowest RTTOV level. We need instead to extrapolate optical depths
! so as replicate behaviour of rttov_transmit.F90 when RTTOV levels are used
! throughout.
      IF (interp_weightfn) THEN
! When interpolating d(opdep)/dp we want to extend the gradient at a constant
! value giving a constant change in optical depth wrt pressure.
        WHERE (zlnpo2(:) > zlnpi2(kni-1))
          zpz(kni-1, 1:kno-1, prof)   = 1.0_jprb
          istart(1:kno-1, prof)       = kni - 1
          iend(1:kno-1, prof)         = kni - 1
        ENDWHERE
      ELSE
        WHERE (ppo(:) > pvlev(kni))
          zdp(:)                = (ppo(:) - pvlev(kni)) / (pvlev(kni) - pvlev(kni - 1))
          zpz(kni - 1, :, prof) =  - zdp(:)
          zpz(kni, :, prof)     = zdp(:) + 1.0_jprb
          istart(:, prof)       = kni - 1
          iend(:, prof)         = kni
        ENDWHERE
      ENDIF
    ENDIF
  ENDDO

  IF (opts%interpolation%lgradp) THEN
    zpz_ad = 0.0_jprb
    ppo_ad = 0.0_jprb
  ENDIF

  DO kn = 1, nchannels! loop over channels for profile passed in this time
    IF (llevels_different) THEN
      prof = chanprof(kn)%prof
    ELSE
      prof = 1
    ENDIF

    IF (interp_weightfn) THEN

      pvlev(1:kni) = ProfIn(prof)%p(1:kni)
      IF (top_flag(prof)) pvlev(1) = ProfOut(prof)%p(1)
      ppo(1:kno)   = ProfOut(prof)%p(1:kno)

      IF (thermal(kn)) THEN

        IF (opts%interpolation%lgradp) THEN
          doddpi(1:kni-1) = ((OpdepIn%atm_level(1:kni-1, kn)) - (OpdepIn%atm_level(2:kni, kn))) / &
                            (pvlev(1:kni-1) - pvlev(2:kni))

          DO jo = 1, kno-1! o/p levels
      ! for viewing path of radiation generated by atmosphere
            doddpo(jo) = &
              & SUM(zpz(istart(jo, prof):iend(jo, prof), jo, prof) * doddpi(istart(jo, prof):iend(jo, prof)))
          ENDDO
        ENDIF

        doddpi_ad(:) = 0._jprb
        doddpo_ad(:) = 0._jprb

        DO jo = kno, 1, -1! o/p levels
    ! for viewing path of radiation generated by atmosphere

          IF (jo == 1) THEN
            IF (pvlev(1) >= ppo(1)) THEN
              OpdepIn_ad%atm_level(1, kn) = OpdepIn_ad%atm_level(1, kn) + OpdepOut_ad%atm_level(jo, kn)
            ELSE
              ! Find the input levels straddling the top-most output level
              DO ji = 2, kni
                IF (pvlev(ji) > ppo(1)) EXIT
              ENDDO
              ! Log-linear interpolation of OpdepIn
              OpdepIn_ad%atm_level(ji-1, kn) = OpdepIn_ad%atm_level(ji-1, kn) + &
                                               OpdepOut_ad%atm_level(jo, kn) * LOG(pvlev(ji)/ppo(1)) / &
                                               LOG(pvlev(ji)/pvlev(ji-1))
              OpdepIn_ad%atm_level(ji, kn) = OpdepIn_ad%atm_level(ji, kn) + &
                                             OpdepOut_ad%atm_level(jo, kn) * LOG(ppo(1)/pvlev(ji-1)) / &
                                             LOG(pvlev(ji)/pvlev(ji-1))
              IF (opts%interpolation%lgradp) THEN
                ppo_ad(1,prof) = ppo_ad(1,prof) + (OpdepIn%atm_level(ji, kn) - OpdepIn%atm_level(ji-1, kn)) * &
                                                  (OpdepOut_ad%atm_level(jo, kn) / ppo(1)) / &
                                                  LOG(pvlev(ji)/pvlev(ji-1))
              ENDIF
            ENDIF
          ELSE
            IF (opts%interpolation%lgradp) THEN
              ppo_ad(jo-1,prof) = ppo_ad(jo-1,prof) - doddpo(jo-1) * OpdepOut_ad%atm_level(jo, kn)
              ppo_ad(jo,prof)   = ppo_ad(jo,prof)   + doddpo(jo-1) * OpdepOut_ad%atm_level(jo, kn)
            ENDIF

            OpdepOut_ad%atm_level(jo-1, kn) = OpdepOut_ad%atm_level(jo-1, kn) + OpdepOut_ad%atm_level(jo, kn)
            doddpo_ad(jo-1) = doddpo_ad(jo-1) - OpdepOut_ad%atm_level(jo, kn) * (ppo(jo-1) - ppo(jo))
          ENDIF

          IF (jo < kno) THEN
            IF (opts%interpolation%lgradp) THEN
              zpz_ad(istart(jo, prof):iend(jo, prof), jo, prof) = &
                zpz_ad(istart(jo, prof):iend(jo, prof), jo, prof) + &
                doddpo_ad(jo) * doddpi(istart(jo, prof):iend(jo, prof))
            ENDIF

            doddpi_ad(istart(jo, prof):iend(jo, prof)) = doddpi_ad(istart(jo, prof):iend(jo, prof)) + &
              & zpz(istart(jo, prof):iend(jo, prof), jo, prof) * doddpo_ad(jo)
          ENDIF

        ENDDO

        OpdepIn_ad%atm_level(1:kni-1, kn) = OpdepIn_ad%atm_level(1:kni-1, kn) + &
                                            doddpi_ad(1:kni-1) / (pvlev(1:kni-1) - pvlev(2:kni))
        OpdepIn_ad%atm_level(2:kni, kn)   = OpdepIn_ad%atm_level(2:kni, kn) - &
                                            doddpi_ad(1:kni-1) / (pvlev(1:kni-1) - pvlev(2:kni))
      ENDIF


    ! for solar beam
    ! weights for any channel of sun_level same as for any channel of atm_level
    ! nb there could be a problem with variable incidence angles from raytracing
      IF (solar(kn)) THEN

        IF (opts%interpolation%lgradp) THEN
          doddpi(1:kni-1) = ((OpdepIn%sun_level_path2(1:kni-1, kn)) - (OpdepIn%sun_level_path2(2:kni, kn))) / &
                            (pvlev(1:kni-1) - pvlev(2:kni))

          DO jo = 1, kno-1! o/p levels
            doddpo(jo) = &
              & SUM(zpz(istart(jo, prof):iend(jo, prof), jo, prof) * doddpi(istart(jo, prof):iend(jo, prof)))
          ENDDO
        ENDIF

        doddpi_ad(:) = 0._jprb
        doddpo_ad(:) = 0._jprb

        DO jo = kno, 1, -1! o/p levels

          IF (jo == 1) THEN
            IF (pvlev(1) >= ppo(1)) THEN
              OpdepIn_ad%sun_level_path2(1, kn) = OpdepIn_ad%sun_level_path2(1, kn) + &
                                                  OpdepOut_ad%sun_level_path2(jo, kn)
            ELSE
              ! Find the input levels straddling the top-most output level
              DO ji = 2, kni
                IF (pvlev(ji) > ppo(1)) EXIT
              ENDDO
              ! Log-linear interpolation of OpdepIn
              OpdepIn_ad%sun_level_path2(ji-1, kn) = OpdepIn_ad%sun_level_path2(ji-1, kn) + &
                                               OpdepOut_ad%sun_level_path2(jo, kn) * LOG(pvlev(ji)/ppo(1)) / &
                                               LOG(pvlev(ji)/pvlev(ji-1))
              OpdepIn_ad%sun_level_path2(ji, kn) = OpdepIn_ad%sun_level_path2(ji, kn) + &
                                             OpdepOut_ad%sun_level_path2(jo, kn) * LOG(ppo(1)/pvlev(ji-1)) / &
                                             LOG(pvlev(ji)/pvlev(ji-1))
              IF (opts%interpolation%lgradp) THEN
                ppo_ad(1,prof) = ppo_ad(1,prof) + &
                                 (OpdepIn%sun_level_path2(ji, kn) - OpdepIn%sun_level_path2(ji-1, kn)) * &
                                 (OpdepOut_ad%sun_level_path2(jo, kn) / ppo(1)) / &
                                 LOG(pvlev(ji)/pvlev(ji-1))
              ENDIF
            ENDIF
          ELSE

            IF (opts%interpolation%lgradp) THEN
              ppo_ad(jo-1,prof) = ppo_ad(jo-1,prof) - doddpo(jo-1) * OpdepOut_ad%sun_level_path2(jo, kn)
              ppo_ad(jo,prof)   = ppo_ad(jo,prof)   + doddpo(jo-1) * OpdepOut_ad%sun_level_path2(jo, kn)
            ENDIF

            OpdepOut_ad%sun_level_path2(jo-1, kn) = OpdepOut_ad%sun_level_path2(jo-1, kn) + &
                                                    OpdepOut_ad%sun_level_path2(jo, kn)
            doddpo_ad(jo-1) = doddpo_ad(jo-1) - OpdepOut_ad%sun_level_path2(jo, kn) * (ppo(jo-1) - ppo(jo))
          ENDIF

          IF (jo < kno) THEN
            IF (opts%interpolation%lgradp) THEN
              zpz_ad(istart(jo, prof):iend(jo, prof), jo, prof) = &
                zpz_ad(istart(jo, prof):iend(jo, prof), jo, prof) + &
                doddpo_ad(jo) * doddpi(istart(jo, prof):iend(jo, prof))
            ENDIF

            doddpi_ad(istart(jo, prof):iend(jo, prof)) = doddpi_ad(istart(jo, prof):iend(jo, prof)) + &
              & zpz(istart(jo, prof):iend(jo, prof), jo, prof) * doddpo_ad(jo)
          ENDIF

        ENDDO

        OpdepIn_ad%sun_level_path2(1:kni-1, kn) = OpdepIn_ad%sun_level_path2(1:kni-1, kn) + &
                                            doddpi_ad(1:kni-1) / (pvlev(1:kni-1) - pvlev(2:kni))
        OpdepIn_ad%sun_level_path2(2:kni, kn)   = OpdepIn_ad%sun_level_path2(2:kni, kn) - &
                                            doddpi_ad(1:kni-1) / (pvlev(1:kni-1) - pvlev(2:kni))
      ENDIF

    ELSE

      IF (opts%interpolation%lgradp) THEN
        DO jo = 1, kno! o/p levels
  ! for viewing path of radiation generated by atmosphere
  ! nb for solar beam
  ! weights for any channel of sun_level same as for any channel of atm_level
  ! nb there could be a problem with variable incidence angles from raytracing
          IF (thermal(kn)) zpz_ad(istart(jo, prof):iend(jo, prof), jo, prof) = &
            & zpz_ad(istart(jo, prof):iend(jo, prof), jo, prof) +      &
            & OpdepIn%atm_level(istart(jo, prof):iend(jo, prof), kn) * OpdepOut_ad%atm_level(jo, kn)
          IF (solar(kn)) THEN
            zpz_ad(istart(jo, prof):iend(jo, prof), jo, prof) =   &
              zpz_ad(istart(jo, prof):iend(jo, prof), jo, prof) + &
              OpdepIn%sun_level_path2(istart(jo, prof):iend(jo, prof), kn) * OpdepOut_ad%sun_level_path2(jo, kn)
          ENDIF
        ENDDO  ! target levels
      ENDIF

      DO jo = 1, kno
  ! for fwd model weights are given by zpz (ignore zps and zpzps)
        IF (thermal(kn)) OpdepIn_ad%atm_level(istart(jo, prof):iend(jo, prof), kn) =      &
          & OpdepIn_ad%atm_level(istart(jo, prof):iend(jo, prof), kn) +  &
          & zpz(istart(jo, prof):iend(jo, prof), jo, prof) * OpdepOut_ad%atm_level(jo, kn)
        IF (solar(kn)) THEN
          OpdepIn_ad%sun_level_path2(istart(jo, prof):iend(jo, prof), kn) =   &
            OpdepIn_ad%sun_level_path2(istart(jo, prof):iend(jo, prof), kn) + &
            zpz(istart(jo, prof):iend(jo, prof), jo, prof) * OpdepOut_ad%sun_level_path2(jo, kn)
        ENDIF
      ENDDO
  ! target levels
    ENDIF
  ENDDO
! channels
  IF (opts%interpolation%lgradp) THEN
    DO prof = 1, nprofiles
! source levels
      pvlev(1:kni) = ProfIn(prof)%p(1:kni)
      IF (interp_weightfn) THEN
        IF (top_flag(prof)) pvlev(1) = ProfOut(prof)%p(1)
      ENDIF
      zlnpi(1:kni) = LOG(pvlev(1:kni))
! target levels
      ppo(1:kno)   = ProfOut(prof)%p(1:kno)
      zlnpo(1:kno) = LOG(ppo(1:kno))

      IF (interp_weightfn) THEN
        zlnpi2(1:kni-1) = 0.5_jprb * (zlnpi(1:kni-1) + zlnpi(2:kni))
        zlnpo2(1:kno-1) = 0.5_jprb * (zlnpo(1:kno-1) + zlnpo(2:kno))
        WHERE (zlnpo2(:) > zlnpi2(kni-1))
          zpz_ad(kni-1, :, prof) = 0._jprb
        ENDWHERE
      ELSE
        WHERE (ppo(:) > pvlev(kni))
          zdp_ad(:)                = zpz_ad(kni, :, prof)
          zpz_ad(kni, :, prof)     = 0.0_jprb
          zdp_ad(:)                = zdp_ad(:) - zpz_ad(kni - 1, :, prof)
          zpz_ad(kni - 1, :, prof) = 0.0_jprb
          ppo_ad(:,prof)           = ppo_ad(:,prof) + zdp_ad(:) / (pvlev(kni) - pvlev(kni - 1))
        ENDWHERE
      ENDIF

      zlnpi_ad = 0.0_jprb
      zlnpo_ad = 0.0_jprb
      IF (interp_weightfn) THEN
        zlnpi2_ad = 0.0_jprb
        zlnpo2_ad = 0.0_jprb
        CALL rttov_layeravg_ad(        &
              & zlnpo2,                &
              & zlnpo2_ad,             &
              & zlnpi2,                &
              & zlnpi2_ad,             &
              & kno-1_jpim,            &
              & kni-1_jpim,            &
              & zpz(:, :, prof),       &
              & zpz_ad(:, :, prof),    &
              & istart(1:kno-1, prof), &
              & iend(1:kno-1, prof),   &
              & interp_mode)
        zlnpo_ad(1:kno-1) = zlnpo_ad(1:kno-1) + 0.5_jprb * zlnpo2_ad(1:kno-1)
        zlnpo_ad(2:kno)   = zlnpo_ad(2:kno)   + 0.5_jprb * zlnpo2_ad(1:kno-1)
      ELSE
        CALL rttov_layeravg_ad( &
              & zlnpo,              &
              & zlnpo_ad,           &
              & zlnpi,              &
              & zlnpi_ad,           &
              & kno,                &
              & kni,                &
              & zpz(:, :, prof),    &
              & zpz_ad(:, :, prof), &
              & istart(:, prof),    &
              & iend(:, prof),      &
              & interp_mode)
      ENDIF
      ppo_ad(1:kno,prof)        = ppo_ad(1:kno,prof) + 1._jprb / ppo(1:kno) * zlnpo_ad(1:kno)
      ProfOut_ad(prof)%p(1:kno) = ProfOut_ad(prof)%p(1:kno) + ppo_ad(1:kno,prof)
    ENDDO
  ENDIF
  DEALLOCATE(zpz)
  IF (opts%interpolation%lgradp) DEALLOCATE(zpz_ad)
!
  IF (LHOOK) CALL DR_HOOK('RTTOV_INTAVG_CHAN_AD', 1_jpim, ZHOOK_HANDLE)
!
END SUBROUTINE rttov_intavg_chan_ad
