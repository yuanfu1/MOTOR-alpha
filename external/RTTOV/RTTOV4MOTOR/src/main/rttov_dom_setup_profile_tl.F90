! Description:
!> @file
!!   Set up input profiles of optical depths and SSAs for DOM algorithm - TL.
!
!> @brief
!!   Set up input profiles of optical depths and SSAs for DOM algorithm - TL.
!!
!! @details
!!   See direct model for details. The members of the TL profiles_dom array
!!   are already allocated in the associated allocation subroutine.
!!
!! @param[out]    err                       status on exit
!! @param[in]     opts                      options to configure the simulations
!! @param[in]     coef                      optical depth coefficients structure for instrument to simulate
!! @param[in]     chanprof                  specifies channels and profiles to simulate
!! @param[in]     dosolar                   flag to indicate if this allocation is for solar-enabled channels
!! @param[in]     chanflag                  flags to indicate which channels should be allocated (either channels with
!!                                          significant thermal or solar contributions)
!! @param[in]     do_rayleigh_dom           flag to indicate DOM Rayleigh simulation
!! @param[in]     nlayers                   number of input profile layers (nlevels-1)
!! @param[in]     aux_prof                  auxiliary profile variables
!! @param[in]     aux_prof_tl               auxiliary profile variable perturbations
!! @param[in]     opdp_path                 gas absorption optical depths as calculated by RTTOV
!! @param[in]     opdp_path_tl              gas absorption optical depth perturbations
!! @param[in]     angles                    information on simulation geometry
!! @param[in]     ircld                     information on cloud columns
!! @param[in]     transmission_scatt_ir     aerosol/cloud optical depths and total layer optical depths and SSAs
!! @param[in,out] transmission_scatt_ir_tl  perturbations in aerosol/cloud optical depths and total layer optical depths and SSAs
!! @param[in]     profiles_dom              array of rttov_profile_dom structures
!! @param[in,out] profiles_dom_tl           array of rttov_profile_dom structures containing perturbations
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
!    Copyright 2016, EUMETSAT, All Rights Reserved.
!
SUBROUTINE rttov_dom_setup_profile_tl(  &
              err,                      &
              opts,                     &
              coef,                     &
              chanprof,                 &
              dosolar,                  &
              chanflag,                 &
              do_rayleigh_dom,          &
              nlayers,                  &
              aux_prof,                 &
              aux_prof_tl,              &
              opdp_path,                &
              opdp_path_tl,             &
              angles,                   &
              ircld,                    &
              transmission_scatt_ir,    &
              transmission_scatt_ir_tl, &
              profiles_dom,             &
              profiles_dom_tl)
!INTF_OFF
#include "throw.h"
!INTF_ON
  USE parkind1, ONLY : jplm, jpim, jprb

  USE rttov_types, ONLY :          &
      rttov_options,               &
      rttov_coef,                  &
      rttov_chanprof,              &
      rttov_profile_aux,           &
      rttov_geometry,              &
      rttov_ircld,                 &
      rttov_transmission_scatt_ir, &
      rttov_profile_dom
!INTF_OFF
  USE rttov_const, ONLY : min_od
!INTF_ON
  IMPLICIT NONE

  INTEGER(jpim),                     INTENT(OUT)   :: err
  TYPE(rttov_options),               INTENT(IN)    :: opts
  TYPE(rttov_coef),                  INTENT(IN)    :: coef
  TYPE(rttov_chanprof),              INTENT(IN)    :: chanprof(:)
  LOGICAL(jplm),                     INTENT(IN)    :: dosolar
  LOGICAL(jplm),                     INTENT(IN)    :: chanflag(SIZE(chanprof))
  LOGICAL(jplm),                     INTENT(IN)    :: do_rayleigh_dom
  INTEGER(jpim),                     INTENT(IN)    :: nlayers
  TYPE(rttov_profile_aux),           INTENT(IN)    :: aux_prof
  TYPE(rttov_profile_aux),           INTENT(IN)    :: aux_prof_tl
  REAL(jprb),                        INTENT(IN)    :: opdp_path(:,:)
  REAL(jprb),                        INTENT(IN)    :: opdp_path_tl(:,:)
  TYPE(rttov_geometry),              INTENT(IN)    :: angles(:)
  TYPE(rttov_ircld),                 INTENT(IN)    :: ircld
  TYPE(rttov_transmission_scatt_ir), INTENT(IN)    :: transmission_scatt_ir
  TYPE(rttov_transmission_scatt_ir), INTENT(INOUT) :: transmission_scatt_ir_tl
  TYPE(rttov_profile_dom),           INTENT(IN)    :: profiles_dom(0:,:)
  TYPE(rttov_profile_dom),           INTENT(INOUT) :: profiles_dom_tl(0:,:)
!INTF_END

#include "rttov_errorreport.interface"

  INTEGER(jpim) :: i, col, coli, lay, su
  INTEGER(jpim) :: nchanprof
  INTEGER(jpim) :: prof, chan
  INTEGER(jpim) :: domnlayers                 ! Number of layers for DOM profile
  INTEGER(jpim) :: maxnlayers                 ! Lowest full layer included from input profile
  INTEGER(jpim) :: laysurf                    ! Layer containing surface
  REAL(jprb)    :: patheff_r
  REAL(jprb)    :: opdep                      ! Accumulated optical depth from TOA
  INTEGER(jpim) :: laymaptodom(nlayers)       ! Maps input layers to DOM layers

  REAL(jprb)    :: clropdep(nlayers)
  REAL(jprb)    :: clropdep_tl(nlayers)
  REAL(jprb)    :: absopdep(nlayers)
  REAL(jprb)    :: absopdep_tl(nlayers)
  REAL(jprb)    :: scaopdep(nlayers)
  REAL(jprb)    :: scaopdep_tl(nlayers)
  REAL(jprb)    :: totopdep_tl(nlayers)
  REAL(jprb)    :: surfclropdep_tl
  REAL(jprb)    :: surfabsopdep_tl
  REAL(jprb)    :: surfscaopdep
  REAL(jprb)    :: surfscaopdep_tl
  REAL(jprb)    :: surftotopdep_tl
  REAL(jprb)    :: totopdep_init              ! Nadir optical depth at top input level
  REAL(jprb)    :: totopdep_init_tl
  REAL(jprb)    :: layerod_tl(0:1,nlayers)

  REAL(jprb), POINTER :: layerod(:,:,:)
  REAL(jprb), POINTER :: ssa_tl(:,:,:)

  ! --------------------------------------------------------------------------------------

  TRY

  nchanprof = SIZE(chanprof)

  IF (dosolar) THEN
    ssa_tl => transmission_scatt_ir_tl%ssa_solar
    layerod => transmission_scatt_ir%layerod_solar
  ELSE
    ssa_tl => transmission_scatt_ir_tl%ssa_thermal
    layerod => transmission_scatt_ir%layerod_thermal
  ENDIF

  DO i = 1, nchanprof

    IF (.NOT. chanflag(i)) CYCLE
    prof = chanprof(i)%prof
    chan = chanprof(i)%chan

    ! ------------------------------------------------------------------------
    ! Calculate the clear-sky optical depths
    ! ------------------------------------------------------------------------

    IF (dosolar) THEN
      ! Remember that opdp_path contains optical depths on sun-surface-satellite path ("path2") so
      ! we need to divide by patheff = 1 / mu_eff to get the nadir optical depths. We also apply the
      ! gamma correction which is applied elsewhere (opdpscattir and transmit) to other optical depths.

      patheff_r = 1._jprb / (1._jprb / angles(prof)%coszen + 1._jprb / angles(prof)%coszen_sun)

      ! Nadir opdep at first level (no scattering particles above this level) - zero if spacetop is true.
      totopdep_init = - opdp_path(1,i) * coef%ff_gam(chan) * patheff_r
      totopdep_init_tl = - opdp_path_tl(1,i) * coef%ff_gam(chan) * patheff_r

      ! Nadir clear-sky optical depths of input layers
      clropdep(:) = -(opdp_path(2:nlayers+1,i) - opdp_path(1:nlayers,i)) * coef%ff_gam(chan) * patheff_r
      clropdep_tl(:) = -(opdp_path_tl(2:nlayers+1,i) - opdp_path_tl(1:nlayers,i)) * coef%ff_gam(chan) * patheff_r
    ELSE
      ! Nadir opdep at first level (no scattering particles above this level) - zero if spacetop is true.
      totopdep_init = - opdp_path(1,i) * coef%ff_gam(chan) * angles(prof)%coszen
      totopdep_init_tl = - opdp_path_tl(1,i) * coef%ff_gam(chan) * angles(prof)%coszen

      ! Nadir clear-sky optical depths of input layers
      clropdep(:) = -(opdp_path(2:nlayers+1,i) - opdp_path(1:nlayers,i)) * &
                     coef%ff_gam(chan) * angles(prof)%coszen
      clropdep_tl(:) = -(opdp_path_tl(2:nlayers+1,i) - opdp_path_tl(1:nlayers,i)) * &
                        coef%ff_gam(chan) * angles(prof)%coszen
    ENDIF

    ! ------------------------------------------------------------------------
    ! Calculate the non-cloudy/cloudy total layer optical depths and SSAs
    ! ------------------------------------------------------------------------

    IF (opts%rt_ir%addclouds) THEN
      su = 1
    ELSE
      su = 0
    ENDIF

    DO coli = 0, su
      layerod_tl(coli,:) = clropdep_tl(:) + &
                           transmission_scatt_ir_tl%opdpabs(coli,:,i) + &
                           transmission_scatt_ir_tl%opdpsca(coli,:,i)
    ENDDO

    ssa_tl(:,:,i) = 0._jprb
    IF (opts%rt_ir%addaerosl .OR. do_rayleigh_dom) THEN
      WHERE (layerod(0,:,i) > 0._jprb)
        ssa_tl(0,:,i) = &
            transmission_scatt_ir_tl%opdpsca(0,:,i) / layerod(0,:,i) - &
            transmission_scatt_ir%opdpsca(0,:,i) * layerod_tl(0,:) / layerod(0,:,i) ** 2_jpim
      ENDWHERE
      WHERE (transmission_scatt_ir%opdpsca(0,:,i) > layerod(0,:,i)) ssa_tl(0,:,i) = 0._jprb
    ENDIF
    IF (opts%rt_ir%addclouds) THEN
      WHERE (layerod(1,:,i) > 0._jprb)
        ssa_tl(1,:,i) = &
            transmission_scatt_ir_tl%opdpsca(1,:,i) / layerod(1,:,i) - &
            transmission_scatt_ir%opdpsca(1,:,i) * layerod_tl(1,:) / layerod(1,:,i) ** 2_jpim
      ENDWHERE
      WHERE (transmission_scatt_ir%opdpsca(1,:,i) > layerod(1,:,i)) ssa_tl(1,:,i) = 0._jprb
    ENDIF


    ! ------------------------------------------------------------------------
    ! Loop over cloud columns to determine the number of layers to pass into
    ! DOM for each column independently
    ! ------------------------------------------------------------------------
    DO col = 0, ircld%ncolumn(prof)

      ! ------------------------------------------------------------------------
      ! Calculate the combined optical parameters for all scattering particles
      ! ------------------------------------------------------------------------
      ! These are the nadir values for the user input layers

      DO lay = 1, nlayers
        coli = ircld%icldarr(col,lay,prof)
        absopdep(lay) = transmission_scatt_ir%opdpabs(coli,lay,i)
        scaopdep(lay) = transmission_scatt_ir%opdpsca(coli,lay,i)
        absopdep_tl(lay) = transmission_scatt_ir_tl%opdpabs(coli,lay,i)
        scaopdep_tl(lay) = transmission_scatt_ir_tl%opdpsca(coli,lay,i)
      ENDDO

      totopdep_tl(:) = clropdep_tl(:) + absopdep_tl(:) + scaopdep_tl(:)

      ! ------------------------------------------------------------------------
      ! Calculate number of layers
      ! ------------------------------------------------------------------------

      ! maxnlayers is the layer above the one immediately above/containing the surface
      ! (recall nearestlev_surf is the nearest *level* below or equal to surface).
      ! If 2m pressure is exactly on a pressure level, maxnlayers *excludes*
      ! the layer immediately above the surface (this is dealt with below).
      ! If 2m pressure is below the input profile, maxnlayers = nlayers - 1
      IF (aux_prof%s(prof)%pfraction_surf < 0._jprb) THEN
        ! 2m pressure below bottom of input profile
        maxnlayers = nlayers - 1
      ELSE
        maxnlayers = aux_prof%s(prof)%nearestlev_surf - 2
      ENDIF

      ! Find the lowest layer to consider based on *absorption* optical depth
      ! (assume clropdep is due to absorption, but for v9 predictors it includes
      ! Rayleigh scattering extinction for visible channels).
      ! If the profile is optically thick set surface to .FALSE.
      opdep = totopdep_init
      DO lay = 1, maxnlayers
        opdep = opdep + clropdep(lay) + absopdep(lay)

        IF (opts%rt_ir%dom_opdep_threshold > 0._jprb) THEN
          IF (opdep > opts%rt_ir%dom_opdep_threshold) THEN
            maxnlayers = lay
            EXIT
          ENDIF
        ENDIF
      ENDDO

      ! Examine optical parameters to see if we can combine layers together
      ! This combines layers with the same properties including the "spacetop" layer
      ! and any partial near-surface layer.
      ! If spacetop is true then lay=1 is associated with domnlayers=1.
      ! If spacetop is false then the "spacetop" layer is associated with domnlayers=1
      ! and lay=1 is associated with domnlayers=1 if the SSA is zero or domnlayers=2 otherwise.
      ! Then domnlayers is incremented for each additional layer down to maxnlayers
      ! only if one or both of the adjacent layer SSAs are zero.
      ! Finally the near-surface layer is considered.

      ! Input layers (lay) are mapped to DOM layers (domnlayers) via laymaptodom(:)
      DO lay = 1, maxnlayers
        IF (lay == 1) THEN
          IF (opts%interpolation%spacetop .OR. .NOT. scaopdep(lay) > 0._jprb) THEN
            ! There is no layer above first input level or top input layer has no scattering particles
            domnlayers = 1
          ELSE
            ! There is a top layer above first input level
            domnlayers = 2
          ENDIF
        ELSE

          IF (dosolar) THEN
            ! -------------------------------------------------------------------
            ! Combine consecutive clear layers: increment domnlayers if
            ! either this layer's or previous layer's SSA is non-zero
            IF (scaopdep(lay) > 0._jprb .OR. scaopdep(lay-1) > 0._jprb) THEN
              domnlayers = domnlayers + 1
            ENDIF
            ! -------------------------------------------------------------------
          ELSE
            ! -------------------------------------------------------------------
            ! Don't combine layers at all
            domnlayers = domnlayers + 1
            ! -------------------------------------------------------------------
          ENDIF
        ENDIF
        laymaptodom(lay) = domnlayers
      ENDDO ! input layers

      IF (profiles_dom(col,i)%surface) THEN
        ! Deal with the partial surface level (if there is one)
        ! This is either combined with the previous layer or appears as an additional layer
        ! If it lies below the bottom of the input profile, it is assumed clear

        laysurf = aux_prof%s(prof)%nearestlev_surf - 1

        surfscaopdep = (1._jprb - aux_prof%s(prof)%pfraction_surf) * scaopdep(laysurf)

        surfclropdep_tl = (1._jprb - aux_prof%s(prof)%pfraction_surf) * clropdep_tl(laysurf) - &
                          aux_prof_tl%s(prof)%pfraction_surf * clropdep(laysurf)
        surfabsopdep_tl = (1._jprb - aux_prof%s(prof)%pfraction_surf) * absopdep_tl(laysurf) - &
                          aux_prof_tl%s(prof)%pfraction_surf * absopdep(laysurf)
        surfscaopdep_tl = (1._jprb - aux_prof%s(prof)%pfraction_surf) * scaopdep_tl(laysurf) - &
                          aux_prof_tl%s(prof)%pfraction_surf * scaopdep(laysurf)

        surftotopdep_tl = surfclropdep_tl + surfabsopdep_tl + surfscaopdep_tl

        ! ** This currently only combines clear layers **
        ! Increase domnlayers only if near-surf layer differs to the one above
        IF ((.NOT. dosolar) .OR. surfscaopdep > 0._jprb .OR. scaopdep(maxnlayers) > 0._jprb) THEN
          domnlayers = domnlayers + 1
        ENDIF
      ENDIF

      ! ------------------------------------------------------------------------
      ! Populate arrays
      ! ------------------------------------------------------------------------

      ! This should never happen
      IF (domnlayers /= profiles_dom(col,i)%nlayers) THEN
        err = errorstatus_fatal
        THROWM(err.NE.0, 'TL - domnlayers mismatch')
      ENDIF

      ! Deal with all full layers above surface
      profiles_dom_tl(col,i)%layerod(:) = 0._jprb

      ! If there is non-zero opdep above user top level, add this to the first dom layer opdep
      ! (this layer has no scattering particles)
      profiles_dom_tl(col,i)%layerod(1) = totopdep_init_tl

      DO lay = 1, maxnlayers
        profiles_dom_tl(col,i)%layerod(laymaptodom(lay)) = profiles_dom_tl(col,i)%layerod(laymaptodom(lay)) + &
                                                           totopdep_tl(lay)
      ENDDO

      IF (profiles_dom(col,i)%surface) THEN
        ! Handle partial layer above surface - either added to last layer from above or treated as new (last) layer
        profiles_dom_tl(col,i)%layerod(domnlayers) = profiles_dom_tl(col,i)%layerod(domnlayers) + &
                                                     surftotopdep_tl
      ENDIF

      IF (.NOT. dosolar) THEN
        WHERE (profiles_dom(col,i)%layerod(:) <= min_od) profiles_dom_tl(col,i)%layerod(:) = 0._jprb
      ENDIF

    ENDDO ! cloud columns

  ENDDO ! chanprof

  CATCH
END SUBROUTINE rttov_dom_setup_profile_tl
