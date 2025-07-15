! Description:
!> @file
!!   Set up input profiles of optical depths and SSAs for DOM algorithm - AD.
!
!> @brief
!!   Set up input profiles of optical depths and SSAs for DOM algorithm - AD.
!!
!! @details
!!   See direct model for details. The members of the AD profiles_dom array
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
!! @param[in,out] aux_prof_ad               auxiliary profile variable increments
!! @param[in]     opdp_path                 gas absorption optical depths as calculated by RTTOV
!! @param[in,out] opdp_path_ad              gas absorption optical depth increments
!! @param[in]     angles                    information on simulation geometry
!! @param[in]     ircld                     information on cloud columns
!! @param[in]     transmission_scatt_ir     aerosol/cloud optical depths and total layer optical depths and SSAs
!! @param[in,out] transmission_scatt_ir_ad  increments in aerosol/cloud optical depths and total layer optical depths and SSAs
!! @param[in]     profiles_dom              array of rttov_profile_dom structures
!! @param[in,out] profiles_dom_ad           array of rttov_profile_dom structures containing increments
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
SUBROUTINE rttov_dom_setup_profile_ad(  &
              err,                      &
              opts,                     &
              coef,                     &
              chanprof,                 &
              dosolar,                  &
              chanflag,                 &
              do_rayleigh_dom,          &
              nlayers,                  &
              aux_prof,                 &
              aux_prof_ad,              &
              opdp_path,                &
              opdp_path_ad,             &
              angles,                   &
              ircld,                    &
              transmission_scatt_ir,    &
              transmission_scatt_ir_ad, &
              profiles_dom,             &
              profiles_dom_ad)
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
  TYPE(rttov_profile_aux),           INTENT(INOUT) :: aux_prof_ad
  REAL(jprb),                        INTENT(IN)    :: opdp_path(:,:)
  REAL(jprb),                        INTENT(INOUT) :: opdp_path_ad(:,:)
  TYPE(rttov_geometry),              INTENT(IN)    :: angles(:)
  TYPE(rttov_ircld),                 INTENT(IN)    :: ircld
  TYPE(rttov_transmission_scatt_ir), INTENT(IN)    :: transmission_scatt_ir
  TYPE(rttov_transmission_scatt_ir), INTENT(INOUT) :: transmission_scatt_ir_ad
  TYPE(rttov_profile_dom),           INTENT(IN)    :: profiles_dom(0:,:)
  TYPE(rttov_profile_dom),           INTENT(INOUT) :: profiles_dom_ad(0:,:)
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
  REAL(jprb)    :: clropdep_ad(nlayers)
  REAL(jprb)    :: absopdep(nlayers)
  REAL(jprb)    :: absopdep_ad(nlayers)
  REAL(jprb)    :: scaopdep(nlayers)
  REAL(jprb)    :: scaopdep_ad(nlayers)
  REAL(jprb)    :: totopdep_ad(nlayers)
  REAL(jprb)    :: surfclropdep_ad
  REAL(jprb)    :: surfabsopdep_ad
  REAL(jprb)    :: surfscaopdep
  REAL(jprb)    :: surfscaopdep_ad
  REAL(jprb)    :: surftotopdep_ad
  REAL(jprb)    :: totopdep_init              ! Nadir optical depth at top input level
  REAL(jprb)    :: totopdep_init_ad
  REAL(jprb)    :: layerod_ad(0:1,nlayers)

  REAL(jprb), POINTER :: layerod(:,:,:)
  REAL(jprb), POINTER :: ssa_ad(:,:,:)

  ! --------------------------------------------------------------------------------------

  TRY

  nchanprof = SIZE(chanprof)

  IF (dosolar) THEN
    ssa_ad => transmission_scatt_ir_ad%ssa_solar
    layerod => transmission_scatt_ir%layerod_solar
  ELSE
    ssa_ad => transmission_scatt_ir_ad%ssa_thermal
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

      ! Nadir clear-sky optical depths of input layers
      clropdep(:) = -(opdp_path(2:nlayers+1,i) - opdp_path(1:nlayers,i)) * coef%ff_gam(chan) * patheff_r
    ELSE
      ! Nadir opdep at first level (no scattering particles above this level) - zero if spacetop is true.
      totopdep_init = - opdp_path(1,i) * coef%ff_gam(chan) * angles(prof)%coszen

      ! Nadir clear-sky optical depths of input layers
      clropdep(:) = -(opdp_path(2:nlayers+1,i) - opdp_path(1:nlayers,i)) * &
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

    layerod_ad  = 0._jprb
    clropdep_ad = 0._jprb

    ! ------------------------------------------------------------------------
    ! Loop over cloud columns to determine the number of layers to pass into
    ! DOM for each column independently
    ! ------------------------------------------------------------------------
    DO col = ircld%ncolumn(prof), 0, -1

      ! ------------------------------------------------------------------------
      ! Calculate the combined optical parameters for all scattering particles
      ! ------------------------------------------------------------------------
      ! These are the nadir values for the user input layers

      DO lay = 1, nlayers
        coli = ircld%icldarr(col,lay,prof)
        absopdep(lay) = transmission_scatt_ir%opdpabs(coli,lay,i)
        scaopdep(lay) = transmission_scatt_ir%opdpsca(coli,lay,i)
      ENDDO

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
        THROWM(err.NE.0, 'AD - domnlayers mismatch')
      ENDIF


      ! Adjoint calculations

      totopdep_init_ad = 0._jprb
      absopdep_ad      = 0._jprb
      scaopdep_ad      = 0._jprb
      totopdep_ad      = 0._jprb
      surfclropdep_ad  = 0._jprb
      surfabsopdep_ad  = 0._jprb
      surfscaopdep_ad  = 0._jprb
      surftotopdep_ad  = 0._jprb

      IF (.NOT. dosolar) THEN
        WHERE (profiles_dom(col,i)%layerod(:) <= min_od) profiles_dom_ad(col,i)%layerod(:) = 0._jprb
      ENDIF

      IF (profiles_dom(col,i)%surface) THEN
        surftotopdep_ad = surftotopdep_ad + profiles_dom_ad(col,i)%layerod(domnlayers)
      ENDIF

      DO lay = maxnlayers, 1, -1
        totopdep_ad(lay) = totopdep_ad(lay) + profiles_dom_ad(col,i)%layerod(laymaptodom(lay))
      ENDDO

      totopdep_init_ad = totopdep_init_ad + profiles_dom_ad(col,i)%layerod(1)

!       profiles_dom_ad(col,i)%layerod(:) = 0._jprb


      IF (profiles_dom(col,i)%surface) THEN
        surfclropdep_ad = surfclropdep_ad + surftotopdep_ad
        surfabsopdep_ad = surfabsopdep_ad + surftotopdep_ad
        surfscaopdep_ad = surfscaopdep_ad + surftotopdep_ad

        aux_prof_ad%s(prof)%pfraction_surf = aux_prof_ad%s(prof)%pfraction_surf - &
            surfscaopdep_ad * scaopdep(laysurf)
        aux_prof_ad%s(prof)%pfraction_surf = aux_prof_ad%s(prof)%pfraction_surf - &
            surfabsopdep_ad * absopdep(laysurf)
        aux_prof_ad%s(prof)%pfraction_surf = aux_prof_ad%s(prof)%pfraction_surf - &
            surfclropdep_ad * clropdep(laysurf)

        scaopdep_ad(laysurf) = scaopdep_ad(laysurf) + &
            surfscaopdep_ad * (1._jprb - aux_prof%s(prof)%pfraction_surf)
        absopdep_ad(laysurf) = absopdep_ad(laysurf) + &
            surfabsopdep_ad * (1._jprb - aux_prof%s(prof)%pfraction_surf)
        clropdep_ad(laysurf) = clropdep_ad(laysurf) + &
            surfclropdep_ad * (1._jprb - aux_prof%s(prof)%pfraction_surf)
      ENDIF

      clropdep_ad = clropdep_ad + totopdep_ad
      absopdep_ad = absopdep_ad + totopdep_ad
      scaopdep_ad = scaopdep_ad + totopdep_ad

      DO lay = nlayers, 1, -1
        coli = ircld%icldarr(col,lay,prof)

        transmission_scatt_ir_ad%opdpsca(coli,lay,i) = &
            transmission_scatt_ir_ad%opdpsca(coli,lay,i) + scaopdep_ad(lay)
        transmission_scatt_ir_ad%opdpabs(coli,lay,i) = &
            transmission_scatt_ir_ad%opdpabs(coli,lay,i) + absopdep_ad(lay)
      ENDDO

    ENDDO ! cloud columns


    ! ------------------------------------------------------------------------
    ! Calculate the non-cloudy/cloudy total layer optical depths and SSAs
    ! ------------------------------------------------------------------------

    IF (opts%rt_ir%addaerosl .OR. do_rayleigh_dom) THEN
      WHERE (transmission_scatt_ir%opdpsca(0,:,i) > layerod(0,:,i)) ssa_ad(0,:,i) = 0._jprb
      WHERE (layerod(0,:,i) > 0._jprb)
        layerod_ad(0,:) = layerod_ad(0,:) - &
            ssa_ad(0,:,i) * transmission_scatt_ir%opdpsca(0,:,i) / layerod(0,:,i) ** 2_jpim
        transmission_scatt_ir_ad%opdpsca(0,:,i) = &
            transmission_scatt_ir_ad%opdpsca(0,:,i) + ssa_ad(0,:,i) / layerod(0,:,i)
      ENDWHERE
    ENDIF
    IF (opts%rt_ir%addclouds) THEN
      WHERE (transmission_scatt_ir%opdpsca(1,:,i) > layerod(1,:,i)) ssa_ad(1,:,i) = 0._jprb
      WHERE (layerod(1,:,i) > 0._jprb)
        layerod_ad(1,:) = layerod_ad(1,:) - &
            ssa_ad(1,:,i) * transmission_scatt_ir%opdpsca(1,:,i) / layerod(1,:,i) ** 2_jpim
        transmission_scatt_ir_ad%opdpsca(1,:,i) = &
            transmission_scatt_ir_ad%opdpsca(1,:,i) + ssa_ad(1,:,i) / layerod(1,:,i)
      ENDWHERE
    ENDIF
!     ssa_ad(:,:,i) = 0._jprb

    DO coli = 0, su
      transmission_scatt_ir_ad%opdpsca(coli,:,i) = &
          transmission_scatt_ir_ad%opdpsca(coli,:,i) + layerod_ad(coli,:)
      transmission_scatt_ir_ad%opdpabs(coli,:,i) = &
          transmission_scatt_ir_ad%opdpabs(coli,:,i) + layerod_ad(coli,:)
      clropdep_ad(:) = clropdep_ad(:) + layerod_ad(coli,:)
    ENDDO


    ! ------------------------------------------------------------------------
    ! Calculate the clear-sky optical depths
    ! ------------------------------------------------------------------------

    IF (dosolar) THEN
      opdp_path_ad(2:nlayers+1,i) = opdp_path_ad(2:nlayers+1,i) - &
          clropdep_ad(:) * coef%ff_gam(chan) * patheff_r
      opdp_path_ad(1:nlayers,i) = opdp_path_ad(1:nlayers,i) + &
          clropdep_ad(:) * coef%ff_gam(chan) * patheff_r

      opdp_path_ad(1,i) = opdp_path_ad(1,i) - &
          totopdep_init_ad * coef%ff_gam(chan) * patheff_r
    ELSE
      opdp_path_ad(2:nlayers+1,i) = opdp_path_ad(2:nlayers+1,i) - &
          clropdep_ad(:) * coef%ff_gam(chan) * angles(prof)%coszen
      opdp_path_ad(1:nlayers,i) = opdp_path_ad(1:nlayers,i) + &
          clropdep_ad(:) * coef%ff_gam(chan) * angles(prof)%coszen

      opdp_path_ad(1,i) = opdp_path_ad(1,i) - &
          totopdep_init_ad * coef%ff_gam(chan) * angles(prof)%coszen
    ENDIF

  ENDDO ! chanprof

  CATCH
END SUBROUTINE rttov_dom_setup_profile_ad
