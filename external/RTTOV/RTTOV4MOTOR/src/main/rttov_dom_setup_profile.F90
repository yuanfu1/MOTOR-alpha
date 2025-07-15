! Description:
!> @file
!!   Set up input profiles of optical depths and SSAs for DOM algorithm.
!
!> @brief
!!   Set up input profiles of optical depths and SSAs for DOM algorithm.
!!
!! @details
!!   The optical depths for the DOM algorithm are stored in the profiles_dom
!!   structure while the single-scattering albedos (SSAs) are stored in the
!!   trans_scatt_ir structure. The total layer optical depths (including
!!   aerosols and clouds) on the input layers are also stored in the
!!   trans_scatt_ir array for re-use in the TL/AD/K.
!!
!!   The DOM algorithm run-time is proportional to the number of layers on
!!   which it is run. The algorithm is called separately for to treat the 
!!   solar source term and the emission source term. For solar simulations
!!   any adjacent clear layers may be combined together in the profiles_dom
!!   structure to decrease run-time without any impact on radiances. If DOM
!!   Rayleigh scattering is enabled, then affected layers cannot be combined
!!   and the run-time will increase significantly.
!!
!!   For thermal simulations, combining layers is not possible because the
!!   emission source term depends on the layer temperature. However for some IR
!!   channels if the absorption optical depth reaches some user-defined
!!   threshold then the remainder of the profile below this can be ignored
!!   which can reduce run-time. (This is applied to solar simulations also
!!   but rarely if ever would be applied in practice as the channels are
!!   relatively transparent).
!!
!!   The number of layers for each DOM profile is determined dynamically
!!   here for each cloud column and channel so the members of the
!!   profiles_dom array are allocated here rather than in the profiles_dom
!!   allocation subroutine (for the direct model).
!!
!! @param[in]     opts                   options to configure the simulations
!! @param[in]     coef                   optical depth coefficients structure for instrument to simulate
!! @param[in]     chanprof               specifies channels and profiles to simulate
!! @param[in]     dosolar                flag to indicate if this allocation is for solar-enabled channels
!! @param[in]     chanflag               flags to indicate which channels should be allocated (either channels with
!!                                       significant thermal or solar contributions)
!! @param[in]     do_rayleigh_dom        flag to indicate DOM Rayleigh simulation
!! @param[in]     nlayers                number of input profile layers (nlevels-1)
!! @param[in]     aux_prof               auxiliary profile variables
!! @param[in]     opdp_path              gas absorption optical depths as calculated by RTTOV
!! @param[in]     angles                 information on simulation geometry
!! @param[in]     ircld                  information on cloud columns
!! @param[in,out] transmission_scatt_ir  input aerosol/cloud optical depths and output total layer optical depths and SSAs
!! @param[in,out] profiles_dom           array of rttov_profile_dom structures
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
SUBROUTINE rttov_dom_setup_profile(  &
              opts,                  &
              coef,                  &
              chanprof,              &
              dosolar,               &
              chanflag,              &
              do_rayleigh_dom,       &
              nlayers,               &
              aux_prof,              &
              opdp_path,             &
              angles,                &
              ircld,                 &
              transmission_scatt_ir, &
              profiles_dom)

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

  TYPE(rttov_options),               INTENT(IN)    :: opts
  TYPE(rttov_coef),                  INTENT(IN)    :: coef
  TYPE(rttov_chanprof),              INTENT(IN)    :: chanprof(:)
  LOGICAL(jplm),                     INTENT(IN)    :: dosolar
  LOGICAL(jplm),                     INTENT(IN)    :: chanflag(SIZE(chanprof))
  LOGICAL(jplm),                     INTENT(IN)    :: do_rayleigh_dom
  INTEGER(jpim),                     INTENT(IN)    :: nlayers
  TYPE(rttov_profile_aux),           INTENT(IN)    :: aux_prof
  REAL(jprb),                        INTENT(IN)    :: opdp_path(:,:)
  TYPE(rttov_geometry),              INTENT(IN)    :: angles(:)
  TYPE(rttov_ircld),                 INTENT(IN)    :: ircld
  TYPE(rttov_transmission_scatt_ir), INTENT(INOUT) :: transmission_scatt_ir
  TYPE(rttov_profile_dom),           INTENT(INOUT) :: profiles_dom(0:,:)
!INTF_END

  INTEGER(jpim) :: i, col, coli, lay, su
  INTEGER(jpim) :: nchanprof
  INTEGER(jpim) :: prof, chan
  INTEGER(jpim) :: domnlayers                 ! Number of layers for DOM profile
  INTEGER(jpim) :: maxnlayers                 ! Lowest full layer included from input profile
  INTEGER(jpim) :: laysurf                    ! Layer containing surface
  REAL(jprb)    :: patheff_r
  REAL(jprb)    :: opdep                      ! Accumulated optical depth from TOA
  INTEGER(jpim) :: laymaptodom(nlayers)       ! Maps input layers to DOM layers
  INTEGER(jpim) :: laymaptouser(nlayers+1)    ! Maps DOM layers to input layers

  REAL(jprb)    :: clropdep(nlayers)          ! Layer nadir gas absorption optical depth
  REAL(jprb)    :: absopdep(nlayers)          ! Layer nadir aer/cld abs optical depth
  REAL(jprb)    :: scaopdep(nlayers)          ! Layer nadir aer/cld sca optical depth
  REAL(jprb)    :: totopdep(nlayers)          ! Total nadir layer optical depth
  REAL(jprb)    :: surfclropdep               ! As above but for near-surface partial layer
  REAL(jprb)    :: surfabsopdep               !
  REAL(jprb)    :: surfscaopdep               !
  REAL(jprb)    :: surftotopdep               !
  REAL(jprb)    :: totopdep_init              ! Nadir optical depth at top input level

  REAL(jprb), POINTER :: ssa(:,:,:), layerod(:,:,:)

  ! --------------------------------------------------------------------------------------

! This is the structure which needs populating:
! profiles_dom(maxncolumns,nchanprof):  (maxncolumns refers to cloud columns)
!   %nlayers                      # number of layers for this channel/profile/cloud column
!   %surface                      # flag to indicate if surface reflectance should be included
!   %layerod(nlayers)             # single layer nadir total opdeps
!   %laymap(nlayers)              # index of user layer for each DOM layer containing scattering material

  nchanprof = SIZE(chanprof)
  profiles_dom(:,:)%nlayers = 0

  IF (dosolar) THEN
    ssa => transmission_scatt_ir%ssa_solar
    layerod => transmission_scatt_ir%layerod_solar
  ELSE
    ssa => transmission_scatt_ir%ssa_thermal
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

    DO coli = 0, su
      layerod(coli,:,i) = clropdep(:) + &
                          transmission_scatt_ir%opdpabs(coli,:,i) + &
                          transmission_scatt_ir%opdpsca(coli,:,i)
    ENDDO

    ssa(:,:,i) = 0._jprb
    IF (opts%rt_ir%addaerosl .OR. do_rayleigh_dom) THEN
      WHERE (layerod(0,:,i) > 0._jprb)
        ssa(0,:,i) = transmission_scatt_ir%opdpsca(0,:,i) / layerod(0,:,i)
      ENDWHERE
    ENDIF
    IF (opts%rt_ir%addclouds) THEN
      WHERE (layerod(1,:,i) > 0._jprb)
        ssa(1,:,i) = transmission_scatt_ir%opdpsca(1,:,i) / layerod(1,:,i)
      ENDWHERE
    ENDIF
    WHERE (ssa(:,:,i) > 1._jprb) ssa(:,:,i) = 1._jprb

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
      ENDDO

      totopdep(:) = clropdep(:) + absopdep(:) + scaopdep(:)

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
      profiles_dom(col,i)%surface = .TRUE.
      DO lay = 1, maxnlayers
        opdep = opdep + clropdep(lay) + absopdep(lay)

        IF (opts%rt_ir%dom_opdep_threshold > 0._jprb) THEN
          IF (opdep > opts%rt_ir%dom_opdep_threshold) THEN
            profiles_dom(col,i)%surface = .FALSE.
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

      ! This array holds the mapping from DOM layers to user layers: only scattering layers
      ! are considered; clear layers are set to zero because they are combined.
      laymaptouser(:) = 0

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
        IF (scaopdep(lay) > 0._jprb) laymaptouser(domnlayers) = lay
      ENDDO ! input layers

      IF (profiles_dom(col,i)%surface) THEN
        ! Deal with the layer immediately above the surface:
        ! take the fraction of the optical depths in the layer according to pfraction_surf.
        ! If the surface lies below the bottom of the input profile this extends the lowest
        ! layer down to the surface.
        ! If the surface lies exactly on a pressure level this is OK because we excluded
        ! the layer above the surface when computing maxnlayers above.

        laysurf = aux_prof%s(prof)%nearestlev_surf - 1

        surfclropdep = (1._jprb - aux_prof%s(prof)%pfraction_surf) * clropdep(laysurf)
        surfabsopdep = (1._jprb - aux_prof%s(prof)%pfraction_surf) * absopdep(laysurf)
        surfscaopdep = (1._jprb - aux_prof%s(prof)%pfraction_surf) * scaopdep(laysurf)

        surftotopdep = surfclropdep + surfabsopdep + surfscaopdep

        ! ** This currently only combines clear layers **
        ! Increase domnlayers only if near-surf layer differs to the one above
        IF ((.NOT. dosolar) .OR. surfscaopdep > 0._jprb .OR. scaopdep(maxnlayers) > 0._jprb) THEN
          domnlayers = domnlayers + 1
          IF (surfscaopdep > 0._jprb) laymaptouser(domnlayers) = laysurf
        ENDIF
      ENDIF

      ! domnlayers is the number of layers for the DOM calculation for this column/channel
      profiles_dom(col,i)%nlayers = domnlayers

      ! ------------------------------------------------------------------------
      ! Allocate DOM profile arrays
      ! ------------------------------------------------------------------------

      ALLOCATE(profiles_dom(col,i)%layerod(domnlayers), &
               profiles_dom(col,i)%laymap(domnlayers))

      ! ------------------------------------------------------------------------
      ! Populate arrays
      ! ------------------------------------------------------------------------

      profiles_dom(col,i)%laymap(:) = laymaptouser(1:domnlayers)

      ! Deal with all full layers above surface
      profiles_dom(col,i)%layerod(:) = 0._jprb

      ! If there is non-zero opdep above user top level, add this to the first dom layer opdep
      ! (this layer has no scattering particles)
      profiles_dom(col,i)%layerod(1) = totopdep_init

      DO lay = 1, maxnlayers
        profiles_dom(col,i)%layerod(laymaptodom(lay)) = profiles_dom(col,i)%layerod(laymaptodom(lay)) + &
                                                        totopdep(lay)
      ENDDO

      IF (profiles_dom(col,i)%surface) THEN
        ! Handle partial layer above surface - either added to last layer from above or treated as new (last) layer
        profiles_dom(col,i)%layerod(domnlayers) = profiles_dom(col,i)%layerod(domnlayers) + &
                                                  surftotopdep
      ENDIF

      IF (.NOT. dosolar) THEN
        WHERE (profiles_dom(col,i)%layerod(:) < min_od) profiles_dom(col,i)%layerod(:) = min_od
      ENDIF

    ENDDO ! cloud columns

  ENDDO ! chanprof

END SUBROUTINE rttov_dom_setup_profile
