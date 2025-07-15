! Description:
!> @file
!!   Compute cloud columns using selected cloud overlap scheme for 
!!   visible/IR scattering simulations
!
!> @brief
!!   Compute cloud columns using selected cloud overlap scheme for 
!!   visible/IR scattering simulations
!!
!! @details
!!   For visible/IR cloud simulations the vertical profile is divided into a
!!   number of cloud columns, each with a unique distribution of
!!   cloud among the layers. Each column has an associated weight. The
!!   radiative transfer equation is solved for each column and the TOA
!!   radiances are linearly combined using their respective column weights
!!   to obtain the total cloudy radiance.
!!
!!   Within RTTOV cloud column zero is always the clear-sky column. The
!!   overlap calculation and total number of columns varies from profile to
!!   profile.
!!
!!   Currently two cloud overlap options are available, selected via the
!!   opts\%rt_ir\%cloud_overlap option:
!!
!!   1 => maximum/random overlap, as described in:
!!   Matricardi, M., 2005 The inclusion of aerosols and clouds in RTIASI, the
!!   ECMWF fast radiative transfer model for the Infrared Atmospheric Sounding
!!   Interferometer. ECMWF Technical Memorandum 474.
!!
!!   This can result in a large number of cloud columns which increases the
!!   memory requirement and computational burden of the simulation. The
!!   cldcol_threshold parameter in the options structure can be used to
!!   ignore cloud columns with weights below the specified threshold (they are
!!   excluded and the clear column weight is increased in compensation). See
!!   the user guide for more details.
!!
!!   2 => simple cloud overlap
!!   This generates just two columns: one clear and one cloudy. The single cloud
!!   fraction is calculated as the maximum cloud fraction in the input profile 
!!   above the boundary layer. This is much more efficient in CPU and memory, 
!!   but is only appropriate for higher-peaking channels and is only recommended
!!   for advanced users who understand what they are doing.
!!
!!   The pressure level which defines the top of boundary layer can be specified
!!   in the opts\%rt_ir\%cc_low_cloud_top option (default: 750 hPa).
!!
!!   After computing the information related to the cloud overlap scheme
!!   this routine also scales the cloud concentrations by the relevant
!!   cloud fraction(s) if opts\%rt_ir\%grid_box_avg_cloud is true meaning the
!!   input cloud concentrations are grid box layer averages.
!!
!! @param[in]     opts_rt_ir          visible/IR-specific options structure
!! @param[in]     profiles            input atmospheric profiles and surface variables
!! @param[in,out] profiles_int        input atmospheric profiles converted to RTTOV internal units
!! @param[in,out] ircld               cloud column data
!! @param[out]    ncolumns            largest number of cloud columns across all profiles
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
SUBROUTINE rttov_cloud_overlap( &
              opts_rt_ir,   &
              profiles,     &
              profiles_int, &
              ircld,        &
              ncolumns)

  USE rttov_types, ONLY : rttov_profile, rttov_ircld, rttov_opts_rt_ir
  USE parkind1, ONLY : jpim
!INTF_OFF
  USE parkind1, ONLY : jprb
  USE rttov_const, ONLY : realtol, cloud_overlap_simple, cfrac_min
  USE yomhook, ONLY : LHOOK, DR_HOOK
!INTF_ON
  IMPLICIT NONE
  TYPE(rttov_opts_rt_ir),  INTENT(IN)    :: opts_rt_ir
  TYPE(rttov_profile),     INTENT(IN)    :: profiles(:)
  TYPE(rttov_profile),     INTENT(INOUT) :: profiles_int(SIZE(profiles))
  TYPE(rttov_ircld),       INTENT(INOUT) :: ircld
  INTEGER(KIND=jpim),      INTENT(OUT)   :: ncolumns
!INTF_END

  INTEGER(KIND=jpim) :: i, j, icol, ijcol, ilay
  REAL   (KIND=jprb) :: delta_cfrac
  INTEGER(KIND=jpim) :: ibdy_layer(1), imax_cfrac(1)
  REAL   (KIND=jprb) :: cfrac_max, cloud_tot(profiles(1)%nlayers)
  INTEGER(KIND=jpim) :: nprofiles
  REAL   (KIND=jprb) :: ntot(profiles(1)%nlevels)
  REAL   (KIND=jprb) :: ZHOOK_HANDLE
!-----End of header-------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_CLOUD_OVERLAP', 0_jpim, ZHOOK_HANDLE)
  nprofiles = SIZE(profiles)

  IF (opts_rt_ir%cloud_overlap == cloud_overlap_simple) THEN

    ! A simpler, faster, but quite approximate "Cmax" single-column approach. 
    ! Main benefit is that it is much more memory efficient. Intended mainly
    ! for mid- and upper-tropospheric channels.
    DO j = 1, nprofiles

      ircld%icldarr(:, :, j) = 0_jpim

      ! Find maximum cloud fraction and cloud amount above the boundary layer
      ! NB pressure is on levels, cloud in layers, so p(1:nlevels-1) is the top of the layer

      ! NAG compiler with optimisation gives wrong answers with this...
!       ibdy_layer = MINLOC(ABS(profiles(j)%p(1:profiles(j)%nlevels-1) - opts_rt_ir%cc_low_cloud_top))
      ! ...so instead use this
      ibdy_layer = MINLOC(ABS(profiles(j)%p(:) - opts_rt_ir%cc_low_cloud_top))
      ibdy_layer(1) = MIN(ibdy_layer(1), profiles(j)%nlayers)

      imax_cfrac = MAXLOC(profiles(j)%cfrac(1:ibdy_layer(1)))
      cfrac_max  = profiles(j)%cfrac(imax_cfrac(1))
      cloud_tot(:) = SUM(profiles_int(j)%cloud(:,:), DIM=1)

      ! Ignoring trivial amounts of cloud and precip
      IF (ANY(cloud_tot(1:ibdy_layer(1)) > 1E-6_jprb) .AND. cfrac_max > 1E-3_jprb) THEN

        ! Cloudy column required
        ircld%ncolumn(j) = 1
        ircld%xcol(1,j)  = 0._jprb
        ircld%xcol(2,j)  = cfrac_max
        ircld%xcolclr(j) = 1._jprb - cfrac_max
        WHERE (cloud_tot > 1E-6_jprb) ircld%icldarr(1,:,j) = 1

        IF (opts_rt_ir%grid_box_avg_cloud .AND. .NOT. opts_rt_ir%user_cld_opt_param) THEN
          profiles_int(j)%cloud(:,:) = profiles_int(j)%cloud(:,:) / cfrac_max
        ENDIF

      ELSE

        ! Pure clear-sky
        ircld%ncolumn(j)  = 0
        ircld%xcol(1:2,j) = 0._jprb
        ircld%xcolclr(j)  = 1._jprb

      ENDIF

    ENDDO
    ncolumns = MAXVAL(ircld%ncolumn(1:nprofiles))

  ELSE

    ncolumns  = 0_jpim
    delta_cfrac = 10._jprb * EPSILON(1._jprb)
    DO j = 1, nprofiles
      ircld%iflag(:, j)      = 1_jpim
      ircld%icldarr(:, :, j) = 0_jpim
      ircld%cldcfr(:, j)     = 0._jprb
      ircld%xcol(:, j)       = 0._jprb
      ircld%xcolmin(:, j)    = 0._jprb
      ircld%xcolmax(:, j)    = 0._jprb
      ircld%flag(:, j)       = .FALSE.
      ntot = 0._jprb
      DO ilay = 1, profiles(1)%nlayers
        IF (profiles(j)%cfrac(ilay) > 0._jprb) THEN
          ircld%cldcfr(ilay, j) = profiles(j)%cfrac(ilay)
        ENDIF
      ENDDO

      ! Check for overcast layers and identical cfrac values on consecutive layers
      ! These need adjusting to ensure TL/AD/K models are correct
      DO ilay = 2, profiles(1)%nlayers
        IF (ircld%cldcfr(ilay, j) > 0.) THEN

          ! Check for overcast layers
          IF (ircld%cldcfr(ilay, j) >= 1._jprb) THEN
            ircld%cldcfr(ilay, j) = 1._jprb - ilay*delta_cfrac
          ENDIF

          ! Check for identical adjacent cfrac (note that this won't always work if cldcol_threshold is +ve)
          IF (ircld%cldcfr(ilay, j) == ircld%cldcfr(ilay-1, j)) THEN
            ircld%cldcfr(ilay, j) = ircld%cldcfr(ilay, j) - SIGN(delta_cfrac, ircld%cldcfr(ilay, j)-0.5_jprb)
          ENDIF

        ENDIF
      ENDDO

  !---------Compute the cumulative cloud coverage using the maximum-random---------
  !         overlap assumption
      ntot(1)             = 1._jprb - ircld%cldcfr(1, j)
      ircld%ntotref(1, j) = ntot(1)
      DO ilay = 2, profiles(1)%nlayers
        ircld%maxcov(ilay, j) = (1._jprb - MAX(ircld%cldcfr(ilay - 1, j), ircld%cldcfr(ilay, j)))
        IF (ntot(ilay - 1) == 0._jprb .AND. ABS(ircld%cldcfr(ilay - 1, j) - 1._jprb) < realtol) THEN
          ntot(ilay) = ircld%maxcov(ilay, j)
        ELSEIF (ABS(ircld%cldcfr(ilay - 1, j) - 1._jprb) < realtol) THEN
          ntot(ilay) = ntot(ilay - 1)
        ELSE
          ntot(ilay) = ntot(ilay - 1) * ircld%maxcov(ilay, j) / (1._jprb - ircld%cldcfr(ilay - 1, j))
        ENDIF
        ircld%ntotref(ilay, j) = ntot(ilay)
      ENDDO
      DO ilay = 1, profiles(1)%nlayers
        IF (ircld%cldcfr(ilay, j) > 0._jprb) THEN
          ntot(ilay) = 1._jprb - ntot(ilay)
        ELSE
          ntot(ilay) = 0._jprb
        ENDIF
      ENDDO
  !---------Determine the limits of each column----------------------------------
      ircld%icount(j) = 1_jpim
      DO ilay = 1, profiles(1)%nlayers
        ircld%xcolmax(ilay, j)    = ntot(ilay)
        ircld%xcolmin(ilay, j)    = ntot(ilay) - ircld%cldcfr(ilay, j)
        ircld%xcolminref(ilay, j) = ircld%xcolmin(ilay, j)
        IF (ircld%xcolmin(ilay, j) < 0._jprb) ircld%xcolmin(ilay, j) = 0._jprb
        IF (ircld%xcolmax(ilay, j) > 0._jprb) THEN
          ircld%xcol(ircld%icount(j), j) = ircld%xcolmin(ilay, j)
          ircld%icount(j)                = ircld%icount(j) + 1
          ircld%xcol(ircld%icount(j), j) = ircld%xcolmax(ilay, j)
          ircld%icount(j)                = ircld%icount(j) + 1
        ENDIF
      ENDDO
      ircld%xcolref(:, j) = ircld%xcol(:, j)
  !---------Re-arrange the limits of each column in ascending order---------------
      loop1 : DO icol = 2, ircld%icount(j) - 1
        ircld%a(icol, j) = ircld%xcol(icol, j)
        DO i = icol - 1, 1, -1
          IF (ircld%xcol(i, j) <= ircld%a(icol, j)) THEN
            ircld%xcol(i + 1, j) = ircld%a(icol, j)
            ircld%iflag(icol, j) = i
            ircld%flag(icol, j)  = .TRUE.
            CYCLE loop1
          ELSE
            ircld%xcol(i + 1, j) = ircld%xcol(i, j)
          ENDIF
        ENDDO
        ircld%xcol(1, j) = ircld%a(icol, j)
      ENDDO loop1
      ircld%icount1(j)    = ircld%icount(j) - 1
      ircld%iloop(j)      = 0_jpim
      ircld%iloopin(:, j) = 0_jpim
      outer : DO
        ircld%iloop(j)                       = ircld%iloop(j) + 1
        ircld%icount1ref(ircld%iloop(j), j)  = ircld%icount1(j)
        ircld%xcolref2(ircld%iloop(j), :, j) = ircld%xcol(:, j)
        inner : DO icol = 1, ircld%icount1(j) - 1
          ircld%iloopin(ircld%iloop(j), j) = ircld%iloopin(ircld%iloop(j), j) + 1
          IF (ircld%xcol(icol, j) == ircld%xcol(icol + 1, j)) THEN
            IF (ircld%xcol(icol, j) < 1._jprb) THEN
              ircld%icount1(j) = ircld%icount1(j) - 1
              DO ijcol = icol, ircld%icount1(j)
                ircld%xcol(ijcol, j) = ircld%xcol(ijcol + 1, j)
              ENDDO
              CYCLE outer
            ELSE
              EXIT outer
            ENDIF
          ENDIF
        ENDDO inner
        EXIT outer
      ENDDO outer
      ircld%ncolumn(j) = ircld%icount1(j) - 1
  !---------Compute the weight of the clear column------------------------------
      DO icol = 1, ircld%ncolumn(j)
        DO ilay = 1, profiles(1)%nlayers
          IF (ircld%xcolmin(ilay, j) <= ircld%xcol(icol, j) .AND. ircld%xcolmax(ilay, j) >= ircld%xcol(icol + 1, j)) THEN
            ircld%icldarr(icol, ilay, j) = 1_jpim
          ENDIF
        ENDDO
      ENDDO
      IF (ircld%ncolumn(j) == -1_jpim) THEN
        ircld%ncolumn(j) = 0_jpim
      ENDIF
      ircld%xcolclr(j)    = 1._jprb - (ircld%xcol(ircld%ncolumn(j) + 1, j) - ircld%xcol(1, j))
      ircld%ncolumnref(j) = ircld%ncolumn(j)
  !---------Consider only the columns whose weight is greater than cldcol_threshold-------
      IF (ircld%ncolumn(j) /= 0_jpim) THEN
        ircld%icouncol(j) = 0_jpim
        DO icol = 1, ircld%ncolumn(j)
          IF ((ircld%xcol(icol + 1, j) - ircld%xcol(icol, j)) >= opts_rt_ir%cldcol_threshold) THEN
            ircld%icouncol(j)                    = ircld%icouncol(j) + 1
            ircld%indexcol(ircld%icouncol(j), j) = icol
          ENDIF
        ENDDO
        IF (ircld%icouncol(j) /= 0_jpim) THEN
          DO icol = 1, ircld%icouncol(j)
            DO ilay = 1, profiles(1)%nlayers
              ircld%icldarr(icol, ilay, j) = ircld%icldarr(ircld%indexcol(icol, j), ilay, j)
            ENDDO
          ENDDO
          DO icol = 1, ircld%icouncol(j) + 1
            IF (icol == 1_jpim) THEN
              ircld%xcol(icol, j) = ircld%xcol(ircld%indexcol(icol, j), j)
            ELSE
              ircld%xcol(icol, j) = ircld%xcol(icol - 1, j) +      &
                  (ircld%xcol(ircld%indexcol(icol - 1, j) + 1, j) - ircld%xcol(ircld%indexcol(icol - 1, j), j))
            ENDIF
          ENDDO
          ircld%xcolclr(j) = 1._jprb - (ircld%xcol(ircld%icouncol(j) + 1, j) - ircld%xcol(1, j))
          ircld%ncolumn(j) = ircld%icouncol(j)
        ELSE
          ircld%xcolclr(j) = 1._jprb
          ircld%ncolumn(j) = 0_jpim
        ENDIF
      ELSE
        ircld%ncolumn(j) = 0_jpim
      ENDIF
      ncolumns = MAX(ncolumns, ircld%ncolumn(j))

      IF (opts_rt_ir%grid_box_avg_cloud .AND. .NOT. opts_rt_ir%user_cld_opt_param) THEN
        DO ilay = 1, profiles(1)%nlayers
          IF (profiles(j)%cfrac(ilay) > cfrac_min) THEN
            profiles_int(j)%cloud(:,ilay) = profiles_int(j)%cloud(:,ilay) / profiles(j)%cfrac(ilay)
          ELSE
            profiles_int(j)%cloud(:,ilay) = 0._jprb
          ENDIF
        ENDDO
      ENDIF

    ENDDO ! profiles

  ENDIF

  IF (LHOOK) CALL DR_HOOK('RTTOV_CLOUD_OVERLAP', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_cloud_overlap
