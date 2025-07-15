! Description:
!> @file
!!   AD of cloud overlap schemes.
!
!> @brief
!!   AD of cloud overlap schemes.
!!
!!
!! @param[in]     opts_rt_ir          visible/IR-specific options structure
!! @param[in]     profiles            input atmospheric profiles and surface variables
!! @param[in,out] profiles_ad         atmospheric profile increments
!! @param[in]     profiles_int        input atmospheric profiles converted to RTTOV internal units
!! @param[in,out] profiles_int_ad     atmospheric profile increments in RTTOV internal units
!! @param[in,out] ircld               cloud column data
!! @param[in,out] ircld_ad            cloud column increments
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
SUBROUTINE rttov_cloud_overlap_ad( &
              opts_rt_ir,      &
              profiles,        &
              profiles_ad,     &
              profiles_int,    &
              profiles_int_ad, &
              ircld,           &
              ircld_ad)

  USE rttov_types, ONLY : rttov_profile, rttov_ircld, rttov_opts_rt_ir
!INTF_OFF
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jpim, jprb
  USE rttov_const, ONLY : realtol, cloud_overlap_simple, cfrac_min
!INTF_ON
  IMPLICIT NONE

  TYPE(rttov_opts_rt_ir),  INTENT(IN)    :: opts_rt_ir
  TYPE(rttov_profile),     INTENT(IN)    :: profiles(:)
  TYPE(rttov_profile),     INTENT(INOUT) :: profiles_ad(SIZE(profiles))
  TYPE(rttov_profile),     INTENT(IN)    :: profiles_int(SIZE(profiles))
  TYPE(rttov_profile),     INTENT(INOUT) :: profiles_int_ad(SIZE(profiles))
  TYPE(rttov_ircld),       INTENT(INOUT) :: ircld
  TYPE(rttov_ircld),       INTENT(INOUT) :: ircld_ad
!INTF_END

  INTEGER(KIND=jpim) :: i, j, icol, ijcol, ilay, ic
  INTEGER(KIND=jpim) :: nprofiles
  INTEGER(KIND=jpim) :: ibdy_layer(1), imax_cfrac(1)
  REAL   (KIND=jprb) :: cfrac_max, cloud_tot(profiles(1)%nlayers) 
  REAL   (KIND=jprb) :: cfrac_max_ad 
  REAL   (KIND=jprb) :: ntot(profiles(1)%nlevels)
  REAL   (KIND=jprb) :: ZHOOK_HANDLE
!-----End of header-------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_CLOUD_OVERLAP_AD', 0_jpim, ZHOOK_HANDLE)
  nprofiles = SIZE(profiles)

  IF (opts_rt_ir%cloud_overlap == cloud_overlap_simple) THEN

    ! A simpler, faster, but quite approximate "Cmax" single-column approach. 
    ! Main benefit is that it is much more memory efficient. Intended mainly
    ! for mid- and upper-tropospheric channels.
    DO j = 1, nprofiles

      ! Find maximum cloud fraction and cloud amount above the boundary layer
      ! NB pressure is on levels, cloud in layers, so p(1:nlevels-1) is the top of the layer
      ibdy_layer = MINLOC(ABS(profiles(j)%p(:) - opts_rt_ir%cc_low_cloud_top))
      ibdy_layer(1) = MIN(ibdy_layer(1), profiles(j)%nlayers)
      imax_cfrac = MAXLOC(profiles(j)%cfrac(1:ibdy_layer(1)))
      cfrac_max  = profiles(j)%cfrac(imax_cfrac(1))
      cloud_tot(:) = SUM(profiles_int(j)%cloud(:,:), DIM=1)

      cfrac_max_ad = 0._jprb

      ! Ignoring trivial amounts of cloud and precip
      IF (ANY(cloud_tot(1:ibdy_layer(1)) > 1E-6_jprb) .AND. cfrac_max > 1E-3_jprb) THEN

        IF (opts_rt_ir%grid_box_avg_cloud .AND. .NOT. opts_rt_ir%user_cld_opt_param) THEN
          ! The value of profiles_int(iprof)%cloud required is that before the cfrac_max
          ! scaling hence it is multiplied by cfrac_max
          cfrac_max_ad = cfrac_max_ad - &
              SUM(profiles_int_ad(j)%cloud(:,:) * profiles_int(j)%cloud(:,:)) / cfrac_max

          ! This must be done after the above
          profiles_int_ad(j)%cloud(:,:) = profiles_int_ad(j)%cloud(:,:) / cfrac_max
        ENDIF

        ! Cloudy column required
        cfrac_max_ad = cfrac_max_ad + ircld_ad%xcol(2,j)
        cfrac_max_ad = cfrac_max_ad - ircld_ad%xcolclr(j)

      ENDIF

      profiles_ad(j)%cfrac(imax_cfrac(1)) = profiles_ad(j)%cfrac(imax_cfrac(1)) + cfrac_max_ad

    ENDDO

  ELSE

    DO j = 1, nprofiles

      IF (opts_rt_ir%grid_box_avg_cloud .AND. .NOT. opts_rt_ir%user_cld_opt_param) THEN
        DO ilay = 1, profiles(1)%nlayers
          IF (profiles(j)%cfrac(ilay) > cfrac_min) THEN
            ! The value of profiles_int(iprof)%cloud required is that before the cfrac scaling hence it is multiplied by cfrac
            profiles_ad(j)%cfrac(ilay) = profiles_ad(j)%cfrac(ilay) - &
                SUM(profiles_int_ad(j)%cloud(:,ilay) * profiles_int(j)%cloud(:,ilay)) / profiles(j)%cfrac(ilay)

            ! This must be done after the above
            profiles_int_ad(j)%cloud(:,ilay) = profiles_int_ad(j)%cloud(:,ilay) / profiles(j)%cfrac(ilay)
          ELSE
            profiles_int_ad(j)%cloud(:,ilay) = 0._jprb
          ENDIF
        ENDDO
      ENDIF

      loop2 : DO icol = 2, ircld%icount(j) - 1
        DO i = icol - 1, 1, -1
          ircld%xcolref1(icol, i, j) = ircld%xcolref(i, j)
          IF (ircld%xcolref(i, j) <= ircld%a(icol, j)) THEN
            ircld%xcolref(i + 1, j) = ircld%a(icol, j)
            CYCLE loop2
          ELSE
            ircld%xcolref(i + 1, j) = ircld%xcolref(i, j)
          ENDIF
        ENDDO
        i = 0
        ircld%xcolref(i + 1, j) = ircld%a(icol, j)
      ENDDO loop2
      ircld_ad%a(:, j)       = 0._jprb
      ircld_ad%maxcov(:, j)  = 0._jprb
      ircld_ad%cldcfr(:, j)  = 0._jprb
      ircld_ad%xcolmin(:, j) = 0._jprb
      ircld_ad%xcolmax(:, j) = 0._jprb
      ntot(:)                = 0._jprb
  !---------Consider only the columns whose weight is greater than cldcol_threshold-------
      IF (ircld%ncolumnref(j) /= 0_jpim) THEN
        IF (ircld%icouncol(j) /= 0_jpim) THEN
          ircld_ad%xcol(ircld%icouncol(j) + 1, j) = ircld_ad%xcol(ircld%icouncol(j) + 1, j) - ircld_ad%xcolclr(j)
          ircld_ad%xcol(1, j)                     = ircld_ad%xcol(1, j) + ircld_ad%xcolclr(j)
          ircld_ad%xcolclr(j)                     = 0._jprb
          DO icol = ircld%icouncol(j) + 1, 1, -1
            IF (icol == 1) THEN
              IF (ircld%indexcol(icol, j) /= icol) THEN
                ircld_ad%xcol(ircld%indexcol(icol, j), j) =      &
                    ircld_ad%xcol(ircld%indexcol(icol, j), j) + ircld_ad%xcol(icol, j)
                ircld_ad%xcol(icol, j)                    = 0._jprb
              ENDIF
            ELSE
              ircld_ad%xcol(icol - 1, j) = ircld_ad%xcol(icol - 1, j) + ircld_ad%xcol(icol, j)
              IF ((ircld%indexcol(icol - 1, j) + 1) /= icol) THEN
                ircld_ad%xcol(ircld%indexcol(icol - 1, j) + 1, j) =      &
                    ircld_ad%xcol(ircld%indexcol(icol - 1, j) + 1, j) + ircld_ad%xcol(icol, j)
              ELSE
                ircld_ad%xcol(ircld%indexcol(icol - 1, j) + 1, j) = ircld_ad%xcol(icol, j)
              ENDIF
              IF (ircld%indexcol(icol - 1, j) /= icol) THEN
                ircld_ad%xcol(ircld%indexcol(icol - 1, j), j) =      &
                    ircld_ad%xcol(ircld%indexcol(icol - 1, j), j) - ircld_ad%xcol(icol, j)
              ELSE
                ircld_ad%xcol(ircld%indexcol(icol - 1, j), j) =  - ircld_ad%xcol(icol, j)
              ENDIF
              IF (((ircld%indexcol(icol - 1, j) + 1) /= icol) .AND. (ircld%indexcol(icol - 1, j) /= icol)) THEN
                ircld_ad%xcol(icol, j) = 0._jprb
              ENDIF
            ENDIF
          ENDDO
        ELSE
          ircld_ad%xcolclr(j) = 0._jprb
        ENDIF
      ENDIF
  !---------Compute the weight of the clear column--------------------------------
      ircld_ad%xcol(ircld%ncolumnref(j) + 1, j) = ircld_ad%xcol(ircld%ncolumnref(j) + 1, j) - ircld_ad%xcolclr(j)
      ircld_ad%xcol(1, j)                       = ircld_ad%xcol(1, j) + ircld_ad%xcolclr(j)
  !---------Re-arrange the limits of each column in ascending order---------------
      outer : DO i = ircld%iloop(j), 1, -1
        inner : DO icol = ircld%iloopin(i, j), 1, -1
          IF (ircld%xcolref2(i, icol, j) == ircld%xcolref2(i, icol + 1, j)) THEN
            IF (ircld%xcolref2(i, icol, j) < 1._jprb) THEN
              ircld%icount1ref(i, j) = ircld%icount1ref(i, j) - 1
              DO ijcol = ircld%icount1ref(i, j), icol, -1
                ircld_ad%xcol(ijcol + 1, j) = ircld_ad%xcol(ijcol + 1, j) + ircld_ad%xcol(ijcol, j)
                ircld_ad%xcol(ijcol, j)     = 0._jprb
              ENDDO
            ENDIF
          ENDIF
        ENDDO inner
      ENDDO outer
      loop1 : DO icol = ircld%icount(j) - 1, 2, -1
        IF (.NOT. ircld%flag(icol, j)) THEN
          ircld_ad%a(icol, j) = ircld_ad%a(icol, j) + ircld_ad%xcol(1, j)
          ircld_ad%xcol(1, j) = 0._jprb                                      !
        ENDIF
        DO i = ircld%iflag(icol, j), icol - 1
          IF (ircld%xcolref1(icol, i, j) <= ircld%a(icol, j)) THEN
            ircld_ad%a(icol, j)     = ircld_ad%a(icol, j) + ircld_ad%xcol(i + 1, j)
            ircld_ad%xcol(i + 1, j) = 0._jprb
          ELSE
            ircld_ad%xcol(i, j)     = ircld_ad%xcol(i, j) + ircld_ad%xcol(i + 1, j)
            ircld_ad%xcol(i + 1, j) = 0._jprb
          ENDIF
        ENDDO
        ircld_ad%xcol(icol, j) = ircld_ad%xcol(icol, j) + ircld_ad%a(icol, j)
        ircld_ad%a(icol, j)    = 0._jprb
      ENDDO loop1
  !---------Determine the limits of each column----------------------------------
      ic = ircld%icount(j)
      DO ilay = profiles(1)%nlayers, 1, -1
        IF (ircld%xcolmax(ilay, j) > 0._jprb) THEN
          ic = ic - 1
          ircld_ad%xcolmax(ilay, j) = ircld_ad%xcolmax(ilay, j) + ircld_ad%xcol(ic, j)
          ic = ic - 1
          ircld_ad%xcolmin(ilay, j) = ircld_ad%xcolmin(ilay, j) + ircld_ad%xcol(ic, j)
        ENDIF
        IF (ircld%xcolminref(ilay, j) < 0._jprb) THEN
          ircld_ad%xcolmin(ilay, j) = 0._jprb
        ENDIF
        ntot(ilay)               = ntot(ilay) + ircld_ad%xcolmin(ilay, j)
        ircld_ad%cldcfr(ilay, j) = ircld_ad%cldcfr(ilay, j) - ircld_ad%xcolmin(ilay, j)
        ntot(ilay)               = ntot(ilay) + ircld_ad%xcolmax(ilay, j)
      ENDDO
  !---------Compute the cumulative cloud coverage using the maximum-random---------
  !         overlap assumption
      DO ilay = profiles(1)%nlayers, 1, -1
        IF (ircld%cldcfr(ilay, j) > 0._jprb) THEN
          ntot(ilay) =  - ntot(ilay)
        ELSE
          ntot(ilay) = 0._jprb
        ENDIF
      ENDDO
      DO ilay = profiles(1)%nlayers, 2, -1
        IF (ircld%ntotref(ilay - 1, j) == 0._jprb .AND. ABS(ircld%cldcfr(ilay - 1, j) - 1._jprb) < realtol) THEN
          ircld_ad%maxcov(ilay, j) = ircld_ad%maxcov(ilay, j) + ntot(ilay)
        ELSEIF (ABS(ircld%cldcfr(ilay - 1, j) - 1._jprb) < realtol) THEN
          ntot(ilay - 1)               =      &
              ntot(ilay - 1) + ntot(ilay)
        ELSE
          ntot(ilay - 1)               =      &
              ntot(ilay - 1) + ntot(ilay) * ircld%maxcov(ilay, j) / (1._jprb - ircld%cldcfr(ilay - 1, j))
          ircld_ad%cldcfr(ilay - 1, j) = ircld_ad%cldcfr(ilay - 1, j) +      &
              ntot(ilay) * ircld%ntotref(ilay - 1, j) * ircld%maxcov(ilay, j) / (1._jprb - ircld%cldcfr(ilay - 1, j)) ** 2
          ircld_ad%maxcov(ilay, j)     =      &
              ircld_ad%maxcov(ilay, j) + ntot(ilay) * ircld%ntotref(ilay - 1, j) / (1._jprb - ircld%cldcfr(ilay - 1, j))
        ENDIF
        IF ((ircld%cldcfr(ilay - 1, j) > ircld%cldcfr(ilay, j))) THEN
          ircld_ad%cldcfr(ilay - 1, j) = ircld_ad%cldcfr(ilay - 1, j) - ircld_ad%maxcov(ilay, j)
          ircld_ad%maxcov(ilay, j)     = 0._jprb
        ELSE
          ircld_ad%cldcfr(ilay, j) = ircld_ad%cldcfr(ilay, j) - ircld_ad%maxcov(ilay, j)
          ircld_ad%maxcov(ilay, j) = 0._jprb
        ENDIF
      ENDDO
      ircld_ad%cldcfr(1, j) = ircld_ad%cldcfr(1, j) - ntot(1)

      DO ilay = profiles(1)%nlayers, 1, -1
        IF (profiles(j)%cfrac(ilay) > 0._jprb) THEN
          profiles_ad(j)%cfrac(ilay) = profiles_ad(j)%cfrac(ilay) + ircld_ad%cldcfr(ilay, j)
        ENDIF
      ENDDO
      ircld_ad%cldcfr(:, j)  = 0._jprb
      ircld_ad%xcolmin(:, j) = 0._jprb
      ircld_ad%xcolmax(:, j) = 0._jprb
      ircld_ad%xcol(:, j) = 0._jprb
      ntot = 0._jprb
    ENDDO

  ENDIF

  IF (LHOOK) CALL DR_HOOK('RTTOV_CLOUD_OVERLAP_AD', 1_jpim, ZHOOK_HANDLE)
!      ircld_ad(:)%xcolclr=0._jprb
END SUBROUTINE rttov_cloud_overlap_ad
