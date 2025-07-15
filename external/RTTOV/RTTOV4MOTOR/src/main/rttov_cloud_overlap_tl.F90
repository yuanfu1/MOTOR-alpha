! Description:
!> @file
!!   TL of cloud overlap schemes.
!
!> @brief
!!   TL of cloud overlap schemes.
!!
!!
!! @param[in]     opts_rt_ir          visible/IR-specific options structure
!! @param[in]     profiles            input atmospheric profiles and surface variables
!! @param[in]     profiles_tl         atmospheric profile perturbations
!! @param[in]     profiles_int        input atmospheric profiles converted to RTTOV internal units
!! @param[in,out] profiles_int_tl     atmospheric profile perturbations in RTTOV internal units
!! @param[in,out] ircld               cloud column data
!! @param[in,out] ircld_tl            cloud column perturbations
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
SUBROUTINE rttov_cloud_overlap_tl( &
              opts_rt_ir,      &
              profiles,        &
              profiles_tl,     &
              profiles_int,    &
              profiles_int_tl, &
              ircld,           &
              ircld_tl)

  USE rttov_types, ONLY : rttov_profile, rttov_ircld, rttov_opts_rt_ir
!INTF_OFF
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jpim, jprb
  USE rttov_const, ONLY : realtol, cloud_overlap_simple, cfrac_min
!INTF_ON
  IMPLICIT NONE

  TYPE(rttov_opts_rt_ir),  INTENT(IN)    :: opts_rt_ir
  TYPE(rttov_profile),     INTENT(IN)    :: profiles(:)
  TYPE(rttov_profile),     INTENT(IN)    :: profiles_tl(SIZE(profiles))
  TYPE(rttov_profile),     INTENT(IN)    :: profiles_int(SIZE(profiles))
  TYPE(rttov_profile),     INTENT(INOUT) :: profiles_int_tl(SIZE(profiles))
  TYPE(rttov_ircld),       INTENT(INOUT) :: ircld
  TYPE(rttov_ircld),       INTENT(INOUT) :: ircld_tl
!INTF_END

  INTEGER(KIND=jpim) :: i, j, icol, ijcol, ilay
  INTEGER(KIND=jpim) :: nprofiles
  INTEGER(KIND=jpim) :: ibdy_layer(1), imax_cfrac(1)
  REAL   (KIND=jprb) :: cfrac_max, cloud_tot(profiles(1)%nlayers)
  REAL   (KIND=jprb) :: cfrac_max_tl 
  REAL   (KIND=jprb) :: ntot(profiles(1)%nlevels)
  REAL   (KIND=jprb) :: ZHOOK_HANDLE
!-----End of header-------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_CLOUD_OVERLAP_TL', 0_jpim, ZHOOK_HANDLE)
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
      cfrac_max_tl = profiles_tl(j)%cfrac(imax_cfrac(1))
      cloud_tot(:) = SUM(profiles_int(j)%cloud(:,:), DIM=1)

      ! Ignoring trivial amounts of cloud and precip
      IF (ANY(cloud_tot(1:ibdy_layer(1)) > 1E-6_jprb) .AND. cfrac_max > 1E-3_jprb) THEN

        ! Cloudy column required
        ircld_tl%xcol(1,j)  = 0._jprb
        ircld_tl%xcol(2,j)  = cfrac_max_tl
        ircld_tl%xcolclr(j) = -1._jprb*cfrac_max_tl

        IF (opts_rt_ir%grid_box_avg_cloud .AND. .NOT. opts_rt_ir%user_cld_opt_param) THEN
          ! The value of profiles_int(iprof)%cloud required is that before the cfrac_max
          ! scaling hence it is multiplied by cfrac_max
          profiles_int_tl(j)%cloud(:,:) = profiles_int_tl(j)%cloud(:,:) / cfrac_max - &
                                          cfrac_max_tl * profiles_int(j)%cloud(:,:) / cfrac_max
        ENDIF

      ELSE

        ! Pure clear-sky
        ircld_tl%xcol(1:2,j) = 0._jprb
        ircld_tl%xcolclr(j)  = 0._jprb

      ENDIF

    ENDDO

  ELSE

    DO j = 1, nprofiles
      ircld_tl%cldcfr(:, j)  = 0._jprb
      ircld_tl%xcolmin(:, j) = 0._jprb
      ircld_tl%xcolmax(:, j) = 0._jprb
      ircld_tl%xcol(:, j) = 0._jprb ! important when there is no cloud, as it would be undefined otherwise
      ntot = 0._jprb
      DO ilay = 1, profiles(1)%nlayers
        IF (profiles(j)%cfrac(ilay) > 0._jprb) THEN
          ircld_tl%cldcfr(ilay, j) = profiles_tl(j)%cfrac(ilay)
        ENDIF
      ENDDO
      ircld%icount(j) = 0_jpim
  !---------Compute the cumulative cloud coverage using the maximum-random---------
  !         overlap assumption
      ntot(1)         =  - ircld_tl%cldcfr(1, j)
      DO ilay = 2, profiles(1)%nlayers
        IF ((ircld%cldcfr(ilay - 1, j) > ircld%cldcfr(ilay, j))) THEN
          ircld_tl%maxcov(ilay, j) =  - ircld_tl%cldcfr(ilay - 1, j)
        ELSE
          ircld_tl%maxcov(ilay, j) =  - ircld_tl%cldcfr(ilay, j)
        ENDIF
        IF (ircld%ntotref(ilay - 1, j) == 0._jprb .AND. ABS(ircld%cldcfr(ilay - 1, j) - 1._jprb) < realtol) THEN
          ntot(ilay) = ircld_tl%maxcov(ilay, j)
        ELSEIF (ABS(ircld%cldcfr(ilay - 1, j) - 1._jprb) < realtol) THEN
          ntot(ilay) = ntot(ilay - 1)
        ELSE
          ntot(ilay) = ntot(ilay - 1) * ircld%maxcov(ilay, j) / (1._jprb - ircld%cldcfr(ilay - 1, j)) +      &
              ircld_tl%cldcfr(ilay - 1, j) * ircld%ntotref(ilay - 1, j) * ircld%maxcov(ilay, j) /            &
              (1._jprb - ircld%cldcfr(ilay - 1, j)) ** 2 +                                                   &
              ircld_tl%maxcov(ilay, j) * ircld%ntotref(ilay - 1, j) / (1._jprb - ircld%cldcfr(ilay - 1, j))
        ENDIF
      ENDDO
      DO ilay = 1, profiles(1)%nlayers
        IF (ircld%cldcfr(ilay, j) > 0._jprb) THEN
          ntot(ilay)   =  - ntot(ilay)
        ELSE
          ntot(ilay)   = 0._jprb
        ENDIF
      ENDDO
  !---------Determine the limits of each column----------------------------------
      ircld%icount(j) = 1_jpim
      DO ilay = 1, profiles(1)%nlayers
        ircld_tl%xcolmax(ilay, j) = ntot(ilay)
        ircld_tl%xcolmin(ilay, j) = ntot(ilay) - ircld_tl%cldcfr(ilay, j)
        IF (ircld%xcolminref(ilay, j) < 0._jprb) ircld_tl%xcolmin(ilay, j) = 0._jprb
        IF (ircld%xcolmax(ilay, j) > 0._jprb) THEN
          ircld_tl%xcol(ircld%icount(j), j) = ircld_tl%xcolmin(ilay, j)
          ircld%icount(j)                   = ircld%icount(j) + 1
          ircld_tl%xcol(ircld%icount(j), j) = ircld_tl%xcolmax(ilay, j)
          ircld%icount(j)                   = ircld%icount(j) + 1
        ENDIF
      ENDDO
  !---------Re-arrange the limits of each column in ascending order---------------
      loop1 : DO icol = 2, ircld%icount(j) - 1
        ircld%a(icol, j)    = ircld%xcolref(icol, j)
        ircld_tl%a(icol, j) = ircld_tl%xcol(icol, j)
        DO i = icol - 1, 1, -1
          IF (ircld%xcolref(i, j) <= ircld%a(icol, j)) THEN
            ircld%xcolref(i + 1, j) = ircld%a(icol, j)
            ircld_tl%xcol(i + 1, j) = ircld_tl%a(icol, j)
            CYCLE loop1
          ELSE
            ircld%xcolref(i + 1, j) = ircld%xcolref(i, j)
            ircld_tl%xcol(i + 1, j) = ircld_tl%xcol(i, j)
          ENDIF
        ENDDO
        ircld%xcolref(1, j) = ircld%a(icol, j)
        ircld_tl%xcol(1, j) = ircld_tl%a(icol, j)
      ENDDO loop1
      outer : DO i = 1, ircld%iloop(j)
        inner : DO icol = 1, ircld%iloopin(i, j)
          IF (ircld%xcolref2(i, icol, j) == ircld%xcolref2(i, icol + 1, j)) THEN
            IF (ircld%xcolref2(i, icol, j) < 1._jprb) THEN
              ircld%icount1ref(i, j) = ircld%icount1ref(i, j) - 1
              DO ijcol = icol, ircld%icount1ref(i, j)
                ircld_tl%xcol(ijcol, j) = ircld_tl%xcol(ijcol + 1, j)
              ENDDO
            ENDIF
          ENDIF
        ENDDO inner
      ENDDO outer
  !---------Compute the weight of the clear column------------------------------
      ircld_tl%xcolclr(j) =  - ircld_tl%xcol(ircld%ncolumnref(j) + 1, j) + ircld_tl%xcol(1, j)
  !---------Consider only the columns whose weight is greater than cldcol_threshold-------
      IF (ircld%ncolumnref(j) /= 0_jpim) THEN
        IF (ircld%icouncol(j) /= 0_jpim) THEN
          DO icol = 1, ircld%icouncol(j) + 1
            IF (icol == 1) THEN
              ircld_tl%xcol(icol, j) = ircld_tl%xcol(ircld%indexcol(icol, j), j)
            ELSE
              ircld_tl%xcol(icol, j) = ircld_tl%xcol(icol - 1, j) +      &
                  (ircld_tl%xcol(ircld%indexcol(icol - 1, j) + 1, j) - ircld_tl%xcol(ircld%indexcol(icol - 1, j), j))
            ENDIF
          ENDDO
          ircld_tl%xcolclr(j) =  - ircld_tl%xcol(ircld%icouncol(j) + 1, j) + ircld_tl%xcol(1, j)
        ELSE
          ircld_tl%xcolclr(j) = 0._jprb
        ENDIF
      ENDIF

      IF (opts_rt_ir%grid_box_avg_cloud .AND. .NOT. opts_rt_ir%user_cld_opt_param) THEN
        DO ilay = 1, profiles(1)%nlayers
          IF (profiles(j)%cfrac(ilay) > cfrac_min) THEN
            ! The value of profiles_int(iprof)%cloud required is that before the cfrac scaling hence it is multiplied by cfrac
            profiles_int_tl(j)%cloud(:,ilay) = profiles_int_tl(j)%cloud(:,ilay) / profiles(j)%cfrac(ilay) - &
                profiles_tl(j)%cfrac(ilay) * profiles_int(j)%cloud(:,ilay) / profiles(j)%cfrac(ilay)
          ELSE
            profiles_int_tl(j)%cloud(:,ilay) = 0._jprb
          ENDIF
        ENDDO
      ENDIF

    ENDDO ! profiles

  ENDIF

  IF (LHOOK) CALL DR_HOOK('RTTOV_CLOUD_OVERLAP_TL', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_cloud_overlap_tl
