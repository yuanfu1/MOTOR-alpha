! Description:
!> @file
!!   K of cloud overlap schemes.
!
!> @brief
!!   K of cloud overlap schemes.
!!
!!
!! @param[in]     opts_rt_ir          visible/IR-specific options structure
!! @param[in]     chanprof            specifies channels and profiles to simulate
!! @param[in]     profiles            input atmospheric profiles and surface variables
!! @param[in,out] profiles_k          atmospheric profile increments
!! @param[in]     profiles_int        input atmospheric profiles converted to RTTOV internal units
!! @param[in,out] profiles_int_k      atmospheric profile increments in RTTOV internal units
!! @param[in,out] ircld               cloud column data
!! @param[in,out] ircld_k             cloud column increments
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
SUBROUTINE rttov_cloud_overlap_k( &
              opts_rt_ir,     &
              chanprof,       &
              profiles,       &
              profiles_k,     &
              profiles_int,   &
              profiles_int_k, &
              ircld,          &
              ircld_k)

  USE rttov_types, ONLY : rttov_profile, rttov_ircld, rttov_opts_rt_ir, &
                          rttov_chanprof
!INTF_OFF
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jpim, jprb
  USE rttov_const, ONLY : realtol, cloud_overlap_simple, cfrac_min
!INTF_ON
  IMPLICIT NONE
  TYPE(rttov_opts_rt_ir),  INTENT(IN)    :: opts_rt_ir
  TYPE(rttov_chanprof),    INTENT(IN)    :: chanprof(:)
  TYPE(rttov_profile),     INTENT(IN)    :: profiles(:)
  TYPE(rttov_profile),     INTENT(INOUT) :: profiles_k(SIZE(chanprof))
  TYPE(rttov_profile),     INTENT(IN)    :: profiles_int(SIZE(profiles))
  TYPE(rttov_profile),     INTENT(INOUT) :: profiles_int_k(SIZE(chanprof))
  TYPE(rttov_ircld),       INTENT(INOUT) :: ircld
  TYPE(rttov_ircld),       INTENT(INOUT) :: ircld_k
!INTF_END

  INTEGER(KIND=jpim) :: i, j, icol, ijcol, ilay, ic, jpk
  INTEGER(KIND=jpim) :: nchannels
  INTEGER(KIND=jpim) :: ibdy_layer(1), imax_cfrac(1)
  REAL   (KIND=jprb) :: cfrac_max, cloud_tot(profiles(1)%nlayers) 
  REAL   (KIND=jprb) :: cfrac_max_k 
  REAL   (KIND=jprb) :: ntot(profiles(1)%nlevels)
  INTEGER(KIND=jpim) :: icount1ref(2*profiles(1)%nlayers)
  REAL   (KIND=jprb) :: ZHOOK_HANDLE
!-----End of header-------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_CLOUD_OVERLAP_K', 0_jpim, ZHOOK_HANDLE)
  nchannels = SIZE(chanprof)

  IF (opts_rt_ir%cloud_overlap == cloud_overlap_simple) THEN

    ! A simpler, faster, but quite approximate "Cmax" single-column approach. 
    ! Main benefit is that it is much more memory efficient. Intended mainly
    ! for mid- and upper-tropospheric channels.
    DO j = 1, nchannels
      jpk = chanprof(j)%prof

      ! Find maximum cloud fraction and cloud amount above the boundary layer
      ! NB pressure is on levels, cloud in layers, so p(1:nlevels-1) is the top of the layer
      ibdy_layer = MINLOC(ABS(profiles(jpk)%p(:) - opts_rt_ir%cc_low_cloud_top))
      ibdy_layer(1) = MIN(ibdy_layer(1), profiles(jpk)%nlayers)
      imax_cfrac = MAXLOC(profiles(jpk)%cfrac(1:ibdy_layer(1)))
      cfrac_max  = profiles(jpk)%cfrac(imax_cfrac(1))
      cloud_tot(:) = SUM(profiles_int(jpk)%cloud(:,:), DIM=1)

      cfrac_max_k = 0._jprb

      ! Ignoring trivial amounts of cloud and precip
      IF (ANY(cloud_tot(1:ibdy_layer(1)) > 1E-6_jprb) .AND. cfrac_max > 1E-3_jprb) THEN

        IF (opts_rt_ir%grid_box_avg_cloud .AND. .NOT. opts_rt_ir%user_cld_opt_param) THEN
          ! The value of profiles_int(iprof)%cloud required is that before the cfrac_max
          ! scaling hence it is multiplied by cfrac_max
          cfrac_max_k = cfrac_max_k - &
              SUM(profiles_int_k(j)%cloud(:,:) * profiles_int(jpk)%cloud(:,:)) / cfrac_max

          ! This must be done after the above
          profiles_int_k(j)%cloud(:,:) = profiles_int_k(j)%cloud(:,:) / cfrac_max
        ENDIF

        ! Cloudy column required
        cfrac_max_k = cfrac_max_k + ircld_k%xcol(2,j)
        cfrac_max_k = cfrac_max_k - ircld_k%xcolclr(j)

      ENDIF

      profiles_k(j)%cfrac(imax_cfrac(1)) = profiles_k(j)%cfrac(imax_cfrac(1)) + cfrac_max_k

    ENDDO

  ELSE

    DO j = 1, nchannels
      jpk = chanprof(j)%prof

      IF (opts_rt_ir%grid_box_avg_cloud .AND. .NOT. opts_rt_ir%user_cld_opt_param) THEN
        DO ilay = 1, profiles(1)%nlayers
          IF (profiles(jpk)%cfrac(ilay) > cfrac_min) THEN
            ! The value of profiles_int(prof)%cloud required is that before the cfrac scaling hence it is multiplied by cfrac
            profiles_k(j)%cfrac(ilay) = profiles_k(j)%cfrac(ilay) - &
                SUM(profiles_int_k(j)%cloud(:,ilay) * profiles_int(jpk)%cloud(:,ilay)) / profiles(jpk)%cfrac(ilay)

            ! This must be done after the above
            profiles_int_k(j)%cloud(:,ilay) = profiles_int_k(j)%cloud(:,ilay) / profiles(jpk)%cfrac(ilay)
          ELSE
            profiles_int_k(j)%cloud(:,ilay) = 0._jprb
          ENDIF
        ENDDO
      ENDIF

      loop2 : DO icol = 2, ircld%icount(jpk) - 1
        DO i = icol - 1, 1, -1
          ircld%xcolref1(icol, i, jpk) = ircld%xcolref(i, jpk)
          IF (ircld%xcolref(i, jpk) <= ircld%a(icol, jpk)) THEN
            ircld%xcolref(i + 1, jpk) = ircld%a(icol, jpk)
            CYCLE loop2
          ELSE
            ircld%xcolref(i + 1, jpk) = ircld%xcolref(i, jpk)
          ENDIF
        ENDDO
        i = 0
        ircld%xcolref(i + 1, jpk) = ircld%a(icol, jpk)
      ENDDO loop2
      ircld_k%a(:, j)       = 0._jprb
      ircld_k%maxcov(:, j)  = 0._jprb
      ircld_k%cldcfr(:, j)  = 0._jprb
      ircld_k%xcolmin(:, j) = 0._jprb
      ircld_k%xcolmax(:, j) = 0._jprb
      ntot(:)               = 0._jprb
  !---------Consider only the columns whose weight is greater than cldcol_threshold-------
      IF (ircld%ncolumnref(jpk) /= 0_jpim) THEN
        IF (ircld%icouncol(jpk) /= 0_jpim) THEN
          ircld_k%xcol(ircld%icouncol(jpk) + 1, j) = ircld_k%xcol(ircld%icouncol(jpk) + 1, j) - ircld_k%xcolclr(j)
          ircld_k%xcol(1, j)                       = ircld_k%xcol(1, j) + ircld_k%xcolclr(j)
          ircld_k%xcolclr(j)                       = 0._jprb
          DO icol = ircld%icouncol(jpk) + 1, 1, -1
            IF (icol == 1) THEN
              IF (ircld%indexcol(icol, jpk) /= icol) THEN
                ircld_k%xcol(ircld%indexcol(icol, jpk), j) =      &
                    ircld_k%xcol(ircld%indexcol(icol, jpk), j) + ircld_k%xcol(icol, j)
                ircld_k%xcol(icol, j)                      = 0._jprb
              ENDIF
            ELSE
              ircld_k%xcol(icol - 1, j) = ircld_k%xcol(icol - 1, j) + ircld_k%xcol(icol, j)
              IF ((ircld%indexcol(icol - 1, jpk) + 1) /= icol) THEN
                ircld_k%xcol(ircld%indexcol(icol - 1, jpk) + 1, j) =      &
                    ircld_k%xcol(ircld%indexcol(icol - 1, jpk) + 1, j) + ircld_k%xcol(icol, j)
              ELSE
                ircld_k%xcol(ircld%indexcol(icol - 1, jpk) + 1, j) = ircld_k%xcol(icol, j)
              ENDIF
              IF (ircld%indexcol(icol - 1, jpk) /= icol) THEN
                ircld_k%xcol(ircld%indexcol(icol - 1, jpk), j) =      &
                    ircld_k%xcol(ircld%indexcol(icol - 1, jpk), j) - ircld_k%xcol(icol, j)
              ELSE
                ircld_k%xcol(ircld%indexcol(icol - 1, jpk), j) =  - ircld_k%xcol(icol, j)
              ENDIF
              IF (((ircld%indexcol(icol - 1, jpk) + 1) /= icol) .AND. (ircld%indexcol(icol - 1, jpk) /= icol)) THEN
                ircld_k%xcol(icol, j) = 0._jprb
              ENDIF
            ENDIF
          ENDDO
        ELSE
          ircld_k%xcolclr(j) = 0._jprb
        ENDIF
      ENDIF
  !---------Compute the weight of the clear column--------------------------------
      ircld_k%xcol(ircld%ncolumnref(jpk) + 1, j) = ircld_k%xcol(ircld%ncolumnref(jpk) + 1, j) - ircld_k%xcolclr(j)
      ircld_k%xcol(1, j)                         = ircld_k%xcol(1, j) + ircld_k%xcolclr(j)
  !---------Re-arrange the limits of each column in ascending order---------------
      outer : DO i = ircld%iloop(jpk), 1, -1
        icount1ref(i) = ircld%icount1ref(i, jpk)
        inner : DO icol = ircld%iloopin(i, jpk), 1, -1
          IF (ircld%xcolref2(i, icol, jpk) == ircld%xcolref2(i, icol + 1, jpk)) THEN
            IF (ircld%xcolref2(i, icol, jpk) < 1._jprb) THEN
              icount1ref(i) = icount1ref(i) - 1
              DO ijcol = icount1ref(i), icol, -1
                ircld_k%xcol(ijcol + 1, j) = ircld_k%xcol(ijcol + 1, j) + ircld_k%xcol(ijcol, j)
                ircld_k%xcol(ijcol, j)     = 0._jprb
              ENDDO
            ENDIF
          ENDIF
        ENDDO inner
      ENDDO outer
      loop1 : DO icol = ircld%icount(jpk) - 1, 2, -1
        IF (.NOT. ircld%flag(icol, jpk)) THEN
          ircld_k%a(icol, j) = ircld_k%a(icol, j) + ircld_k%xcol(1, j)
          ircld_k%xcol(1, j) = 0._jprb                                    !
        ENDIF
        DO i = ircld%iflag(icol, jpk), icol - 1
          IF (ircld%xcolref1(icol, i, jpk) <= ircld%a(icol, jpk)) THEN
            ircld_k%a(icol, j)     = ircld_k%a(icol, j) + ircld_k%xcol(i + 1, j)
            ircld_k%xcol(i + 1, j) = 0._jprb
          ELSE
            ircld_k%xcol(i, j)     = ircld_k%xcol(i, j) + ircld_k%xcol(i + 1, j)
            ircld_k%xcol(i + 1, j) = 0._jprb
          ENDIF
        ENDDO
        ircld_k%xcol(icol, j) = ircld_k%xcol(icol, j) + ircld_k%a(icol, j)
        ircld_k%a(icol, j)    = 0._jprb
      ENDDO loop1
  !---------Determine the limits of each column----------------------------------
      ic = ircld%icount(jpk)
      DO ilay = profiles(1)%nlayers, 1, -1
        IF (ircld%xcolmax(ilay, jpk) > 0._jprb) THEN
          ic = ic - 1
          ircld_k%xcolmax(ilay, j) = ircld_k%xcolmax(ilay, j) + ircld_k%xcol(ic, j)
          ic = ic - 1
          ircld_k%xcolmin(ilay, j) = ircld_k%xcolmin(ilay, j) + ircld_k%xcol(ic, j)
        ENDIF
        IF (ircld%xcolminref(ilay, jpk) < 0._jprb) THEN
          ircld_k%xcolmin(ilay, j) = 0._jprb
        ENDIF
        ntot(ilay)              = ntot(ilay) + ircld_k%xcolmin(ilay, j)
        ircld_k%cldcfr(ilay, j) = ircld_k%cldcfr(ilay, j) - ircld_k%xcolmin(ilay, j)
        ntot(ilay)              = ntot(ilay) + ircld_k%xcolmax(ilay, j)
      ENDDO
  !---------Compute the cumulative cloud coverage using the maximum-random---------
  !         overlap assumption
      DO ilay = profiles(1)%nlayers, 1, -1
        IF (ircld%cldcfr(ilay, jpk) > 0._jprb) THEN
          ntot(ilay) =  - ntot(ilay)
        ELSE
          ntot(ilay) = 0._jprb
        ENDIF
      ENDDO
      DO ilay = profiles(1)%nlayers, 2, -1
        IF (ircld%ntotref(ilay - 1, jpk) == 0._jprb .AND. ABS(ircld%cldcfr(ilay - 1, jpk) - 1._jprb) < realtol) THEN
          ircld_k%maxcov(ilay, j) = ircld_k%maxcov(ilay, j) + ntot(ilay)
        ELSEIF (ABS(ircld%cldcfr(ilay - 1, jpk) - 1._jprb) < realtol) THEN
          ntot(ilay - 1)              =      &
              ntot(ilay - 1) + ntot(ilay)
        ELSE
          ntot(ilay - 1)              =      &
              ntot(ilay - 1) + ntot(ilay) * ircld%maxcov(ilay, jpk) / (1._jprb - ircld%cldcfr(ilay - 1, jpk))
          ircld_k%cldcfr(ilay - 1, j) = ircld_k%cldcfr(ilay - 1, j) +                &
              ntot(ilay) * ircld%ntotref(ilay - 1, jpk) * ircld%maxcov(ilay, jpk) /  &
              (1._jprb - ircld%cldcfr(ilay - 1, jpk)) ** 2
          ircld_k%maxcov(ilay, j)     =      &
              ircld_k%maxcov(ilay, j) + ntot(ilay) * ircld%ntotref(ilay - 1, jpk) / (1._jprb - ircld%cldcfr(ilay - 1, jpk))
        ENDIF
        IF ((ircld%cldcfr(ilay - 1, jpk) > ircld%cldcfr(ilay, jpk))) THEN
          ircld_k%cldcfr(ilay - 1, j) = ircld_k%cldcfr(ilay - 1, j) - ircld_k%maxcov(ilay, j)
          ircld_k%maxcov(ilay, j)     = 0._jprb
        ELSE
          ircld_k%cldcfr(ilay, j) = ircld_k%cldcfr(ilay, j) - ircld_k%maxcov(ilay, j)
          ircld_k%maxcov(ilay, j) = 0._jprb
        ENDIF
      ENDDO
      ircld_k%cldcfr(1, j) = ircld_k%cldcfr(1, j) - ntot(1)

      DO ilay = profiles(1)%nlayers, 1, -1
        IF (profiles(jpk)%cfrac(ilay) > 0._jprb) THEN
          profiles_k(j)%cfrac(ilay) = profiles_k(j)%cfrac(ilay) + ircld_k%cldcfr(ilay, j)
        ENDIF
      ENDDO
      ircld_k%cldcfr(:, j)  = 0._jprb
      ircld_k%xcolmin(:, j) = 0._jprb
      ircld_k%xcolmax(:, j) = 0._jprb
      ircld_k%xcol(:, j) = 0._jprb
      ntot = 0._jprb
    ENDDO
  ENDIF
  IF (LHOOK) CALL DR_HOOK('RTTOV_CLOUD_OVERLAP_K', 1_jpim, ZHOOK_HANDLE)
!      ircld_k(:)%xcolclr=0._jprb
END SUBROUTINE rttov_cloud_overlap_k
