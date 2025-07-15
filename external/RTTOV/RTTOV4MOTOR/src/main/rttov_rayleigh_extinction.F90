! Description:
!> @file
!!   Calculates Rayleigh scattering extinction.
!!
!> @brief
!!   Calculates Rayleigh scattering extinction.
!!
!! @details
!!   This calculation can only be used with v13 predictors (the LBLRTM
!!   calculations for v9 predictors included Rayleigh extinction in the mixed
!!   gases).
!!
!!   The parameterisation is applied on user levels and follows Bucholtz (1995).
!!   The channel-averaged values of the Bucholtz parameterised volume extinction
!!   coefficients are stored in the coefficient file. These are used to compute
!!   the Rayleigh extinction optical depth in each layer. No account is taken
!!   of the polychromatic nature of the transmittances beyond the calculation of
!!   the channel-averaged extinction coefficients: in practice the errors due to
!!   this approximation are small and the calculation is fast.
!!
!!   The extinction coefficients b_s in the rtcoef files are valid at a given
!!   temperature and pressure (ray_ts = t_s, ray_ps = p_s in rttov_const). From
!!   the paper, the extinction coefficient b at general temperature t and
!!   pressure p is:
!!
!!     b = b_s * (p / t) * (t_s / p_s)
!!
!!   The layer optical depth is given by integrating the extinction coefficient
!!   wrt height over the layer. Applying the ideal gas law and using the
!!   hydrostatic equation we obtain the following expression for the layer
!!   optical depth tau:
!!
!!     tau = b_s * (t_s / p_s) * R / (g * Mair) * delta-p * sec(theta)
!!
!!   where
!!      R is the gas constant
!!      g is gravitational acceleration (assumed constant)
!!      Mair is molar mass of dry air (ignore variable water vapour as impact
!!          is small)
!!      delta-p is the pressure difference across the layer
!!      theta is the local zenith angle of the view path in the layer
!!
!!   In the general case, the Rayleigh optical depths are added to the level-to-
!!   space solar gas optical depths, but if Rayleigh DOM calculations are
!!   enabled then the nadir Rayleigh extinction is instead stored in a separate
!!   array for use later on.
!!
!!   Bucholtz, A., 1995: Rayleigh-scattering calculations for the terrestrial
!!   atmosphere. Applied Optics, 34, 15, 2765-2773.
!!
!! @param[in]     opts                  RTTOV options
!! @param[in]     chanprof              specifies channels and profiles to simulate
!! @param[in]     do_rayleigh_dom       flag to indicate DOM Rayleigh simulation
!! @param[in]     solar                 per-channel flag to indicate if solar simulations are being performed
!! @param[in]     profiles              input profiles
!! @param[in]     coef                  rttov_coef structure
!! @param[in]     raytracing            RTTOV raytracing structure
!! @param[in,out] od_sun_path2          Total level-space solar optical depths on path2
!! @param[in,out] trans_scatt_ir        VIS/IR scattering transmittance structure (Rayleigh scattering nadir
!!                                          optical depth for DOM Rayleigh)
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
!    Copyright 2019, EUMETSAT, All Rights Reserved.
!
SUBROUTINE rttov_rayleigh_extinction(opts,            &
                                     chanprof,        &
                                     do_rayleigh_dom, &
                                     solar,           &
                                     profiles,        &
                                     coef,            &
                                     raytracing,      &
                                     od_sun_path2,    &
                                     trans_scatt_ir)

  USE parkind1, ONLY : jprb, jplm
  USE rttov_types, ONLY : &
    rttov_options,    &
    rttov_chanprof,   &
    rttov_profile,    &
    rttov_coef,       &
    rttov_raytracing, &
    rttov_transmission_scatt_ir
!INTF_OFF
  USE parkind1, ONLY : jpim
  USE rttov_const, ONLY : ray_ps, ray_ts, rgc, gravity, mair
!INTF_ON
  IMPLICIT NONE

  TYPE(rttov_options),               INTENT(IN)    :: opts
  TYPE(rttov_chanprof),              INTENT(IN)    :: chanprof(:)
  LOGICAL(jplm),                     INTENT(IN)    :: do_rayleigh_dom
  LOGICAL(jplm),                     INTENT(IN)    :: solar(:)
  TYPE(rttov_profile),               INTENT(IN)    :: profiles(:)
  TYPE(rttov_coef),                  INTENT(IN)    :: coef
  TYPE(rttov_raytracing),            INTENT(IN)    :: raytracing
  REAL(jprb),                        INTENT(INOUT) :: od_sun_path2(:,:)
  TYPE(rttov_transmission_scatt_ir), INTENT(INOUT) :: trans_scatt_ir
!INTF_END

  INTEGER(jpim) :: j, lay, lev, chan, prof, nchanprof
  REAL(jprb)    :: fac, od_lay, od_ray(profiles(1)%nlevels)

  nchanprof = SIZE(chanprof)

  DO j = 1, nchanprof
    prof = chanprof(j)%prof
    chan = chanprof(j)%chan

    ! Skip calculation for channels above max wavelength
    IF (.NOT. solar(j) .OR. 10000._jprb / coef%ff_cwn(chan) > opts%rt_ir%rayleigh_max_wavelength) CYCLE

    od_ray(1) = 0._jprb
    fac = coef%ss_rayleigh_ext(chan) * (ray_ts / ray_ps) * rgc / (gravity * mair)
    DO lev = 2, profiles(1)%nlevels
      ! Skip calculation for whole layers above min pressure
      IF (profiles(prof)%p(lev) < opts%rt_ir%rayleigh_min_pressure) THEN
        od_ray(lev) = 0._jprb
        CYCLE
      ENDIF

      lay = lev - 1

      ! Calculation of nadir layer optical depth
      od_lay = fac * (profiles(prof)%p(lev) - profiles(prof)%p(lev-1))

      IF (do_rayleigh_dom) THEN
        ! Store nadir Rayleigh scattering extinction for DOM Rayleigh
        trans_scatt_ir%ray_sca(lay,j) = od_lay
      ELSE
        ! Accumulate level-to-space optical depths (on path2)
        od_ray(lev) = od_ray(lev-1) + od_lay * raytracing%patheff(lay,prof)
      ENDIF
    ENDDO

    IF (.NOT. do_rayleigh_dom) THEN
      ! Accumulate Rayleigh extinction in total optical depths (NB these are negative)
      od_sun_path2(:,j) = od_sun_path2(:,j) - od_ray(:)
    ENDIF
  ENDDO ! chanprof

END SUBROUTINE rttov_rayleigh_extinction
