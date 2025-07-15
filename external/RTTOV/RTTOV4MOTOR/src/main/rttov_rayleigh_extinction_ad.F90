! Description:
!> @file
!!   AD of Rayleigh scattering extinction calculation.
!!
!> @brief
!!   AD of Rayleigh scattering extinction calculation.
!!
!! @param[in]     opts                  RTTOV options
!! @param[in]     chanprof              specifies channels and profiles to simulate
!! @param[in]     do_rayleigh_dom       flag to indicate DOM Rayleigh simulation
!! @param[in]     solar                 per-channel flag to indicate if solar simulations are being performed
!! @param[in]     profiles              input profiles
!! @param[in,out] profiles_ad           input profile increments
!! @param[in]     coef                  rttov_coef structure
!! @param[in]     raytracing            raytracing structure
!! @param[in,out] raytracing_ad         raytracing increments
!! @param[in]     od_sun_path2_ad       total level-space solar optical depth increments on path2
!! @param[in]     trans_scatt_ir_ad     VIS/IR scattering transmittance increments
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
SUBROUTINE rttov_rayleigh_extinction_ad(opts,            &
                                        chanprof,        &
                                        do_rayleigh_dom, &
                                        solar,           &
                                        profiles,        &
                                        profiles_ad,     &
                                        coef,            &
                                        raytracing,      &
                                        raytracing_ad,   &
                                        od_sun_path2_ad, &
                                        trans_scatt_ir_ad)

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
  TYPE(rttov_profile),               INTENT(INOUT) :: profiles_ad(:)
  TYPE(rttov_coef),                  INTENT(IN)    :: coef
  TYPE(rttov_raytracing),            INTENT(IN)    :: raytracing
  TYPE(rttov_raytracing),            INTENT(INOUT) :: raytracing_ad
  REAL(jprb),                        INTENT(IN)    :: od_sun_path2_ad(:,:)
  TYPE(rttov_transmission_scatt_ir), INTENT(IN)    :: trans_scatt_ir_ad
!INTF_END

  INTEGER(jpim) :: j, lay, lev, chan, prof, nchanprof
  REAL(jprb)    :: fac, od_lay, od_lay_ad, od_ray_ad(profiles(1)%nlevels)

  nchanprof = SIZE(chanprof)

  DO j = 1, nchanprof
    prof = chanprof(j)%prof
    chan = chanprof(j)%chan

    ! Skip calculation for channels above max wavelength
    IF (.NOT. solar(j) .OR. 10000._jprb / coef%ff_cwn(chan) > opts%rt_ir%rayleigh_max_wavelength) CYCLE

    fac = coef%ss_rayleigh_ext(chan) * (ray_ts / ray_ps) * rgc / (gravity * mair)

    IF (.NOT. do_rayleigh_dom) od_ray_ad(:) = -od_sun_path2_ad(:,j)

    DO lev = profiles(1)%nlevels, 2, -1
      ! Skip calculation for whole layers above min pressure
      IF (profiles(prof)%p(lev) < opts%rt_ir%rayleigh_min_pressure) CYCLE

      lay = lev - 1

      IF (do_rayleigh_dom) THEN
        od_lay_ad = trans_scatt_ir_ad%ray_sca(lay,j)
      ELSE
        od_lay = fac * (profiles(prof)%p(lev) - profiles(prof)%p(lev-1))

        od_lay_ad = od_ray_ad(lev) * raytracing%patheff(lay,prof)
        raytracing_ad%patheff(lay,prof) = raytracing_ad%patheff(lay,prof) + &
                                          od_ray_ad(lev) * od_lay
        od_ray_ad(lev-1) = od_ray_ad(lev-1) + od_ray_ad(lev)
      ENDIF

      IF (opts%interpolation%lgradp) THEN
        profiles_ad(prof)%p(lev) = profiles_ad(prof)%p(lev) + fac * od_lay_ad
        profiles_ad(prof)%p(lev-1) = profiles_ad(prof)%p(lev-1) - fac * od_lay_ad
      ENDIF
    ENDDO
    od_ray_ad(1) = 0._jprb

  ENDDO ! chanprof

END SUBROUTINE rttov_rayleigh_extinction_ad
