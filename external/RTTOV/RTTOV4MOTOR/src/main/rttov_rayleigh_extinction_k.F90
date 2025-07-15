! Description:
!> @file
!!   K of Rayleigh scattering extinction calculation.
!!
!> @brief
!!   K of Rayleigh scattering extinction calculation.
!!
!! @param[in]     opts                  RTTOV options
!! @param[in]     chanprof              specifies channels and profiles to simulate
!! @param[in]     do_rayleigh_dom       flag to indicate DOM Rayleigh simulation
!! @param[in]     solar                 per-channel flag to indicate if solar simulations are being performed
!! @param[in]     profiles              input profiles
!! @param[in,out] profiles_k            input profile increments
!! @param[in]     coef                  rttov_coef structure
!! @param[in]     raytracing            raytracing structure
!! @param[in,out] raytracing_k          raytracing increments
!! @param[in]     od_sun_path2_k        total level-space solar optical depth increments on path2
!! @param[in]     trans_scatt_ir_k      VIS/IR scattering transmittance increments
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
SUBROUTINE rttov_rayleigh_extinction_k(opts,            &
                                       chanprof,        &
                                       do_rayleigh_dom, &
                                       solar,           &
                                       profiles,        &
                                       profiles_k,      &
                                       coef,            &
                                       raytracing,      &
                                       raytracing_k,    &
                                       od_sun_path2_k,  &
                                       trans_scatt_ir_k)

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
  TYPE(rttov_profile),               INTENT(INOUT) :: profiles_k(:)
  TYPE(rttov_coef),                  INTENT(IN)    :: coef
  TYPE(rttov_raytracing),            INTENT(IN)    :: raytracing
  TYPE(rttov_raytracing),            INTENT(INOUT) :: raytracing_k
  REAL(jprb),                        INTENT(IN)    :: od_sun_path2_k(:,:)
  TYPE(rttov_transmission_scatt_ir), INTENT(IN)    :: trans_scatt_ir_k
!INTF_END

  INTEGER(jpim) :: j, lay, lev, chan, prof, nchanprof
  REAL(jprb)    :: fac, od_lay, od_lay_k, od_ray_k(profiles(1)%nlevels)

  nchanprof = SIZE(chanprof)

  DO j = 1, nchanprof
    prof = chanprof(j)%prof
    chan = chanprof(j)%chan

    ! Skip calculation for channels above max wavelength
    IF (.NOT. solar(j) .OR. 10000._jprb / coef%ff_cwn(chan) > opts%rt_ir%rayleigh_max_wavelength) CYCLE

    fac = coef%ss_rayleigh_ext(chan) * (ray_ts / ray_ps) * rgc / (gravity * mair)

    IF (.NOT. do_rayleigh_dom) od_ray_k(:) = -od_sun_path2_k(:,j)

    DO lev = profiles(1)%nlevels, 2, -1
      ! Skip calculation for whole layers above min pressure
      IF (profiles(prof)%p(lev) < opts%rt_ir%rayleigh_min_pressure) CYCLE

      lay = lev - 1

      IF (do_rayleigh_dom) THEN
        od_lay_k = trans_scatt_ir_k%ray_sca(lay,j)
      ELSE
        od_lay = fac * (profiles(prof)%p(lev) - profiles(prof)%p(lev-1))

        od_lay_k = od_ray_k(lev) * raytracing%patheff(lay,prof)
        raytracing_k%patheff(lay,j) = raytracing_k%patheff(lay,j) + &
                                      od_ray_k(lev) * od_lay
        od_ray_k(lev-1) = od_ray_k(lev-1) + od_ray_k(lev)
      ENDIF

      IF (opts%interpolation%lgradp) THEN
        profiles_k(j)%p(lev) = profiles_k(j)%p(lev) + fac * od_lay_k
        profiles_k(j)%p(lev-1) = profiles_k(j)%p(lev-1) - fac * od_lay_k
      ENDIF
    ENDDO
    od_ray_k(1) = 0._jprb

  ENDDO ! chanprof

END SUBROUTINE rttov_rayleigh_extinction_k
