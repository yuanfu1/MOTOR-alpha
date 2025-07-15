! Description:
!> @file
!!   TL of Rayleigh scattering extinction calculation.
!!
!> @brief
!!   TL of Rayleigh scattering extinction calculation.
!!
!! @param[in]     opts                  RTTOV options
!! @param[in]     chanprof              specifies channels and profiles to simulate
!! @param[in]     do_rayleigh_dom       flag to indicate DOM Rayleigh simulation
!! @param[in]     solar                 per-channel flag to indicate if solar simulations are being performed
!! @param[in]     profiles              input profiles
!! @param[in]     profiles_tl           input profile perturbations
!! @param[in]     coef                  rttov_coef structure
!! @param[in]     raytracing            raytracing structure
!! @param[in]     raytracing_tl         raytracing structure containing perturbations
!! @param[in,out] od_sun_path2_tl       total level-space solar optical depth perturbations on path2
!! @param[in,out] trans_scatt_ir_tl     VIS/IR scattering transmittance perturbations
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
SUBROUTINE rttov_rayleigh_extinction_tl(opts,            &
                                        chanprof,        &
                                        do_rayleigh_dom, &
                                        solar,           &
                                        profiles,        &
                                        profiles_tl,     &
                                        coef,            &
                                        raytracing,      &
                                        raytracing_tl,   &
                                        od_sun_path2_tl, &
                                        trans_scatt_ir_tl)

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
  TYPE(rttov_profile),               INTENT(IN)    :: profiles_tl(:)
  TYPE(rttov_coef),                  INTENT(IN)    :: coef
  TYPE(rttov_raytracing),            INTENT(IN)    :: raytracing
  TYPE(rttov_raytracing),            INTENT(IN)    :: raytracing_tl
  REAL(jprb),                        INTENT(INOUT) :: od_sun_path2_tl(:,:)
  TYPE(rttov_transmission_scatt_ir), INTENT(INOUT) :: trans_scatt_ir_tl
!INTF_END

  INTEGER(jpim) :: j, lay, lev, chan, prof, nchanprof
  REAL(jprb)    :: fac, od_lay, od_lay_tl, od_ray_tl(profiles(1)%nlevels)

  nchanprof = SIZE(chanprof)

  DO j = 1, nchanprof
    prof = chanprof(j)%prof
    chan = chanprof(j)%chan

    ! Skip calculation for channels above max wavelength
    IF (.NOT. solar(j) .OR. 10000._jprb / coef%ff_cwn(chan) > opts%rt_ir%rayleigh_max_wavelength) CYCLE

    od_ray_tl(1) = 0._jprb
    fac = coef%ss_rayleigh_ext(chan) * (ray_ts / ray_ps) * rgc / (gravity * mair)
    DO lev = 2, profiles(1)%nlevels
      ! Skip calculation for whole layers above min pressure
      IF (profiles(prof)%p(lev) < opts%rt_ir%rayleigh_min_pressure) THEN
        od_ray_tl(lev) = 0._jprb
        CYCLE
      ENDIF

      lay = lev - 1

      ! Calculation of nadir layer optical depth
      od_lay_tl = 0._jprb

      IF (opts%interpolation%lgradp) THEN
        od_lay_tl = fac * (profiles_tl(prof)%p(lev) - profiles_tl(prof)%p(lev-1))
      ENDIF

      IF (do_rayleigh_dom) THEN
        ! Store nadir Rayleigh scattering extinction for DOM Rayleigh
        trans_scatt_ir_tl%ray_sca(lay,j) = od_lay_tl
      ELSE
        ! Accumulate level-to-space optical depths (on path2)
        od_lay = fac * (profiles(prof)%p(lev) - profiles(prof)%p(lev-1))
        od_ray_tl(lev) = od_ray_tl(lev-1) + &
                         od_lay_tl * raytracing%patheff(lay,prof) + &
                         od_lay * raytracing_tl%patheff(lay,prof)
      ENDIF
    ENDDO

    IF (.NOT. do_rayleigh_dom) THEN
      ! Accumulate Rayleigh extinction in total optical depths (NB these are negative)
      od_sun_path2_tl(:,j) = od_sun_path2_tl(:,j) - od_ray_tl(:)
    ENDIF
  ENDDO ! chanprof

END SUBROUTINE rttov_rayleigh_extinction_tl
