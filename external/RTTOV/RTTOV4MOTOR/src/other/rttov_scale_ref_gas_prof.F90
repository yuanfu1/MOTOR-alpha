! Description:
!> @file
!!   Populate optional trace gas profiles with scaled copies of the RTTOV
!!   reference (background) profiles.
!
!> @brief
!!   Populate optional trace gas profiles with scaled copies of the RTTOV
!!   reference (background) profiles.
!!
!! @details
!!   This subroutine makes it simple to populate an array of profile
!!   structures with optional trace gas profiles that are scaled
!!   copies of the RTTOV reference (background) profiles. This makes it
!!   easy, for example, to provide CO2 profiles with a lower maximum
!!   value representing an atmosphere from earlier in the satellite era.
!!
!!   The scaling is done either by specifying a column-integrated
!!   gas amount in kg/m^2 (or in Dobson units for O3) or by specifying
!!   the maximum value in ppmv over dry air for each gas that is required.
!!
!!   The inputs are: a coefficients structure which has been populated
!!   from a coefficients file (this is used to obtain the reference
!!   gas profiles - these are independent of the instrument); and an
!!   array of profile structures which should be populated with all
!!   profile data other than the trace gases to be provided by this
!!   subroutine.
!!
!!   The reference profile(s) for each gas for which an input quantity
!!   is specified are interpolated onto each profile's pressure levels,
!!   and the profile is scaled and the units converted to match those
!!   specified in gas_units. On output the profiles array contains 
!!   profiles for all specified trace gases.
!!
!!   Note that you must have set the relevant option (e.g. ozone_data or
!!   co2_data) in the options structure to indicate that you are passing
!!   profiles for the relevant gas into RTTOV: this ensures that the
!!   gas array(s) are allocated in the profiles structure.
!!
!!
!! @param[out]    err            status on exit
!! @param[in]     coefs          RTTOV coefficients structure
!! @param[inout]  profiles       array of profile structures
!! @param[in]     o3_col_int     total column-integrated O3 (kg/m^2), optional
!! @param[in]     o3_col_int_du  total column-integrated O3 (Dobson units), optional
!! @param[in]     o3_max_ppmv    maximum O3 value (ppmv over dry air), optional
!! @param[in]     co2_col_int    total column-integrated CO2 (kg/m^2), optional
!! @param[in]     co2_max_ppmv   maximum CO2 value (ppmv over dry air), optional
!! @param[in]     n2o_col_int    total column-integrated N2O (kg/m^2), optional
!! @param[in]     n2o_max_ppmv   maximum N2O value (ppmv over dry air), optional
!! @param[in]     co_col_int     total column-integrated CO (kg/m^2), optional
!! @param[in]     co_max_ppmv    maximum CO value (ppmv over dry air), optional
!! @param[in]     ch4_col_int    total column-integrated CH4 (kg/m^2), optional
!! @param[in]     ch4_max_ppmv   maximum CH4 value (ppmv over dry air), optional
!! @param[in]     so2_col_int    total column-integrated SO2 (kg/m^2), optional
!! @param[in]     so2_max_ppmv   maximum SO2 value (ppmv over dry air), optional
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
!    Copyright 2018, EUMETSAT, All Rights Reserved.
!
SUBROUTINE rttov_scale_ref_gas_prof( &
               err,           &
               coefs,         &
               profiles,      &
               o3_col_int,    &
               o3_col_int_du, &
               o3_max_ppmv,   &
               co2_col_int,   &
               co2_max_ppmv,  &
               n2o_col_int,   &
               n2o_max_ppmv,  &
               co_col_int,    &
               co_max_ppmv,   &
               ch4_col_int,   &
               ch4_max_ppmv,  &
               so2_col_int,   &
               so2_max_ppmv)

  USE parkind1, ONLY : jprb, jpim
  USE rttov_types, ONLY : rttov_coefs, rttov_profile
!INTF_OFF
#include "throw.h"
  USE parkind1, ONLY : jplm
  USE rttov_const, ONLY : &
    gas_id_watervapour, &
    gas_id_ozone,       &
    gas_id_co2,         &
    gas_id_n2o,         &
    gas_id_co,          &
    gas_id_ch4,         &
    gas_id_so2,         &
    mair, gas_mass,     &
    gas_unit_ppmv,      &
    gas_unit_specconc,  &
    gas_unit_ppmvdry,   &
    rho_ozone_stp,      &
    gravity
!INTF_ON

  IMPLICIT NONE

  INTEGER(jpim),       INTENT(OUT)          :: err
  TYPE(rttov_coefs),   INTENT(IN)           :: coefs
  TYPE(rttov_profile), INTENT(INOUT)        :: profiles(:)
  REAL(jprb),          INTENT(IN), OPTIONAL :: o3_col_int
  REAL(jprb),          INTENT(IN), OPTIONAL :: o3_col_int_du
  REAL(jprb),          INTENT(IN), OPTIONAL :: o3_max_ppmv
  REAL(jprb),          INTENT(IN), OPTIONAL :: co2_col_int
  REAL(jprb),          INTENT(IN), OPTIONAL :: co2_max_ppmv
  REAL(jprb),          INTENT(IN), OPTIONAL :: n2o_col_int
  REAL(jprb),          INTENT(IN), OPTIONAL :: n2o_max_ppmv
  REAL(jprb),          INTENT(IN), OPTIONAL :: co_col_int
  REAL(jprb),          INTENT(IN), OPTIONAL :: co_max_ppmv
  REAL(jprb),          INTENT(IN), OPTIONAL :: ch4_col_int
  REAL(jprb),          INTENT(IN), OPTIONAL :: ch4_max_ppmv
  REAL(jprb),          INTENT(IN), OPTIONAL :: so2_col_int
  REAL(jprb),          INTENT(IN), OPTIONAL :: so2_max_ppmv
!INTF_END

#include "rttov_errorreport.interface"

  INTEGER(jpim) :: prof, nlev
  LOGICAL(jplm) :: htfrtc
  REAL(jprb), ALLOCATABLE :: pref(:), gasref(:)

  TRY

  htfrtc = ASSOCIATED(coefs%coef_htfrtc%p)

  ! Check inputs against profile and coefficients

  IF ((PRESENT(o3_col_int) .OR. PRESENT(o3_col_int_du) .OR. PRESENT(o3_max_ppmv))) THEN
    IF (.NOT. ASSOCIATED(profiles(1)%o3)) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0,'O3 array not allocated in profiles structure, ensure ozone_data is true')
    ENDIF
    IF (coefs%coef%nozone == 0 .AND. .NOT. htfrtc) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0,'Coefficient file does not support variable O3')
    ENDIF
  ENDIF

  IF ((PRESENT(co2_col_int) .OR. PRESENT(co2_max_ppmv))) THEN
    IF (.NOT. ASSOCIATED(profiles(1)%co2)) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0,'CO2 array not allocated in profiles structure, ensure co2_data is true')
    ENDIF
    IF (coefs%coef%nco2 == 0 .AND. .NOT. htfrtc) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0,'Coefficient file does not support variable CO2')
    ENDIF
  ENDIF

  IF ((PRESENT(n2o_col_int) .OR. PRESENT(n2o_max_ppmv))) THEN
    IF (.NOT. ASSOCIATED(profiles(1)%n2o)) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0,'N2O array not allocated in profiles structure, ensure n2o_data is true')
    ENDIF
    IF (coefs%coef%nn2o == 0 .AND. .NOT. htfrtc) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0,'Coefficient file does not support variable N2O')
    ENDIF
  ENDIF

  IF ((PRESENT(co_col_int) .OR. PRESENT(co_max_ppmv))) THEN
    IF (.NOT. ASSOCIATED(profiles(1)%co)) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0,'CO array not allocated in profiles structure, ensure co_data is true')
    ENDIF
    IF (coefs%coef%nco == 0 .AND. .NOT. htfrtc) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0,'Coefficient file does not support variable CO')
    ENDIF
  ENDIF

  IF ((PRESENT(ch4_col_int) .OR. PRESENT(ch4_max_ppmv))) THEN
    IF (.NOT. ASSOCIATED(profiles(1)%ch4)) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0,'CH4 array not allocated in profiles structure, ensure ch4_data is true')
    ENDIF
    IF (coefs%coef%nch4 == 0 .AND. .NOT. htfrtc) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0,'Coefficient file does not support variable CH4')
    ENDIF
  ENDIF

  IF ((PRESENT(so2_col_int) .OR. PRESENT(so2_max_ppmv))) THEN
    IF (.NOT. ASSOCIATED(profiles(1)%so2)) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0,'SO2 array not allocated in profiles structure, ensure so2_data is true')
    ENDIF
    IF (coefs%coef%nso2 == 0 .AND. .NOT. htfrtc) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0,'Coefficient file does not support variable SO2')
    ENDIF
  ENDIF

  ! Allocate pressure and gas arrays

  IF (htfrtc) THEN
    nlev = coefs%coef_htfrtc%n_p + 1
    ALLOCATE(pref(nlev), gasref(nlev))
    pref(1) = 0.005_jprb
    pref(2:nlev) = coefs%coef_htfrtc%p
  ELSE
    nlev = SIZE(coefs%coef%ref_prfl_p)
    ALLOCATE(pref(nlev), gasref(nlev))
    pref = coefs%coef%ref_prfl_p
  ENDIF

  ! Calculate scaled profiles

  DO prof = 1, SIZE(profiles)

    IF ((PRESENT(o3_col_int) .OR. PRESENT(o3_col_int_du) .OR. PRESENT(o3_max_ppmv))) THEN
      IF (htfrtc) THEN
        gasref(2:nlev) = coefs%coef_htfrtc%mixed_ref_frac(:,1)
        gasref(1) = gasref(2)
      ELSE
        gasref = coefs%coef%bkg_prfl_mr(:,coefs%coef%fmv_gas_pos(gas_id_ozone))
      ENDIF
      CALL process_gas(err, gas_id_ozone, pref, gasref, &
                       profiles(prof)%gas_units, profiles(prof)%p, profiles(prof)%q, &
                       profiles(prof)%o3, o3_max_ppmv, o3_col_int, o3_col_int_du)
      THROW(err.NE.0)
    ENDIF

    IF ((PRESENT(co2_col_int) .OR. PRESENT(co2_max_ppmv))) THEN
      IF (htfrtc) THEN
        gasref(2:nlev) = coefs%coef_htfrtc%mixed_ref_frac(:,2)
        gasref(1) = gasref(2)
      ELSE
        gasref = coefs%coef%bkg_prfl_mr(:,coefs%coef%fmv_gas_pos(gas_id_co2))
      ENDIF
      CALL process_gas(err, gas_id_co2, pref, gasref, &
                       profiles(prof)%gas_units, profiles(prof)%p, profiles(prof)%q, &
                       profiles(prof)%co2, co2_max_ppmv, co2_col_int)
      THROW(err.NE.0)
    ENDIF

    IF ((PRESENT(n2o_col_int) .OR. PRESENT(n2o_max_ppmv))) THEN
      IF (htfrtc) THEN
        gasref(2:nlev) = coefs%coef_htfrtc%mixed_ref_frac(:,3)
        gasref(1) = gasref(2)
      ELSE
        gasref = coefs%coef%bkg_prfl_mr(:,coefs%coef%fmv_gas_pos(gas_id_n2o))
      ENDIF
      CALL process_gas(err, gas_id_n2o, pref, gasref, &
                       profiles(prof)%gas_units, profiles(prof)%p, profiles(prof)%q, &
                       profiles(prof)%n2o, n2o_max_ppmv, n2o_col_int)
      THROW(err.NE.0)
    ENDIF

    IF ((PRESENT(co_col_int) .OR. PRESENT(co_max_ppmv))) THEN
      IF (htfrtc) THEN
        gasref(2:nlev) = coefs%coef_htfrtc%mixed_ref_frac(:,4)
        gasref(1) = gasref(2)
      ELSE
        gasref = coefs%coef%bkg_prfl_mr(:,coefs%coef%fmv_gas_pos(gas_id_co))
      ENDIF
      CALL process_gas(err, gas_id_co, pref, gasref, &
                       profiles(prof)%gas_units, profiles(prof)%p, profiles(prof)%q, &
                       profiles(prof)%co, co_max_ppmv, co_col_int)
      THROW(err.NE.0)
    ENDIF

    IF ((PRESENT(ch4_col_int) .OR. PRESENT(ch4_max_ppmv))) THEN
      IF (htfrtc) THEN
        gasref(2:nlev) = coefs%coef_htfrtc%mixed_ref_frac(:,5)
        gasref(1) = gasref(2)
      ELSE
        gasref = coefs%coef%bkg_prfl_mr(:,coefs%coef%fmv_gas_pos(gas_id_ch4))
      ENDIF
      CALL process_gas(err, gas_id_ch4, pref, gasref, &
                       profiles(prof)%gas_units, profiles(prof)%p, profiles(prof)%q, &
                       profiles(prof)%ch4, ch4_max_ppmv, ch4_col_int)
      THROW(err.NE.0)
    ENDIF

    IF ((PRESENT(so2_col_int) .OR. PRESENT(so2_max_ppmv))) THEN
      IF (htfrtc) THEN
        gasref(2:nlev) = coefs%coef_htfrtc%mixed_ref_frac(:,6)
        gasref(1) = gasref(2)
      ELSE
        gasref = coefs%coef%bkg_prfl_mr(:,coefs%coef%fmv_gas_pos(gas_id_so2))
      ENDIF
      CALL process_gas(err, gas_id_so2, pref, gasref, &
                       profiles(prof)%gas_units, profiles(prof)%p, profiles(prof)%q, &
                       profiles(prof)%so2, so2_max_ppmv, so2_col_int)
      THROW(err.NE.0)
    ENDIF

  ENDDO

  DEALLOCATE(pref, gasref)

  CATCH

CONTAINS

  SUBROUTINE interp(p_ref, p_out, gas_ref, gas_out)

    REAL(jprb),       INTENT(IN)           :: p_ref(:)
    REAL(jprb),       INTENT(IN)           :: p_out(:)
    REAL(jprb),       INTENT(IN)           :: gas_ref(:)
    REAL(jprb),       INTENT(INOUT)        :: gas_out(:)

    INTEGER(jpim) :: j, k, nlevels_coef, nlevels

    nlevels = profiles(1)%nlevels
    nlevels_coef = SIZE(p_ref)

    k = 2
    DO j = 1, nlevels
      IF (p_out(j) <= p_ref(1)) THEN
        gas_out(j) = gas_ref(1)
      ELSE IF (p_out(j) >= p_ref(nlevels_coef)) THEN
        gas_out(j) = gas_ref(nlevels_coef)
      ELSE
        DO
          IF (p_out(j) <= p_ref(k)) EXIT
          k = k + 1
        ENDDO
        gas_out(j) = ((p_ref(k) - p_out(j)) * gas_ref(k-1) + &
                      (p_out(j) - p_ref(k-1)) * gas_ref(k)) / (p_ref(k) - p_ref(k-1))
      ENDIF
    ENDDO

  END SUBROUTINE

  SUBROUTINE process_gas(err, gas_id, p_ref, gas_ref, gas_units, &
                         p_out, q_out, gas_out, max_ppmv, col_int, col_int_du)

    INTEGER(jpim),    INTENT(OUT)          :: err
    INTEGER(jpim),    INTENT(IN)           :: gas_id
    REAL(jprb),       INTENT(IN)           :: p_ref(:)
    REAL(jprb),       INTENT(IN)           :: gas_ref(:)
    INTEGER(jpim),    INTENT(IN)           :: gas_units
    REAL(jprb),       INTENT(IN)           :: p_out(:)
    REAL(jprb),       INTENT(IN)           :: q_out(:)
    REAL(jprb),       INTENT(INOUT)        :: gas_out(:)
    REAL(jprb),       INTENT(IN), OPTIONAL :: max_ppmv
    REAL(jprb),       INTENT(IN), OPTIONAL :: col_int
    REAL(jprb),       INTENT(IN), OPTIONAL :: col_int_du

    INTEGER(jpim) :: i, nlevels
    REAL(jprb) :: scale_factor, col_int_ref

    TRY

    ! gas_ref is in units of ppmv dry
    ! q_out has units specified in gas_units
    ! gas_out must also be in gas_units

    IF (PRESENT(col_int) .OR. PRESENT(col_int_du)) THEN
      ! ---------------------------------
      ! Scale to column-integrated value
      ! ---------------------------------

      ! Interpolate reference profile onto profile levels

      CALL interp(p_ref, p_out, gas_ref, gas_out)

      ! Convert reference profile from ppmv dry to kg/kg moist

      IF (gas_units == gas_unit_ppmv) THEN
        ! q_out is ppmv moist
        gas_out(:) = gas_mass(gas_id) * gas_out(:) * (1._jprb - 1.E-06_jprb * q_out(:)) / &
                     (mair * 1.E06_jprb + (gas_mass(gas_id_watervapour) - mair) * q_out(:))

      ELSE IF (gas_units == gas_unit_specconc) THEN
        ! q_out is kg/kg moist
        gas_out(:) = gas_mass(gas_id) * gas_out(:) * (1._jprb - q_out(:)) / (mair * 1.E06_jprb)

      ELSE IF (gas_units <= gas_unit_ppmvdry) THEN
        ! q_out is ppmv dry
        gas_out(:) = gas_mass(gas_id) * gas_out(:) / &
                     (mair * 1.E06_jprb + gas_mass(gas_id_watervapour) * q_out(:))

      ELSE
        err = errorstatus_fatal
        THROWM(err.NE.0,'Unknown gas_units')
      ENDIF

      ! Calculate column-integrated value in kg/m^2 (follows OPS code)
      !   SUM(delta_pressure * layer_mean_gas)

      nlevels = profiles(1)%nlevels
      col_int_ref = 0._jprb
      DO i = 1, nlevels - 1
        col_int_ref = col_int_ref + 100._jprb * (p_out(i+1) - p_out(i)) * &
                                    0.5_jprb * (gas_out(i) + gas_out(i+1))
      END DO
      col_int_ref = col_int_ref / gravity

      ! Compute scale factor

      IF (PRESENT(col_int_du)) THEN
        ! Convert Dobson units to kg/m^2
        scale_factor = col_int_du * rho_ozone_stp / col_int_ref
      ELSE
        scale_factor = col_int / col_int_ref
      ENDIF

      ! Scale profile and convert units if required (gas_out is kg/kg moist)

      IF (gas_units == gas_unit_ppmv) THEN
        ! kg/kg moist to ppmv moist (q_out is ppmv moist)
        gas_out(:) = scale_factor * gas_out(:) * &
                     (mair * 1.E06_jprb + (gas_mass(gas_id_watervapour) - mair) * q_out(:)) / gas_mass(gas_id)

      ELSE IF (gas_units <= gas_unit_ppmvdry) THEN
        ! kg/kg moist to ppmv dry (q_out is ppmv dry)
        gas_out(:) = scale_factor * gas_out(:) * &
                     (mair * 1.E06_jprb + gas_mass(gas_id_watervapour) * q_out(:)) / gas_mass(gas_id)

      ELSE IF (gas_units == gas_unit_specconc) THEN
        ! No conversion required
        gas_out(:) = scale_factor * gas_out(:)

      ELSE
        err = errorstatus_fatal
        THROWM(err.NE.0,'Unknown gas_units')
      ENDIF

    ELSE IF (PRESENT(max_ppmv)) THEN
      ! ---------------------------------
      ! Scale to maximum ppmv value
      ! ---------------------------------

      ! Compute scale factor

      scale_factor = max_ppmv / MAXVAL(gas_ref(:))

      ! Interpolate reference profile onto profile levels

      CALL interp(p_ref, p_out, gas_ref, gas_out)

      ! Scale profile and convert units if required (gas_out is ppmv dry)

      IF (gas_units == gas_unit_ppmv) THEN
        ! ppmv dry to ppmv moist (q_out is ppmv moist)
        gas_out(:) = scale_factor * gas_out(:) * (1._jprb - 1.E-06_jprb * q_out(:))

      ELSE IF (gas_units == gas_unit_specconc) THEN
        ! ppmv dry to kg/kg moist (q_out is kg/kg moist)
        gas_out(:) = scale_factor * gas_mass(gas_id) * gas_out(:) * (1._jprb - q_out(:)) / (mair * 1.E06_jprb)

      ELSE IF (gas_units <= gas_unit_ppmvdry) THEN
        ! No conversion required
        gas_out(:) = scale_factor * gas_out(:)

      ELSE
        err = errorstatus_fatal
        THROWM(err.NE.0,'Unknown gas_units')
      ENDIF

    ENDIF

    CATCH
  END SUBROUTINE

END SUBROUTINE rttov_scale_ref_gas_prof
