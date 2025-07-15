! Description:
!> @file
!!   Modify coefficient regression limits.
!
!> @brief
!!   Modify coefficient regression limits.
!!
!> @details
!!   The coefficient files contain the min/max envelope of the profile set
!!   which is used to train RTTOV. The regression limits which are used to
!!   trigger warnings about exceeding regression limits or to which input
!!   profiles are clipped if apply_reg_limits is TRUE are derived from the
!!   envelope.
!!
!!   By default the temperature minimum/maximum limits are 90%/110% of the
!!   envelope min/max respectively. For gases the stretch is 20% each way
!!   (i.e. 80%/120% for min/max).
!!
!!   This subroutine allows the profile limits to be modified by applying
!!   different scaling factors. Separate scaling factors can be specified for
!!   temperature and each gas (only variable gases supported by the coefficient
!!   file have any effect).
!!
!!   Each factor supplied is applied as follows:
!!      limit_minimum = envelope_minimum * (1. - factor)
!!      limit_maximum = envelope_maximum * (1. + factor)
!!
!!   In the default case factor is 0.1 for temperature and 0.2 for gases.
!!
!!   Only limits for which factors are supplied are modified: other profile
!!   limits remain at their default (or previously modified) values.
!!
!! @param[out]    err             status on exit
!! @param[in,out] coefs           RTTOV coefficients structure (already loaded)
!! @param[in]     tfac            factor for calculating temperature profile limits, optional
!! @param[in]     qfac            factor for calculating water vapour profile limits, optional
!! @param[in]     o3fac           factor for calculating ozone profile limits, optional
!! @param[in]     co2fac          factor for calculating CO2 profile limits, optional
!! @param[in]     n2ofac          factor for calculating N2O profile limits, optional
!! @param[in]     cofac           factor for calculating CO profile limits, optional
!! @param[in]     ch4fac          factor for calculating CH4 profile limits, optional
!! @param[in]     so2fac          factor for calculating SO2 profile limits, optional
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
SUBROUTINE rttov_set_coef_limits(err, coefs, tfac, qfac, o3fac, co2fac, n2ofac, cofac, ch4fac, so2fac)

#include "throw.h"

  USE rttov_types, ONLY : rttov_coefs
  USE parkind1, ONLY : jpim, jprb
!INTF_OFF
  USE rttov_const, ONLY :  &
        version_compatible_min, &
        version_compatible_max, &
        gas_id_watervapour,     &
        gas_id_ozone,           &
        gas_id_co2,             &
        gas_id_n2o,             &
        gas_id_co,              &
        gas_id_ch4,             &
        gas_id_so2
!INTF_ON
  IMPLICIT NONE

  INTEGER(jpim),     INTENT(OUT)            :: err
  TYPE(rttov_coefs), INTENT(INOUT)          :: coefs
  REAL(jprb),        INTENT(IN),   OPTIONAL :: tfac
  REAL(jprb),        INTENT(IN),   OPTIONAL :: qfac
  REAL(jprb),        INTENT(IN),   OPTIONAL :: o3fac
  REAL(jprb),        INTENT(IN),   OPTIONAL :: co2fac
  REAL(jprb),        INTENT(IN),   OPTIONAL :: n2ofac
  REAL(jprb),        INTENT(IN),   OPTIONAL :: cofac
  REAL(jprb),        INTENT(IN),   OPTIONAL :: ch4fac
  REAL(jprb),        INTENT(IN),   OPTIONAL :: so2fac
!INTF_END

  INTEGER(jpim) :: i, j

#include "rttov_errorreport.interface"

  TRY

  IF (coefs%coef%id_comp_lvl < version_compatible_min .OR. &
      coefs%coef%id_comp_lvl > version_compatible_max) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, "Version of coefficient file is incompatible with RTTOV library")
  ENDIF

  IF (PRESENT(tfac)) THEN
    coefs%coef%lim_prfl_tmin = coefs%coef%env_prfl_tmin * (1._jprb - tfac)
    coefs%coef%lim_prfl_tmax = coefs%coef%env_prfl_tmax * (1._jprb + tfac)
  ENDIF

  DO i = 1, coefs%coef%fmv_gas
    SELECT CASE (coefs%coef%fmv_gas_id(i))
    CASE (gas_id_watervapour)
      IF (PRESENT(qfac)) THEN
        j = coefs%coef%fmv_gas_pos(gas_id_watervapour)
        coefs%coef%lim_prfl_gmin(:,j) = coefs%coef%env_prfl_gmin(:,j) * (1._jprb - qfac)
        coefs%coef%lim_prfl_gmax(:,j) = coefs%coef%env_prfl_gmax(:,j) * (1._jprb + qfac)
      ENDIF
    CASE (gas_id_ozone)
      IF (PRESENT(o3fac)) THEN
        j = coefs%coef%fmv_gas_pos(gas_id_ozone)
        coefs%coef%lim_prfl_gmin(:,j) = coefs%coef%env_prfl_gmin(:,j) * (1._jprb - o3fac)
        coefs%coef%lim_prfl_gmax(:,j) = coefs%coef%env_prfl_gmax(:,j) * (1._jprb + o3fac)
      ENDIF
    CASE (gas_id_co2)
      IF (PRESENT(co2fac)) THEN
        j = coefs%coef%fmv_gas_pos(gas_id_co2)
        coefs%coef%lim_prfl_gmin(:,j) = coefs%coef%env_prfl_gmin(:,j) * (1._jprb - co2fac)
        coefs%coef%lim_prfl_gmax(:,j) = coefs%coef%env_prfl_gmax(:,j) * (1._jprb + co2fac)
      ENDIF
    CASE (gas_id_n2o)
      IF (PRESENT(n2ofac)) THEN
        j = coefs%coef%fmv_gas_pos(gas_id_n2o)
        coefs%coef%lim_prfl_gmin(:,j) = coefs%coef%env_prfl_gmin(:,j) * (1._jprb - n2ofac)
        coefs%coef%lim_prfl_gmax(:,j) = coefs%coef%env_prfl_gmax(:,j) * (1._jprb + n2ofac)
      ENDIF
    CASE (gas_id_co)
      IF (PRESENT(cofac)) THEN
        j = coefs%coef%fmv_gas_pos(gas_id_co)
        coefs%coef%lim_prfl_gmin(:,j) = coefs%coef%env_prfl_gmin(:,j) * (1._jprb - cofac)
        coefs%coef%lim_prfl_gmax(:,j) = coefs%coef%env_prfl_gmax(:,j) * (1._jprb + cofac)
      ENDIF
    CASE (gas_id_ch4)
      IF (PRESENT(ch4fac)) THEN
        j = coefs%coef%fmv_gas_pos(gas_id_ch4)
        coefs%coef%lim_prfl_gmin(:,j) = coefs%coef%env_prfl_gmin(:,j) * (1._jprb - ch4fac)
        coefs%coef%lim_prfl_gmax(:,j) = coefs%coef%env_prfl_gmax(:,j) * (1._jprb + ch4fac)
      ENDIF
    CASE (gas_id_so2)
      IF (PRESENT(so2fac)) THEN
        j = coefs%coef%fmv_gas_pos(gas_id_so2)
        coefs%coef%lim_prfl_gmin(:,j) = coefs%coef%env_prfl_gmin(:,j) * (1._jprb - so2fac)
        coefs%coef%lim_prfl_gmax(:,j) = coefs%coef%env_prfl_gmax(:,j) * (1._jprb + so2fac)
      ENDIF
    END SELECT
  ENDDO

  CATCH
END SUBROUTINE rttov_set_coef_limits
