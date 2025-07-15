! Description:
!> @file
!!   Internal module defining helper subroutines for handling coefficients.
!
!> @brief
!!   Internal module defining helper subroutines for handling coefficients.
!!
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
MODULE rttov_fast_coef_utils_mod

  USE rttov_types, ONLY : &
    rttov_coef, & 
    rttov_fast_coef

  USE parkind1, ONLY : jpim, jprb, jplm

  USE rttov_const, ONLY :  &
    gas_id_mixed,          &
    gas_id_watervapour,    &
    gas_id_ozone,          &
    gas_id_wvcont,         &
    gas_id_co2,            &
    gas_id_n2o,            &
    gas_id_co,             &
    gas_id_ch4,            &
    gas_id_so2,            &
    upper_layer_bound_not_set

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: set_pointers, nullify_pointers, set_fastcoef_level_bounds

CONTAINS

  !> Set the fast_coef%[gas] pointers to the fast_coef%gas_array(:) elements
  !! @param[inout]    fast_coef     thermal or solar coefficients (within rttov_coef structure)
  !! @param[in]       gas_pos       index of gas in gasarray
  !! @param[in]       gas_id        ID of gas (from rttov_const)
  SUBROUTINE set_pointers(fast_coef, gas_pos, gas_id)
    TYPE(rttov_fast_coef), INTENT(INOUT), TARGET :: fast_coef
    INTEGER(KIND=jpim), INTENT(IN) :: gas_pos
    INTEGER(KIND=jpim), INTENT(IN) :: gas_id

    SELECT CASE(gas_id)
    CASE(gas_id_mixed)
      fast_coef%mixedgas => fast_coef%gasarray(gas_pos)%coef
    CASE(gas_id_watervapour)
      fast_coef%watervapour => fast_coef%gasarray(gas_pos)%coef
    CASE(gas_id_ozone)
      fast_coef%ozone => fast_coef%gasarray(gas_pos)%coef
    CASE(gas_id_wvcont)
      fast_coef%wvcont => fast_coef%gasarray(gas_pos)%coef
    CASE(gas_id_co2)
      fast_coef%co2 => fast_coef%gasarray(gas_pos)%coef
    CASE(gas_id_n2o)
      fast_coef%n2o => fast_coef%gasarray(gas_pos)%coef
    CASE(gas_id_co)
      fast_coef%co => fast_coef%gasarray(gas_pos)%coef
    CASE(gas_id_ch4)
      fast_coef%ch4 => fast_coef%gasarray(gas_pos)%coef
    CASE(gas_id_so2)
      fast_coef%so2 => fast_coef%gasarray(gas_pos)%coef
    END SELECT
  END SUBROUTINE set_pointers

  !> Nullify the fast_coef%[gas] pointers
  !! @param[inout]    fast_coef     thermal or solar coefficients (within rttov_coef structure)
  !! @param[in]       gas_id        ID of gas (from rttov_const)
  SUBROUTINE nullify_pointers(fast_coef, gas_id)
    TYPE(rttov_fast_coef), INTENT(INOUT), TARGET :: fast_coef
    INTEGER(KIND=jpim), INTENT(IN) :: gas_id

    SELECT CASE(gas_id)
    CASE(gas_id_mixed)
      NULLIFY(fast_coef%mixedgas)
    CASE(gas_id_watervapour)
      NULLIFY(fast_coef%watervapour)
    CASE(gas_id_ozone)
      NULLIFY(fast_coef%ozone)
    CASE(gas_id_wvcont)
      NULLIFY(fast_coef%wvcont)
    CASE(gas_id_co2)
      NULLIFY(fast_coef%co2)
    CASE(gas_id_n2o)
      NULLIFY(fast_coef%n2o)
    CASE(gas_id_co)
      NULLIFY(fast_coef%co)
    CASE(gas_id_ch4)
      NULLIFY(fast_coef%ch4)
    CASE(gas_id_so2)
      NULLIFY(fast_coef%so2)
    END SELECT
  END SUBROUTINE nullify_pointers

  !> Fine upper/lower layer bounds for which optical depth calculation is
  !! required. This allows the option of faster simulations by omitting
  !! calculations on layers to which the sensor is insensitive on a channel-
  !! by-channel basis.
  !! @param[inout]    coef          optical depth coefficient structure
  !! @param[in]       fastcoef      thermal or solar coefficients (within rttov_coef structure)
  !! @param[in]       thermal       true for thermal coefs, false for solar
  SUBROUTINE set_fastcoef_level_bounds(coef, fastcoef, thermal)

    TYPE(rttov_coef),      INTENT(INOUT) :: coef
    TYPE(rttov_fast_coef), INTENT(IN)    :: fastcoef(:)
    LOGICAL(KIND=jplm),    INTENT(IN)    :: thermal

    INTEGER(jpim) :: i, j, k, ind, gas_id, gas_pos

    ind = 1
    IF (.NOT. thermal) ind = 2

    DO k = 1, coef%fmv_gas ! ngases
      gas_id = coef%fmv_gas_id(k)
      gas_pos = coef%fmv_gas_pos(gas_id)

      DO i = 1, coef%fmv_chn
        IF (ASSOCIATED(fastcoef(i)%gasarray(gas_pos)%coef)) THEN
          DO j = 1, coef%fmv_lvl(1)-1
            IF (ANY(fastcoef(i)%gasarray(gas_pos)%coef(:,j) /= 0._jprb)) THEN
              coef%bounds(1,gas_pos,i,ind) = j
              EXIT
            ENDIF
          ENDDO
          DO j = coef%fmv_lvl(1)-1, 1, -1
            IF (ANY(fastcoef(i)%gasarray(gas_pos)%coef(:,j) /= 0._jprb)) THEN
              coef%bounds(2,gas_pos,i,ind) = j
              EXIT
            ENDIF
          ENDDO
        ELSE ! coefs not allocated at all
          coef%bounds(1,gas_pos,i,ind) = upper_layer_bound_not_set
          coef%bounds(2,gas_pos,i,ind) = 0
        ENDIF
      ENDDO
    ENDDO
  END SUBROUTINE set_fastcoef_level_bounds

END MODULE rttov_fast_coef_utils_mod
