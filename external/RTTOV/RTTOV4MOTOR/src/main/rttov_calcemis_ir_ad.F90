! Description:
!> @file
!!   AD of IR surface emissivity calculation.
!
!> @brief
!!   AD of IR surface emissivity calculation.
!!
!! @details
!!   For land and sea-ice surface types and for sea
!!   surfaces when ISEM is used the AD is not modified as
!!   there is no dependence on any input profile variables.
!!
!!   For sea surfaces where the physically-based emissivity
!!   model is used the wind speed and, where applicable, Tskin
!!   AD values are updated.
!!
!! @param[out]    err            status on exit
!! @param[in]     opts           options to configure the simulations
!! @param[in]     chanprof       specifies channels and profiles to simulate
!! @param[in]     profiles       input atmospheric profiles and surface variables
!! @param[in,out] profiles_ad    gradient wrt atmospheric profile and surface variables
!! @param[in]     coefs          coefficients structure for instrument to simulate
!! @param[in]     thermal        flag to indicate channels with thermal emission
!! @param[in]     calcemis       flags for internal RTTOV surface emissivity calculation
!! @param[in,out] emissivity_ad  gradient wrt surface emissivities
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
SUBROUTINE rttov_calcemis_ir_ad( &
              err,         &
              opts,        &
              chanprof,    &
              profiles,    &
              profiles_ad, &
              coefs,       &
              thermal,     &
              calcemis,    &
              emissivity_ad)
!INTF_OFF
#include "throw.h"
!INTF_ON
  USE rttov_types, ONLY :  &
       rttov_options,      &
       rttov_chanprof,     &
       rttov_coefs,        &
       rttov_profile
  USE parkind1, ONLY : jprb, jpim, jplm
!INTF_OFF
  USE rttov_const, ONLY :   &
       realtol,             &
       surftype_sea,        &
       surftype_land,       &
       surftype_seaice,     &
       errorstatus_success, &
       errorstatus_fatal
  USE yomhook, ONLY : LHOOK, DR_HOOK
!INTF_ON
  IMPLICIT NONE
  INTEGER(KIND=jpim),   INTENT(OUT)   :: err
  TYPE(rttov_options),  INTENT(IN)    :: opts
  TYPE(rttov_chanprof), INTENT(IN)    :: chanprof(:)
  TYPE(rttov_profile),  INTENT(IN)    :: profiles(:)
  TYPE(rttov_profile),  INTENT(INOUT) :: profiles_ad(SIZE(profiles))
  TYPE(rttov_coefs),    INTENT(IN)    :: coefs
  LOGICAL(KIND=jplm),   INTENT(IN)    :: thermal(SIZE(chanprof))
  LOGICAL(KIND=jplm),   INTENT(IN)    :: calcemis(SIZE(chanprof))
  REAL(KIND=jprb),      INTENT(INOUT) :: emissivity_ad(SIZE(chanprof))
!INTF_END

#include "rttov_errorreport.interface"

  INTEGER(KIND=jpim) :: j, chan, iprof
  REAL(KIND=jprb)    :: windsp, aems, bems, cems
  REAL(KIND=jprb)    :: windsp_ad, aems_ad, bems_ad, cems_ad
  REAL(KIND=jprb)    :: expf, fac
  INTEGER(KIND=jpim) :: nchannels
  REAL(KIND=jprb)    :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------

  TRY
  IF (LHOOK) CALL DR_HOOK('RTTOV_CALCEMIS_IR_AD', 0_jpim, ZHOOK_HANDLE)
  nchannels = SIZE(chanprof)
  aems_ad = 0._jprb
  bems_ad = 0._jprb
  cems_ad = 0._jprb
  windsp_ad = 0._jprb
  DO j = 1, nchannels

    IF (.NOT. thermal(j)) CYCLE
    IF (.NOT. calcemis(j)) CYCLE

    chan  = chanprof(j)%chan
    iprof = chanprof(j)%prof

    IF (opts%rt_ir%pc%addpc .AND. profiles(iprof)%skin%surftype == surftype_sea) THEN

      windsp  = SQRT(profiles(iprof)%s2m%u ** 2_jpim + profiles(iprof)%s2m%v ** 2_jpim)
      aems    = coefs%coef_pccomp%emiss_c1(chan) + coefs%coef_pccomp%emiss_c2(chan) * windsp + &
                coefs%coef_pccomp%emiss_c3(chan) * windsp ** 2_jpim
      bems    = coefs%coef_pccomp%emiss_c4(chan) + coefs%coef_pccomp%emiss_c5(chan) * windsp + &
                coefs%coef_pccomp%emiss_c6(chan) * windsp ** 2_jpim
      cems    = coefs%coef_pccomp%emiss_c7(chan) + coefs%coef_pccomp%emiss_c8(chan) * windsp
      expf    = EXP( ((coefs%coef_pccomp%emiss_c9(chan) - 60._jprb) ** 2_jpim - &
                      (profiles(iprof)%zenangle - coefs%coef_pccomp%emiss_c9(chan)) ** 2_jpim) / cems )
      fac     = - ((coefs%coef_pccomp%emiss_c9(chan) - 60._jprb) ** 2_jpim -   &
                   (profiles(iprof)%zenangle - coefs%coef_pccomp%emiss_c9(chan)) ** 2_jpim) / cems ** 2_jpim
      aems_ad = aems_ad + emissivity_ad(j)
      bems_ad = bems_ad + emissivity_ad(j) * expf
      cems_ad = cems_ad + bems * expf * emissivity_ad(j) * fac
      aems_ad = aems_ad - emissivity_ad(j) * expf
      cems_ad = cems_ad - aems * expf * emissivity_ad(j) * fac
      emissivity_ad(j) = 0._jprb
      windsp_ad  = windsp_ad + coefs%coef_pccomp%emiss_c8(chan) * cems_ad
      windsp_ad  = windsp_ad + (coefs%coef_pccomp%emiss_c5(chan) + &
                                2._jprb * coefs%coef_pccomp%emiss_c6(chan) * windsp) * bems_ad
      windsp_ad  = windsp_ad + (coefs%coef_pccomp%emiss_c2(chan) + &
                                2._jprb * coefs%coef_pccomp%emiss_c3(chan) * windsp) * aems_ad
      aems_ad = 0._jprb
      bems_ad = 0._jprb
      cems_ad = 0._jprb
      IF (windsp > realtol) THEN
        profiles_ad(iprof)%s2m%u = profiles_ad(iprof)%s2m%u + profiles(iprof)%s2m%u * windsp_ad / windsp
        profiles_ad(iprof)%s2m%v = profiles_ad(iprof)%s2m%v + profiles(iprof)%s2m%v * windsp_ad / windsp
      ELSE
        profiles_ad(iprof)%s2m%u = profiles_ad(iprof)%s2m%u + windsp_ad / SQRT(2._jprb)
        profiles_ad(iprof)%s2m%v = profiles_ad(iprof)%s2m%v + windsp_ad / SQRT(2._jprb)
      ENDIF
      windsp_ad = 0._jprb

    ELSE ! Not using PC-RTTOV sea surface emissivity model

      IF (profiles(iprof)%skin%surftype == surftype_land .OR. &
          profiles(iprof)%skin%surftype == surftype_seaice) THEN

        ! Nothing to do

      ELSE

        IF (opts%rt_ir%ir_sea_emis_model == 2) THEN

          windsp  = SQRT(profiles(iprof)%s2m%u ** 2_jpim + profiles(iprof)%s2m%v ** 2_jpim)
          aems    = coefs%coef%iremis_coef(1,chan) + coefs%coef%iremis_coef(2,chan) * windsp + &
                    coefs%coef%iremis_coef(3,chan) * windsp ** 2_jpim
          bems    = coefs%coef%iremis_coef(4,chan) + coefs%coef%iremis_coef(5,chan) * windsp + &
                    coefs%coef%iremis_coef(6,chan) * windsp ** 2_jpim
          cems    = coefs%coef%iremis_coef(7,chan) + coefs%coef%iremis_coef(8,chan) * windsp
          expf    = EXP( ((coefs%coef%iremis_coef(9,chan) - coefs%coef%iremis_angle0) ** 2_jpim - &
                          (profiles(iprof)%zenangle - coefs%coef%iremis_coef(9,chan)) ** 2_jpim) / cems )
          fac     = - ((coefs%coef%iremis_coef(9,chan) - coefs%coef%iremis_angle0) ** 2_jpim - &
                       (profiles(iprof)%zenangle - coefs%coef%iremis_coef(9,chan)) ** 2_jpim) / cems ** 2_jpim

          IF (coefs%coef%iremis_ncoef > 9) THEN
            IF (ABS(coefs%coef%iremis_coef(10,chan)) > 0._jprb .AND. &
                ABS(coefs%coef%iremis_coef(11,chan)) > 0._jprb) THEN
              profiles_ad(iprof)%skin%t = profiles_ad(iprof)%skin%t + &
                  emissivity_ad(j) * &
                  (coefs%coef%iremis_coef(10,chan) + &
                    EXP(coefs%coef%iremis_coef(11,chan) * profiles(iprof)%zenangle**2 / coefs%coef%iremis_angle0**2))
            ENDIF
          ENDIF

          aems_ad = aems_ad + emissivity_ad(j)
          bems_ad = bems_ad + emissivity_ad(j) * expf
          cems_ad = cems_ad + bems * expf * emissivity_ad(j) * fac
          aems_ad = aems_ad - emissivity_ad(j) * expf
          cems_ad = cems_ad - aems * expf * emissivity_ad(j) * fac
          emissivity_ad(j) = 0._jprb
          windsp_ad  = windsp_ad + coefs%coef%iremis_coef(8,chan) * cems_ad
          windsp_ad  = windsp_ad + (coefs%coef%iremis_coef(5,chan) + &
                                    2._jprb * coefs%coef%iremis_coef(6,chan) * windsp) * bems_ad
          windsp_ad  = windsp_ad + (coefs%coef%iremis_coef(2,chan) + &
                                    2._jprb * coefs%coef%iremis_coef(3,chan) * windsp) * aems_ad
          aems_ad = 0._jprb
          bems_ad = 0._jprb
          cems_ad = 0._jprb
          IF (windsp > realtol) THEN
            profiles_ad(iprof)%s2m%u = profiles_ad(iprof)%s2m%u + profiles(iprof)%s2m%u * windsp_ad / windsp
            profiles_ad(iprof)%s2m%v = profiles_ad(iprof)%s2m%v + profiles(iprof)%s2m%v * windsp_ad / windsp
          ELSE
            profiles_ad(iprof)%s2m%u = profiles_ad(iprof)%s2m%u + windsp_ad / SQRT(2._jprb)
            profiles_ad(iprof)%s2m%v = profiles_ad(iprof)%s2m%v + windsp_ad / SQRT(2._jprb)
          ENDIF
          windsp_ad = 0._jprb

        ELSE IF (opts%rt_ir%ir_sea_emis_model == 1) THEN

          ! ISEM6
          ! No dependence on active profile variables

        ELSE
          err = errorstatus_fatal
          THROWM(err.NE.0,'Unknown IR sea surface emissivity model')
        ENDIF

      ENDIF

    ENDIF
  ENDDO
  IF (LHOOK) CALL DR_HOOK('RTTOV_CALCEMIS_IR_AD', 1_jpim, ZHOOK_HANDLE)
  CATCH
  IF (LHOOK) CALL DR_HOOK('RTTOV_CALCEMIS_IR_AD', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_calcemis_ir_ad
