! Description:
!> @file
!!   Computes IR surface emissivities.
!
!> @brief
!!   Computes IR surface emissivities.
!!
!! @details
!!   This subroutine provides emissivity values for each
!!   channel/profile where calcemis is set to TRUE.
!!
!!   For land and sea-ice surfaces constant values are
!!   assigned.
!!
!!   For sea surfaces two models are available, selected
!!   via the options structure:
!!   ISEM is the older model with only a zenith angle dependence.
!!   The newer model also includes a wind speed dependence and
!!   a Tskin depedence in some parts of the spectrum.
!!
!!   PC-RTTOV currently uses it's own sea surface emissivity
!!   model, but this will be changed in the future.
!!
!!   RTTOV-6 IR surface emissivity report, V. Sherlock at:
!!   http://www.metoffice.com/research/interproj/nwpsaf/rtm/papers/isem6.pdf
!!
!! @param[out]    err            status on exit
!! @param[in]     opts           options to configure the simulations
!! @param[in]     chanprof       specifies channels and profiles to simulate
!! @param[in]     profiles       input atmospheric profiles and surface variables
!! @param[in]     geometry       internal geometry structure
!! @param[in]     coefs          coefficients structure for instrument to simulate
!! @param[in]     thermal        flag to indicate channels with thermal emission
!! @param[in]     calcemis       flags for internal RTTOV surface emissivity calculation
!! @param[in,out] emissivity     updated with surface emissivities on exit
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
SUBROUTINE rttov_calcemis_ir( &
              err,         &
              opts,        &
              chanprof,    &
              profiles,    &
              geometry,    &
              coefs,       &
              thermal,     &
              calcemis,    &
              emissivity)
!INTF_OFF
#include "throw.h"
!INTF_ON
  USE rttov_types, ONLY :  &
       rttov_options,      &
       rttov_chanprof,     &
       rttov_coefs,        &
       rttov_profile,      &
       rttov_geometry
  USE parkind1, ONLY : jprb, jpim, jplm
!INTF_OFF
  USE rttov_const, ONLY :   &
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
  TYPE(rttov_geometry), INTENT(IN)    :: geometry(SIZE(profiles))
  TYPE(rttov_coefs),    INTENT(IN)    :: coefs
  LOGICAL(KIND=jplm),   INTENT(IN)    :: thermal(SIZE(chanprof))
  LOGICAL(KIND=jplm),   INTENT(IN)    :: calcemis(SIZE(chanprof))
  REAL(KIND=jprb),      INTENT(INOUT) :: emissivity(SIZE(chanprof))
!INTF_END

#include "rttov_errorreport.interface"

  INTEGER(KIND=jpim) :: j, chan, iprof
  REAL(KIND=jprb)    :: windsp, aems, bems, cems
  INTEGER(KIND=jpim) :: nchannels
  REAL(KIND=jprb)    :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------

  TRY
  IF (LHOOK) CALL DR_HOOK('RTTOV_CALCEMIS_IR', 0_jpim, ZHOOK_HANDLE)
  nchannels = SIZE(chanprof)
  DO j = 1, nchannels

    IF (.NOT. thermal(j)) CYCLE
    IF (.NOT. calcemis(j)) CYCLE

    chan  = chanprof(j)%chan
    iprof = chanprof(j)%prof

    IF (opts%rt_ir%pc%addpc .AND. profiles(iprof)%skin%surftype == surftype_sea) THEN

      windsp  = SQRT(profiles(iprof)%s2m%u ** 2 + profiles(iprof)%s2m%v ** 2)
      aems    = coefs%coef_pccomp%emiss_c1(chan) + &
                coefs%coef_pccomp%emiss_c2(chan) * windsp + &
                coefs%coef_pccomp%emiss_c3(chan) * windsp ** 2_jpim
      bems    = coefs%coef_pccomp%emiss_c4(chan) + &
                coefs%coef_pccomp%emiss_c5(chan) * windsp + &
                coefs%coef_pccomp%emiss_c6(chan) * windsp ** 2_jpim
      cems    = coefs%coef_pccomp%emiss_c7(chan) + coefs%coef_pccomp%emiss_c8(chan) * windsp
      emissivity(j) = aems + (bems - aems) * &
                EXP(((coefs%coef_pccomp%emiss_c9(chan) - 60._jprb) ** 2_jpim - &
                     (profiles(iprof)%zenangle - coefs%coef_pccomp%emiss_c9(chan)) ** 2_jpim) / cems)

    ELSE ! Not using PC-RTTOV sea surface emissivity model

      ! Use a fixed value over land and seaice

      IF (profiles(iprof)%skin%surftype == surftype_land) THEN
        emissivity(j) = 0.98_jprb
      ELSE IF (profiles(iprof)%skin%surftype == surftype_seaice) THEN
        emissivity(j) = 0.99_jprb
      ELSE

        IF (opts%rt_ir%ir_sea_emis_model == 2) THEN

          ! New IR sea surface emissivity model

          windsp  = SQRT(profiles(iprof)%s2m%u ** 2 + profiles(iprof)%s2m%v ** 2)

          aems    = coefs%coef%iremis_coef(1,chan) + &
                    coefs%coef%iremis_coef(2,chan) * windsp + &
                    coefs%coef%iremis_coef(3,chan) * windsp ** 2_jpim
          bems    = coefs%coef%iremis_coef(4,chan) + &
                    coefs%coef%iremis_coef(5,chan) * windsp + &
                    coefs%coef%iremis_coef(6,chan) * windsp ** 2_jpim
          cems    = coefs%coef%iremis_coef(7,chan) + coefs%coef%iremis_coef(8,chan) * windsp
          emissivity(j) = aems + (bems - aems) * &
                    EXP(((coefs%coef%iremis_coef(9,chan) - coefs%coef%iremis_angle0) ** 2_jpim - &
                         (profiles(iprof)%zenangle - coefs%coef%iremis_coef(9,chan)) ** 2_jpim) / cems)

          ! Tskin-dependent parameterisation: not applied in all parts of the spectrum
          IF (coefs%coef%iremis_ncoef > 9) THEN
            IF (ABS(coefs%coef%iremis_coef(10,chan)) > 0._jprb .AND. &
                ABS(coefs%coef%iremis_coef(11,chan)) > 0._jprb) THEN
              emissivity(j) = emissivity(j) + &
                  (profiles(iprof)%skin%t - coefs%coef%iremis_tskin0) * &
                  (coefs%coef%iremis_coef(10,chan) + &
                    EXP(coefs%coef%iremis_coef(11,chan) * profiles(iprof)%zenangle**2 / coefs%coef%iremis_angle0**2))
            ENDIF
          ENDIF

        ELSE IF (opts%rt_ir%ir_sea_emis_model == 1) THEN

          ! ISEM6: emissivity is a polynomial in normalized zenith angle
          emissivity(j) = coefs%coef%ssirem_a0(chan) - &
                          coefs%coef%ssirem_a1(chan) * geometry(iprof)%normzen ** coefs%coef%ssirem_xzn1(chan) -  &
                          coefs%coef%ssirem_a2(chan) * geometry(iprof)%normzen ** coefs%coef%ssirem_xzn2(chan)

        ELSE
          err = errorstatus_fatal
          THROWM(err.NE.0,'Unknown IR sea surface emissivity model')
        ENDIF

      ENDIF

    ENDIF
  ENDDO
  IF (LHOOK) CALL DR_HOOK('RTTOV_CALCEMIS_IR', 1_jpim, ZHOOK_HANDLE)
  CATCH
  IF (LHOOK) CALL DR_HOOK('RTTOV_CALCEMIS_IR', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_calcemis_ir
