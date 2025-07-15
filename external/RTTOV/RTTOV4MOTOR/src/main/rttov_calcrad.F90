! Description:
!> @file
!!   Calculate Planck radiances from atmospheric temperatures
!
!> @brief
!!   Calculate Planck radiances from atmospheric temperatures
!!
!! @details
!!   The Planck function is calculated using the channel central wavenumbers
!!   as specified in the coefficient file. To adjust for the finite spectral
!!   bandwidth the temperatures are modified using band correction coefficients
!!   also stored in the coefficient file.
!!
!!   Radiances units are mW/cm-1/ster/m2 and temperature units are Kelvin.
!!   Radiances are not calculated for channels at wavelengths below 3um  (i.e.
!!   with insignificant thermal emission).
!!
!!   Reference:
!!   Sharp, J.C, 1983: A comparison of approximate methods for converting
!!   between radiance and equivalent black body temperature for a radiometer
!!   channel, Met Office technical report (see docs/ directory)
!!
!! @param[in]     addcosmic        flag for calculation of cosmic MW background radiance
!! @param[in]     chanprof         specifies channels and profiles to simulate
!! @param[in]     profiles         input atmospheric profiles and surface variables
!! @param[in]     coef             optical depth coefficient structure
!! @param[in]     thermal          flag to indicate channels with thermal emission
!! @param[in,out] auxrad           auxiliary radiance structure
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
SUBROUTINE rttov_calcrad( &
              addcosmic, &
              chanprof,  &
              profiles,  &
              coef,      &
              thermal,   &
              auxrad)

  USE rttov_types, ONLY : rttov_chanprof, rttov_coef, rttov_profile, rttov_radiance_aux
  USE parkind1, ONLY : jplm
!INTF_OFF
  USE parkind1, ONLY : jprb, jpim
  USE rttov_const, ONLY : tcosmic
  USE rttov_math_mod, ONLY : planck
!INTF_ON
  IMPLICIT NONE

  LOGICAL(KIND=jplm),       INTENT(IN)    :: addcosmic
  TYPE(rttov_chanprof),     INTENT(IN)    :: chanprof(:)
  TYPE(rttov_profile),      INTENT(IN)    :: profiles(:)
  TYPE(rttov_coef),         INTENT(IN)    :: coef
  LOGICAL(KIND=jplm),       INTENT(IN)    :: thermal(SIZE(chanprof))
  TYPE(rttov_radiance_aux), INTENT(INOUT) :: auxrad
!INTF_END

  INTEGER(KIND=jpim) :: chan, prof, ichan
  INTEGER(KIND=jpim) :: nchanprof
!- End of header ------------------------------------------------------

  nchanprof = SIZE(chanprof)

  DO ichan = 1, nchanprof
    IF (.NOT. thermal(ichan)) CYCLE
    chan = chanprof(ichan)%chan
    prof = chanprof(ichan)%prof

! Calculate effective temperatures, store in auxrad
    IF (coef%ff_val_bc) THEN
      auxrad%skin_t_eff(ichan)  = coef%ff_bco(chan) + coef%ff_bcs(chan) * &
                                  profiles(prof)%skin%t
      auxrad%surf_t_eff(ichan)  = coef%ff_bco(chan) + coef%ff_bcs(chan) * &
                                  profiles(prof)%s2m%t
      auxrad%air_t_eff(:,ichan) = coef%ff_bco(chan) + coef%ff_bcs(chan) * &
                                  profiles(prof)%t(:)
    ELSE
      auxrad%skin_t_eff(ichan)  = profiles(prof)%skin%t
      auxrad%surf_t_eff(ichan)  = profiles(prof)%s2m%t
      auxrad%air_t_eff(:,ichan) = profiles(prof)%t(:)
    ENDIF

    IF (addcosmic) THEN
      ! addcosmic is true for MW only
      IF (coef%ff_val_bc) THEN
        auxrad%cosmic_t_eff(ichan) = coef%ff_bco(chan) + coef%ff_bcs(chan) * tcosmic
      ELSE
        auxrad%cosmic_t_eff(ichan) = tcosmic
      ENDIF
      CALL planck(coef%planck1(chan), coef%planck2(chan), &
                  auxrad%cosmic_t_eff(ichan), auxrad%cosmic(ichan))
    ELSE
      auxrad%cosmic(ichan) = 0.0_jprb
    ENDIF

    CALL planck(coef%planck1(chan), coef%planck2(chan), &
                auxrad%skin_t_eff(ichan), auxrad%skin(ichan))

    CALL planck(coef%planck1(chan), coef%planck2(chan), &
                auxrad%surf_t_eff(ichan), auxrad%surfair(ichan))

    CALL planck(coef%planck1(chan), coef%planck2(chan), &
                auxrad%air_t_eff(:,ichan), auxrad%air(:,ichan))
  ENDDO
END SUBROUTINE rttov_calcrad
