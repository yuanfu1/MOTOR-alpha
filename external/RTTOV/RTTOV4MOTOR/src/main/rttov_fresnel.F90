! Description:
!> @file
!!   Computes Fresnel reflectance for flat water surface.
!
!> @brief
!!   Computes Fresnel reflectance for flat water surface.
!!
!! @param[in]     chanprof       specifies channels and profiles to simulate
!! @param[in]     calcrefl       flags for internal RTTOV surface BRDF calculation
!! @param[in]     profiles       input atmospheric profiles and surface variables
!! @param[in]     solar          flag to indicate channels with solar contribution
!! @param[in]     coef           optical depth coefficient structure
!! @param[in]     sunglint       internal structure for sea surface BRDF model variables
!! @param[out]    fresnrefl      calculated Fresnel reflectance
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
SUBROUTINE rttov_fresnel( &
              chanprof,  &
              calcrefl,  &
              profiles,  &
              solar,     &
              coef,      &
              sunglint,  &
              fresnrefl)

  USE rttov_types, ONLY :  &
         rttov_chanprof, &
         rttov_profile,  &
         rttov_coef,     &
         rttov_sunglint
  USE parkind1, ONLY : jprb, jplm
!INTF_OFF
  USE rttov_const, ONLY :       &
         surftype_sea,          &
         watertype_fresh_water, &
         watertype_ocean_water
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jpim
!INTF_ON
  IMPLICIT NONE

  TYPE(rttov_chanprof),  INTENT(IN)  :: chanprof(:)
  LOGICAL(KIND=jplm),    INTENT(IN)  :: calcrefl(SIZE(chanprof))
  TYPE(rttov_profile),   INTENT(IN)  :: profiles(:)
  LOGICAL(KIND=jplm),    INTENT(IN)  :: solar(SIZE(chanprof))
  TYPE(rttov_coef),      INTENT(IN)  :: coef
  TYPE(rttov_sunglint),  INTENT(IN)  :: sunglint
  REAL(KIND=jprb),       INTENT(OUT) :: fresnrefl(SIZE(chanprof))
!INTF_END

  INTEGER(KIND=jpim) :: i, chan, prof
  REAL   (KIND=jprb) :: sincsi
  REAL   (KIND=jprb) :: coscsi
  COMPLEX(KIND=jprb) :: crefpar
  COMPLEX(KIND=jprb) :: crefparj
  COMPLEX(KIND=jprb) :: crefperp
  COMPLEX(KIND=jprb) :: crefperpj
  COMPLEX(KIND=jprb) :: sincsi1
  COMPLEX(KIND=jprb) :: coscsi1
  INTEGER(KIND=jpim) :: nchannels
  COMPLEX(KIND=jprb) :: waopc(SIZE(chanprof))
  REAL   (KIND=jprb) :: ZHOOK_HANDLE
!-----End of header-------------------------------------------------------------

  IF (LHOOK) CALL DR_HOOK('RTTOV_FRESNEL', 0_jpim, ZHOOK_HANDLE)
  nchannels = SIZE(chanprof)
  DO i = 1, nchannels
    chan = chanprof(i)%chan
    prof = chanprof(i)%prof
    IF (solar(i) .AND. calcrefl(i)) THEN
      IF (profiles(prof)%skin%surftype == surftype_sea) THEN
        IF (profiles(prof)%skin%watertype == watertype_fresh_water) THEN
          waopc(i) = coef%woc_waopc_fw(chan)
        ELSE IF (profiles(prof)%skin%watertype == watertype_ocean_water) THEN
          waopc(i) = coef%woc_waopc_ow(chan)
        ENDIF
        coscsi       = COS(sunglint%s(prof)%omega)
        sincsi       = SIN(sunglint%s(prof)%omega)
        sincsi1      = sincsi / waopc(i)
        coscsi1      = SQRT(1._jprb - sincsi1 ** 2_jpim)
        crefpar      =  - (waopc(i) * coscsi - coscsi1) / (waopc(i) * coscsi + coscsi1)
        crefparj     = CONJG(crefpar)
        crefperp     = (coscsi - waopc(i) * coscsi1) / (coscsi + waopc(i) * coscsi1)
        crefperpj    = CONJG(crefperp)
        fresnrefl(i) = REAL (0.5_jprb * (crefpar * crefparj + crefperp * crefperpj), jprb)
      ENDIF
    ENDIF
  ENDDO
  IF (LHOOK) CALL DR_HOOK('RTTOV_FRESNEL', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_fresnel
