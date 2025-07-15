! Description:
!> @file
!!   AD of Fresnel reflectance calculation for flat water surface.
!
!> @brief
!!   AD of Fresnel reflectance calculation for flat water surface.
!!
!! @param[in]     chanprof       specifies channels and profiles to simulate
!! @param[in]     calcrefl       flags for internal RTTOV surface BRDF calculation
!! @param[in]     profiles       input atmospheric profiles and surface variables
!! @param[in]     solar          flag to indicate channels with solar contribution
!! @param[in]     coef           optical depth coefficient structure
!! @param[in]     sunglint       internal structure for sea surface BRDF model variables
!! @param[in,out] sunglint_ad    gradient wrt sea surface BRDF model variables
!! @param[in]     fresnrefl_ad   gradient wrt Fresnel reflectance
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
SUBROUTINE rttov_fresnel_ad( &
              chanprof,     &
              calcrefl,     &
              profiles,     &
              solar,        &
              coef,         &
              sunglint,     &
              sunglint_ad,  &
              fresnrefl_ad)

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

  TYPE(rttov_chanprof),  INTENT(IN)    :: chanprof(:)
  LOGICAL(KIND=jplm),    INTENT(IN)    :: calcrefl(SIZE(chanprof))
  TYPE(rttov_profile),   INTENT(IN)    :: profiles(:)
  LOGICAL(KIND=jplm),    INTENT(IN)    :: solar(SIZE(chanprof))
  TYPE(rttov_coef),      INTENT(IN)    :: coef
  TYPE(rttov_sunglint),  INTENT(IN)    :: sunglint
  TYPE(rttov_sunglint),  INTENT(INOUT) :: sunglint_ad
  REAL(KIND=jprb),       INTENT(IN)    :: fresnrefl_ad(SIZE(chanprof))
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
  REAL   (KIND=jprb) :: sincsi_ad
  REAL   (KIND=jprb) :: coscsi_ad
  COMPLEX(KIND=jprb) :: crefpar_ad
  COMPLEX(KIND=jprb) :: crefparj_ad
  COMPLEX(KIND=jprb) :: crefperp_ad
  COMPLEX(KIND=jprb) :: crefperpj_ad
  COMPLEX(KIND=jprb) :: sincsi1_ad
  COMPLEX(KIND=jprb) :: coscsi1_ad
  INTEGER(KIND=jpim) :: nchannels
  COMPLEX(KIND=jprb) :: waopc(SIZE(chanprof))
  REAL   (KIND=jprb) :: ZHOOK_HANDLE
!-----End of header-------------------------------------------------------------

  IF (LHOOK) CALL DR_HOOK('RTTOV_FRESNEL_AD', 0_jpim, ZHOOK_HANDLE)
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

        coscsi                    = COS(sunglint%s(prof)%omega)
        sincsi                    = SIN(sunglint%s(prof)%omega)
        sincsi1                   = sincsi / waopc(i)
        coscsi1                   = SQRT(1._jprb - sincsi1 ** 2_jpim)
        crefpar                   =  - (waopc(i) * coscsi - coscsi1) / (waopc(i) * coscsi + coscsi1)
        crefparj                  = CONJG(crefpar)
        crefperp                  = (coscsi - waopc(i) * coscsi1) / (coscsi + waopc(i) * coscsi1)
        crefperpj                 = CONJG(crefperp)

        crefpar_ad                = fresnrefl_ad(i) * 0.5_jprb * crefparj
        crefparj_ad               = fresnrefl_ad(i) * 0.5_jprb * crefpar
        crefperp_ad               = fresnrefl_ad(i) * 0.5_jprb * crefperpj
        crefperpj_ad              = fresnrefl_ad(i) * 0.5_jprb * crefperp
        crefperp_ad               = crefperp_ad + CONJG(crefperpj_ad)
        crefperpj_ad              = 0
        coscsi_ad                 = REAL(crefperp_ad / (coscsi + waopc(i) * coscsi1), jprb)
        coscsi1_ad                =  - crefperp_ad * waopc(i) / (coscsi + waopc(i) * coscsi1)
        coscsi_ad                 = REAL(coscsi_ad - crefperp_ad * (coscsi - waopc(i) * coscsi1) / &
                                         (coscsi + waopc(i) * coscsi1) ** 2_jpim, jprb)
        coscsi1_ad                = coscsi1_ad - crefperp_ad * waopc(i) * (coscsi - waopc(i) * coscsi1) / &
                                    (coscsi + waopc(i) * coscsi1) ** 2_jpim
        crefperp_ad               = 0
        crefpar_ad                = crefpar_ad + CONJG(crefparj_ad)
        crefparj_ad               = 0
        coscsi_ad                 = REAL(coscsi_ad - crefpar_ad * waopc(i) / (waopc(i) * coscsi + coscsi1), jprb)
        coscsi1_ad                = coscsi1_ad + crefpar_ad / (waopc(i) * coscsi + coscsi1)
        coscsi_ad                 = REAL(coscsi_ad + crefpar_ad * waopc(i) * ( - coscsi1 + waopc(i) * coscsi) / &
                                         (waopc(i) * coscsi + coscsi1) ** 2_jpim, jprb)
        coscsi1_ad                = coscsi1_ad + crefpar_ad * ( - coscsi1 + waopc(i) * coscsi) / &
                                    (waopc(i) * coscsi + coscsi1) ** 2_jpim
        crefpar_ad                = 0
        sincsi1_ad                =  - coscsi1_ad * sincsi1 / coscsi1
        coscsi1_ad                = 0
        sincsi_ad                 = REAL(sincsi1_ad / waopc(i), jprb)
        sincsi1_ad                = 0
        sunglint_ad%s(prof)%omega = sunglint_ad%s(prof)%omega + sincsi_ad * coscsi
        sunglint_ad%s(prof)%omega = sunglint_ad%s(prof)%omega - coscsi_ad * sincsi
        sincsi_ad                 = 0
        coscsi_ad                 = 0
      ENDIF
    ENDIF
  ENDDO
  IF (LHOOK) CALL DR_HOOK('RTTOV_FRESNEL_AD', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_fresnel_ad
