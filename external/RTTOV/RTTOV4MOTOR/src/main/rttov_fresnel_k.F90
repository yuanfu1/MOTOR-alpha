! Description:
!> @file
!!   K of Fresnel reflectance calculation for flat water surface.
!
!> @brief
!!   K of Fresnel reflectance calculation for flat water surface.
!!
!! @param[in]     chanprof       specifies channels and profiles to simulate
!! @param[in]     calcrefl       flags for internal RTTOV surface BRDF calculation
!! @param[in]     profiles       input atmospheric profiles and surface variables
!! @param[in]     solar          flag to indicate channels with solar contribution
!! @param[in]     coef           optical depth coefficient structure
!! @param[in]     sunglint       internal structure for sea surface BRDF model variables
!! @param[in,out] sunglint_k     Jacobian wrt sea surface BRDF model variables
!! @param[in]     fresnrefl_k    Jacobian wrt Fresnel reflectance
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
SUBROUTINE rttov_fresnel_k( &
              chanprof,    &
              calcrefl,    &
              profiles,    &
              solar,       &
              coef,        &
              sunglint,    &
              sunglint_k,  &
              fresnrefl_k)


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
  TYPE(rttov_sunglint),  INTENT(INOUT) :: sunglint_k
  REAL(KIND=jprb),       INTENT(IN)    :: fresnrefl_k(SIZE(chanprof))
!INTF_END

  INTEGER(KIND=jpim) :: i, prof, chan
  REAL   (KIND=jprb) :: sincsi
  REAL   (KIND=jprb) :: coscsi
  COMPLEX(KIND=jprb) :: crefpar
  COMPLEX(KIND=jprb) :: crefparj
  COMPLEX(KIND=jprb) :: crefperp
  COMPLEX(KIND=jprb) :: crefperpj
  COMPLEX(KIND=jprb) :: sincsi1
  COMPLEX(KIND=jprb) :: coscsi1
  REAL   (KIND=jprb) :: sincsi_k
  REAL   (KIND=jprb) :: coscsi_k
  COMPLEX(KIND=jprb) :: crefpar_k
  COMPLEX(KIND=jprb) :: crefparj_k
  COMPLEX(KIND=jprb) :: crefperp_k
  COMPLEX(KIND=jprb) :: crefperpj_k
  COMPLEX(KIND=jprb) :: sincsi1_k
  COMPLEX(KIND=jprb) :: coscsi1_k
  INTEGER(KIND=jpim) :: nchannels
  COMPLEX(KIND=jprb) :: waopc(SIZE(chanprof))
  REAL   (KIND=jprb) :: ZHOOK_HANDLE
!-----End of header-------------------------------------------------------------

  IF (LHOOK) CALL DR_HOOK('RTTOV_FRESNEL_K', 0_jpim, ZHOOK_HANDLE)
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
        coscsi                = COS(sunglint%s(prof)%omega)
        sincsi                = SIN(sunglint%s(prof)%omega)
        sincsi1               = sincsi / waopc(i)
        coscsi1               = SQRT(1._jprb - sincsi1 ** 2_jpim)
        crefpar               =  - (waopc(i) * coscsi - coscsi1) / (waopc(i) * coscsi + coscsi1)
        crefparj              = CONJG(crefpar)
        crefperp              = (coscsi - waopc(i) * coscsi1) / (coscsi + waopc(i) * coscsi1)
        crefperpj             = CONJG(crefperp)

        crefpar_k             = fresnrefl_k(i) * 0.5_jprb * crefparj
        crefparj_k            = fresnrefl_k(i) * 0.5_jprb * crefpar
        crefperp_k            = fresnrefl_k(i) * 0.5_jprb * crefperpj
        crefperpj_k           = fresnrefl_k(i) * 0.5_jprb * crefperp
        crefperp_k            = crefperp_k + CONJG(crefperpj_k)
        crefperpj_k           = 0
        coscsi_k              = REAL(crefperp_k / (coscsi + waopc(i) * coscsi1), jprb)
        coscsi1_k             =  - crefperp_k * waopc(i) / (coscsi + waopc(i) * coscsi1)
        coscsi_k              = REAL(coscsi_k - crefperp_k * (coscsi - waopc(i) * coscsi1) / &
                                     (coscsi + waopc(i) * coscsi1) ** 2_jpim, jprb)
        coscsi1_k             = coscsi1_k - crefperp_k * waopc(i) * (coscsi - waopc(i) * coscsi1) / &
                                (coscsi + waopc(i) * coscsi1) ** 2_jpim
        crefperp_k            = 0
        crefpar_k             = crefpar_k + CONJG(crefparj_k)
        crefparj_k            = 0
        coscsi_k              = REAL(coscsi_k - crefpar_k * waopc(i) / (waopc(i) * coscsi + coscsi1), jprb)
        coscsi1_k             = coscsi1_k + crefpar_k / (waopc(i) * coscsi + coscsi1)
        coscsi_k              = REAL(coscsi_k + crefpar_k * waopc(i) * ( - coscsi1 + waopc(i) * coscsi) / &
                                     (waopc(i) * coscsi + coscsi1) ** 2_jpim, jprb)
        coscsi1_k             = coscsi1_k + crefpar_k * ( - coscsi1 + waopc(i) * coscsi) / &
                                (waopc(i) * coscsi + coscsi1) ** 2_jpim
        crefpar_k             = 0
        sincsi1_k             =  - coscsi1_k * sincsi1 / coscsi1
        coscsi1_k             = 0
        sincsi_k              = REAL(sincsi1_k / waopc(i), jprb)
        sincsi1_k             = 0
        sunglint_k%s(i)%omega = sunglint_k%s(i)%omega + sincsi_k * coscsi
        sunglint_k%s(i)%omega = sunglint_k%s(i)%omega - coscsi_k * sincsi
        sincsi_k              = 0
        coscsi_k              = 0
      ENDIF
    ENDIF
  ENDDO
  IF (LHOOK) CALL DR_HOOK('RTTOV_FRESNEL_K', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_fresnel_k
