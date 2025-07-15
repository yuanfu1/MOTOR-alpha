! Description:
!> @file
!!   TL of Fresnel reflectance calculation for flat water surface.
!
!> @brief
!!   TL of Fresnel reflectance calculation for flat water surface.
!!
!! @param[in]     chanprof       specifies channels and profiles to simulate
!! @param[in]     calcrefl       flags for internal RTTOV surface BRDF calculation
!! @param[in]     profiles       input atmospheric profiles and surface variables
!! @param[in]     solar          flag to indicate channels with solar contribution
!! @param[in]     coef           optical depth coefficient structure
!! @param[in]     sunglint       internal structure for sea surface BRDF model variables
!! @param[in]     sunglint_tl    sea surface BRDF model variable perturbations
!! @param[out]    fresnrefl_tl   calculated Fresnel reflectance perturbations
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
SUBROUTINE rttov_fresnel_tl( &
              chanprof,     &
              calcrefl,     &
              profiles,     &
              solar,        &
              coef,         &
              sunglint,     &
              sunglint_tl,  &
              fresnrefl_tl)

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
  TYPE(rttov_sunglint),  INTENT(IN)  :: sunglint_tl
  REAL(KIND=jprb),       INTENT(OUT) :: fresnrefl_tl(SIZE(chanprof))
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
  REAL   (KIND=jprb) :: sincsi_tl
  REAL   (KIND=jprb) :: coscsi_tl
  COMPLEX(KIND=jprb) :: crefpar_tl
  COMPLEX(KIND=jprb) :: crefparj_tl
  COMPLEX(KIND=jprb) :: crefperp_tl
  COMPLEX(KIND=jprb) :: crefperpj_tl
  COMPLEX(KIND=jprb) :: sincsi1_tl
  COMPLEX(KIND=jprb) :: coscsi1_tl
  INTEGER(KIND=jpim) :: nchannels
  COMPLEX(KIND=jprb) :: waopc(SIZE(chanprof))
  REAL   (KIND=jprb) :: ZHOOK_HANDLE
!-----End of header-------------------------------------------------------------

  IF (LHOOK) CALL DR_HOOK('RTTOV_FRESNEL_TL', 0_jpim, ZHOOK_HANDLE)
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
        coscsi          = COS(sunglint%s(prof)%omega)
        sincsi          = SIN(sunglint%s(prof)%omega)
        sincsi1         = sincsi / waopc(i)
        coscsi1         = SQRT(1._jprb - sincsi1 ** 2_jpim)
        crefpar         =  - (waopc(i) * coscsi - coscsi1) / (waopc(i) * coscsi + coscsi1)
        crefparj        = CONJG(crefpar)
        crefperp        = (coscsi - waopc(i) * coscsi1) / (coscsi + waopc(i) * coscsi1)
        crefperpj       = CONJG(crefperp)
        coscsi_tl       =  - sincsi * sunglint_tl%s(prof)%omega
        sincsi_tl       = coscsi * sunglint_tl%s(prof)%omega
        sincsi1_tl      = sincsi_tl / waopc(i)
        coscsi1_tl      =  - sincsi1_tl * sincsi1 / coscsi1
        crefpar_tl      =                                                             &
          &  - (waopc(i) * coscsi_tl - coscsi1_tl) / (waopc(i) * coscsi + coscsi1) +  &
          & (waopc(i) * coscsi_tl + coscsi1_tl) *                                     &
          & (waopc(i) * coscsi - coscsi1) / (waopc(i) * coscsi + coscsi1) ** 2_jpim
        crefparj_tl     = CONJG(crefpar_tl)
        crefperp_tl     =                                                          &
          & (coscsi_tl - waopc(i) * coscsi1_tl) / (coscsi + waopc(i) * coscsi1) -  &
          & (coscsi_tl + waopc(i) * coscsi1_tl) *                                  &
          & (coscsi - waopc(i) * coscsi1) / (coscsi + waopc(i) * coscsi1) ** 2_jpim
        crefperpj_tl    = CONJG(crefperp_tl)
        fresnrefl_tl(i) = REAL (crefpar_tl * 0.5_jprb * crefparj + crefparj_tl * 0.5_jprb * crefpar + &
                                crefperp_tl * 0.5_jprb * crefperpj + crefperpj_tl * 0.5_jprb * crefperp, jprb)
      ENDIF
    ENDIF
  ENDDO
  IF (LHOOK) CALL DR_HOOK('RTTOV_FRESNEL_TL', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_fresnel_tl
