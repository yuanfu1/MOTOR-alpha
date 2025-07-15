! Description:
!> @file
!!   K of calculation of variables required for cloud and aerosol simulations.
!
!> @brief
!!   K of calculation of variables required for cloud and aerosol simulations.
!!
!! @param[in]     opts                 options to configure the simulations
!! @param[in]     chanprof             channels and profiles to simulate
!! @param[in]     profiles             input profiles
!! @param[in,out] profiles_k           input profile increments
!! @param[in]     profiles_int         profiles in internal units
!! @param[in,out] profiles_int_k       profile increments in internal units
!! @param[in]     aux                  auxiliary profile data structure
!! @param[in,out] aux_k                auxiliary profile data increments
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
SUBROUTINE rttov_profaux_cldaer_k( &
              opts,           &
              chanprof,       &
              profiles,       &
              profiles_k,     &
              profiles_int,   &
              profiles_int_k, &
              aux,            &
              aux_k)

  USE rttov_types, ONLY :  &
        rttov_chanprof,    &
        rttov_options,     &
        rttov_profile,     &
        rttov_profile_aux
!INTF_OFF
  USE rttov_const, ONLY : dgmin_clw, dgmax_clw, dgmin_baum, dgmax_baum, &
                          clw_scheme_deff, ice_scheme_baum, nwcl_max, &
                          martin_clwde_min
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jpim, jprb
  USE rttov_scattering_mod, ONLY: calc_rel_hum_k
!INTF_ON
  IMPLICIT NONE

  TYPE(rttov_options),     INTENT(IN)    :: opts
  TYPE(rttov_chanprof),    INTENT(IN)    :: chanprof(:)
  TYPE(rttov_profile),     INTENT(IN)    :: profiles(:)
  TYPE(rttov_profile),     INTENT(INOUT) :: profiles_k(SIZE(chanprof))
  TYPE(rttov_profile),     INTENT(IN)    :: profiles_int(SIZE(profiles))
  TYPE(rttov_profile),     INTENT(INOUT) :: profiles_int_k(SIZE(chanprof))
  TYPE(rttov_profile_aux), INTENT(IN)    :: aux
  TYPE(rttov_profile_aux), INTENT(INOUT) :: aux_k
!INTF_END

  INTEGER(jpim) :: nchanprof, nlayers
  INTEGER(jpim) :: i, prof, lev, lay
  REAL(jprb) :: cloud, cloud_k, clwde_min
  REAL(jprb) :: zradipou_upp, zradipou_low
  REAL(jprb) :: zmcfarq, ztempc, bwyser
  REAL(jprb) :: zmcfarq_k, ztempc_k, bwyser_k
  REAL(jprb), PARAMETER :: rtt      = 273.15_jprb
  REAL(jprb), PARAMETER :: rtou_upp = -20._jprb
  REAL(jprb), PARAMETER :: rtou_low = -60._jprb
  REAL(jprb), PARAMETER :: nft = (SQRT(3._jprb) + 4._jprb) / (3._jprb * SQRT(3._jprb))
  REAL(jprb), PARAMETER :: amcfarq = 1.78449_jprb
  REAL(jprb), PARAMETER :: bmcfarq = 0.281301_jprb
  REAL(jprb), PARAMETER :: cmcfarq = 0.0177166_jprb
  REAL(jprb)    :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------

  IF (LHOOK) CALL DR_HOOK('RTTOV_PROFAUX_CLDAER_K', 0_jpim, ZHOOK_HANDLE)

  nchanprof = SIZE(chanprof)
  nlayers = profiles(1)%nlayers

!-----------------------------------------------------------------------------------------
! K of relative humidity calculation
!-----------------------------------------------------------------------------------------
  IF (opts%rt_ir%addaerosl .AND. .NOT. opts%rt_ir%user_aer_opt_param) THEN
    CALL calc_rel_hum_k(opts, chanprof, profiles, profiles_k, profiles_int_k, aux, aux_k)
  ENDIF


  IF (opts%rt_ir%addclouds .AND. .NOT. opts%rt_ir%user_cld_opt_param) THEN

!-----------------------------------------------------------------------------------------
! Deff CLW scheme
!-----------------------------------------------------------------------------------------
    DO i = 1, nchanprof
      prof = chanprof(i)%prof
      IF (profiles(prof)%clw_scheme == clw_scheme_deff) THEN
        DO lay = 1, nlayers
          cloud = SUM(profiles_int(prof)%cloud(1:nwcl_max,lay)) * 1.E-3_jprb ! kg/m^3
          IF (cloud > 0._jprb) THEN
            IF (profiles(prof)%clwde(lay) > 0._jprb) THEN
              clwde_min = dgmin_clw
            ELSE
              IF (profiles(prof)%clwde_param == 1) THEN
                clwde_min = martin_clwde_min
              ENDIF
            ENDIF

            IF (ABS(aux%clw_dg_ref(lay,prof)) > 0._jprb) THEN
              IF (aux%clw_dg_ref(lay,prof) < clwde_min) aux_k%clw_dg(lay,i) = 0._jprb
              IF (aux%clw_dg_ref(lay,prof) > dgmax_clw) aux_k%clw_dg(lay,i) = 0._jprb
            ENDIF

            IF (profiles(prof)%clwde(lay) > 0._jprb) THEN
              profiles_k(i)%clwde(lay) = profiles_k(i)%clwde(lay) + aux_k%clw_dg(lay,i)
            ELSE
              IF (profiles(prof)%clwde_param == 1) THEN
                cloud_k = (1._jprb / 3._jprb) * aux_k%clw_dg(lay,i) * &
                          aux%clw_dg(lay,prof) / cloud
                profiles_int_k(i)%cloud(1:nwcl_max,lay) = &
                  profiles_int_k(i)%cloud(1:nwcl_max,lay) + cloud_k * 1.E-3_jprb
              ENDIF
            ENDIF
          ENDIF
          aux_k%clw_dg(lay,i) = 0._jprb
        ENDDO
      ENDIF
    ENDDO

!------------------------------------------------------------
! K for ice water content into effective generalized diameter
!------------------------------------------------------------
    zmcfarq_k = 0._jprb
    ztempc_k  = 0._jprb
    bwyser_k  = 0._jprb

    ! Calculate upper and lower limits for Ou-Liou effective size
    !
    zradipou_upp = 326.3_jprb + rtou_upp * (12.42_jprb + rtou_upp * (0.197_jprb + rtou_upp * 0.0012_jprb))
    zradipou_low = 326.3_jprb + rtou_low * (12.42_jprb + rtou_low * (0.197_jprb + rtou_low * 0.0012_jprb))
    !
    ! and convert these to the "generalized" effective size used here (using McFarquhar et al 2003 equation),
    ! not forgetting the factor of 2 to convert from McFarquhar's radius to a diameter
    !
    zradipou_upp =  - 1.56_jprb + zradipou_upp * (0.388_jprb + zradipou_upp * 0.00051_jprb)
    zradipou_upp = 2.0_jprb * zradipou_upp
    zradipou_low =  - 1.56_jprb + zradipou_low * (0.388_jprb + zradipou_low * 0.00051_jprb)
    zradipou_low = 2.0_jprb * zradipou_low

    DO i = 1, nchanprof
      prof = chanprof(i)%prof
      IF (.NOT. profiles(prof)%ice_scheme == ice_scheme_baum) CYCLE
      DO lay = nlayers, 1, -1
        lev = lay + 1
        IF (profiles_int(prof)%cloud(6,lay) > 0._jprb) THEN
          IF (ABS(aux%ice_dg_ref(lay,prof)) > 0._jprb) THEN
            IF (aux%ice_dg_ref(lay,prof) < dgmin_baum) aux_k%ice_dg(lay,i) = 0._jprb
            IF (aux%ice_dg_ref(lay,prof) > dgmax_baum) aux_k%ice_dg(lay,i) = 0._jprb
          ENDIF
          ! Use effective diameter from input profile if specified
          IF (profiles(prof)%icede(lay) > 0._jprb) THEN
            profiles_k(i)%icede(lay) = profiles_k(i)%icede(lay) + aux_k%ice_dg(lay,i)
          ELSE
            IF (profiles(prof)%icede_param == 1) THEN
              ! Scheme by Ou and Liou, 1995, Atmos. Res., 35, 127-138.
              !
              ! Take Ou-Liou scheme as being valid only between 20C and 60C
              !
              IF (aux%ice_dg_ref(lay,prof) < zradipou_low .OR. aux%ice_dg_ref(lay,prof) > zradipou_upp) THEN
                aux_k%ice_dg(lay,i)      = 0._jprb
              ELSE
                ztempc = profiles(prof)%t(lev) - rtt
                aux_k%fac1_ice_dg(lay,i) = aux_k%fac1_ice_dg(lay,i) + &
                    2._jprb * aux_k%ice_dg(lay,i) * &
                    (0.388_jprb + 0.00051_jprb * 2 * aux%fac1_ice_dg(lay,prof))
                ztempc_k                 = ztempc_k + aux_k%fac1_ice_dg(lay,i) * &
                    (12.42_jprb + 2._jprb * 0.197_jprb * ztempc + 3._jprb * 0.0012_jprb * ztempc * ztempc)
                profiles_k(i)%t(lev)     = profiles_k(i)%t(lev) + ztempc_k
                ztempc_k                 = 0._jprb
              ENDIF
            ELSE IF (profiles(prof)%icede_param == 2) THEN
              ! Scheme by Wyser et al. (see McFarquhar et al. (2003))
              IF (profiles(prof)%t(lev) < 273._jprb) THEN
                bwyser = -2.0_jprb + (0.001_jprb * ((273._jprb - profiles(prof)%t(lev)) ** 1.5_jprb) * &
                                      LOG10(profiles_int(prof)%cloud(6,lay) / 50._jprb))

                aux_k%fac1_ice_dg(lay,i) = aux_k%fac1_ice_dg(lay,i) + &
                    (aux_k%ice_dg(lay,i) / nft) * 2._jprb * 4._jprb * SQRT(3._jprb) / 9._jprb
                bwyser_k                 = bwyser_k + aux_k%fac1_ice_dg(lay,i) * &
                    (203.3_jprb + 2 * bwyser * 37.91_jprb + 3 * bwyser * bwyser * 2.3696_jprb)

                profiles_k(i)%t(lev)           = profiles_k(i)%t(lev) - &
                    bwyser_k * 0.001_jprb * 1.5_jprb * ((273._jprb - profiles(prof)%t(lev)) ** 0.5_jprb) * &
                    LOG10(profiles_int(prof)%cloud(6,lay) / 50._jprb)
                profiles_int_k(i)%cloud(6,lay) = profiles_int_k(i)%cloud(6,lay) + &
                    bwyser_k * 0.001_jprb * ((273._jprb - profiles(prof)%t(lev)) ** 1.5_jprb) / &
                    profiles_int(prof)%cloud(6,lay) / LOG(10._jprb)
                bwyser_k = 0._jprb
              ENDIF
            ELSE IF (profiles(prof)%icede_param == 3) THEN
              ! Scheme by Boudala et al., 2002, Int. J. Climatol., 22, 1267-1284.
              ztempc                         = profiles(prof)%t(lev) - rtt
              profiles_int_k(i)%cloud(6,lay) = profiles_int_k(i)%cloud(6,lay) + &
                  aux_k%ice_dg(lay,i) * 53.005_jprb * (profiles_int(prof)%cloud(6,lay) ** (0.06_jprb - 1)) * &
                  EXP(0.013_jprb * ztempc) * 0.06_jprb
              ztempc_k                       = ztempc_k + &
                  aux_k%ice_dg(lay,i) * 53.005_jprb * ((profiles_int(prof)%cloud(6,lay)) ** 0.06_jprb) * 0.013_jprb * &
                  EXP(0.013_jprb * ztempc)
              profiles_k(i)%t(lev)           = profiles_k(i)%t(lev) + ztempc_k
              ztempc_k                       = 0._jprb
            ELSE IF (profiles(prof)%icede_param == 4) THEN
              ! Scheme by McFarquhar et al. (2003)
              zmcfarq                        = profiles_int(prof)%cloud(6,lay)
              aux_k%fac1_ice_dg(lay,i)       = aux_k%fac1_ice_dg(lay,i) + aux_k%ice_dg(lay,i) * 2.0_jprb
              zmcfarq_k                      = zmcfarq_k + aux_k%fac1_ice_dg(lay,i) * &
                  10.0_jprb ** (amcfarq + bmcfarq * LOG10(zmcfarq) + cmcfarq * LOG10(zmcfarq) * LOG10(zmcfarq)) * &
                  (bmcfarq + 2._jprb * cmcfarq * LOG10(zmcfarq)) / zmcfarq
              profiles_int_k(i)%cloud(6,lay) = profiles_int_k(i)%cloud(6,lay) + zmcfarq_k
              zmcfarq_k                      = 0._jprb
            ENDIF
          ENDIF
        ENDIF
        aux_k%ice_dg(lay,i) = 0._jprb
      ENDDO
    ENDDO
  ENDIF

  IF (LHOOK) CALL DR_HOOK('RTTOV_PROFAUX_CLDAER_K', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_profaux_cldaer_k
