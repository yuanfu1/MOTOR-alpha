! Description:
!> @file
!!   TL of calculation of variables required for cloud and aerosol simulations.
!
!> @brief
!!   TL of calculation of variables required for cloud and aerosol simulations.
!!
!! @param[in]     opts                 options to configure the simulations
!! @param[in]     profiles             input profiles
!! @param[in]     profiles_tl          input profile perturbations
!! @param[in]     profiles_int         profiles in internal units
!! @param[in]     profiles_int_tl      profile perturbations in internal units
!! @param[in]     aux                  auxiliary profile data structure
!! @param[in,out] aux_tl               auxiliary profile data perturbations
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
SUBROUTINE rttov_profaux_cldaer_tl( &
              opts,            &
              profiles,        &
              profiles_tl,     &
              profiles_int,    &
              profiles_int_tl, &
              aux,             &
              aux_tl)

  USE rttov_types, ONLY :  &
        rttov_options,     &
        rttov_profile,     &
        rttov_profile_aux
!INTF_OFF
  USE rttov_const, ONLY : dgmin_clw, dgmax_clw, dgmin_baum, dgmax_baum, &
                          clw_scheme_deff, ice_scheme_baum, nwcl_max, &
                          martin_clwde_min
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jpim, jprb
  USE rttov_scattering_mod, ONLY: calc_rel_hum_tl
!INTF_ON

  IMPLICIT NONE

  TYPE(rttov_options),     INTENT(IN)    :: opts
  TYPE(rttov_profile),     INTENT(IN)    :: profiles(:)
  TYPE(rttov_profile),     INTENT(IN)    :: profiles_tl(:)
  TYPE(rttov_profile),     INTENT(IN)    :: profiles_int(SIZE(profiles))
  TYPE(rttov_profile),     INTENT(IN)    :: profiles_int_tl(SIZE(profiles))
  TYPE(rttov_profile_aux), INTENT(IN)    :: aux
  TYPE(rttov_profile_aux), INTENT(INOUT) :: aux_tl
!INTF_END

  INTEGER(jpim) :: nprofiles, nlayers
  INTEGER(jpim) :: iprof, lev, lay
  REAL(jprb) :: cloud, cloud_tl, clwde_min
  REAL(jprb) :: zradipou_upp, zradipou_low
  REAL(jprb) :: zmcfarq, ztempc, bwyser
  REAL(jprb) :: zmcfarq_tl, ztempc_tl, bwyser_tl
  REAL(jprb), PARAMETER :: rtt      = 273.15_jprb
  REAL(jprb), PARAMETER :: rtou_upp = -20._jprb
  REAL(jprb), PARAMETER :: rtou_low = -60._jprb
  REAL(jprb), PARAMETER :: nft = (SQRT(3._jprb) + 4._jprb) / (3._jprb * SQRT(3._jprb))
  REAL(jprb), PARAMETER :: amcfarq = 1.78449_jprb
  REAL(jprb), PARAMETER :: bmcfarq = 0.281301_jprb
  REAL(jprb), PARAMETER :: cmcfarq = 0.0177166_jprb
  REAL(jprb) :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------

  IF (LHOOK) CALL DR_HOOK('RTTOV_PROFAUX_CLDAER_TL', 0_jpim, ZHOOK_HANDLE)

  nprofiles = SIZE(profiles)
  nlayers = profiles(1)%nlayers

  IF (opts%rt_ir%addclouds .AND. .NOT. opts%rt_ir%user_cld_opt_param) THEN
    DO iprof = 1, nprofiles

!-----------------------------------------------------------------------------------------
! Deff CLW scheme
!-----------------------------------------------------------------------------------------
      IF (profiles(iprof)%clw_scheme == clw_scheme_deff) THEN
        DO lay = 1, nlayers
          aux_tl%clw_dg(lay,iprof) = 0._jprb
          cloud = SUM(profiles_int(iprof)%cloud(1:nwcl_max,lay)) * 1.E-3_jprb ! kg/m^3
          IF (cloud > 0._jprb) THEN
            IF (profiles(iprof)%clwde(lay) > 0._jprb) THEN
              aux_tl%clw_dg(lay,iprof) = profiles_tl(iprof)%clwde(lay)
              clwde_min = dgmin_clw
            ELSE
              IF (profiles(iprof)%clwde_param == 1) THEN
                cloud_tl = SUM(profiles_int_tl(iprof)%cloud(1:nwcl_max,lay)) * 1.E-3_jprb
                aux_tl%clw_dg(lay,iprof) = &
                  (1._jprb / 3._jprb) * cloud_tl * aux%clw_dg(lay,iprof) / cloud
                clwde_min = martin_clwde_min
              ENDIF
            ENDIF

            IF (ABS(aux%clw_dg_ref(lay,iprof)) > 0._jprb) THEN
              IF (aux%clw_dg_ref(lay,iprof) < clwde_min) aux_tl%clw_dg(lay,iprof) = 0._jprb
              IF (aux%clw_dg_ref(lay,iprof) > dgmax_clw) aux_tl%clw_dg(lay,iprof) = 0._jprb
            ENDIF
          ENDIF
        ENDDO
      ENDIF

!--------------------------------------------------------------
! TL for ice water content into effective generalized diameter
!--------------------------------------------------------------
      IF (profiles(iprof)%ice_scheme == ice_scheme_baum) THEN
        IF (profiles(iprof)%icede_param == 1) THEN
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
        ENDIF

        DO lay = 1, nlayers
          lev = lay + 1
          aux_tl%ice_dg(lay,iprof) = 0._jprb
          IF (profiles_int(iprof)%cloud(6,lay) > 0._jprb) THEN
            ! Use effective diameter from input profile if specified
            IF (profiles(iprof)%icede(lay) > 0._jprb) THEN
              aux_tl%ice_dg(lay,iprof) = profiles_tl(iprof)%icede(lay)
            ELSE
              IF (profiles(iprof)%icede_param == 1) THEN
                !Scheme by Ou and Liou, 1995, Atmos. Res., 35, 127-138.
                !
                ! Take Ou-Liou scheme as being valid only between 20C and 60C
                !
                IF (aux%ice_dg_ref(lay,iprof) < zradipou_low .OR. aux%ice_dg_ref(lay,iprof) > zradipou_upp) THEN
                  aux_tl%ice_dg(lay,iprof)      = 0._jprb
                ELSE
                  ztempc                        = profiles(iprof)%t(lev) - rtt
                  ztempc_tl                     = profiles_tl(iprof)%t(lev)
                  aux_tl%fac1_ice_dg(lay,iprof) = &
                      ztempc_tl * (12.42_jprb + 2._jprb * 0.197_jprb * ztempc + 3._jprb * 0.0012_jprb * ztempc * ztempc)
                  aux_tl%ice_dg(lay,iprof)      = 2._jprb * aux_tl%fac1_ice_dg(lay,iprof) * &
                      (0.388_jprb + 0.00051_jprb * 2 * aux%fac1_ice_dg(lay,iprof))
                ENDIF
              ELSE IF (profiles(iprof)%icede_param == 2) THEN
                !Scheme by Wyser et al. (see McFarquhar et al. (2003))
                IF (profiles(iprof)%t(lev) < 273._jprb) THEN
                  bwyser = -2.0_jprb + (0.001_jprb * ((273._jprb - profiles(iprof)%t(lev)) ** 1.5_jprb) * &
                                        LOG10(profiles_int(iprof)%cloud(6,lay) / 50._jprb))
                  bwyser_tl = &
                      -0.001_jprb * 1.5_jprb * ((273._jprb - profiles(iprof)%t(lev)) ** 0.5_jprb) * &
                      profiles_tl(iprof)%t(lev) * LOG10(profiles_int(iprof)%cloud(6,lay) / 50._jprb) + &
                      0.001_jprb * ((273._jprb - profiles(iprof)%t(lev)) ** 1.5_jprb) * &
                      profiles_int_tl(iprof)%cloud(6,lay) / profiles_int(iprof)%cloud(6,lay) / LOG(10._jprb)
                  aux_tl%fac1_ice_dg(lay,iprof) = &
                      (203.3_jprb + 2 * bwyser * 37.91_jprb + 3 * bwyser * bwyser * 2.3696_jprb) * bwyser_tl
                  aux_tl%ice_dg(lay,iprof) = 2._jprb * 4._jprb * (aux_tl%fac1_ice_dg(lay,iprof) / nft) * SQRT(3._jprb) / 9._jprb
                ELSE
                  aux_tl%ice_dg(lay,iprof) = 0._jprb
                ENDIF
              ELSE IF (profiles(iprof)%icede_param == 3) THEN
                !Scheme by Boudala et al., 2002, Int. J. Climatol., 22, 1267-1284.
                ztempc                   = profiles(iprof)%t(lev) - rtt
                ztempc_tl                = profiles_tl(iprof)%t(lev)
                aux_tl%ice_dg(lay,iprof) = 53.005_jprb * 0.06_jprb * &
                  (profiles_int(iprof)%cloud(6,lay) ** (0.06_jprb - 1)) * profiles_int_tl(iprof)%cloud(6,lay) * &
                  EXP(0.013_jprb * ztempc) + 53.005_jprb * &
                  ((profiles_int(iprof)%cloud(6,lay)) ** 0.06_jprb) * 0.013_jprb * ztempc_tl * EXP(0.013_jprb * ztempc)
              ELSE IF (profiles(iprof)%icede_param == 4) THEN
                ! Scheme by McFarquhar et al. (2003)
                zmcfarq                       = profiles_int(iprof)%cloud(6,lay)
                zmcfarq_tl                    = profiles_int_tl(iprof)%cloud(6,lay)
                aux_tl%fac1_ice_dg(lay,iprof) = &
                    10.0_jprb ** (amcfarq + bmcfarq * LOG10(zmcfarq) + cmcfarq * LOG10(zmcfarq) * LOG10(zmcfarq)) * &
                    (bmcfarq + 2._jprb * cmcfarq * LOG10(zmcfarq)) / zmcfarq * zmcfarq_tl
                aux_tl%ice_dg(lay,iprof)      = 2.0_jprb * aux_tl%fac1_ice_dg(lay,iprof)
              ENDIF
            ENDIF

            IF (ABS(aux%ice_dg_ref(lay,iprof)) > 0._jprb) THEN
              IF (aux%ice_dg_ref(lay,iprof) < dgmin_baum) aux_tl%ice_dg(lay,iprof) = 0._jprb
              IF (aux%ice_dg_ref(lay,iprof) > dgmax_baum) aux_tl%ice_dg(lay,iprof) = 0._jprb
            ENDIF
          ENDIF
        ENDDO
      ENDIF
    ENDDO
  ENDIF

!-----------------------------------------------------------------------------------------
! TL of relative humidity calculation
!-----------------------------------------------------------------------------------------
  IF (opts%rt_ir%addaerosl .AND. .NOT. opts%rt_ir%user_aer_opt_param) THEN
    CALL calc_rel_hum_tl(opts, profiles, profiles_tl, profiles_int_tl, aux, aux_tl)
  ENDIF

  IF (LHOOK) CALL DR_HOOK('RTTOV_PROFAUX_CLDAER_TL', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_profaux_cldaer_tl
