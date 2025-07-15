! Description:
!> @file
!!   AD of calculation of variables required for cloud and aerosol simulations.
!
!> @brief
!!   AD of calculation of variables required for cloud and aerosol simulations.
!!
!! @param[in]     opts                 options to configure the simulations
!! @param[in]     profiles             input profiles
!! @param[in,out] profiles_ad          input profile increments
!! @param[in]     profiles_int         profiles in internal units
!! @param[in,out] profiles_int_ad      profile increments in internal units
!! @param[in]     aux                  auxiliary profile data structure
!! @param[in,out] aux_ad               auxiliary profile data increments
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
SUBROUTINE rttov_profaux_cldaer_ad( &
              opts,            &
              profiles,        &
              profiles_ad,     &
              profiles_int,    &
              profiles_int_ad, &
              aux,             &
              aux_ad)

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
  USE rttov_scattering_mod, ONLY: calc_rel_hum_ad
!INTF_ON
  IMPLICIT NONE

  TYPE(rttov_options),     INTENT(IN)    :: opts
  TYPE(rttov_profile),     INTENT(IN)    :: profiles(:)
  TYPE(rttov_profile),     INTENT(INOUT) :: profiles_ad(:)
  TYPE(rttov_profile),     INTENT(IN)    :: profiles_int(SIZE(profiles))
  TYPE(rttov_profile),     INTENT(INOUT) :: profiles_int_ad(SIZE(profiles))
  TYPE(rttov_profile_aux), INTENT(IN)    :: aux
  TYPE(rttov_profile_aux), INTENT(INOUT) :: aux_ad
!INTF_END

  INTEGER(jpim) :: nprofiles, nlayers
  INTEGER(jpim) :: iprof, lev, lay
  REAL(jprb) :: cloud, cloud_ad, clwde_min
  REAL(jprb) :: zradipou_upp, zradipou_low
  REAL(jprb) :: zmcfarq, ztempc, bwyser
  REAL(jprb) :: zmcfarq_ad, ztempc_ad, bwyser_ad
  REAL(jprb), PARAMETER :: rtt      = 273.15_jprb
  REAL(jprb), PARAMETER :: rtou_upp = -20._jprb
  REAL(jprb), PARAMETER :: rtou_low = -60._jprb
  REAL(jprb), PARAMETER :: nft = (SQRT(3._jprb) + 4._jprb) / (3._jprb * SQRT(3._jprb))
  REAL(jprb), PARAMETER :: amcfarq = 1.78449_jprb
  REAL(jprb), PARAMETER :: bmcfarq = 0.281301_jprb
  REAL(jprb), PARAMETER :: cmcfarq = 0.0177166_jprb
  REAL(jprb) :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------

  IF (LHOOK) CALL DR_HOOK('RTTOV_PROFAUX_CLDAER_AD', 0_jpim, ZHOOK_HANDLE)

  nprofiles = SIZE(profiles)
  nlayers = profiles(1)%nlayers

!-----------------------------------------------------------------------------------------
! AD of relative humidity calculation
!-----------------------------------------------------------------------------------------
  IF (opts%rt_ir%addaerosl .AND. .NOT. opts%rt_ir%user_aer_opt_param) THEN
    CALL calc_rel_hum_ad(opts, profiles, profiles_ad, profiles_int_ad, aux, aux_ad)
  ENDIF


  IF (opts%rt_ir%addclouds .AND. .NOT. opts%rt_ir%user_cld_opt_param) THEN

    DO iprof = 1, nprofiles

!-----------------------------------------------------------------------------------------
! Deff CLW scheme
!-----------------------------------------------------------------------------------------
      IF (profiles(iprof)%clw_scheme == clw_scheme_deff) THEN
        DO lay = 1, nlayers
          cloud = SUM(profiles_int(iprof)%cloud(1:nwcl_max,lay)) * 1.E-3_jprb ! kg/m^3
          IF (cloud > 0._jprb) THEN
            IF (profiles(iprof)%clwde(lay) > 0._jprb) THEN
              clwde_min = dgmin_clw
            ELSE
              IF (profiles(iprof)%clwde_param == 1) THEN
                clwde_min = martin_clwde_min
              ENDIF
            ENDIF

            IF (ABS(aux%clw_dg_ref(lay,iprof)) > 0._jprb) THEN
              IF (aux%clw_dg_ref(lay,iprof) < clwde_min) aux_ad%clw_dg(lay,iprof) = 0._jprb
              IF (aux%clw_dg_ref(lay,iprof) > dgmax_clw) aux_ad%clw_dg(lay,iprof) = 0._jprb
            ENDIF

            IF (profiles(iprof)%clwde(lay) > 0._jprb) THEN
              profiles_ad(iprof)%clwde(lay) = profiles_ad(iprof)%clwde(lay) + aux_ad%clw_dg(lay,iprof)
            ELSE
              IF (profiles(iprof)%clwde_param == 1) THEN
                cloud_ad = (1._jprb / 3._jprb) * aux_ad%clw_dg(lay,iprof) * &
                           aux%clw_dg(lay,iprof) / cloud
                profiles_int_ad(iprof)%cloud(1:nwcl_max,lay) = &
                  profiles_int_ad(iprof)%cloud(1:nwcl_max,lay) + cloud_ad * 1.E-3_jprb
              ENDIF
            ENDIF
          ENDIF
          aux_ad%clw_dg(lay,iprof) = 0._jprb
        ENDDO
      ENDIF

!-----------------------------------------------------------------------------------------
! AD for ice water content into effective generalized diameter
!-----------------------------------------------------------------------------------------
      IF (profiles(iprof)%ice_scheme == ice_scheme_baum) THEN
        zmcfarq_ad = 0._jprb
        ztempc_ad  = 0._jprb
        bwyser_ad  = 0._jprb

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

        DO lay = nlayers, 1,  - 1
          lev = lay + 1
          IF (profiles_int(iprof)%cloud(6,lay) > 0._jprb) THEN
            IF (ABS(aux%ice_dg_ref(lay,iprof)) > 0._jprb) THEN
              IF (aux%ice_dg_ref(lay,iprof) < dgmin_baum) aux_ad%ice_dg(lay,iprof) = 0._jprb
              IF (aux%ice_dg_ref(lay,iprof) > dgmax_baum) aux_ad%ice_dg(lay,iprof) = 0._jprb
            ENDIF
            ! Use effective diameter from input profile if specified
            IF (profiles(iprof)%icede(lay) > 0._jprb) THEN
              profiles_ad(iprof)%icede(lay) = profiles_ad(iprof)%icede(lay) + aux_ad%ice_dg(lay,iprof)
            ELSE
              IF (profiles(iprof)%icede_param == 1) THEN
                !Scheme by Ou and Liou, 1995, Atmos. Res., 35, 127-138.
                !
                ! Take Ou-Liou scheme as being valid only between 20C and 60C
                !
                IF (aux%ice_dg_ref(lay,iprof) < zradipou_low .OR. aux%ice_dg_ref(lay,iprof) > zradipou_upp) THEN
                  aux_ad%ice_dg(lay,iprof)      = 0._jprb
                ELSE
                  ztempc = profiles(iprof)%t(lev) - rtt
                  aux_ad%fac1_ice_dg(lay,iprof) = aux_ad%fac1_ice_dg(lay,iprof) + &
                      2._jprb * aux_ad%ice_dg(lay,iprof)* &
                      (0.388_jprb + 0.00051_jprb * 2 * aux%fac1_ice_dg(lay,iprof))
                  ztempc_ad                     = ztempc_ad + aux_ad%fac1_ice_dg(lay,iprof) * &
                      (12.42_jprb + 2._jprb * 0.197_jprb * ztempc + 3._jprb * 0.0012_jprb * ztempc * ztempc)
                  profiles_ad(iprof)%t(lev)     = profiles_ad(iprof)%t(lev) + ztempc_ad
                  ztempc_ad                     = 0._jprb
                ENDIF
              ELSE IF (profiles(iprof)%icede_param == 2) THEN
                !Scheme by Wyser et al. (see McFarquhar et al. (2003))
                IF (profiles(iprof)%t(lev) < 273._jprb) THEN
                  bwyser = -2.0_jprb + (0.001_jprb * ((273._jprb - profiles(iprof)%t(lev)) ** 1.5_jprb) * &
                                        LOG10(profiles_int(iprof)%cloud(6,lay) / 50._jprb))

                  aux_ad%fac1_ice_dg(lay,iprof) = aux_ad%fac1_ice_dg(lay,iprof) + &
                      (aux_ad%ice_dg(lay,iprof) / nft) * 2._jprb * 4._jprb * SQRT(3._jprb) / 9._jprb
                  bwyser_ad                     = bwyser_ad + aux_ad%fac1_ice_dg(lay,iprof) * &
                      (203.3_jprb + 2 * bwyser * 37.91_jprb + 3 * bwyser * bwyser * 2.3696_jprb)

                  profiles_ad(iprof)%t(lev)           = profiles_ad(iprof)%t(lev) - &
                      bwyser_ad * 0.001_jprb * 1.5_jprb * ((273._jprb - profiles(iprof)%t(lev)) ** 0.5_jprb) * &
                      LOG10(profiles_int(iprof)%cloud(6,lay) / 50._jprb)
                  profiles_int_ad(iprof)%cloud(6,lay) = profiles_int_ad(iprof)%cloud(6,lay) + &
                    bwyser_ad * 0.001_jprb * ((273._jprb - profiles(iprof)%t(lev)) ** 1.5_jprb) / &
                    profiles_int(iprof)%cloud(6,lay) / LOG(10._jprb)
                  bwyser_ad = 0._jprb
                ENDIF
              ELSE IF (profiles(iprof)%icede_param == 3) THEN
                !Scheme by Boudala et al., 2002, Int. J. Climatol., 22, 1267-1284.
                ztempc                              = profiles(iprof)%t(lev) - rtt
                profiles_int_ad(iprof)%cloud(6,lay) = profiles_int_ad(iprof)%cloud(6,lay) + &
                    aux_ad%ice_dg(lay,iprof) * 53.005_jprb * (profiles_int(iprof)%cloud(6,lay) ** (0.06_jprb - 1)) * &
                    EXP(0.013_jprb * ztempc) * 0.06_jprb
                ztempc_ad                           = ztempc_ad + &
                    aux_ad%ice_dg(lay,iprof) * 53.005_jprb * ((profiles_int(iprof)%cloud(6,lay)) ** 0.06_jprb) * 0.013_jprb * &
                    EXP(0.013_jprb * ztempc)
                profiles_ad(iprof)%t(lev)           = profiles_ad(iprof)%t(lev) + ztempc_ad
                ztempc_ad                           = 0._jprb
              ELSE IF (profiles(iprof)%icede_param == 4) THEN
                ! Scheme by McFarquhar et al. (2003)
                zmcfarq                             = profiles_int(iprof)%cloud(6,lay)
                aux_ad%fac1_ice_dg(lay,iprof)       = aux_ad%fac1_ice_dg(lay,iprof) + aux_ad%ice_dg(lay,iprof) * 2.0_jprb
                zmcfarq_ad                          = zmcfarq_ad + aux_ad%fac1_ice_dg(lay,iprof) * &
                    10.0_jprb ** (amcfarq + bmcfarq * LOG10(zmcfarq) + cmcfarq * LOG10(zmcfarq) * LOG10(zmcfarq)) * &
                    (bmcfarq + 2._jprb * cmcfarq * LOG10(zmcfarq)) / zmcfarq
                profiles_int_ad(iprof)%cloud(6,lay) = profiles_int_ad(iprof)%cloud(6,lay) + zmcfarq_ad
                zmcfarq_ad                          = 0._jprb
              ENDIF
            ENDIF
          ENDIF
          aux_ad%ice_dg(lay,iprof) = 0._jprb
        ENDDO
      ENDIF
    ENDDO
  ENDIF

  IF (LHOOK) CALL DR_HOOK('RTTOV_PROFAUX_CLDAER_AD', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_profaux_cldaer_ad
