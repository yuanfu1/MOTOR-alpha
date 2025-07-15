! Description:
!> @file
!!   Calculate some variables required for cloud and aerosol simulations.
!
!> @brief
!!   Calculate some variables required for cloud and aerosol simulations.
!!
!> @details
!!   This subroutine calculates the following information:
!!   - clw cloud Deff parameterisation for VIS/IR scattering
!!   - ice cloud Deff parameterisations for VIS/IR scattering
!!   - relative humidity for aerosol VIS/IR scattering
!!
!! @param[in]     opts                 options to configure the simulations
!! @param[in]     profiles             input profiles
!! @param[in]     profiles_int         profiles in internal units
!! @param[in,out] aux                  auxiliary profile data structure
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
SUBROUTINE rttov_profaux_cldaer(opts, profiles, profiles_int, aux)

  USE rttov_types, ONLY :  &
        rttov_options,     &
        rttov_profile,     &
        rttov_profile_aux
!INTF_OFF
  USE rttov_const, ONLY : dgmin_clw, dgmax_clw, dgmin_baum, dgmax_baum, &
                          clw_scheme_deff, ice_scheme_baum, nwcl_max, pi, &
                          surftype_land, rho_water, martin_clwde_min, &
                          martin_k_land, martin_ntot_land, &
                          martin_k_sea, martin_ntot_sea
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jpim, jprb
  USE rttov_scattering_mod, ONLY: calc_rel_hum
!INTF_ON
  IMPLICIT NONE
  TYPE(rttov_options),     INTENT(IN)    :: opts
  TYPE(rttov_profile),     INTENT(IN)    :: profiles(:)
  TYPE(rttov_profile),     INTENT(IN)    :: profiles_int(:)
  TYPE(rttov_profile_aux), INTENT(INOUT) :: aux
!INTF_END

  INTEGER(jpim) :: nprofiles, nlayers
  INTEGER(jpim) :: iprof, lev, lay
  REAL(jprb) :: cloud, clwde_min, k, ntot
  REAL(jprb) :: zradipou_upp, zradipou_low
  REAL(jprb) :: zmcfarq, ztempc, bwyser
  REAL(jprb), PARAMETER :: rtt      = 273.15_jprb
  REAL(jprb), PARAMETER :: rtou_upp = -20._jprb ! Equivalent to zradipou_upp = 133.1007782 um
  REAL(jprb), PARAMETER :: rtou_low = -60._jprb ! Equivalent to zradipou_low = 22.00015420 um
  REAL(jprb), PARAMETER :: nft = (SQRT(3._jprb) + 4._jprb) / (3._jprb * SQRT(3._jprb))
  REAL(jprb), PARAMETER :: amcfarq = 1.78449_jprb
  REAL(jprb), PARAMETER :: bmcfarq = 0.281301_jprb
  REAL(jprb), PARAMETER :: cmcfarq = 0.0177166_jprb
  REAL(jprb) :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------

  IF (LHOOK) CALL DR_HOOK('RTTOV_PROFAUX_CLDAER', 0_jpim, ZHOOK_HANDLE)

  nprofiles = SIZE(profiles)
  nlayers = profiles(1)%nlayers

  IF (opts%rt_ir%addclouds .AND. .NOT. opts%rt_ir%user_cld_opt_param) THEN

    DO iprof = 1, nprofiles

      !------------------------------------------------------------------------
      ! For Deff CLW scheme store effective diameter
      !------------------------------------------------------------------------
      IF (profiles(iprof)%clw_scheme == clw_scheme_deff) THEN
        DO lay = 1, nlayers
          aux%clw_dg(lay,iprof) = 0._jprb
          cloud = SUM(profiles_int(iprof)%cloud(1:nwcl_max,lay)) * 1.E-3_jprb ! kg/m^3
          IF (cloud > 0._jprb) THEN
            IF (profiles(iprof)%clwde(lay) > 0._jprb) THEN
              ! Use effective diameter from input profile if specified
              aux%clw_dg_ref(lay,iprof) = profiles(iprof)%clwde(lay)
              aux%clw_dg(lay,iprof)     = profiles(iprof)%clwde(lay)
              ! Minimum size determined by optical property limits
              clwde_min = dgmin_clw
            ELSE
              IF (profiles(iprof)%clwde_param == 1) THEN
                ! Martin et al (1994) parameterisation
                IF (profiles(iprof)%skin%surftype == surftype_land) THEN
                  k = martin_k_land
                  ntot = martin_ntot_land
                ELSE
                  k = martin_k_sea
                  ntot = martin_ntot_sea
                ENDIF
                aux%clw_dg_ref(lay,iprof) = 2._jprb * 1.E6_jprb * &
                  (0.75_jprb * cloud / (pi * k * ntot * rho_water))**(1._jprb / 3._jprb)

                aux%clw_dg(lay,iprof) = aux%clw_dg_ref(lay,iprof)
                ! More realistic minimum size imposed to avoid small values at low LWC
                clwde_min = martin_clwde_min
              ENDIF
            ENDIF

            IF (ABS(aux%clw_dg_ref(lay,iprof)) > 0._jprb) THEN
              ! The minimum is specified above, maximum is determined by optical property limits
              IF (aux%clw_dg_ref(lay,iprof) < clwde_min) aux%clw_dg(lay,iprof) = clwde_min
              IF (aux%clw_dg_ref(lay,iprof) > dgmax_clw) aux%clw_dg(lay,iprof) = dgmax_clw
            ENDIF
          ENDIF
        ENDDO
      ENDIF

      !------------------------------------------------------------------------
      ! For ice clouds convert IWC into effective generalized diameter
      !------------------------------------------------------------------------
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
          aux%ice_dg(lay,iprof) = 0._jprb
          IF (profiles_int(iprof)%cloud(6,lay) > 0._jprb) THEN
            ! Use effective diameter from input profile if specified
            IF (profiles(iprof)%icede(lay) > 0._jprb) THEN
              aux%ice_dg_ref(lay,iprof) = profiles(iprof)%icede(lay)
              aux%ice_dg(lay,iprof)     = profiles(iprof)%icede(lay)
            ELSE
              IF (profiles(iprof)%icede_param == 1) THEN
                ! Scheme by Ou and Liou, 1995, Atmos. Res., 35, 127-138.
                ztempc                     = profiles(iprof)%t(lev) - rtt
                ! intermediate factors in calculating the generalized effective diameter
                aux%fac1_ice_dg(lay,iprof) = 326.3_jprb + ztempc * (12.42_jprb + ztempc * (0.197_jprb + ztempc * 0.0012_jprb))
                aux%ice_dg_ref(lay,iprof)  = 2._jprb * &
                  (-1.56_jprb + aux%fac1_ice_dg(lay,iprof) * (0.388_jprb + aux%fac1_ice_dg(lay,iprof) * 0.00051_jprb))
                !
                ! Take Ou-Liou scheme as being valid only between -20C and -60C
                !
                aux%ice_dg(lay,iprof)      = MAX(aux%ice_dg_ref(lay,iprof), zradipou_low)
                aux%ice_dg(lay,iprof)      = MIN(aux%ice_dg(lay,iprof), zradipou_upp)
              ELSE IF (profiles(iprof)%icede_param == 2) THEN
                ! Scheme by Wyser et al. (see McFarquhar et al. (2003))
                bwyser = -2.0_jprb
                IF (profiles(iprof)%t(lev) < 273._jprb) THEN
                  bwyser = bwyser + (0.001_jprb * ((273._jprb - profiles(iprof)%t(lev)) ** 1.5_jprb) * &
                                     LOG10(profiles_int(iprof)%cloud(6,lay) / 50._jprb))
                ENDIF
                aux%fac1_ice_dg(lay,iprof) = 377.4_jprb + bwyser * (203.3_jprb + bwyser * (37.91_jprb + bwyser * 2.3696_jprb))
                aux%ice_dg_ref(lay,iprof)  = 2._jprb * 4._jprb * (aux%fac1_ice_dg(lay,iprof) / nft) * SQRT(3._jprb) / 9._jprb
                aux%ice_dg(lay,iprof)      = aux%ice_dg_ref(lay,iprof)
              ELSE IF (profiles(iprof)%icede_param == 3) THEN
                ! Scheme by Boudala et al., 2002, Int. J. Climatol., 22, 1267-1284.
                ztempc                    = profiles(iprof)%t(lev) - rtt
                aux%ice_dg_ref(lay,iprof) = 53.005_jprb * ((profiles_int(iprof)%cloud(6,lay)) ** 0.06_jprb) * &
                                            EXP(0.013_jprb * ztempc)
                aux%ice_dg(lay,iprof)     = aux%ice_dg_ref(lay,iprof)
              ELSE IF (profiles(iprof)%icede_param == 4) THEN
                ! Scheme by McFarquhar et al. (2003)
                zmcfarq                    = profiles_int(iprof)%cloud(6,lay)
                aux%fac1_ice_dg(lay,iprof) = &
                  10.0_jprb ** (amcfarq + (bmcfarq * LOG10(zmcfarq)) + (cmcfarq * LOG10(zmcfarq) * LOG10(zmcfarq)))
                aux%ice_dg_ref(lay,iprof)  = 2.0_jprb * aux%fac1_ice_dg(lay,iprof)
                aux%ice_dg(lay,iprof)      = aux%ice_dg_ref(lay,iprof)
              ENDIF
            ENDIF

            IF (ABS(aux%ice_dg_ref(lay,iprof)) > 0._jprb) THEN
              IF (aux%ice_dg_ref(lay,iprof) < dgmin_baum) aux%ice_dg(lay,iprof) = dgmin_baum
              IF (aux%ice_dg_ref(lay,iprof) > dgmax_baum) aux%ice_dg(lay,iprof) = dgmax_baum
            ENDIF
          ENDIF
        ENDDO
      ENDIF ! Baum/SSEC
    ENDDO ! iprof
  ENDIF ! addclouds

  !----------------------------------------------------------------------------
  ! If using pre-defined aerosol types, do relative humidity calculation
  !----------------------------------------------------------------------------
  IF (opts%rt_ir%addaerosl .AND. .NOT. opts%rt_ir%user_aer_opt_param) THEN
    CALL calc_rel_hum(profiles, profiles_int, aux)
  ENDIF

  IF (LHOOK) CALL DR_HOOK('RTTOV_PROFAUX_CLDAER', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_profaux_cldaer
