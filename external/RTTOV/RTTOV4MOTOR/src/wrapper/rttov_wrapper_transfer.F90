! Description:
!> @file
!!   Defines subroutines which transfer data between the wrapper
!!   input/output arrays and the RTTOV structures.
!
!> @brief
!!   Defines subroutines which transfer data between the wrapper
!!   input/output arrays and the RTTOV structures.
!!
!! @details
!!   Defines subroutines which transfer data between the wrapper
!!   input/output arrays and the RTTOV structures:
!!
!!   rttov_copy_to_profiles         - copy from arrays to profiles
!!   rttov_copy_to_opt_param        - copy from arrays to cld/aer optical parameter structures
!!   rttov_scatt_copy_to_profiles   - copy from arrays to profiles for RTTOV-SCATT
!!   rttov_copy_from_profiles_k     - copy from profiles_k to arrays
!!   rttov_copy_from_cld_profiles_k - copy from cld_profiles_k to arrays
!!   rttov_copy_to_radiance_k       - copy from arrays to radiance_k
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
MODULE rttov_wrapper_transfer

#include "throw.h"

IMPLICIT NONE

PRIVATE
PUBLIC :: rttov_copy_to_profiles,         &
          rttov_copy_to_opt_param,        &
          rttov_scatt_copy_to_profiles,   &
          rttov_copy_from_profiles_k,     &
          rttov_copy_from_cld_profiles_k, &
          rttov_copy_to_radiance_k

CONTAINS

!> Copy profile data from input arrays to RTTOV profile structure. The arrays contain
!! all input profile data, but profiles may be smaller than this so only a subset of
!! data are copied starting at profile number iprof1.
!! @param         rth               instrument structure
!! @param[in,out] profiles          array of profile structures to populate
!! @param[in]     iprof1            starting profile number in array data
!! @param[in]     datetimes         profile dates/times: year, month, day, hour, min, sec
!! @param[in]     angles            profile angles: zenang, aziang, sunzenang, sunaziang
!! @param[in]     surfgeom          profile surface geometry: lat, lon, elevation
!! @param[in]     surftype          profile surface type: surftype, watertype
!! @param[in]     skin              profile skin data: skin T, salinity, snow_frac, foam_frac, fastem_coefsx5
!! @param[in]     s2m               profile s2m data: 2m p, 2m t, 2m q, 10m wind u, v, wind-fetch
!! @param[in]     p                 pressure profiles
!! @param[in]     t                 temperature profiles
!! @param[in]     gas_units         profile gas units
!! @param[in]     gas_id            list of gas IDs specified in gases argument
!! @param[in]     gases             gas, aerosol and cloud profiles
!! @param[in]     mmr_cldaer        profile mmr_cldaer flag; optional
!! @param[in]     simplecloud       profile simple cloud data: ctp, cfraction; optional
!! @param[in]     clwscheme         profile clw scheme: clw_scheme, clwde_param; optional
!! @param[in]     icecloud          profile ice scheme: ice_scheme, icede_param; optional
!! @param[in]     zeeman            profile B-field data: Be, cosbk; optional
SUBROUTINE rttov_copy_to_profiles( &
    rth,                      &
    profiles,                 &
    iprof1,                   &
    datetimes,                &
    angles, surfgeom,         &
    surftype, skin, s2m,      &
    p, t,                     &
    gas_units, gas_id, gases, &
    mmr_cldaer,               &
    simplecloud,              &
    clwscheme, icecloud, zeeman)

  USE parkind1, ONLY : jpim, jprb, jplm
  USE rttov_types, ONLY : rttov_profile
  USE rttov_const, ONLY : nwcl_max, naer_opac, naer_cams, aer_id_opac, aer_id_cams
  USE rttov_wrapper_handle

  IMPLICIT NONE

  TYPE(rttovwrapperhandle_type), POINTER              :: rth
  TYPE(rttov_profile),           INTENT(INOUT)        :: profiles(:)
  INTEGER(jpim),                 INTENT(IN)           :: iprof1
  INTEGER(jpim),                 INTENT(IN)           :: datetimes(:,:)   !(6,nprofiles)
  REAL(jprb),                    INTENT(IN)           :: angles(:,:)      !(4,nprofiles)
  REAL(jprb),                    INTENT(IN)           :: surfgeom(:,:)    !(3,nprofiles)
  INTEGER(jpim),                 INTENT(IN)           :: surftype(:,:)    !(2,nprofiles)
  REAL(jprb),                    INTENT(IN)           :: skin(:,:)        !(9,nprofiles)
  REAL(jprb),                    INTENT(IN)           :: s2m(:,:)         !(6,nprofiles)
  REAL(jprb),                    INTENT(IN)           :: p(:,:)           !(nlevels,nprofiles)
  REAL(jprb),                    INTENT(IN)           :: t(:,:)           !(nlevels,nprofiles)
  INTEGER(jpim),                 INTENT(IN)           :: gas_units
  INTEGER(jpim),                 INTENT(IN)           :: gas_id(:)        !(ngases)
  REAL(jprb),                    INTENT(IN)           :: gases(:,:,:)     !(nlevels,nprofiles,ngases)
  LOGICAL(jplm),                 INTENT(IN), OPTIONAL :: mmr_cldaer
  REAL(jprb),                    INTENT(IN), OPTIONAL :: simplecloud(:,:) !(2,nprofiles)
  INTEGER(jpim),                 INTENT(IN), OPTIONAL :: clwscheme(:,:)   !(2,nprofiles)
  INTEGER(jpim),                 INTENT(IN), OPTIONAL :: icecloud(:,:)    !(2,nprofiles)
  REAL(jprb),                    INTENT(IN), OPTIONAL :: zeeman(:,:)      !(2,nprofiles)

  INTEGER(jpim) :: i, j, k, g, nprofiles, ngases
!------------------------------------------------------------------------------

  nprofiles = SIZE(profiles)
  ngases = SIZE(gas_id)

  DO i = 1, nprofiles
    j = iprof1 + i - 1

    profiles(i)%p(:) = p(:,j)
    profiles(i)%t(:) = t(:,j)

    profiles(i)%gas_units = gas_units
    IF (PRESENT(mmr_cldaer)) profiles(i)%mmr_cldaer = mmr_cldaer

    DO g = 1, ngases
      IF (gas_id(g) == gas_id_q) THEN
        profiles(i)%q(:) = gases(:,j,g)
      ELSEIF (gas_id(g) == gas_id_o3 .AND. rth%opts%rt_all%ozone_data) THEN
        profiles(i)%o3(:) = gases(:,j,g)
      ELSEIF (gas_id(g) == gas_id_co2 .AND. rth%opts%rt_all%co2_data) THEN
        profiles(i)%co2(:) = gases(:,j,g)
      ELSEIF (gas_id(g) == gas_id_n2o .AND. rth%opts%rt_all%n2o_data) THEN
        profiles(i)%n2o(:) = gases(:,j,g)
      ELSEIF (gas_id(g) == gas_id_co .AND. rth%opts%rt_all%co_data) THEN
        profiles(i)%co(:) = gases(:,j,g)
      ELSEIF (gas_id(g) == gas_id_ch4 .AND. rth%opts%rt_all%ch4_data) THEN
        profiles(i)%ch4(:) = gases(:,j,g)
      ELSEIF (gas_id(g) == gas_id_so2 .AND. rth%opts%rt_all%so2_data) THEN
        profiles(i)%so2(:) = gases(:,j,g)
      ELSEIF (gas_id(g) == gas_id_clw .AND. rth%opts%rt_mw%clw_data) THEN
        profiles(i)%clw(:) = gases(:,j,g)
      ENDIF

      IF (rth%opts%rt_ir%addclouds) THEN
        IF (.NOT. rth%opts%rt_ir%user_cld_opt_param) THEN
          DO k = 1, nwcl_max
            IF (gas_id(g) == gas_id_lwc(k)) THEN
              profiles(i)%cloud(k,:) = gases(1:profiles(i)%nlayers,j,g)
            ENDIF
          ENDDO
          IF (gas_id(g) == gas_id_iwc) THEN
            profiles(i)%cloud(nwcl_max+1,:) = gases(1:profiles(i)%nlayers,j,g)
          ELSEIF (gas_id(g) == gas_id_icede) THEN
            profiles(i)%icede(:) = gases(1:profiles(i)%nlayers,j,g)
          ELSEIF (gas_id(g) == gas_id_clwde) THEN
            profiles(i)%clwde(:) = gases(1:profiles(i)%nlayers,j,g)
          ENDIF
        ENDIF
        IF (gas_id(g) == gas_id_cfrac) THEN
          profiles(i)%cfrac(:) = gases(1:profiles(i)%nlayers,j,g)
        ENDIF
      ENDIF

      IF (rth%opts%rt_ir%addaerosl .AND. .NOT. rth%opts%rt_ir%user_aer_opt_param) THEN
        DO k = 1, SIZE(profiles(i)%aerosols, DIM=1)
          IF (rth%coefs%coef_scatt%optp_aer%id == aer_id_opac .AND. &
              gas_id(g) >= gas_id_aer_opac(1) .AND. gas_id(g) <= gas_id_aer_opac(naer_opac)) THEN
            IF (gas_id(g) == gas_id_aer_opac(k)) THEN
              profiles(i)%aerosols(k,:) = gases(1:profiles(i)%nlayers,j,g)
            ENDIF
          ELSEIF (rth%coefs%coef_scatt%optp_aer%id == aer_id_cams .AND. &
                  gas_id(g) >= gas_id_aer_cams(1) .AND. gas_id(g) <= gas_id_aer_cams(naer_cams)) THEN
            IF (gas_id(g) == gas_id_aer_cams(k)) THEN
              profiles(i)%aerosols(k,:) = gases(1:profiles(i)%nlayers,j,g)
            ENDIF
          ELSEIF (gas_id(g) >= gas_id_aer_user_min .AND. gas_id(g) <= gas_id_aer_user_max) THEN
            IF (k == gas_id(g) - gas_id_aer_user_min + 1) THEN
              profiles(i)%aerosols(k,:) = gases(1:profiles(i)%nlayers,j,g)
            ENDIF
          ENDIF
        ENDDO
      ENDIF
    ENDDO

    profiles(i)%s2m%p              = s2m(1,j)
    profiles(i)%s2m%t              = s2m(2,j)
    profiles(i)%s2m%q              = s2m(3,j)
    profiles(i)%s2m%u              = s2m(4,j)
    profiles(i)%s2m%v              = s2m(5,j)
    profiles(i)%s2m%wfetc          = s2m(6,j)

    profiles(i)%skin%surftype      = surftype(1,j)
    profiles(i)%skin%watertype     = surftype(2,j)
    profiles(i)%skin%t             = skin(1,j)
    profiles(i)%skin%salinity      = skin(2,j)
    profiles(i)%skin%snow_fraction = skin(3,j)
    profiles(i)%skin%foam_fraction = skin(4,j)
    profiles(i)%skin%fastem(:)     = skin(5:9,j)

    profiles(i)%zenangle           = angles(1,j)
    profiles(i)%azangle            = angles(2,j)
    profiles(i)%sunzenangle        = angles(3,j)
    profiles(i)%sunazangle         = angles(4,j)

    profiles(i)%latitude           = surfgeom(1,j)
    profiles(i)%longitude          = surfgeom(2,j)
    profiles(i)%elevation          = surfgeom(3,j)

    profiles(i)%date(:)            = datetimes(1:3,j)
    profiles(i)%time(:)            = datetimes(4:6,j)

    IF (PRESENT(simplecloud)) THEN
      profiles(i)%ctp                = simplecloud(1,j)
      profiles(i)%cfraction          = simplecloud(2,j)
    ELSE
      profiles(i)%ctp                = 500._jprb
      profiles(i)%cfraction          = 0._jprb
    ENDIF

    IF (PRESENT(clwscheme)) THEN
      profiles(i)%clw_scheme         = clwscheme(1,j)
      profiles(i)%clwde_param        = clwscheme(2,j)
    ELSE
      profiles(i)%clw_scheme         = 1
      profiles(i)%clwde_param        = 1
    ENDIF

    IF (PRESENT(icecloud)) THEN
      profiles(i)%ice_scheme         = icecloud(1,j)
      profiles(i)%icede_param        = icecloud(2,j)
    ELSE
      profiles(i)%ice_scheme         = 1
      profiles(i)%icede_param        = 2 ! recommended parameterisation
    ENDIF

    IF (PRESENT(zeeman)) THEN
      profiles(i)%be                 = zeeman(1,j)
      profiles(i)%cosbk              = zeeman(2,j)
    ELSE
      profiles(i)%be                 = 0._jprb
      profiles(i)%cosbk              = 0._jprb
    ENDIF

  ENDDO

END SUBROUTINE rttov_copy_to_profiles


!> Copy aerosol and/or cloud optical property data from input arrays to rttov_opt_param
!! structures. The arrays contain all input profile data, but RTTOV may be called on a
!! subset of this so data are copied starting at profile number iprof1.
!! @param         rth               instrument structure
!! @param[in,out] aer_opt_param     aerosol optical property structure to populate
!! @param[in,out] cld_opt_param     cloud optical property structure to populate
!! @param[in]     iprof1            starting profile number in array data
!! @param[in]     totnprofiles      total number of profiles in input data
!! @param[in]     nprofiles         number of profiles to copy
!! @param[in]     nchannels         number of channels being simulated per profile
!! @param[in]     nlayers           number of layers in profiles
!! @param[in]     aer_nmom          number of aerosol Legendre moments (excluding zeroth moment)
!! @param[in]     aer_nphangle      size of aerosol phase function arrays
!! @param[in]     aer_phangle       aerosol phase function angle grid
!! @param[in]     aer_asb           aerosol abs, sca and bpr values
!! @param[in]     aer_legcoef       aerosol Legendre coefficients
!! @param[in]     aer_pha           aerosol phase functions
!! @param[in]     cld_nmom          number of cloud Legendre moments (excluding zeroth moment)
!! @param[in]     cld_nphangle      size of cloud phase function arrays
!! @param[in]     cld_phangle       cloud phase function angle grid
!! @param[in]     cld_asb           cloud abs, sca and bpr values
!! @param[in]     cld_legcoef       cloud Legendre coefficients
!! @param[in]     cld_pha           cloud phase functions
SUBROUTINE rttov_copy_to_opt_param( &
    rth,                           &
    aer_opt_param, cld_opt_param,  &
    iprof1, totnprofiles,          &
    nprofiles, nchannels, nlayers, &
    aer_nmom, aer_nphangle,        &
    aer_phangle,                   &
    aer_asb, aer_legcoef, aer_pha, &
    cld_nmom, cld_nphangle,        &
    cld_phangle,                   &
    cld_asb, cld_legcoef, cld_pha)

  USE parkind1, ONLY : jpim, jprb
  USE rttov_types, ONLY : rttov_opt_param
  USE rttov_const, ONLY : vis_scatt_dom, ir_scatt_dom
  USE rttov_wrapper_handle

  IMPLICIT NONE

  TYPE(rttovwrapperhandle_type), POINTER       :: rth
  TYPE(rttov_opt_param),         INTENT(INOUT) :: aer_opt_param
  TYPE(rttov_opt_param),         INTENT(INOUT) :: cld_opt_param
  INTEGER(jpim),                 INTENT(IN)    :: iprof1
  INTEGER(jpim),                 INTENT(IN)    :: totnprofiles
  INTEGER(jpim),                 INTENT(IN)    :: nprofiles
  INTEGER(jpim),                 INTENT(IN)    :: nchannels
  INTEGER(jpim),                 INTENT(IN)    :: nlayers
  INTEGER(jpim),                 INTENT(IN)    :: aer_nmom
  INTEGER(jpim),                 INTENT(IN)    :: aer_nphangle
  REAL(jprb),                    INTENT(IN)    :: aer_phangle(aer_nphangle)
  REAL(jprb),                    INTENT(IN)    :: aer_asb(nlayers,nchannels,totnprofiles,3)
  REAL(jprb),                    INTENT(IN)    :: aer_legcoef(aer_nmom+1,nlayers,nchannels,totnprofiles)
  REAL(jprb),                    INTENT(IN)    :: aer_pha(aer_nphangle,nlayers,nchannels,totnprofiles)
  INTEGER(jpim),                 INTENT(IN)    :: cld_nmom
  INTEGER(jpim),                 INTENT(IN)    :: cld_nphangle
  REAL(jprb),                    INTENT(IN)    :: cld_phangle(cld_nphangle)
  REAL(jprb),                    INTENT(IN)    :: cld_asb(nlayers,nchannels,totnprofiles,3)
  REAL(jprb),                    INTENT(IN)    :: cld_legcoef(cld_nmom+1,nlayers,nchannels,totnprofiles)
  REAL(jprb),                    INTENT(IN)    :: cld_pha(cld_nphangle,nlayers,nchannels,totnprofiles)

  INTEGER(jpim) :: i, j, lo, hi
!------------------------------------------------------------------------------

  DO i = 1, nprofiles
    j = iprof1 + i - 1

    lo = (i - 1) * nchannels + 1
    hi = lo + nchannels - 1

    IF (rth%opts%rt_ir%addaerosl) THEN
      aer_opt_param%abs(:,lo:hi) = aer_asb(:,:,j,1)
      aer_opt_param%sca(:,lo:hi) = aer_asb(:,:,j,2)
      aer_opt_param%bpr(:,lo:hi) = aer_asb(:,:,j,3)

      IF ((rth%opts%rt_ir%vis_scatt_model == vis_scatt_dom .AND. rth%opts%rt_ir%addsolar) .OR. &
           rth%opts%rt_ir%ir_scatt_model == ir_scatt_dom) THEN
        aer_opt_param%legcoef(:,:,lo:hi) = aer_legcoef(:,:,:,j)
      ENDIF

      IF (rth%opts%rt_ir%addsolar) THEN
        aer_opt_param%phangle = aer_phangle
        aer_opt_param%pha(:,:,lo:hi) = aer_pha(:,:,:,j)
      ENDIF
    ENDIF

    IF (rth%opts%rt_ir%addclouds) THEN
      cld_opt_param%abs(:,lo:hi) = cld_asb(:,:,j,1)
      cld_opt_param%sca(:,lo:hi) = cld_asb(:,:,j,2)
      cld_opt_param%bpr(:,lo:hi) = cld_asb(:,:,j,3)

      IF ((rth%opts%rt_ir%vis_scatt_model == vis_scatt_dom .AND. rth%opts%rt_ir%addsolar) .OR. &
           rth%opts%rt_ir%ir_scatt_model == ir_scatt_dom) THEN
        cld_opt_param%legcoef(:,:,lo:hi) = cld_legcoef(:,:,:,j)
      ENDIF

      IF (rth%opts%rt_ir%addsolar) THEN
        cld_opt_param%phangle = cld_phangle
        cld_opt_param%pha(:,:,lo:hi) = cld_pha(:,:,:,j)
      ENDIF
    ENDIF

  ENDDO

END SUBROUTINE rttov_copy_to_opt_param


!> Copy profile data from input arrays to RTTOV profile and cld_profile structures for
!! RTTOV-SCATT simulations. The arrays contain all input profile data, but profiles may
!! be smaller than this so only a subset of data are copied starting at profile number
!! iprof1.
!! @param         rth               instrument structure
!! @param[in,out] profiles          array of profile structures to populate
!! @param[in,out] cld_profiles      array of RTTOV-SCATT cld_profile structures to populate
!! @param[in]     iprof1            starting profile number in array data
!! @param[in]     datetimes         profile dates/times: year, month, day, hour, min, sec
!! @param[in]     angles            profile angles: zenang, aziang
!! @param[in]     surfgeom          profile surface geometry: lat, lon, elevation
!! @param[in]     surftype          profile surface type
!! @param[in]     skin              profile skin data: skin T, salinity, foam_frac, fastem_coefsx5
!! @param[in]     s2m               profile s2m data: 2m p, 2m t, 2m q, 10m wind u, v
!! @param[in]     zeeman            profile B-field data: Be, cosbk
!! @param[in]     p                 pressure profiles
!! @param[in]     t                 temperature profiles
!! @param[in]     gas_units         profile gas units
!! @param[in]     gas_id            list of gas IDs specified in gases argument
!! @param[in]     gases             gas, cloud and hydrometeor profiles
!! @param[in]     ph                pressure half-level profiles
!! @param[in]     cfrac             user cloud fraction
SUBROUTINE rttov_scatt_copy_to_profiles( &
    rth,                      &
    profiles,                 &
    cld_profiles,             &
    iprof1,                   &
    datetimes,                &
    angles, surfgeom,         &
    surftype, skin, s2m,      &
    zeeman,                   &
    p, t, gas_units,          &
    gas_id, gases, ph, cfrac)

  USE parkind1, ONLY : jpim, jprb
  USE rttov_types, ONLY : rttov_profile, rttov_profile_cloud
  USE rttov_wrapper_handle

  IMPLICIT NONE

  TYPE(rttovwrapperhandle_type), POINTER       :: rth
  TYPE(rttov_profile),           INTENT(INOUT) :: profiles(:)
  TYPE(rttov_profile_cloud),     INTENT(INOUT) :: cld_profiles(:)
  INTEGER(jpim),                 INTENT(IN)    :: iprof1
  INTEGER(jpim),                 INTENT(IN)    :: datetimes(:,:)   !(6,nprofiles)
  REAL(jprb),                    INTENT(IN)    :: angles(:,:)      !(2,nprofiles)
  REAL(jprb),                    INTENT(IN)    :: surfgeom(:,:)    !(3,nprofiles)
  INTEGER(jpim),                 INTENT(IN)    :: surftype(:)      !(nprofiles)
  REAL(jprb),                    INTENT(IN)    :: skin(:,:)        !(8,nprofiles)
  REAL(jprb),                    INTENT(IN)    :: s2m(:,:)         !(5,nprofiles)
  REAL(jprb),                    INTENT(IN)    :: zeeman(:,:)      !(2,nprofiles)
  REAL(jprb),                    INTENT(IN)    :: p(:,:)           !(nlevels,nprofiles)
  REAL(jprb),                    INTENT(IN)    :: t(:,:)           !(nlevels,nprofiles)
  INTEGER(jpim),                 INTENT(IN)    :: gas_units
  INTEGER(jpim),                 INTENT(IN)    :: gas_id(:)        !(ngases)
  REAL(jprb),                    INTENT(IN)    :: gases(:,:,:)     !(nlevels,nprofiles,ngases)
  REAL(jprb),                    INTENT(IN)    :: ph(:,:)          !(nlevels+1,nprofiles)
  REAL(jprb),                    INTENT(IN)    :: cfrac(:)         !(nprofiles)

  INTEGER(jpim) :: i, j, k, g, nprofiles, ngases
!------------------------------------------------------------------------------

  nprofiles = SIZE(profiles)
  ngases = SIZE(gas_id)

  DO i = 1, nprofiles
    j = iprof1 + i - 1

    profiles(i)%p(:) = p(:,j)
    profiles(i)%t(:) = t(:,j)

    profiles(i)%gas_units = gas_units

    DO g = 1, ngases
      IF (gas_id(g) == gas_id_q) THEN
        profiles(i)%q(:) = gases(:,j,g)
      ELSEIF (gas_id(g) == gas_id_o3 .AND. rth%opts%rt_all%ozone_data) THEN
        profiles(i)%o3(:) = gases(:,j,g)
      ELSEIF (gas_id(g) == gas_id_scatt_hydro_frac) THEN
        cld_profiles(i)%hydro_frac(:,1) = gases(:,j,g)
      ELSEIF (gas_id(g) == gas_id_scatt_rain) THEN
        cld_profiles(i)%hydro(:,1) = gases(:,j,g)
      ELSEIF (gas_id(g) == gas_id_scatt_snow) THEN
        cld_profiles(i)%hydro(:,2) = gases(:,j,g)
      ELSEIF (gas_id(g) == gas_id_scatt_graupel) THEN
        cld_profiles(i)%hydro(:,3) = gases(:,j,g)
      ELSEIF (gas_id(g) == gas_id_scatt_clw) THEN
        cld_profiles(i)%hydro(:,4) = gases(:,j,g)
      ELSEIF (gas_id(g) == gas_id_scatt_ciw) THEN
        cld_profiles(i)%hydro(:,5) = gases(:,j,g)
      ENDIF

      DO k = 1, SIZE(cld_profiles(i)%hydro, DIM=2)
        IF (gas_id(g) >= gas_id_scatt_hydro_min .AND. &
            gas_id(g) <= gas_id_scatt_hydro_max) THEN
          IF (k == gas_id(g) - gas_id_scatt_hydro_min + 1) THEN
            cld_profiles(i)%hydro(:,k) = gases(:,j,g)
          ENDIF
        ENDIF
      ENDDO

      DO k = 1, SIZE(cld_profiles(i)%hydro_frac, DIM=2)
        IF (gas_id(g) >= gas_id_scatt_hydro_frac_min .AND. &
            gas_id(g) <= gas_id_scatt_hydro_frac_max) THEN
          IF (k == gas_id(g) - gas_id_scatt_hydro_frac_min + 1) THEN
            cld_profiles(i)%hydro_frac(:,k) = gases(:,j,g)
          ENDIF
        ENDIF
      ENDDO
    ENDDO

    cld_profiles(i)%ph(:) = ph(:,j)
    IF (rth%opts_scatt%lusercfrac) cld_profiles(i)%cfrac = cfrac(j)

    profiles(i)%s2m%p              = s2m(1,j)
    profiles(i)%s2m%t              = s2m(2,j)
    profiles(i)%s2m%q              = s2m(3,j)
    profiles(i)%s2m%u              = s2m(4,j)
    profiles(i)%s2m%v              = s2m(5,j)

    profiles(i)%skin%surftype      = surftype(j)
    profiles(i)%skin%watertype     = 1 ! Not used by FASTEM/MW simulations
    profiles(i)%skin%t             = skin(1,j)
    profiles(i)%skin%salinity      = skin(2,j)
    profiles(i)%skin%foam_fraction = skin(3,j)
    profiles(i)%skin%fastem(:)     = skin(4:8,j)

    profiles(i)%zenangle           = angles(1,j)
    profiles(i)%azangle            = angles(2,j)

    profiles(i)%latitude           = surfgeom(1,j)
    profiles(i)%longitude          = surfgeom(2,j)
    profiles(i)%elevation          = surfgeom(3,j)

    profiles(i)%date(:)            = datetimes(1:3,j)
    profiles(i)%time(:)            = datetimes(4:6,j)

    profiles(i)%be                 = zeeman(1,j)
    profiles(i)%cosbk              = zeeman(2,j)
  ENDDO

END SUBROUTINE rttov_scatt_copy_to_profiles


!> Copy Jacobian data from RTTOV profiles_k structures to output arrays.
!! The data are copied to profile positions starting at iprof1 in the output arrays.
!! @param         rth               instrument structure
!! @param[in]     profiles_k        input array of profile Jacobians
!! @param[in]     iprof1            starting profile number in array data
!! @param[in,out] skin_k            output skin Jacobians
!! @param[in,out] s2m_k             output s2m Jacobians
!! @param[in,out] t_k               output temperature Jacobians
!! @param[in]     gas_id            list of gas IDs specified in gases_k argument
!! @param[in,out] gases_k           output gas, aerosol and cloud Jacobians
!! @param[in,out] simplecloud_k     output simple cloud Jacobians, optional
!! @param[in,out] p_k               output pressure Jacobians, optional
SUBROUTINE rttov_copy_from_profiles_k( &
    rth,             &
    profiles_k,      &
    iprof1,          &
    skin_k, s2m_k,   &
    t_k,             &
    gas_id, gases_k, &
    simplecloud_k,   &
    p_k)

  USE parkind1, ONLY : jpim, jprb
  USE rttov_types, ONLY : rttov_profile
  USE rttov_const, ONLY : nwcl_max, naer_opac, naer_cams, aer_id_opac, aer_id_cams
  USE rttov_wrapper_handle

  IMPLICIT NONE

  TYPE(rttovwrapperhandle_type), POINTER       :: rth
  TYPE(rttov_profile),           INTENT(IN)    :: profiles_k(:)
  INTEGER(jpim),                 INTENT(IN)    :: iprof1
  REAL(jprb),                    INTENT(INOUT) :: skin_k(:,:,:)        !(9/8,nchannels,nprofiles)
  REAL(jprb),                    INTENT(INOUT) :: s2m_k(:,:,:)         !(6,nchannels,nprofiles)
  REAL(jprb),                    INTENT(INOUT) :: t_k(:,:,:)           !(nlevels,nchannels,nprofiles)
  INTEGER(jpim),                 INTENT(IN)    :: gas_id(:)            !(ngases)
  REAL(jprb),                    INTENT(INOUT) :: gases_k(:,:,:,:)     !(nlevels,nchannels,nprofiles,ngases)
  REAL(jprb),          OPTIONAL, INTENT(INOUT) :: simplecloud_k(:,:,:) !(2,nchannels,nprofiles)
  REAL(jprb),          OPTIONAL, INTENT(INOUT) :: p_k(:,:,:)           !(nlevels,nchannels,nprofiles)

  INTEGER(jpim) :: i, j, k, t, g, iprof, nchannels, nprofiles, ngases, nlayers
!------------------------------------------------------------------------------

  nchannels = SIZE(t_k, 2)
  nprofiles = SIZE(profiles_k) / nchannels
  ngases = SIZE(gas_id)
  nlayers = SIZE(t_k, 1) - 1

  DO i = 1, nprofiles
    iprof = iprof1 + i - 1

    DO j = 1, nchannels
      k = (i - 1) * nchannels + j

      IF (PRESENT(p_k)) p_k(:,j,iprof) = profiles_k(k)%p(:)
      t_k(:,j,iprof) = profiles_k(k)%t(:)

      DO g = 1, ngases
        IF (gas_id(g) == gas_id_q) THEN
          gases_k(:,j,iprof,g) = profiles_k(k)%q(:)
        ELSEIF (gas_id(g) == gas_id_o3 .AND. rth%opts%rt_all%ozone_data) THEN
          gases_k(:,j,iprof,g) = profiles_k(k)%o3(:)
        ELSEIF (gas_id(g) == gas_id_co2 .AND. rth%opts%rt_all%co2_data) THEN
          gases_k(:,j,iprof,g) = profiles_k(k)%co2(:)
        ELSEIF (gas_id(g) == gas_id_n2o .AND. rth%opts%rt_all%n2o_data) THEN
          gases_k(:,j,iprof,g) = profiles_k(k)%n2o(:)
        ELSEIF (gas_id(g) == gas_id_co .AND. rth%opts%rt_all%co_data) THEN
          gases_k(:,j,iprof,g) = profiles_k(k)%co(:)
        ELSEIF (gas_id(g) == gas_id_ch4 .AND. rth%opts%rt_all%ch4_data) THEN
          gases_k(:,j,iprof,g) = profiles_k(k)%ch4(:)
        ELSEIF (gas_id(g) == gas_id_so2 .AND. rth%opts%rt_all%so2_data) THEN
          gases_k(:,j,iprof,g) = profiles_k(k)%so2(:)
        ELSEIF (gas_id(g) == gas_id_clw .AND. rth%opts%rt_mw%clw_data) THEN
          gases_k(:,j,iprof,g) = profiles_k(k)%clw(:)
        ENDIF

        IF (rth%opts%rt_ir%addclouds) THEN
          IF (.NOT. rth%opts%rt_ir%user_cld_opt_param) THEN
            DO t = 1, nwcl_max
              IF (gas_id(g) == gas_id_lwc(t)) THEN
                gases_k(1:nlayers,j,iprof,g) = profiles_k(k)%cloud(t,:)
              ENDIF
            ENDDO
            IF (gas_id(g) == gas_id_iwc) THEN
              gases_k(1:nlayers,j,iprof,g) = profiles_k(k)%cloud(nwcl_max+1,:)
            ELSEIF (gas_id(g) == gas_id_icede) THEN
              gases_k(1:nlayers,j,iprof,g) = profiles_k(k)%icede(:)
            ELSEIF (gas_id(g) == gas_id_clwde) THEN
              gases_k(1:nlayers,j,iprof,g) = profiles_k(k)%clwde(:)
            ENDIF
          ENDIF
          IF (gas_id(g) == gas_id_cfrac) THEN
            gases_k(1:nlayers,j,iprof,g) = profiles_k(k)%cfrac(:)
          ENDIF
        ENDIF

        IF (rth%opts%rt_ir%addaerosl .AND. .NOT. rth%opts%rt_ir%user_aer_opt_param) THEN
          DO t = 1, SIZE(profiles_k(k)%aerosols, DIM=1)
            IF (rth%coefs%coef_scatt%optp_aer%id == aer_id_opac .AND. &
                gas_id(g) >= gas_id_aer_opac(1) .AND. gas_id(g) <= gas_id_aer_opac(naer_opac)) THEN
              IF (gas_id(g) == gas_id_aer_opac(t)) THEN
                gases_k(1:nlayers,j,iprof,g) = profiles_k(k)%aerosols(t,:)
              ENDIF
            ELSEIF (rth%coefs%coef_scatt%optp_aer%id == aer_id_cams .AND. &
                    gas_id(g) >= gas_id_aer_cams(1) .AND. gas_id(g) <= gas_id_aer_cams(naer_cams)) THEN
              IF (gas_id(g) == gas_id_aer_cams(t)) THEN
                gases_k(1:nlayers,j,iprof,g) = profiles_k(k)%aerosols(t,:)
              ENDIF
            ELSEIF (gas_id(g) >= gas_id_aer_user_min .AND. gas_id(g) <= gas_id_aer_user_max) THEN
              IF (t == gas_id(g) - gas_id_aer_user_min + 1) THEN
                gases_k(1:nlayers,j,iprof,g) = profiles_k(k)%aerosols(t,:)
              ENDIF
            ENDIF
          ENDDO
        ENDIF
      ENDDO

      s2m_k(1,j,iprof) = profiles_k(k)%s2m%p
      s2m_k(2,j,iprof) = profiles_k(k)%s2m%t
      s2m_k(3,j,iprof) = profiles_k(k)%s2m%q
      s2m_k(4,j,iprof) = profiles_k(k)%s2m%u
      s2m_k(5,j,iprof) = profiles_k(k)%s2m%v
      IF (SIZE(s2m_k,1) == 6) s2m_k(6,j,iprof) = profiles_k(k)%s2m%wfetc

      skin_k(1,j,iprof) = profiles_k(k)%skin%t
      skin_k(2,j,iprof) = profiles_k(k)%skin%salinity
      IF (SIZE(skin_k,1) == 9) THEN
        skin_k(3,j,iprof) = profiles_k(k)%skin%snow_fraction
        skin_k(4,j,iprof) = profiles_k(k)%skin%foam_fraction
        skin_k(5:9,j,iprof) = profiles_k(k)%skin%fastem(:)
      ELSE
        skin_k(3,j,iprof) = profiles_k(k)%skin%foam_fraction
        skin_k(4:8,j,iprof) = profiles_k(k)%skin%fastem(:)
      ENDIF

      IF (PRESENT(simplecloud_k)) THEN
        simplecloud_k(1,j,iprof) = profiles_k(k)%ctp
        simplecloud_k(2,j,iprof) = profiles_k(k)%cfraction
      ENDIF

    ENDDO ! channels
  ENDDO ! profiles

END SUBROUTINE rttov_copy_from_profiles_k


!> Copy Jacobian data from RTTOV-SCATT cld_profiles_k structures to output arrays.
!! The data are copied to profile positions starting at iprof1 in the output arrays.
!! @param[in]     cld_profiles_k    input array of cld_profile Jacobians
!! @param[in]     iprof1            starting profile number in array data
!! @param[in]     gas_id            list of gas IDs specified in gases_k argument
!! @param[in,out] gases_k           output gas, cloud and hydrometeor Jacobians
!! @param[in,out] ph_k              output pressure half-level Jacobians
!! @param[in,out] cfrac_k           output user cloud fraction Jacobians
SUBROUTINE rttov_copy_from_cld_profiles_k( &
    cld_profiles_k,  &
    iprof1,          &
    gas_id, gases_k, &
    ph_k, cfrac_k)

  USE parkind1, ONLY : jpim, jprb
  USE rttov_types, ONLY : rttov_profile_cloud
  USE rttov_wrapper_handle

  IMPLICIT NONE

  TYPE(rttov_profile_cloud),     INTENT(IN)    :: cld_profiles_k(:)
  INTEGER(jpim),                 INTENT(IN)    :: iprof1
  INTEGER(jpim),                 INTENT(IN)    :: gas_id(:)            !(ngases)
  REAL(jprb),                    INTENT(INOUT) :: gases_k(:,:,:,:)     !(nlevels,nchannels,nprofiles,ngases)
  REAL(jprb),                    INTENT(INOUT) :: ph_k(:,:,:)          !(nlevels+1,nchannels,nprofiles)
  REAL(jprb),                    INTENT(INOUT) :: cfrac_k(:,:)         !(nchannels,nprofiles)

  INTEGER(jpim) :: i, j, k, t, g, iprof, nchannels, nprofiles, ngases
!------------------------------------------------------------------------------

  nchannels = SIZE(gases_k, 2)
  nprofiles = SIZE(cld_profiles_k) / nchannels
  ngases = SIZE(gas_id)

  DO i = 1, nprofiles
    iprof = iprof1 + i - 1

    DO j = 1, nchannels
      k = (i - 1) * nchannels + j

      ph_k(:,j,iprof) = cld_profiles_k(k)%ph(:)
      cfrac_k(j,iprof) = cld_profiles_k(k)%cfrac

      DO g = 1, ngases
        IF (gas_id(g) == gas_id_scatt_hydro_frac) THEN
          gases_k(:,j,iprof,g) = cld_profiles_k(k)%hydro_frac(:,1)
        ELSEIF (gas_id(g) == gas_id_scatt_rain) THEN
          gases_k(:,j,iprof,g) = cld_profiles_k(k)%hydro(:,1)
        ELSEIF (gas_id(g) == gas_id_scatt_snow) THEN
          gases_k(:,j,iprof,g) = cld_profiles_k(k)%hydro(:,2)
        ELSEIF (gas_id(g) == gas_id_scatt_graupel) THEN
          gases_k(:,j,iprof,g) = cld_profiles_k(k)%hydro(:,3)
        ELSEIF (gas_id(g) == gas_id_scatt_clw) THEN
          gases_k(:,j,iprof,g) = cld_profiles_k(k)%hydro(:,4)
        ELSEIF (gas_id(g) == gas_id_scatt_ciw) THEN
          gases_k(:,j,iprof,g) = cld_profiles_k(k)%hydro(:,5)
        ENDIF

        DO t = 1, SIZE(cld_profiles_k(k)%hydro, DIM=2)
          IF (gas_id(g) >= gas_id_scatt_hydro_min .AND. &
              gas_id(g) <= gas_id_scatt_hydro_max) THEN
            IF (t == gas_id(g) - gas_id_scatt_hydro_min + 1) THEN
              gases_k(:,j,iprof,g) = cld_profiles_k(k)%hydro(:,t)
            ENDIF
          ENDIF
        ENDDO

        DO t = 1, SIZE(cld_profiles_k(k)%hydro_frac, DIM=2)
          IF (gas_id(g) >= gas_id_scatt_hydro_frac_min .AND. &
              gas_id(g) <= gas_id_scatt_hydro_frac_max) THEN
            IF (t == gas_id(g) - gas_id_scatt_hydro_frac_min + 1) THEN
              gases_k(:,j,iprof,g) = cld_profiles_k(k)%hydro_frac(:,t)
            ENDIF
          ENDIF
        ENDDO
      ENDDO

    ENDDO ! channels
  ENDDO ! profiles

END SUBROUTINE rttov_copy_from_cld_profiles_k


!> Copy K perturbation data from input arrays to RTTOV radiance_k structure.
!! The arrays contain all input profile data, but profiles may be smaller than
!! this so only a subset of data are copied starting at profile number iprof1.
!! @param[in,out] radiance_k        radiance perturbation structure to populate
!! @param[in]     iprof1            starting profile number in array data
!! @param[in]     nprofiles         number of profiles to copy
!! @param[in]     bt_k              input BT K perturbations
!! @param[in]     rad_k             input radiance K perturbations, optional
!! @param[in]     zef_k             input radar reflectivity K perturbations, optional
!! @param[in,out] reflectivity_k    radar reflectivity perturbation structure to populate, optional
SUBROUTINE rttov_copy_to_radiance_k( &
    radiance_k,         &
    iprof1,             &
    nprofiles,          &
    bt_k, rad_k, zef_k, &
    reflectivity_k)

  USE parkind1, ONLY : jpim, jprb
  USE rttov_types, ONLY : rttov_radiance, rttov_reflectivity

  IMPLICIT NONE

  TYPE(rttov_radiance),     INTENT(INOUT)           :: radiance_k
  INTEGER(jpim),            INTENT(IN)              :: iprof1
  INTEGER(jpim),            INTENT(IN)              :: nprofiles
  REAL(jprb),               INTENT(IN)              :: bt_k(:,:)      !(nchannels,nprofiles)
  REAL(jprb),               INTENT(IN),    OPTIONAL :: rad_k(:,:)     !(nchannels,nprofiles)
  REAL(jprb),               INTENT(IN),    OPTIONAL :: zef_k(:,:,:)   !(nlevels,nchannels,nprofiles)
  TYPE(rttov_reflectivity), INTENT(INOUT), OPTIONAL :: reflectivity_k

  INTEGER(jpim) :: i, j, lo, hi, nchannels
!------------------------------------------------------------------------------

  nchannels = SIZE(bt_k, 1)

  DO i = 1, nprofiles
    j = iprof1 + i - 1
    lo = (i - 1) * nchannels + 1
    hi = i * nchannels

    ! RTTOV K takes care of using the correct perturbation based on switchrad
    ! and channel type (IR/MW vs VIS/NIR) so copy both bt_k and rad_k into radiance_k
    ! For RTTOV-SCATT passive simulations input perturbations are in BT only

    radiance_k%bt(lo:hi) = bt_k(:,j)
    IF (PRESENT(rad_k)) radiance_k%total(lo:hi) = rad_k(:,j)

    ! For RTTOV-SCATT radar simulations, input perturbations are in radar reflectivity
    IF (PRESENT(zef_k) .AND. PRESENT(reflectivity_k)) THEN
      reflectivity_k%zef(:,lo:hi)  = zef_k(:,:,j)
      reflectivity_k%azef(:,lo:hi) = 0.0_jprb
    ENDIF
  ENDDO

END SUBROUTINE rttov_copy_to_radiance_k

END MODULE rttov_wrapper_transfer
