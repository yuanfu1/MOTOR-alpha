! Description:
!> @file
!!   Allocate/deallocate internal transmission_aux structure.
!
!> @brief
!!   Allocate/deallocate internal transmission_aux structure.
!!
!! @details
!!   This structure holds information related to optical depths
!!   and transmittances.
!!
!!   Some arrays are only required by the direct model instance
!!   of this structure: these are allocated if the optional "direct"
!!   argument is TRUE. If omitted, it is assumed to be FALSE.
!!
!! @param[out]    err               status on exit
!! @param[in]     opts              options to configure the simulations
!! @param[in,out] transmission_aux  transmission_aux structure to (de)allocate
!! @param[in]     nlayers           number of input profile layers (nlevels-1)
!! @param[in]     nchanprof         total number of channels being simulated (SIZE(chanprof))
!! @param[in]     asw               1_jpim => allocate; 0_jpim => deallocate
!! @param[in]     ncolumns          number of cloud columns as calculated by rttov_cloud_overlap
!! @param[in]     init              set .TRUE. to initialise newly allocated structures, optional
!! @param[in]     direct            if FALSE then direct-model-only arrays are not allocated, optional
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
SUBROUTINE rttov_alloc_transmission_aux( &
              err,              &
              opts,             &
              transmission_aux, &
              nlayers,          &
              nchanprof,        &
              asw,              &
              ncolumns,         &
              init,             &
              direct)

!INTF_OFF
#include "throw.h"
!INTF_ON

  USE rttov_types, ONLY : rttov_options, rttov_transmission_aux
  USE parkind1, ONLY : jpim, jplm
!INTF_OFF
  USE rttov_const, ONLY : ir_scatt_dom, vis_scatt_dom, vis_scatt_single, vis_scatt_mfasis
!INTF_ON
  IMPLICIT NONE

  INTEGER(KIND=jpim),           INTENT(OUT)          :: err
  TYPE(rttov_options),          INTENT(IN)           :: opts
  TYPE(rttov_transmission_aux), INTENT(INOUT)        :: transmission_aux
  INTEGER(KIND=jpim),           INTENT(IN)           :: nlayers
  INTEGER(KIND=jpim),           INTENT(IN)           :: nchanprof
  INTEGER(KIND=jpim),           INTENT(IN)           :: asw             ! 1=allocate, 0=deallocate
  INTEGER(KIND=jpim),           INTENT(IN)           :: ncolumns
  LOGICAL(KIND=jplm),           INTENT(IN), OPTIONAL :: init
  LOGICAL(KIND=jplm),           INTENT(IN), OPTIONAL :: direct
!INTF_END
#include "rttov_errorreport.interface"
#include "rttov_init_transmission_aux.interface"

  INTEGER(KIND=jpim) :: nlevels, su
  LOGICAL(KIND=jplm) :: init1, direct1, do_scatt
  LOGICAL(KIND=jplm) :: do_dom_ir, do_dom_vis, do_single_scatt, do_mfasis
!- End of header --------------------------------------------------------

  TRY
  nlevels = nlayers + 1
  init1   = .FALSE.
  IF (PRESENT(init)) init1 = init
  direct1   = .FALSE.
  IF (PRESENT(direct)) direct1 = direct

  IF (opts%rt_ir%addclouds) THEN
    su = 1
  ELSE
    su = 0
  ENDIF

  do_scatt = opts%rt_ir%addclouds .OR. opts%rt_ir%addaerosl
  do_dom_ir = do_scatt .AND. opts%rt_ir%ir_scatt_model == ir_scatt_dom
  do_dom_vis = do_scatt .AND. opts%rt_ir%vis_scatt_model == vis_scatt_dom
  do_single_scatt = do_scatt .AND. opts%rt_ir%vis_scatt_model == vis_scatt_single
  do_mfasis = do_scatt .AND. opts%rt_ir%vis_scatt_model == vis_scatt_mfasis

  IF (asw == 1) THEN
    CALL nullify_struct()

    IF (direct1 .AND. .NOT. do_dom_ir) THEN
      ALLOCATE (transmission_aux%fac1(nlevels, 0:ncolumns, nchanprof), STAT=err)
      THROWM(err.NE.0, "allocation of transmission_aux % fac1")
      ALLOCATE (transmission_aux%surf_fac(0:ncolumns, nchanprof), STAT=err)
      THROWM(err.NE.0, "allocation of transmission_aux % surf_fac")
    ELSE
      NULLIFY(transmission_aux%fac1)
      NULLIFY(transmission_aux%surf_fac)
    ENDIF

    ! Thermal path1
    ALLOCATE (transmission_aux%thermal_path1)

    IF (do_dom_ir) THEN
      ALLOCATE (transmission_aux%thermal_path1%tau_level(nlevels, 0:0, nchanprof), STAT=err)
      THROWM(err.NE.0, "allocation of transmission_aux % thermal_path1 % tau_level")
      ALLOCATE (transmission_aux%thermal_path1%od_singlelayer(0:0, nlayers, nchanprof), STAT=err)
      THROWM(err.NE.0, "allocation of transmission_aux % thermal_path1 % od_singlelayer")
      NULLIFY(transmission_aux%thermal_path1%od_sfrac)
      NULLIFY(transmission_aux%thermal_path1%tau_surf)
      NULLIFY(transmission_aux%thermal_path1%tau_surf_ac)
    ELSE
      ALLOCATE (transmission_aux%thermal_path1%tau_level(nlevels, 0:ncolumns, nchanprof), STAT=err)
      THROWM(err.NE.0, "allocation of transmission_aux % thermal_path1 % tau_level")
      ALLOCATE (transmission_aux%thermal_path1%od_singlelayer(0:su, nlayers, nchanprof), STAT=err)
      THROWM(err.NE.0, "allocation of transmission_aux % thermal_path1 % od_singlelayer")
      ALLOCATE (transmission_aux%thermal_path1%od_sfrac(0:ncolumns, nchanprof), STAT=err)
      THROWM(err.NE.0, "allocation of transmission_aux % thermal_path1 % od_sfrac")
      ALLOCATE (transmission_aux%thermal_path1%tau_surf(0:ncolumns, nchanprof), STAT=err)
      THROWM(err.NE.0, "allocation of transmission_aux % thermal_path1 % tau_surf")
      ALLOCATE (transmission_aux%thermal_path1%tau_surf_ac(0:ncolumns, nchanprof), STAT=err)
      THROWM(err.NE.0, "allocation of transmission_aux % thermal_path1 % tau_surf_ac")
    ENDIF

    IF (direct1 .AND. .NOT. do_dom_ir) THEN
      ALLOCATE (transmission_aux%thermal_path1%fac2(nlevels, 0:ncolumns, nchanprof), STAT=err)
      THROWM(err.NE.0, "allocation of transmission_aux % thermal_path1 % fac2")
      ALLOCATE (transmission_aux%thermal_path1%tau_level_r(nlevels, 0:ncolumns, nchanprof), STAT=err)
      THROWM(err.NE.0, "allocation of transmission_aux % thermal_path1 % tau_level_r")
      ALLOCATE (transmission_aux%thermal_path1%od_singlelayer_r(0:su, nlayers, nchanprof), STAT=err)
      THROWM(err.NE.0, "allocation of transmission_aux % thermal_path1 % od_singlelayer_r")
      ALLOCATE (transmission_aux%thermal_path1%tau_surf_r(0:ncolumns, nchanprof), STAT=err)
      THROWM(err.NE.0, "allocation of transmission_aux % thermal_path1 % tau_surf_r")
      ALLOCATE (transmission_aux%thermal_path1%od_sfrac_r(0:ncolumns, nchanprof), STAT=err)
      THROWM(err.NE.0, "allocation of transmission_aux % thermal_path1 % od_sfrac_r")
    ELSE
      NULLIFY(transmission_aux%thermal_path1%fac2)
      NULLIFY(transmission_aux%thermal_path1%tau_level_r)
      NULLIFY(transmission_aux%thermal_path1%od_singlelayer_r)
      NULLIFY(transmission_aux%thermal_path1%tau_surf_r)
      NULLIFY(transmission_aux%thermal_path1%od_sfrac_r)
    ENDIF

    NULLIFY(transmission_aux%thermal_path1%tau_level_p_r)
    NULLIFY(transmission_aux%thermal_path1%tau_surf_p_r)
    IF (opts%rt_all%do_lambertian .AND. .NOT. do_dom_ir) THEN
      ALLOCATE (transmission_aux%thermal_path1%tau_level_p(nlevels, 0:ncolumns, nchanprof), STAT=err)
      THROWM(err.NE.0, "allocation of transmission_aux % thermal_path1 % tau_level_p")
      ALLOCATE (transmission_aux%thermal_path1%tau_surf_p(0:ncolumns, nchanprof), STAT=err)
      THROWM(err.NE.0, "allocation of transmission_aux % thermal_path1 % tau_surf_p")
      IF (direct1) THEN
        ALLOCATE (transmission_aux%thermal_path1%tau_level_p_r(nlevels, 0:ncolumns, nchanprof), STAT=err)
        THROWM(err.NE.0, "allocation of transmission_aux % thermal_path1 % tau_level_p_r")
        ALLOCATE (transmission_aux%thermal_path1%tau_surf_p_r(0:ncolumns, nchanprof), STAT=err)
        THROWM(err.NE.0, "allocation of transmission_aux % thermal_path1 % tau_surf_p_r")
      ENDIF
    ELSE
      NULLIFY(transmission_aux%thermal_path1%tau_level_p)
      NULLIFY(transmission_aux%thermal_path1%tau_surf_p)
    ENDIF

    IF (opts%rt_ir%addsolar .AND. .NOT. do_mfasis) THEN
      ! Solar path2
      ALLOCATE (transmission_aux%solar_path2)
      ALLOCATE (transmission_aux%solar_path2%tau_level(nlevels, 0:ncolumns, nchanprof), STAT=err)
      THROWM(err.NE.0, "allocation of transmission_aux % solar_path2 % tau_level")
      ALLOCATE (transmission_aux%solar_path2%tau_surf(0:ncolumns, nchanprof), STAT=err)
      THROWM(err.NE.0, "allocation of transmission_aux % solar_path2 % tau_surf")
      ALLOCATE (transmission_aux%solar_path2%tau_surf_ac(0:ncolumns, nchanprof), STAT=err)
      THROWM(err.NE.0, "allocation of transmission_aux % solar_path2 % tau_surf_ac")

      IF (do_single_scatt) THEN
        ALLOCATE (transmission_aux%solar_path2%od_singlelayer(0:su, nlayers, nchanprof), STAT=err)
        THROWM(err.NE.0, "allocation of transmission_aux % solar_path2 % od_singlelayer")
        ALLOCATE (transmission_aux%solar_path2%od_sfrac(0:ncolumns, nchanprof), STAT=err)
        THROWM(err.NE.0, "allocation of transmission_aux % solar_path2 % od_sfrac")
        ALLOCATE (transmission_aux%solar_path2%od_frac_ac(0:ncolumns, nchanprof), STAT=err)
        THROWM(err.NE.0, "allocation of transmission_aux % solar_path2 % od_frac_ac")
      ELSE
        NULLIFY(transmission_aux%solar_path2%od_singlelayer)
        NULLIFY(transmission_aux%solar_path2%od_sfrac)
        NULLIFY(transmission_aux%solar_path2%od_frac_ac)
      ENDIF

      ! Solar path1
      ALLOCATE (transmission_aux%solar_path1)
      ALLOCATE (transmission_aux%solar_path1%tau_level(nlevels, 0:ncolumns, nchanprof), STAT=err)
      THROWM(err.NE.0, "allocation of transmission_aux % solar_path1 % tau_level")
      ALLOCATE (transmission_aux%solar_path1%tau_surf(0:ncolumns, nchanprof), STAT=err)
      THROWM(err.NE.0, "allocation of transmission_aux % solar_path1 % tau_surf")

      IF (do_scatt) THEN
        ALLOCATE (transmission_aux%solar_path1%tau_surf_ac(0:ncolumns, nchanprof), STAT=err)
        THROWM(err.NE.0, "allocation of transmission_aux % solar_path1 % tau_surf_ac")
      ELSE
        NULLIFY(transmission_aux%solar_path1%tau_surf_ac)
      ENDIF

      IF (do_single_scatt) THEN
        ALLOCATE (transmission_aux%solar_path1%od_singlelayer(0:su, nlayers, nchanprof), STAT=err)
        THROWM(err.NE.0, "allocation of transmission_aux % solar_path1 % od_singlelayer")
        ALLOCATE (transmission_aux%solar_path1%od_sfrac(0:ncolumns, nchanprof), STAT=err)
        THROWM(err.NE.0, "allocation of transmission_aux % solar_path1 % od_sfrac")
      ELSE
        NULLIFY(transmission_aux%solar_path1%od_singlelayer)
        NULLIFY(transmission_aux%solar_path1%od_sfrac)
      ENDIF

      IF (direct1 .AND. do_scatt .AND. .NOT. do_dom_vis) THEN
        ALLOCATE (transmission_aux%solar_path1%tau_level_r(nlevels, 0:ncolumns, nchanprof), STAT=err)
        THROWM(err.NE.0, "allocation of transmission_aux % solar_path1 % tau_level_r")
      ELSE
        NULLIFY(transmission_aux%solar_path1%tau_level_r)
      ENDIF

      IF (direct1 .AND. do_single_scatt) THEN
        ALLOCATE (transmission_aux%solar_path1%fac2(nlevels, 0:ncolumns, nchanprof), STAT=err)
        THROWM(err.NE.0, "allocation of transmission_aux % solar_path1 % fac2")
        ALLOCATE (transmission_aux%solar_path1%tau_surf_r(0:ncolumns, nchanprof), STAT=err)
        THROWM(err.NE.0, "allocation of transmission_aux % solar_path1 % tau_surf_r")
      ELSE
        NULLIFY(transmission_aux%solar_path1%fac2)
        NULLIFY(transmission_aux%solar_path1%tau_surf_r)
      ENDIF
    ENDIF

    IF (init1) CALL rttov_init_transmission_aux(opts, transmission_aux)
  ENDIF

  IF (asw == 0) THEN
    IF (ASSOCIATED(transmission_aux%fac1)) DEALLOCATE (transmission_aux%fac1, STAT=err)
    THROWM(err.NE.0, "deallocation of transmission_aux%fac1")
    IF (ASSOCIATED(transmission_aux%surf_fac)) DEALLOCATE (transmission_aux%surf_fac, STAT=err)
    THROWM(err.NE.0, "deallocation of transmission_aux%surf_fac")

    IF(ASSOCIATED(transmission_aux%thermal_path1)) THEN
      IF (ASSOCIATED(transmission_aux%thermal_path1%tau_surf_p)) &
          DEALLOCATE (transmission_aux%thermal_path1%tau_surf_p, STAT=err)
      THROWM(err.NE.0, "allocation of transmission_aux % thermal_path1 % tau_surf_p")
      IF (ASSOCIATED(transmission_aux%thermal_path1%tau_surf_p_r)) &
          DEALLOCATE (transmission_aux%thermal_path1%tau_surf_p_r, STAT=err)
      THROWM(err.NE.0, "allocation of transmission_aux % thermal_path1 % tau_surf_p_r")
      IF (ASSOCIATED(transmission_aux%thermal_path1%tau_level_p)) &
          DEALLOCATE (transmission_aux%thermal_path1%tau_level_p, STAT=err)
      THROWM(err.NE.0, "allocation of transmission_aux % thermal_path1 % tau_level_p")
      IF (ASSOCIATED(transmission_aux%thermal_path1%tau_level_p_r)) &
          DEALLOCATE (transmission_aux%thermal_path1%tau_level_p_r, STAT=err)
      THROWM(err.NE.0, "allocation of transmission_aux % thermal_path1 % tau_level_p_r")

      IF (ASSOCIATED(transmission_aux%thermal_path1%tau_level)) &
          DEALLOCATE (transmission_aux%thermal_path1%tau_level, STAT=err)
      THROWM(err.NE.0, "deallocation of transmission_aux%thermal_path1%tau_level")
      IF (ASSOCIATED(transmission_aux%thermal_path1%od_singlelayer)) &
          DEALLOCATE (transmission_aux%thermal_path1%od_singlelayer, STAT=err)
      THROWM(err.NE.0, "deallocation of transmission_aux%thermal_path1%od_singlelayer")
      IF (ASSOCIATED(transmission_aux%thermal_path1%od_sfrac)) &
          DEALLOCATE (transmission_aux%thermal_path1%od_sfrac, STAT=err)
      THROWM(err.NE.0, "deallocation of transmission_aux%thermal_path1%od_sfrac")
      IF (ASSOCIATED(transmission_aux%thermal_path1%tau_surf)) &
          DEALLOCATE (transmission_aux%thermal_path1%tau_surf, STAT=err)
      THROWM(err.NE.0, "deallocation of transmission_aux%thermal_path1%tau_surf")
      IF (ASSOCIATED(transmission_aux%thermal_path1%tau_surf_ac)) &
          DEALLOCATE (transmission_aux%thermal_path1%tau_surf_ac, STAT=err)
      THROWM(err.NE.0, "deallocation of transmission_aux%thermal_path1%tau_surf_ac")

      IF (ASSOCIATED(transmission_aux%thermal_path1%fac2)) &
          DEALLOCATE (transmission_aux%thermal_path1%fac2, STAT=err)
      THROWM(err.NE.0, "deallocation of transmission_aux%thermal_path1%fac2")
      IF (ASSOCIATED(transmission_aux%thermal_path1%tau_level_r)) &
          DEALLOCATE (transmission_aux%thermal_path1%tau_level_r, STAT=err)
      THROWM(err.NE.0, "deallocation of transmission_aux%thermal_path1%tau_level_r")
      IF (ASSOCIATED(transmission_aux%thermal_path1%tau_surf_r)) &
          DEALLOCATE (transmission_aux%thermal_path1%tau_surf_r, STAT=err)
      THROWM(err.NE.0, "deallocation of transmission_aux%thermal_path1%tau_surf_r")
      IF (ASSOCIATED(transmission_aux%thermal_path1%od_singlelayer_r)) &
          DEALLOCATE (transmission_aux%thermal_path1%od_singlelayer_r, STAT=err)
      THROWM(err.NE.0, "deallocation of transmission_aux%thermal_path1%od_singlelayer_r")
      IF (ASSOCIATED(transmission_aux%thermal_path1%od_sfrac_r)) &
          DEALLOCATE (transmission_aux%thermal_path1%od_sfrac_r, STAT=err)
      THROWM(err.NE.0, "deallocation of transmission_aux%thermal_path1%od_sfrac_r")

      DEALLOCATE (transmission_aux%thermal_path1, STAT=err)
      THROWM(err.NE.0, "deallocation of transmission_aux%thermal_path1")
    ENDIF

    IF (ASSOCIATED(transmission_aux%solar_path2)) THEN
      IF (ASSOCIATED(transmission_aux%solar_path2%tau_level)) &
          DEALLOCATE (transmission_aux%solar_path2%tau_level, STAT=err)
      THROWM(err.NE.0, "deallocation of transmission_aux%solar_path2%tau_level")
      IF (ASSOCIATED(transmission_aux%solar_path2%tau_surf)) &
          DEALLOCATE (transmission_aux%solar_path2%tau_surf, STAT=err)
      THROWM(err.NE.0, "deallocation of transmission_aux%solar_path2%tau_surf")
      IF (ASSOCIATED(transmission_aux%solar_path2%tau_surf_ac)) &
          DEALLOCATE (transmission_aux%solar_path2%tau_surf_ac, STAT=err)
      THROWM(err.NE.0, "deallocation of transmission_aux%solar_path2%tau_surf_ac")
      IF (ASSOCIATED(transmission_aux%solar_path2%od_singlelayer)) &
          DEALLOCATE (transmission_aux%solar_path2%od_singlelayer, STAT=err)
      THROWM(err.NE.0, "deallocation of transmission_aux%solar_path2%od_singlelayer")
      IF (ASSOCIATED(transmission_aux%solar_path2%od_sfrac)) &
          DEALLOCATE (transmission_aux%solar_path2%od_sfrac, STAT=err)
      THROWM(err.NE.0, "deallocation of transmission_aux%solar_path2%od_sfrac")
      IF (ASSOCIATED(transmission_aux%solar_path2%od_frac_ac)) &
          DEALLOCATE (transmission_aux%solar_path2%od_frac_ac, STAT=err)
      THROWM(err.NE.0, "deallocation of transmission_aux%solar_path2%od_frac_ac ")

      DEALLOCATE (transmission_aux%solar_path2, STAT=err)
      THROWM(err.NE.0, "deallocation of transmission_aux%solar_path2")
    ENDIF

    IF (ASSOCIATED(transmission_aux%solar_path1)) THEN
      IF (ASSOCIATED(transmission_aux%solar_path1%tau_surf_ac)) &
          DEALLOCATE (transmission_aux%solar_path1%tau_surf_ac, STAT=err)
      THROWM(err.NE.0, "deallocation of transmission_aux%solar_path1%tau_surf_ac")
      IF (ASSOCIATED(transmission_aux%solar_path1%tau_level_r)) &
          DEALLOCATE (transmission_aux%solar_path1%tau_level_r, STAT=err)
      THROWM(err.NE.0, "deallocation of transmission_aux%solar_path1%tau_level_r")
      IF (ASSOCIATED(transmission_aux%solar_path1%tau_level)) &
          DEALLOCATE (transmission_aux%solar_path1%tau_level, STAT=err)
      THROWM(err.NE.0, "deallocation of transmission_aux%solar_path1%tau_level")
      IF (ASSOCIATED(transmission_aux%solar_path1%tau_surf)) &
          DEALLOCATE (transmission_aux%solar_path1%tau_surf, STAT=err)
      THROWM(err.NE.0, "deallocation of transmission_aux%solar_path1%tau_surf")

      IF (ASSOCIATED(transmission_aux%solar_path1%od_singlelayer)) &
          DEALLOCATE (transmission_aux%solar_path1%od_singlelayer, STAT=err)
      THROWM(err.NE.0, "deallocation of transmission_aux%solar_path1%od_singlelayer")
      IF (ASSOCIATED(transmission_aux%solar_path1%od_sfrac)) &
          DEALLOCATE (transmission_aux%solar_path1%od_sfrac, STAT=err)
      THROWM(err.NE.0, "deallocation of transmission_aux%solar_path1%od_sfrac")
      IF (ASSOCIATED(transmission_aux%solar_path1%fac2)) &
          DEALLOCATE (transmission_aux%solar_path1%fac2, STAT=err)
      THROWM(err.NE.0, "deallocation of transmission_aux%solar_path1%fac2")
      IF (ASSOCIATED(transmission_aux%solar_path1%tau_surf_r)) &
          DEALLOCATE (transmission_aux%solar_path1%tau_surf_r, STAT=err)
      THROWM(err.NE.0, "deallocation of transmission_aux%solar_path1%tau_surf_r")

      DEALLOCATE (transmission_aux%solar_path1, STAT=err)
      THROWM(err.NE.0, "deallocation of transmission_aux%solar_path1")
    ENDIF

    CALL nullify_struct()
  ENDIF
  CATCH

CONTAINS

  SUBROUTINE nullify_struct()
    NULLIFY (transmission_aux%fac1)
    NULLIFY (transmission_aux%surf_fac)
    NULLIFY (transmission_aux%thermal_path1)
    NULLIFY (transmission_aux%solar_path2)
    NULLIFY (transmission_aux%solar_path1)
  END SUBROUTINE nullify_struct

END SUBROUTINE rttov_alloc_transmission_aux
