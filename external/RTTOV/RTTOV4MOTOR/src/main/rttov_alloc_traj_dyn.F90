! Description:
!> @file
!!   Allocate/deallocate structure for "dynamic" internal RTTOV state.
!
!> @brief
!!   Allocate/deallocate structure for "dynamic" internal RTTOV state.
!!
!! @details
!!   The dynamic trajectory structure contains internal state for
!!   each of the RTTOV direct, TL, AD and K models.
!!
!!   This subroutine is used to allocate memory internally within
!!   each RTTOV model. Note that this subroutine takes ncolumns, the
!!   number of cloudy columns in IR scattering simulations,
!!   as an argument. This means that in general it cannot be allocated
!!   outside of RTTOV as ncolumns is determined dynamically at run-time.
!!   However for non-cloudy-IR simulations it would be possible to allocate
!!   this structure with ncolumns=0 before calling rttov_direct and hence to
!!   gain access to some of the internal RTTOV state which could potentially
!!   be useful for debugging purposes.
!!
!! @param[out]    err              status on exit
!! @param[in,out] traj_dyn         RTTOV dynamic trajectory structure (for direct, TL, AD or K model)
!! @param[in]     opts             options to configure the simulations
!! @param[in]     coefs            coefficients structure for instrument to simulate
!! @param[in]     nchanprof        total number of channels being simulated (SIZE(chanprof))
!! @param[in]     nlayers          number of input profile layers (nlevels-1)
!! @param[in]     ncolumns         number of cloud columns as calculated by rttov_cloud_overlap
!! @param[in]     dom_nstr         number of DOM streams (discrete ordinates)
!! @param[in]     thermal          flag to indicate channels with thermal emission contribution
!! @param[in]     solar            flag to indicate channels with solar contribution
!! @param[in]     dothermal        flag to indicate thermal simulations
!! @param[in]     do_mfasis        flag to indicate MFASIS simulation
!! @param[in]     asw              1_jpim => allocate; 0_jpim => deallocate
!! @param[in]     traj_dyn_direct  direct model traj_dyn structure (if calling from TL/AD/K), optional
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
SUBROUTINE rttov_alloc_traj_dyn( &
              err,               &
              traj_dyn,          &
              opts,              &
              coefs,             &
              nchanprof,         &
              nlayers,           &
              ncolumns,          &
              dom_nstr,          &
              thermal,           &
              solar,             &
              dothermal,         &
              do_mfasis,         &
              asw,               &
              traj_dyn_direct)
!INTF_OFF
#include "throw.h"
!INTF_ON
  USE rttov_types, ONLY : rttov_traj_dyn, rttov_options, rttov_coefs
!INTF_OFF
  USE rttov_const, ONLY : vis_scatt_dom, ir_scatt_dom
!INTF_ON
  USE parkind1, ONLY : jpim, jplm

  IMPLICIT NONE

  INTEGER(KIND=jpim),   INTENT(OUT)          :: err
  TYPE(rttov_traj_dyn), INTENT(INOUT)        :: traj_dyn
  TYPE(rttov_options),  INTENT(IN)           :: opts
  TYPE(rttov_coefs),    INTENT(IN)           :: coefs
  INTEGER(KIND=jpim),   INTENT(IN)           :: nchanprof
  INTEGER(KIND=jpim),   INTENT(IN)           :: nlayers
  INTEGER(KIND=jpim),   INTENT(IN)           :: ncolumns
  INTEGER(KIND=jpim),   INTENT(IN)           :: dom_nstr
  LOGICAL(KIND=jplm),   INTENT(IN)           :: thermal(:)
  LOGICAL(KIND=jplm),   INTENT(IN)           :: solar(:)
  LOGICAL(KIND=jplm),   INTENT(IN)           :: dothermal
  LOGICAL(KIND=jplm),   INTENT(IN)           :: do_mfasis
  INTEGER(KIND=jpim),   INTENT(IN)           :: asw
  TYPE(rttov_traj_dyn), INTENT(IN), OPTIONAL :: traj_dyn_direct
!INTF_END

  LOGICAL(KIND=jplm) :: direct

#include "rttov_errorreport.interface"
#include "rttov_alloc_transmission_aux.interface"
#include "rttov_alloc_auxrad_column.interface"
#include "rttov_alloc_trans_scatt_ir.interface"
#include "rttov_alloc_profiles_dom.interface"
#include "rttov_alloc_dom_state.interface"
#include "rttov_alloc_mfasis_refl.interface"

  TRY

  direct = .NOT. PRESENT(traj_dyn_direct)

  IF (do_mfasis) THEN
    CALL rttov_alloc_mfasis_refl( &
                err,                  &
                traj_dyn%mfasis_refl, &
                nchanprof,            &
                asw,                  &
                ncolumns)
    THROWM(err .NE. 0, "allocation of mfasis_refl")
  ENDIF

  IF (dothermal .OR. .NOT. do_mfasis) THEN
    CALL rttov_alloc_transmission_aux( &
            err,                           &
            opts,                          &
            traj_dyn%transmission_aux,     &
            nlayers,                       &
            nchanprof,                     &
            asw,                           &
            ncolumns,                      &
            direct = direct)
    THROWM(err .NE. 0, "allocation of transmission_aux")

    CALL rttov_alloc_auxrad_column( &
            err,                      &
            traj_dyn%auxrad_column,   &
            opts,                     &
            ncolumns,                 &
            nlayers,                  &
            nchanprof,                &
            asw,                      &
            direct = direct)
    THROWM(err .NE. 0, "allocation of auxrad_column")
  ENDIF

  IF (opts%rt_ir%addaerosl .OR. opts%rt_ir%addclouds) THEN
    CALL rttov_alloc_trans_scatt_ir( &
            err,                                    &
            opts,                                   &
            coefs,                                  &
            traj_dyn%transmission_scatt_ir_dyn,     &
            nchanprof,                              &
            nlayers,                                &
            asw,                                    &
            .TRUE._jplm,                            &
            ncolumns,                               &
            dom_nstr,                               &
            thermal,                                &
            solar,                                  &
            direct = direct)
    THROWM(err .NE. 0, "allocation of trans_scatt_ir")

    IF (direct) THEN
      ! Called from direct model, allocate direct model profiles_dom
      IF (opts%rt_ir%addsolar .AND. opts%rt_ir%vis_scatt_model == vis_scatt_dom) THEN
        CALL rttov_alloc_profiles_dom( &
                err,                         &
                traj_dyn%profiles_dom_solar, &
                nchanprof,                   &
                ncolumns,                    &
                solar,                       &
                asw)
        THROWM(err .NE. 0, "allocation of profiles_dom_solar")
      ENDIF

      IF (opts%rt_ir%ir_scatt_model == ir_scatt_dom) THEN
        CALL rttov_alloc_profiles_dom( &
                err,                           &
                traj_dyn%profiles_dom_thermal, &
                nchanprof,                     &
                ncolumns,                      &
                thermal,                       &
                asw)
        THROWM(err .NE. 0, "allocation of profiles_dom_thermal")
      ENDIF

      IF (asw == 1_jpim) NULLIFY(traj_dyn%dom_state_solar, traj_dyn%dom_state_thermal)
    ELSE
      ! Called from TL/AD/K, allocate using sizes from direct model profiles_dom
      IF (opts%rt_ir%addsolar .AND. opts%rt_ir%vis_scatt_model == vis_scatt_dom) THEN
        CALL rttov_alloc_profiles_dom( &
                err,                         &
                traj_dyn%profiles_dom_solar, &
                nchanprof,                   &
                ncolumns,                    &
                solar,                       &
                asw,                         &
                traj_dyn_direct%profiles_dom_solar)
        THROWM(err .NE. 0, "allocation of profiles_dom_solar")
      ENDIF

      IF (opts%rt_ir%ir_scatt_model == ir_scatt_dom) THEN
        CALL rttov_alloc_profiles_dom( &
                err,                           &
                traj_dyn%profiles_dom_thermal, &
                nchanprof,                     &
                ncolumns,                      &
                thermal,                       &
                asw,                           &
                traj_dyn_direct%profiles_dom_thermal)
        THROWM(err .NE. 0, "allocation of profiles_dom_thermal")
      ENDIF

      IF (asw == 1_jpim) NULLIFY(traj_dyn%dom_state_solar, traj_dyn%dom_state_thermal)
    ENDIF

    IF (traj_dyn%from_tladk .OR. .NOT. direct) THEN
      ! rttov_direct called from TL/AD/K or deallocation
      IF (opts%rt_ir%addsolar .AND. opts%rt_ir%vis_scatt_model == vis_scatt_dom) THEN
        CALL rttov_alloc_dom_state( &
                err,                             &
                traj_dyn%dom_state_solar,        &
                nchanprof,                       &
                nlayers,                         &
                ncolumns,                        &
                dom_nstr,                        &
                .TRUE._jplm,                     & ! dosolar
                solar,                           &
                asw)
        THROWM(err .NE. 0, "allocation of dom_state_solar")
      ENDIF

      IF (opts%rt_ir%ir_scatt_model == ir_scatt_dom) THEN
        CALL rttov_alloc_dom_state( &
                err,                               &
                traj_dyn%dom_state_thermal,        &
                nchanprof,                         &
                nlayers,                           &
                ncolumns,                          &
                dom_nstr,                          &
                .FALSE._jplm,                      & ! dothermal
                thermal,                           &
                asw)
        THROWM(err .NE. 0, "allocation of dom_state_thermal")
      ENDIF
    ENDIF

  ENDIF

  traj_dyn%ncolumns = ncolumns

  CATCH
END SUBROUTINE
