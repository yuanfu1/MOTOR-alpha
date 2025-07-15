! Description:
!> @file
!!   Allocate/deallocate structure(s) for internal RTTOV state.
!
!> @brief
!!   Allocate/deallocate structure(s) for internal RTTOV state.
!!
!! @details
!!   The rttov_traj structure contains internal state for RTTOV.
!!   It is not mandatory to allocate these outside of calls to
!!   the direct, TL, AD or K models, but on some architectures
!!   it can be faster to do this when making repeated calls to
!!   RTTOV as the memory allocations are relatively expensive.
!!   There is no disadvantage to using this subroutine to make
!!   these allocations externally to RTTOV.
!!
!!   This subroutine can be called to (de)allocate the direct model
!!   traj structure and also the traj_tl, traj_ad and/or traj_k
!!   structures in a single call.
!!
!! @param[out]    err            status on exit
!! @param[in]     nprofiles      number of profiles being simulated
!! @param[in]     nchanprof      total number of channels being simulated
!! @param[in]     opts           options to configure the simulations
!! @param[in]     nlevels        number of levels in input profiles
!! @param[in]     coefs          coefficients structure for instrument to simulate
!! @param[in]     asw            1_jpim => allocate; 0_jpim => deallocate
!! @param[in,out] traj           rttov_traj structure for the direct model, optional
!! @param[in,out] traj_tl        rttov_traj structure for the TL model, optional
!! @param[in,out] traj_ad        rttov_traj structure for the AD model, optional
!! @param[in,out] traj_k         rttov_traj structure for the K model, optional
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
SUBROUTINE rttov_alloc_traj( &
              err,       &
              nprofiles, &
              nchanprof, &
              opts,      &
              nlevels,   &
              coefs,     &
              asw,       &
              traj,      &
              traj_tl,   &
              traj_ad,   &
              traj_k)
!INTF_OFF
#include "throw.h"
!INTF_ON
  USE parkind1, ONLY : jpim
  USE rttov_types, ONLY :  &
         rttov_options,  &
         rttov_coefs,    &
         rttov_traj
!INTF_OFF
  USE parkind1, ONLY : jplm
!INTF_ON

  IMPLICIT NONE

  INTEGER(KIND=jpim),  INTENT(OUT)              :: err
  INTEGER(KIND=jpim),  INTENT(IN)               :: nprofiles
  INTEGER(KIND=jpim),  INTENT(IN)               :: nchanprof
  TYPE(rttov_options), INTENT(IN)               :: opts
  INTEGER(KIND=jpim),  INTENT(IN)               :: nlevels
  TYPE(rttov_coefs),   INTENT(IN),    TARGET    :: coefs    ! Target attribute needed
  INTEGER(KIND=jpim),  INTENT(IN)               :: asw
  TYPE(rttov_traj),    INTENT(INOUT), OPTIONAL  :: traj, traj_tl, traj_ad, traj_k
!INTF_END

#include "rttov_errorreport.interface"

  TRY

  IF (PRESENT(traj)) THEN
    CALL alloc_traj( &
         err,         &
         nprofiles,   &
         nchanprof,   &
         opts,        &
         nlevels,     &
         coefs,       &
         asw,         &
         .TRUE._jplm, &
         traj)
    THROW(err.NE.0)
  ENDIF
  IF (PRESENT(traj_tl)) THEN
    CALL alloc_traj( &
         err,          &
         nprofiles,    &
         nchanprof,    &
         opts,         &
         nlevels,      &
         coefs,        &
         asw,          &
         .FALSE._jplm, &
         traj_tl)
    THROW(err.NE.0)
  ENDIF
  IF (PRESENT(traj_ad)) THEN
    CALL alloc_traj( &
         err,          &
         nprofiles,    &
         nchanprof,    &
         opts,         &
         nlevels,      &
         coefs,        &
         asw,          &
         .FALSE._jplm, &
         traj_ad)
    THROW(err.NE.0)
  ENDIF
  IF (PRESENT(traj_k)) THEN
    CALL alloc_traj( &
         err,          &
         nchanprof,    &
         nchanprof,    &
         opts,         &
         nlevels,      &
         coefs,        &
         asw,          &
         .FALSE._jplm, &
         traj_k)
    THROW(err.NE.0)
  ENDIF
  CATCH
CONTAINS
  SUBROUTINE alloc_traj( &
                err,       &
                nprofiles, &
                nchanprof, &
                opts,      &
                nlevels,   &
                coefs,     &
                asw,       &
                direct,    &
                traj)
    INTEGER(KIND=jpim),  INTENT(OUT)        :: err
    INTEGER(KIND=jpim),  INTENT(IN)         :: nprofiles
    INTEGER(KIND=jpim),  INTENT(IN)         :: nchanprof
    TYPE(rttov_options), INTENT(IN)         :: opts
    INTEGER(KIND=jpim),  INTENT(IN)         :: nlevels
    TYPE(rttov_coefs),   INTENT(IN), TARGET :: coefs
    INTEGER(KIND=jpim),  INTENT(IN)         :: asw
    LOGICAL(KIND=jplm),  INTENT(IN)         :: direct
    TYPE(rttov_traj),    INTENT(INOUT)      :: traj
#include "rttov_alloc_prof.interface"
#include "rttov_alloc_prof_internal.interface"
#include "rttov_alloc_aux_prof.interface"
#include "rttov_alloc_aux_prof_coef.interface"
#include "rttov_alloc_raytracing.interface"
#include "rttov_alloc_predictor.interface"
#include "rttov_alloc_opdp_path.interface"
#include "rttov_alloc_ircld.interface"
#include "rttov_alloc_trans_scatt_ir.interface"
#include "rttov_alloc_sunglint.interface"
#include "rttov_alloc_auxrad.interface"
    INTEGER(KIND=jpim)  :: nlayers
    TYPE(rttov_options) :: opts_coef
    TRY
    nlayers   = nlevels - 1
    opts_coef                 = opts
    opts_coef%rt_ir%addclouds = .FALSE.
    opts_coef%rt_ir%addaerosl = .FALSE.
    IF (asw .EQ. 1_jpim) THEN
      ALLOCATE (traj%profiles_coef(nprofiles), STAT = err)
      THROWM(err .NE. 0, "allocation of traj%profiles_coef")
      ALLOCATE (traj%profiles_int(nprofiles), STAT = err)
      THROWM(err .NE. 0, "allocation of traj%profiles_int")
! ON LEVELS
! ---------
      CALL rttov_alloc_prof( &
              err,                &
              nprofiles,          &
              traj%profiles_coef, &
              coefs%coef%nlevels, &
              opts_coef,          &
              1_jpim,             &
              coefs = coefs)
      THROW(err .NE. 0)
      CALL rttov_alloc_prof_internal( &
              err,                &
              traj%profiles_int,  &
              nlevels,            &
              opts,               &
              coefs,              &
              1_jpim)
      THROW(err .NE. 0)
      CALL rttov_alloc_raytracing( &
              err,                  &
              nprofiles,            &
              opts%rt_ir%addsolar,  &
              traj%raytracing_coef, &
              coefs%coef%nlevels,   &
              1_jpim)
      THROW(err .NE. 0)
      CALL rttov_alloc_raytracing( &
              err,                 &
              nprofiles,           &
              opts%rt_ir%addsolar, &
              traj%raytracing,     &
              nlevels,             &
              1_jpim)
      THROW(err .NE. 0)
      CALL rttov_alloc_opdp_path( &
              err,            &
              opts,           &
              traj%opdp_path, &
              nlevels,        &
              nchanprof,      &
              1_jpim,         &
              init=.TRUE._jplm)
      THROW(err .NE. 0)
      CALL rttov_alloc_opdp_path( &
              err,                 &
              opts,                &
              traj%opdp_path_coef, &
              coefs%coef%nlevels,  &
              nchanprof,           &
              1_jpim,              &
              init=.TRUE._jplm)
      THROW(err .NE. 0)

! ON LAYERS
! ---------
!
      CALL rttov_alloc_predictor( &
              err,             &
              nprofiles,       &
              traj%predictors, &
              coefs%coef,      &
              1_jpim,          &
              opts%rt_ir%addsolar)
      THROW(err .NE. 0)
      CALL rttov_alloc_ircld( &
              err,            &
              opts,           &
              nprofiles,      &
              traj%ircld,     &
              nlayers,        &
              1_jpim,         &
              direct = direct)
      THROW(err .NE. 0)
      ALLOCATE (traj%diffuse_refl(nchanprof), STAT = err)
      THROWM(err.NE.0,"Allocation of diffuse_refl failed")
      CALL rttov_alloc_aux_prof( &
              err,                  &
              nprofiles,            &
              nlevels,              &
              traj%aux_prof,        &
              opts,                 &
              coefs%coef,           &
              1_jpim)
      THROW(err.NE.0)
      CALL rttov_alloc_aux_prof_coef( &
              err,                  &
              nprofiles,            &
              coefs%coef%nlevels,   &
              traj%aux_prof_coef,   &
              coefs%coef,           &
              1_jpim)
      THROW(err.NE.0)
      IF (opts%rt_ir%addaerosl .OR. opts%rt_ir%addclouds) THEN
        CALL rttov_alloc_trans_scatt_ir( &
                err,                        &
                opts,                       &
                coefs,                      &
                traj%transmission_scatt_ir, &
                nchanprof,                  &
                nlayers,                    &
                1_jpim,                     &
                direct = direct)
        THROWM(err .NE. 0, "allocation of trans_scatt_ir")
      ENDIF
      CALL rttov_alloc_auxrad( &
              err,             &
              traj%auxrad,     &
              nlevels,         &
              nchanprof,       &
              1_jpim,          &
              direct = direct)
      THROWM(err .NE. 0, "allocation of auxrad")
      IF (opts%rt_ir%addsolar) THEN
        ALLOCATE (traj%fresnrefl(nchanprof), STAT = err)
        THROWM(err.NE.0,"allocation of fresnrefl")
        CALL rttov_alloc_sunglint( &
                err,                  &
                traj%sunglint,        &
                opts,                 &
                nprofiles,            &
                coefs%coef%ws_nomega, &
                1_jpim,               &
                init = .TRUE._jplm,   &
                direct = direct)
        THROW(err.NE.0)
      ENDIF
      traj%coefs => coefs
      traj%nchanprof = nchanprof
      traj%opts      = opts
      traj%nlevels   = nlevels
      traj%nlayers   = traj%nlevels - 1
    ELSE
      IF (opts%rt_ir%addsolar) THEN
        DEALLOCATE (traj%fresnrefl, STAT = err)
        THROWM(err.NE.0,"deallocation of fresnrefl")
        CALL rttov_alloc_sunglint( &
                err,                  &
                traj%sunglint,        &
                opts,                 &
                nprofiles,            &
                coefs%coef%ws_nomega, &
                0_jpim)
        THROW(err.NE.0)
      ENDIF
      CALL rttov_alloc_auxrad( &
              err,             &
              traj%auxrad,     &
              nlevels,         &
              nchanprof,       &
              0_jpim)
      THROWM(err .NE. 0, "deallocation of auxrad")
      IF (opts%rt_ir%addaerosl .OR. opts%rt_ir%addclouds) THEN
        CALL rttov_alloc_trans_scatt_ir( &
                err,                        &
                opts,                       &
                coefs,                      &
                traj%transmission_scatt_ir, &
                nchanprof,                  &
                nlayers,                    &
                0_jpim)
        THROWM(err .NE. 0, "deallocation of trans_scatt_ir")
      ENDIF
      CALL rttov_alloc_aux_prof_coef( &
              err,                  &
              nprofiles,            &
              coefs%coef%nlevels,   &
              traj%aux_prof_coef,   &
              coefs%coef,           &
              0_jpim)
      THROW(err.NE.0)
      CALL rttov_alloc_aux_prof( &
              err,                  &
              nprofiles,            &
              nlevels,              &
              traj%aux_prof,        &
              opts,                 &
              coefs%coef,           &
              0_jpim)
      THROW(err.NE.0)
      DEALLOCATE (traj%diffuse_refl, STAT = err)
      THROWM(err.NE.0,"deallocation of diffuse_refl failed")
! ON LAYERS
! ---------
!
      CALL rttov_alloc_ircld( &
              err,            &
              opts,           &
              nprofiles,      &
              traj%ircld,     &
              nlayers,        &
              0_jpim)
      THROW(err .NE. 0)
      CALL rttov_alloc_opdp_path( &
              err,                 &
              opts,                &
              traj%opdp_path_coef, &
              coefs%coef%nlayers,  &
              nchanprof,           &
              0_jpim)
      THROW(err .NE. 0)
      CALL rttov_alloc_opdp_path( &
              err,            &
              opts,           &
              traj%opdp_path, &
              nlayers,        &
              nchanprof,      &
              0_jpim)
      THROW(err .NE. 0)
      CALL rttov_alloc_predictor( &
              err,             &
              nprofiles,       &
              traj%predictors, &
              coefs%coef,      &
              0_jpim,          &
              opts%rt_ir%addsolar)
      THROW(err .NE. 0)
!
! ON LEVELS
! ---------
!
!
      CALL rttov_alloc_raytracing( &
              err,                 &
              nprofiles,           &
              opts%rt_ir%addsolar, &
              traj%raytracing,     &
              nlevels,             &
              0_jpim)
      THROW(err .NE. 0)
      CALL rttov_alloc_raytracing( &
              err,                  &
              nprofiles,            &
              opts%rt_ir%addsolar,  &
              traj%raytracing_coef, &
              coefs%coef%nlevels,   &
              0_jpim)
      THROW(err .NE. 0)
      CALL rttov_alloc_prof_internal( &
              err,                &
              traj%profiles_int,  &
              nlevels,            &
              opts,               &
              coefs,              &
              0_jpim)
      CALL rttov_alloc_prof( &
              err,                &
              nprofiles,          &
              traj%profiles_coef, &
              coefs%coef%nlevels, &
              opts_coef,          &
              0_jpim)
      THROW(err .NE. 0)
      DEALLOCATE (traj%profiles_int, STAT = err)
      THROWM(err .NE. 0, "deallocation of traj%profiles_int")
      DEALLOCATE (traj%profiles_coef, STAT = err)
      THROWM(err .NE. 0, "deallocation of traj%profiles_coef")
    ENDIF
    CATCH
  END SUBROUTINE
END SUBROUTINE 
