! Description:
!> @file
!!   Allocate/deallocate internal transmission_scatt_ir structure.
!
!> @brief
!!   Allocate/deallocate internal transmission_scatt_ir structure.
!!
!! @details
!!   This structure holds information related to VIS/IR scattering
!!   simulations.
!!
!!   Some arrays are only required by the direct model instance
!!   of this structure: these are allocated if the optional "direct"
!!   argument is TRUE. If omitted, it is assumed to be FALSE.
!!
!! @param[out]    err                    status on exit
!! @param[in]     opts                   options to configure the simulations
!! @param[in]     coefs                  RTTOV coefficients structure
!! @param[in,out] trans_scatt_ir         transmission_scatt_ir structure to (de)allocate
!! @param[in]     nchanprof              total number of channels being simulated (SIZE(chanprof))
!! @param[in]     nlayers                number of input profile layers (nlevels-1)
!! @param[in]     asw                    1_jpim => allocate; 0_jpim => deallocate
!! @param[in]     dynamic                if TRUE then allocate dynamically-sized arrays (i.e. those with
!!                                          an ncolumns - dimension), optional
!! @param[in]     ncolumns               number of cloud columns as calculated by rttov_cloud_overlap, optional
!! @param[in]     dom_nstr               number of DOM streams used if DOM solver is being used, optional
!! @param[in]     thermal                per-channel flag to indicate if thermal (emissive) simulations are being
!!                                          performed, optional
!! @param[in]     solar                  per-channel flag to indicate if any solar simulations are being performed, optional
!! @param[in]     init                   set .TRUE. to initialise newly allocated structures, optional
!! @param[in]     direct                 if FALSE then direct-model-only arrays are not allocated, optional
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
SUBROUTINE rttov_alloc_trans_scatt_ir( &
              err,            &
              opts,           &
              coefs,          &
              trans_scatt_ir, &
              nchanprof,      &
              nlayers,        &
              asw,            &
              dynamic,        &
              ncolumns,       &
              dom_nstr,       &
              thermal,        &
              solar,          &
              init,           &
              direct)
!INTF_OFF
#include "throw.h"
!INTF_ON
  USE rttov_types, ONLY : rttov_transmission_scatt_ir, rttov_options, rttov_coefs
  USE parkind1, ONLY : jpim, jplm
!INTF_OFF
  USE rttov_types, ONLY : rttov_scatt_ir_aercld
  USE rttov_const, ONLY : &
    vis_scatt_dom,    &
    vis_scatt_single, &
    vis_scatt_mfasis, &
    ir_scatt_dom,     &
    ir_scatt_chou,    &
    ncldtyp
!INTF_ON
  IMPLICIT NONE

  INTEGER(KIND=jpim),                INTENT(OUT)          :: err
  TYPE(rttov_options),               INTENT(IN)           :: opts
  TYPE(rttov_coefs),                 INTENT(IN)           :: coefs
  TYPE(rttov_transmission_scatt_ir), INTENT(INOUT)        :: trans_scatt_ir
  INTEGER(KIND=jpim),                INTENT(IN)           :: nchanprof
  INTEGER(KIND=jpim),                INTENT(IN)           :: nlayers
  INTEGER(KIND=jpim),                INTENT(IN)           :: asw
  LOGICAL(KIND=jplm),                INTENT(IN), OPTIONAL :: dynamic
  INTEGER(KIND=jpim),                INTENT(IN), OPTIONAL :: ncolumns
  INTEGER(KIND=jpim),                INTENT(IN), OPTIONAL :: dom_nstr
  LOGICAL(KIND=jplm),                INTENT(IN), OPTIONAL :: thermal(nchanprof)
  LOGICAL(KIND=jplm),                INTENT(IN), OPTIONAL :: solar(nchanprof)
  LOGICAL(KIND=jplm),                INTENT(IN), OPTIONAL :: init
  LOGICAL(KIND=jplm),                INTENT(IN), OPTIONAL :: direct
!INTF_END
#include "rttov_errorreport.interface"
#include "rttov_init_trans_scatt_ir.interface"

  INTEGER(KIND=jpim) :: nlevels, sl, su, i, naer, ncld
  LOGICAL(KIND=jplm) :: init1, dynamic1, direct1
  LOGICAL(KIND=jplm) :: do_dom, do_chou, do_dom_ir, do_dom_vis
  LOGICAL(KIND=jplm) :: do_single_scatt, do_mfasis, do_thermal
!- End of header --------------------------------------------------------
  TRY

  IF (asw == 1) THEN
    CALL nullify_struct()

    nlevels = nlayers + 1
    init1   = .FALSE.
    direct1 = .FALSE.
    dynamic1 = .FALSE.
    IF (PRESENT(init)) init1   = init
    IF (PRESENT(direct)) direct1 = direct
    IF (PRESENT(dynamic)) dynamic1 = dynamic
    IF (dynamic1 .AND. (.NOT. PRESENT(ncolumns) .OR. .NOT. PRESENT(dom_nstr) .OR. &
                        .NOT. PRESENT(thermal) .OR. .NOT. PRESENT(solar))) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0, "ncolumns, dom_nstr, thermal and solar arguments required")
    ENDIF

    do_dom_vis = opts%rt_ir%addsolar .AND. opts%rt_ir%vis_scatt_model == vis_scatt_dom
    do_dom_ir = opts%rt_ir%ir_scatt_model == ir_scatt_dom
    do_dom = do_dom_ir .OR. do_dom_vis
    do_chou = opts%rt_ir%ir_scatt_model == ir_scatt_chou
    do_single_scatt = opts%rt_ir%addsolar .AND. opts%rt_ir%vis_scatt_model == vis_scatt_single
    do_mfasis = opts%rt_ir%addsolar .AND. opts%rt_ir%vis_scatt_model == vis_scatt_mfasis
    do_thermal = .TRUE.
    IF (dynamic1) do_thermal = ANY(thermal)

    IF (opts%rt_ir%addaerosl .OR. opts%rt_ir%dom_rayleigh) THEN
      sl = 0
    ELSE
      sl = 1
    ENDIF
    IF (opts%rt_ir%addclouds) THEN
      su = 1
    ELSE
      su = 0
    ENDIF

    ! Layer quantities take one of two values in each column: non-cloudy and cloudy.
    ! 0:0 if aerosol-only, 1:1 if cloud-only, 0:1 for aer+cld
    ! NB in cloud-only case the clear (non-cloudy) column values are all zero for scattering particles

    IF (dynamic1 .AND. (do_thermal .OR. .NOT. do_mfasis)) THEN
      ALLOCATE (trans_scatt_ir%opdpac(nlevels, 0:ncolumns, nchanprof), STAT=err)
      THROWM(err.NE.0, "mem allocation error")

      IF (opts%rt_ir%addsolar) THEN
        ALLOCATE (trans_scatt_ir%opdpacsun(nlevels, 0:ncolumns, nchanprof), STAT=err)
        THROWM(err.NE.0, "mem allocation error")
      ENDIF

      IF (do_dom) THEN
        ALLOCATE (trans_scatt_ir%phasefn(nchanprof), STAT=err)
        THROWM(err.NE.0, "mem allocation error")
        DO i = 1, nchanprof
          NULLIFY(trans_scatt_ir%phasefn(i)%legcoef)
          IF ((do_dom_ir .AND. thermal(i)) .OR. (do_dom_vis .AND. solar(i))) THEN
            ALLOCATE (trans_scatt_ir%phasefn(i)%legcoef(0:dom_nstr, sl:su, nlayers), STAT=err)
            THROWM(err.NE.0, "mem allocation error")
          ENDIF
        ENDDO
      ENDIF
    ELSE

      naer = 1
      IF (opts%rt_ir%addaerosl) THEN
        IF (.NOT. opts%rt_ir%user_aer_opt_param) THEN
          naer = coefs%coef_scatt%optp_aer%ntypes
        ENDIF
      ENDIF

      ncld = 1
      IF (opts%rt_ir%addclouds) THEN
        IF (.NOT. opts%rt_ir%user_cld_opt_param) THEN
          ncld = ncldtyp
        ENDIF
      ENDIF

      IF (do_mfasis) THEN
        IF (opts%rt_ir%addaerosl) THEN
          ALLOCATE (trans_scatt_ir%opdpext(naer, nlayers, nchanprof), STAT=err)
          THROWM(err.NE.0, "mem allocation error")
        ELSEIF (opts%rt_ir%addclouds) THEN
          ALLOCATE (trans_scatt_ir%opdpext(ncld, nlayers, nchanprof), STAT=err)
          THROWM(err.NE.0, "mem allocation error")
        ENDIF
      ENDIF

      IF (opts%rt_ir%dom_rayleigh .AND. do_dom_vis) THEN
        ALLOCATE (trans_scatt_ir%ray_sca(nlayers, nchanprof), STAT=err)
        THROWM(err.NE.0, "mem allocation error")
      ENDIF

      IF (do_thermal .OR. .NOT. do_mfasis) THEN
        IF (opts%rt_ir%addaerosl .OR. opts%rt_ir%dom_rayleigh) THEN
          CALL alloc_scatt_ir_aercld(err, trans_scatt_ir%aer, naer, opts%rt_ir%user_aer_opt_param, .FALSE._jplm)
          THROW(err.NE.0)
        ENDIF

        IF (opts%rt_ir%addclouds) THEN
          CALL alloc_scatt_ir_aercld(err, trans_scatt_ir%cld, ncld, opts%rt_ir%user_cld_opt_param, .TRUE._jplm)
          THROW(err.NE.0)
        ENDIF

        ALLOCATE (trans_scatt_ir%opdpacl(0:su, nlayers, nchanprof), STAT=err)
        THROWM(err.NE.0, "mem allocation error")

        IF (opts%rt_ir%addsolar .AND. .NOT. do_mfasis) THEN
          ALLOCATE (trans_scatt_ir%opdpaclsun(0:su, nlayers, nchanprof), STAT=err)
          THROWM(err.NE.0, "mem allocation error")
        ENDIF
      ENDIF

      IF (do_dom .OR. do_single_scatt) THEN
        ALLOCATE (trans_scatt_ir%opdpabs(0:su, nlayers, nchanprof), STAT=err)
        THROWM(err.NE.0, "mem allocation error")
        ALLOCATE (trans_scatt_ir%opdpsca(0:su, nlayers, nchanprof), STAT=err)
        THROWM(err.NE.0, "mem allocation error")
      ENDIF

      IF (do_dom) THEN
        IF (do_dom_vis) THEN
          ALLOCATE (trans_scatt_ir%ssa_solar(sl:su, nlayers, nchanprof), STAT=err)
          THROWM(err.NE.0, "mem allocation error")
          ALLOCATE (trans_scatt_ir%phup(0:su, nlayers, nchanprof), STAT=err)
          THROWM(err.NE.0, "mem allocation error")
          IF (direct1) THEN
            ALLOCATE (trans_scatt_ir%layerod_solar(0:su, nlayers, nchanprof), STAT=err)
            THROWM(err.NE.0, "mem allocation error")
          ENDIF
        ENDIF
        IF (do_dom_ir) THEN
          ALLOCATE (trans_scatt_ir%ssa_thermal(sl:su, nlayers, nchanprof), STAT=err)
          THROWM(err.NE.0, "mem allocation error")
          IF (direct1) THEN
            ALLOCATE (trans_scatt_ir%layerod_thermal(0:su, nlayers, nchanprof), STAT=err)
            THROWM(err.NE.0, "mem allocation error")
          ENDIF
        ENDIF
      ENDIF

      IF (do_single_scatt) THEN
        ALLOCATE (trans_scatt_ir%opdpext(0:su, nlayers, nchanprof), STAT=err)
        THROWM(err.NE.0, "mem allocation error")
        ALLOCATE (trans_scatt_ir%ssa_solar(0:su, nlayers, nchanprof), STAT=err)
        THROWM(err.NE.0, "mem allocation error")
        ALLOCATE (trans_scatt_ir%phup(0:su, nlayers, nchanprof), STAT=err)
        THROWM(err.NE.0, "mem allocation error")
        ALLOCATE (trans_scatt_ir%phdo(0:su, nlayers, nchanprof), STAT=err)
        THROWM(err.NE.0, "mem allocation error")
      ENDIF
    ENDIF

    IF (init1) CALL rttov_init_trans_scatt_ir(trans_scatt_ir)
  ENDIF

  IF (asw == 0) THEN
    IF (ASSOCIATED(trans_scatt_ir%phasefn)) THEN
      DO i = 1, nchanprof
        IF (ASSOCIATED(trans_scatt_ir%phasefn(i)%legcoef)) &
            DEALLOCATE(trans_scatt_ir%phasefn(i)%legcoef, STAT=err)
        THROWM(err.NE.0, "mem deallocation error")
      ENDDO
      DEALLOCATE(trans_scatt_ir%phasefn, STAT=err)
      THROWM(err.NE.0, "mem deallocation error")
    ENDIF

    IF (ASSOCIATED(trans_scatt_ir%ssa_solar)) DEALLOCATE(trans_scatt_ir%ssa_solar, STAT=err)
    THROWM(err.NE.0, "mem deallocation error")
    IF (ASSOCIATED(trans_scatt_ir%ssa_thermal)) DEALLOCATE(trans_scatt_ir%ssa_thermal, STAT=err)
    THROWM(err.NE.0, "mem deallocation error")
    IF (ASSOCIATED(trans_scatt_ir%layerod_solar)) DEALLOCATE(trans_scatt_ir%layerod_solar, STAT=err)
    THROWM(err.NE.0, "mem deallocation error")
    IF (ASSOCIATED(trans_scatt_ir%layerod_thermal)) DEALLOCATE(trans_scatt_ir%layerod_thermal, STAT=err)
    THROWM(err.NE.0, "mem deallocation error")
    IF (ASSOCIATED(trans_scatt_ir%opdpabs)) DEALLOCATE(trans_scatt_ir%opdpabs, STAT=err)
    THROWM(err.NE.0, "mem deallocation error")
    IF (ASSOCIATED(trans_scatt_ir%opdpsca)) DEALLOCATE(trans_scatt_ir%opdpsca, STAT=err)
    THROWM(err.NE.0, "mem deallocation error")
    IF (ASSOCIATED(trans_scatt_ir%opdpaclsun)) DEALLOCATE(trans_scatt_ir%opdpaclsun, STAT=err)
    THROWM(err.NE.0, "mem deallocation error")
    IF (ASSOCIATED(trans_scatt_ir%opdpacsun)) DEALLOCATE(trans_scatt_ir%opdpacsun, STAT=err)
    THROWM(err.NE.0, "mem deallocation error")
    IF (ASSOCIATED(trans_scatt_ir%opdpacl)) DEALLOCATE(trans_scatt_ir%opdpacl, STAT=err)
    THROWM(err.NE.0, "mem deallocation error")
    IF (ASSOCIATED(trans_scatt_ir%opdpac)) DEALLOCATE(trans_scatt_ir%opdpac, STAT=err)
    THROWM(err.NE.0, "mem deallocation error")
    IF (ASSOCIATED(trans_scatt_ir%opdpext)) DEALLOCATE(trans_scatt_ir%opdpext, STAT=err)
    THROWM(err.NE.0, "mem deallocation error")
    IF (ASSOCIATED(trans_scatt_ir%phup)) DEALLOCATE(trans_scatt_ir%phup, STAT=err)
    THROWM(err.NE.0, "mem deallocation error")
    IF (ASSOCIATED(trans_scatt_ir%phdo)) DEALLOCATE(trans_scatt_ir%phdo, STAT=err)
    THROWM(err.NE.0, "mem deallocation error")
    IF (ASSOCIATED(trans_scatt_ir%ray_sca)) DEALLOCATE(trans_scatt_ir%ray_sca, STAT=err)
    THROWM(err.NE.0, "mem deallocation error")

    IF (ASSOCIATED(trans_scatt_ir%aer)) CALL dealloc_scatt_ir_aercld(err, trans_scatt_ir%aer)
    THROW(err.NE.0)
    IF (ASSOCIATED(trans_scatt_ir%cld)) CALL dealloc_scatt_ir_aercld(err, trans_scatt_ir%cld)
    THROW(err.NE.0)

    CALL nullify_struct()
  ENDIF
  CATCH

CONTAINS

  SUBROUTINE alloc_scatt_ir_aercld(err, aercld, ntypes, user_opt_param, cld)
    INTEGER(jpim),                        INTENT(OUT)   :: err
    TYPE(rttov_scatt_ir_aercld),          POINTER       :: aercld
    INTEGER(jpim),                        INTENT(IN)    :: ntypes
    LOGICAL(jplm),                        INTENT(IN)    :: user_opt_param
    LOGICAL(jplm),                        INTENT(IN)    :: cld
    TRY
    ALLOCATE(aercld, STAT=err)
    THROWM(err.NE.0, "mem allocation error")

    CALL nullify_scatt_ir_aercld(aercld)

    ALLOCATE(aercld%opdpsca(nlayers,nchanprof), &
             aercld%opdpabs(nlayers,nchanprof), &
             aercld%opdp(nlayers,nchanprof), STAT=err)
    THROWM(err.NE.0, "mem allocation error")

    IF (do_chou) THEN
      ALLOCATE(aercld%opdpscabpr(nlayers,nchanprof), STAT=err)
      THROWM(err.NE.0, "mem allocation error")
    ENDIF

    IF (.NOT. user_opt_param .AND. cld) THEN
      ALLOCATE (aercld%partsca(ntypes,nlayers,nchanprof), STAT=err)
      THROWM(err.NE.0, "mem allocation error")
      IF (do_chou) THEN
        ALLOCATE (aercld%partbpr(ntypes,nlayers,nchanprof), STAT=err)
        THROWM(err.NE.0, "mem allocation error")
      ENDIF
    ENDIF

    IF (do_dom) THEN
      ALLOCATE (aercld%sca(nlayers,nchanprof), STAT=err)
      THROWM(err.NE.0, "mem allocation error")
    ENDIF

    IF (opts%rt_ir%addsolar .AND. .NOT. do_mfasis) THEN
      ALLOCATE(aercld%opdpsun(nlayers,nchanprof), STAT=err)
      THROWM(err.NE.0, "mem allocation error")
      IF (do_dom_vis .OR. do_single_scatt) THEN
        ALLOCATE(aercld%phtotup(nlayers,nchanprof), STAT=err)
        THROWM(err.NE.0, "mem allocation error")
        IF (direct1) THEN
          ALLOCATE(aercld%phintup(ntypes,nlayers,nchanprof), STAT=err)
          THROWM(err.NE.0, "mem allocation error")
        ENDIF
      ENDIF
      IF (do_single_scatt) THEN
        ALLOCATE(aercld%phtotdo(nlayers,nchanprof), STAT=err)
        THROWM(err.NE.0, "mem allocation error")
        IF (direct1) THEN
          ALLOCATE(aercld%phintdo(ntypes,nlayers,nchanprof), STAT=err)
          THROWM(err.NE.0, "mem allocation error")
        ENDIF
      ENDIF
    ENDIF
    CATCH
  END SUBROUTINE alloc_scatt_ir_aercld

  SUBROUTINE dealloc_scatt_ir_aercld(err, aercld)
    INTEGER(jpim),                        INTENT(OUT)   :: err
    TYPE(rttov_scatt_ir_aercld),          POINTER       :: aercld
    TRY
    IF (ASSOCIATED(aercld%opdpsca)) DEALLOCATE(aercld%opdpsca, STAT=err)
    THROWM(err.NE.0, "mem deallocation error")
    IF (ASSOCIATED(aercld%opdpabs)) DEALLOCATE(aercld%opdpabs, STAT=err)
    THROWM(err.NE.0, "mem deallocation error")
    IF (ASSOCIATED(aercld%opdpscabpr)) DEALLOCATE(aercld%opdpscabpr, STAT=err)
    THROWM(err.NE.0, "mem deallocation error")
    IF (ASSOCIATED(aercld%opdp)) DEALLOCATE(aercld%opdp, STAT=err)
    THROWM(err.NE.0, "mem deallocation error")
    IF (ASSOCIATED(aercld%opdpsun)) DEALLOCATE(aercld%opdpsun, STAT=err)
    THROWM(err.NE.0, "mem deallocation error")
    IF (ASSOCIATED(aercld%partsca)) DEALLOCATE(aercld%partsca, STAT=err)
    THROWM(err.NE.0, "mem deallocation error")
    IF (ASSOCIATED(aercld%partbpr)) DEALLOCATE(aercld%partbpr, STAT=err)
    THROWM(err.NE.0, "mem deallocation error")
    IF (ASSOCIATED(aercld%sca)) DEALLOCATE(aercld%sca, STAT=err)
    THROWM(err.NE.0, "mem deallocation error")
    IF (ASSOCIATED(aercld%phintup)) DEALLOCATE(aercld%phintup, STAT=err)
    THROWM(err.NE.0, "mem deallocation error")
    IF (ASSOCIATED(aercld%phintdo)) DEALLOCATE(aercld%phintdo, STAT=err)
    THROWM(err.NE.0, "mem deallocation error")
    IF (ASSOCIATED(aercld%phtotup)) DEALLOCATE(aercld%phtotup, STAT=err)
    THROWM(err.NE.0, "mem deallocation error")
    IF (ASSOCIATED(aercld%phtotdo)) DEALLOCATE(aercld%phtotdo, STAT=err)
    THROWM(err.NE.0, "mem deallocation error")
    DEALLOCATE(aercld)
    THROWM(err.NE.0, "mem deallocation error")
    CATCH
  END SUBROUTINE dealloc_scatt_ir_aercld

  SUBROUTINE nullify_scatt_ir_aercld(aercld)
    TYPE(rttov_scatt_ir_aercld), INTENT(INOUT) :: aercld
    NULLIFY(aercld%opdpsca, aercld%opdpabs, aercld%opdpscabpr, &
            aercld%opdp, aercld%opdpsun, aercld%sca, &
            aercld%partsca, aercld%partbpr, &
            aercld%phintup, aercld%phintdo, &
            aercld%phtotup, aercld%phtotdo)
  END SUBROUTINE nullify_scatt_ir_aercld

  SUBROUTINE nullify_struct()
    NULLIFY(trans_scatt_ir%aer, &
            trans_scatt_ir%cld, &
            trans_scatt_ir%phasefn, &
            trans_scatt_ir%ssa_solar, &
            trans_scatt_ir%ssa_thermal, &
            trans_scatt_ir%layerod_solar, &
            trans_scatt_ir%layerod_thermal, &
            trans_scatt_ir%ray_sca, &
            trans_scatt_ir%phup, &
            trans_scatt_ir%phdo, &
            trans_scatt_ir%opdpabs, &
            trans_scatt_ir%opdpsca, &
            trans_scatt_ir%opdpaclsun, &
            trans_scatt_ir%opdpacsun, &
            trans_scatt_ir%opdpacl, &
            trans_scatt_ir%opdpac, &
            trans_scatt_ir%opdpext)
  END SUBROUTINE nullify_struct

END SUBROUTINE rttov_alloc_trans_scatt_ir
