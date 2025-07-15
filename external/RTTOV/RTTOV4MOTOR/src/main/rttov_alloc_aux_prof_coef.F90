! Description:
!> @file
!!   Allocate/deallocate internal coef-level auxiliary profiles structure.
!
!> @brief
!!   Allocate/deallocate internal coef-level auxiliary profiles structure.
!!
!! @details
!!   This structure contains results of various internal calculations
!!   related to the profiles on coefficient levels:
!!   - information related to the near-surface level/layer
!!   - information related to the near-cloud level/layer (for the simple-
!!     cloud scheme)
!!   - quantities used in predictor calculations
!!
!! @param[out]    err               status on exit
!! @param[in]     nprofiles         number of profiles being simulated
!! @param[in]     nlevels           number of levels in input profiles
!! @param[in,out] aux_prof          coef-level auxiliary profiles structure to (de)allocate
!! @param[in]     coef              optical depth coefficients structure
!! @param[in]     asw               1_jpim => allocate; 0_jpim => deallocate
!! @param[in]     init              set .TRUE. to initialise newly allocated structures, optional
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
SUBROUTINE rttov_alloc_aux_prof_coef( &
              err,       &
              nprofiles, &
              nlevels,   &
              aux_prof,  &
              coef,      &
              asw,       &
              init)
!INTF_OFF
#include "throw.h"
!INTF_ON
  USE rttov_types, ONLY : rttov_profile_aux_coef, rttov_coef
  USE parkind1, ONLY : jpim, jplm
  IMPLICIT NONE

  INTEGER(KIND=jpim),           INTENT(OUT)            :: err
  INTEGER(KIND=jpim),           INTENT(IN)             :: nprofiles
  INTEGER(KIND=jpim),           INTENT(IN)             :: nlevels
  TYPE(rttov_profile_aux_coef), INTENT(INOUT)          :: aux_prof
  TYPE(rttov_coef),             INTENT(IN)             :: coef
  INTEGER(KIND=jpim),           INTENT(IN)             :: asw
  LOGICAL(KIND=jplm),           INTENT(IN),   OPTIONAL :: init
!INTF_END

#include "rttov_errorreport.interface"
#include "rttov_init_aux_prof_coef.interface"

  LOGICAL(KIND=jplm) :: init1
  INTEGER(KIND=jpim) :: nlayers

  TRY
  nlayers = nlevels - 1
  init1 = .FALSE.
  IF (PRESENT(init)) init1 = init

  IF (asw == 1) THEN
    CALL nullify_struct()

    ALLOCATE (aux_prof%s(nprofiles), STAT = err)
    THROWM(err .NE. 0, "allocation of aux_prof%s")

    ALLOCATE(aux_prof%pathsat_rsqrt(nlayers, nprofiles), &
             aux_prof%pathsat_sqrt(nlayers, nprofiles),  &
             aux_prof%pathsat_4rt(nlayers, nprofiles), STAT = err)
    THROWM(err .NE. 0, "allocation of aux_prof path quantities")

    ALLOCATE( &
      aux_prof%t_layer(nlayers, nprofiles), aux_prof%w_layer(nlayers, nprofiles), aux_prof%o3_layer(nlayers, nprofiles), &
      aux_prof%dt(nlayers, nprofiles), aux_prof%dtabsdt(nlayers, nprofiles), aux_prof%dto(nlayers, nprofiles), &
      aux_prof%tr(nlayers, nprofiles), aux_prof%tr2(nlayers, nprofiles), &
      aux_prof%tr_r(nlayers, nprofiles), aux_prof%tr_sqrt(nlayers, nprofiles), &
      aux_prof%wr_sqrt(nlayers, nprofiles), aux_prof%wrw_r(nlayers, nprofiles), &
      aux_prof%wr_4rt(nlayers, nprofiles), aux_prof%wr_rsqrt(nlayers, nprofiles), &
      aux_prof%ww_r(nlayers, nprofiles), &
      aux_prof%ow_rsqrt(nlayers, nprofiles), aux_prof%ow_r(nlayers, nprofiles),& 
      aux_prof%or_sqrt(nlayers, nprofiles), aux_prof%ow_sqrt(nlayers, nprofiles), aux_prof%ow_4rt(nlayers, nprofiles), &
      aux_prof%wr(nlayers, nprofiles), aux_prof%or(nlayers, nprofiles), &
      aux_prof%ww(nlayers, nprofiles), &
      aux_prof%ow(nlayers, nprofiles), aux_prof%co2_cm(nprofiles), STAT = err)
    THROWM(err .NE. 0, "allocation of aux_prof layer quantities")

    IF (coef%fmv_model_ver <= 8) THEN
      ALLOCATE(aux_prof%tw(nlayers, nprofiles), aux_prof%tw_sqrt(nlayers, nprofiles), &
               aux_prof%tw_4rt(nlayers, nprofiles), STAT=err)
      THROWM(err .NE. 0, "allocation of aux_prof layer quantities")
    ENDIF

    IF (coef%fmv_model_ver == 9) THEN
      ALLOCATE(aux_prof%tuw(nlayers, nprofiles), STAT=err)
      THROWM(err .NE. 0, "allocation of aux_prof layer quantities")
    ENDIF

    IF (coef%fmv_model_ver >= 8) THEN
      ALLOCATE( &
        aux_prof%twr(nlayers, nprofiles), aux_prof%wwr(nlayers, nprofiles), &
        aux_prof%wwr_r(nlayers, nprofiles), aux_prof%wrwr_r(nlayers, nprofiles), &
        STAT=err)
      THROWM(err .NE. 0, "allocation of aux_prof layer quantities")
  
      IF (coef%nco2 > 0) THEN
        ALLOCATE( &
          aux_prof%co2_layer(nlayers, nprofiles), aux_prof%co2r(nlayers, nprofiles), &
          aux_prof%co2r_sqrt(nlayers, nprofiles), aux_prof%co2w(nlayers, nprofiles), STAT = err)
        THROWM(err .NE. 0, "allocation of aux_prof layer quantities")
      ENDIF
    ENDIF

    IF (coef%fmv_model_ver >= 9) THEN
      ALLOCATE(aux_prof%patheff_rsqrt(nlayers, nprofiles), &
               aux_prof%patheff_sqrt(nlayers, nprofiles),  &
               aux_prof%patheff_4rt(nlayers, nprofiles), STAT = err)
      THROWM(err .NE. 0, "allocation of aux_prof path quantities")

      ALLOCATE( &
        aux_prof%ww_sqrt(nlayers, nprofiles), aux_prof%ww_rsqrt(nlayers, nprofiles), &
        aux_prof%ww_4rt(nlayers, nprofiles), aux_prof%tro(nlayers, nprofiles), STAT=err)
      THROWM(err .NE. 0, "allocation of aux_prof layer quantities")

      IF (coef%nn2o > 0) THEN
        ALLOCATE( &
          aux_prof%n2o_layer(nlayers, nprofiles), aux_prof%n2or(nlayers, nprofiles), &
          aux_prof%n2or_sqrt(nlayers, nprofiles), aux_prof%n2or_4rt(nlayers, nprofiles), &
          aux_prof%n2ow(nlayers, nprofiles), aux_prof%n2owr(nlayers, nprofiles), &
          aux_prof%n2ow_r(nlayers, nprofiles), STAT = err)
        THROWM(err .NE. 0, "allocation of aux_prof layer quantities")
      ENDIF

      IF (coef%nco > 0) THEN
        ALLOCATE( &
          aux_prof%co_layer(nlayers, nprofiles), aux_prof%cor(nlayers, nprofiles), &
          aux_prof%cor_sqrt(nlayers, nprofiles), aux_prof%cor_4rt(nlayers, nprofiles), &
          aux_prof%corw_rsqrt(nlayers, nprofiles), aux_prof%corw_r(nlayers, nprofiles), &
          aux_prof%corw_r4rt(nlayers, nprofiles), aux_prof%cow(nlayers, nprofiles), &
          aux_prof%cow_rsqrt(nlayers, nprofiles), aux_prof%cow_r4rt(nlayers, nprofiles), &
          aux_prof%cowr(nlayers, nprofiles), aux_prof%cowr_r(nlayers, nprofiles), &
          aux_prof%cowr_4rt(nlayers, nprofiles), STAT = err)
        THROWM(err .NE. 0, "allocation of aux_prof layer quantities")
      ENDIF

      IF (coef%nch4 > 0) THEN
        ALLOCATE( &
          aux_prof%ch4_layer(nlayers, nprofiles), aux_prof%ch4r(nlayers, nprofiles), &
          aux_prof%ch4r_sqrt(nlayers, nprofiles), aux_prof%ch4r_4rt(nlayers, nprofiles), &
          aux_prof%ch4w(nlayers, nprofiles), aux_prof%ch4w_4rt(nlayers, nprofiles), &
          aux_prof%ch4w_r(nlayers, nprofiles), aux_prof%ch4wr(nlayers, nprofiles), &
          aux_prof%ch4rw_r(nlayers, nprofiles), STAT = err)
        THROWM(err .NE. 0, "allocation of aux_prof layer quantities")
      ENDIF

      IF (coef%nso2 > 0) THEN
        ALLOCATE( &
          aux_prof%so2_layer(nlayers, nprofiles), aux_prof%so2r(nlayers, nprofiles), &
          aux_prof%so2r_sqrt(nlayers, nprofiles), aux_prof%so2r_4rt(nlayers, nprofiles), &
          aux_prof%so2w(nlayers, nprofiles), aux_prof%so2w_r(nlayers, nprofiles), &
          aux_prof%so2wr(nlayers, nprofiles), aux_prof%so2wr_r(nlayers, nprofiles), &
          aux_prof%so2w_sqrt(nlayers, nprofiles), aux_prof%so2w_4rt(nlayers, nprofiles), &
          aux_prof%so2rw_r(nlayers, nprofiles), aux_prof%so2rwr_r(nlayers, nprofiles), &
          STAT = err)
        THROWM(err .NE. 0, "allocation of aux_prof layer quantities")
      ENDIF
    ENDIF

    IF (init1) CALL rttov_init_aux_prof_coef(aux_prof)
  ENDIF

  IF (asw == 0) THEN
    DEALLOCATE (aux_prof%s, STAT = err)
    THROWM(err .NE. 0, "deallocation of aux_prof%s")

    DEALLOCATE(aux_prof%pathsat_rsqrt, aux_prof%pathsat_sqrt, aux_prof%pathsat_4rt, STAT = err)
    THROWM(err .NE. 0, "deallocation of aux_prof path quantities")

    DEALLOCATE(&
      aux_prof%t_layer, aux_prof%w_layer, aux_prof%o3_layer, &
      aux_prof%dt, aux_prof%dtabsdt, aux_prof%dto, &
      aux_prof%tr, aux_prof%tr2, aux_prof%tr_r, aux_prof%tr_sqrt, &
      aux_prof%wr_sqrt, aux_prof%wrw_r, &
      aux_prof%wr_4rt, aux_prof%wr_rsqrt, &
      aux_prof%ww_r, &
      aux_prof%ow_rsqrt, aux_prof%ow_r, &
      aux_prof%or_sqrt, aux_prof%ow_sqrt, aux_prof%ow_4rt, &
      aux_prof%wr, aux_prof%or, &
      aux_prof%ww, &
      aux_prof%ow, aux_prof%co2_cm, STAT = err)
    THROWM(err .NE. 0, "deallocation of aux_prof layer quantities")

    IF (coef%fmv_model_ver <= 8) THEN
      DEALLOCATE(aux_prof%tw, aux_prof%tw_sqrt, aux_prof%tw_4rt, STAT=err)
      THROWM(err .NE. 0, "deallocation of aux_prof layer quantities")
    ENDIF

    IF (coef%fmv_model_ver == 9) THEN
      DEALLOCATE(aux_prof%tuw, STAT=err)
      THROWM(err .NE. 0, "deallocation of aux_prof layer quantities")
    ENDIF

    IF (coef%fmv_model_ver >= 8) THEN
      DEALLOCATE(aux_prof%twr, aux_prof%wwr, aux_prof%wwr_r, aux_prof%wrwr_r, STAT = err)
      THROWM(err .NE. 0, "deallocation of aux_prof layer quantities")

      IF (coef%nco2 > 0) THEN
        DEALLOCATE(aux_prof%co2_layer, aux_prof%co2r, aux_prof%co2r_sqrt, aux_prof%co2w, STAT = err)
        THROWM(err .NE. 0, "deallocation of aux_prof layer quantities")
      ENDIF
    ENDIF

    IF (coef%fmv_model_ver >= 9) THEN
      DEALLOCATE(aux_prof%patheff_rsqrt, aux_prof%patheff_sqrt, aux_prof%patheff_4rt, STAT = err)
      THROWM(err .NE. 0, "deallocation of aux_prof path quantities")

      DEALLOCATE(aux_prof%ww_sqrt, aux_prof%ww_rsqrt, aux_prof%ww_4rt, aux_prof%tro, STAT = err)
      THROWM(err .NE. 0, "deallocation of aux_prof layer quantities")

      IF (coef%nn2o > 0) THEN
        DEALLOCATE( &
          aux_prof%n2o_layer, aux_prof%n2or, &
          aux_prof%n2or_sqrt, aux_prof%n2or_4rt, &
          aux_prof%n2ow, aux_prof%n2owr, aux_prof%n2ow_r, STAT = err)
        THROWM(err .NE. 0, "deallocation of aux_prof layer quantities")
      ENDIF

      IF (coef%nco > 0) THEN
        DEALLOCATE( &
          aux_prof%co_layer, aux_prof%cor, &
          aux_prof%cor_sqrt, aux_prof%cor_4rt, &
          aux_prof%corw_rsqrt, aux_prof%corw_r, &
          aux_prof%corw_r4rt, aux_prof%cow, &
          aux_prof%cow_rsqrt, aux_prof%cow_r4rt, &
          aux_prof%cowr, aux_prof%cowr_r, &
          aux_prof%cowr_4rt, STAT = err)
        THROWM(err .NE. 0, "deallocation of aux_prof layer quantities")
      ENDIF

      IF (coef%nch4 > 0) THEN
        DEALLOCATE( &
          aux_prof%ch4_layer, aux_prof%ch4r, &
          aux_prof%ch4r_sqrt, aux_prof%ch4r_4rt, &
          aux_prof%ch4w, aux_prof%ch4w_4rt, &
          aux_prof%ch4w_r, aux_prof%ch4wr, aux_prof%ch4rw_r, STAT = err)
        THROWM(err .NE. 0, "deallocation of aux_prof layer quantities")
      ENDIF

      IF (coef%nso2 > 0) THEN
        DEALLOCATE( &
          aux_prof%so2_layer, aux_prof%so2r, &
          aux_prof%so2r_sqrt, aux_prof%so2r_4rt, &
          aux_prof%so2w, aux_prof%so2w_r, &
          aux_prof%so2wr, aux_prof%so2wr_r, &
          aux_prof%so2w_sqrt, aux_prof%so2w_4rt, &
          aux_prof%so2rw_r, aux_prof%so2rwr_r, &
          STAT = err)
        THROWM(err .NE. 0, "deallocation of aux_prof layer quantities")
      ENDIF

    ENDIF

    CALL nullify_struct()
  ENDIF
  CATCH
CONTAINS
  SUBROUTINE nullify_struct()
    NULLIFY (aux_prof%s)
    NULLIFY(aux_prof%pathsat_rsqrt, aux_prof%pathsat_sqrt, aux_prof%pathsat_4rt, &
            aux_prof%patheff_rsqrt, aux_prof%patheff_sqrt, aux_prof%patheff_4rt)
    NULLIFY(&
      aux_prof%t_layer, aux_prof%w_layer, aux_prof%o3_layer, &
      aux_prof%dt, aux_prof%dtabsdt, aux_prof%dto, &
      aux_prof%tr, aux_prof%tr2, aux_prof%tr_r, aux_prof%tr_sqrt, &
      aux_prof%tw_sqrt, aux_prof%wr_sqrt, aux_prof%wrw_r, &
      aux_prof%wr_4rt, aux_prof%wr_rsqrt, &
      aux_prof%tw_4rt, aux_prof%ww_r, &
      aux_prof%ow_rsqrt, aux_prof%ow_r, &
      aux_prof%or_sqrt, aux_prof%ow_sqrt, aux_prof%ow_4rt, &
      aux_prof%wr, aux_prof%or, &
      aux_prof%tw, aux_prof%ww, &
      aux_prof%ow, aux_prof%twr, aux_prof%tuw, &
      aux_prof%co2_cm, &
      aux_prof%wwr, aux_prof%wwr_r, aux_prof%wrwr_r, aux_prof%ww_sqrt, aux_prof%ww_rsqrt, &
      aux_prof%co2_layer, aux_prof%co2r, aux_prof%co2r_sqrt, aux_prof%co2w, &
      aux_prof%ww_4rt, aux_prof%tro, &
      aux_prof%n2o_layer, aux_prof%n2or, &
      aux_prof%n2or_sqrt, aux_prof%n2or_4rt, &
      aux_prof%n2ow, aux_prof%n2owr, aux_prof%n2ow_r, &
      aux_prof%co_layer, aux_prof%cor, &
      aux_prof%cor_sqrt, aux_prof%cor_4rt, &
      aux_prof%corw_rsqrt, aux_prof%corw_r, &
      aux_prof%corw_r4rt, aux_prof%cow, &
      aux_prof%cow_rsqrt, aux_prof%cow_r4rt, &
      aux_prof%cowr, aux_prof%cowr_r, &
      aux_prof%cowr_4rt, &
      aux_prof%ch4_layer, aux_prof%ch4r, &
      aux_prof%ch4r_sqrt, aux_prof%ch4r_4rt, &
      aux_prof%ch4w, aux_prof%ch4w_4rt, &
      aux_prof%ch4w_r, aux_prof%ch4wr, aux_prof%ch4rw_r, &
      aux_prof%so2_layer, aux_prof%so2r, &
      aux_prof%so2r_sqrt, aux_prof%so2r_4rt, &
      aux_prof%so2w, aux_prof%so2w_r, &
      aux_prof%so2wr, aux_prof%so2wr_r, &
      aux_prof%so2w_sqrt, aux_prof%so2w_4rt, &
      aux_prof%so2rw_r, aux_prof%so2rwr_r)
  END SUBROUTINE nullify_struct
END SUBROUTINE rttov_alloc_aux_prof_coef
