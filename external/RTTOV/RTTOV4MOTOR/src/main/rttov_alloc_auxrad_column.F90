! Description:
!> @file
!!   Allocate/deallocate internal dyanmically-sized auxiliary radiance
!!   structure.
!
!> @brief
!!   Allocate/deallocate internal dyanmically-sized auxiliary radiance
!!   structure.
!!
!! @details
!!   This structure holds calculated radiances where arrays are sized
!!   by the number of cloud columns which is determined at runtime.
!!
!!   Some arrays are only required by the direct model instance
!!   of this structure: these are allocated if the optional "direct"
!!   argument is TRUE. If omitted, it is assumed to be FALSE.
!!
!! @param[out]    err                    status on exit
!! @param[in,out] auxrad_column          auxiliary radiance structure to (de)allocate
!! @param[in]     opts                   options to configure the simulations
!! @param[in]     ncolumns               number of cloud columns as calculated by rttov_cloud_overlap
!! @param[in]     nlayers                number of input profile layers (nlevels-1)
!! @param[in]     nchanprof              total number of channels being simulated (SIZE(chanprof))
!! @param[in]     asw                    1_jpim => allocate; 0_jpim => deallocate
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
SUBROUTINE rttov_alloc_auxrad_column( &
              err,           &
              auxrad_column, &
              opts,          &
              ncolumns,      &
              nlayers,       &
              nchanprof,     &
              asw,           &
              init,          &
              direct)

  USE parkind1, ONLY : jpim, jplm
  USE rttov_types, ONLY : rttov_options, rttov_radiance_aux
!INTF_OFF
  USE rttov_const, ONLY : vis_scatt_single, vis_scatt_mfasis, ir_scatt_dom
#include "throw.h"
!INTF_ON
  IMPLICIT NONE
  INTEGER(KIND=jpim),       INTENT(OUT)          :: err
  TYPE(rttov_radiance_aux), INTENT(INOUT)        :: auxrad_column
  TYPE(rttov_options),      INTENT(IN)           :: opts
  INTEGER(KIND=jpim),       INTENT(IN)           :: ncolumns
  INTEGER(KIND=jpim),       INTENT(IN)           :: nlayers
  INTEGER(KIND=jpim),       INTENT(IN)           :: nchanprof
  INTEGER(KIND=jpim),       INTENT(IN)           :: asw
  LOGICAL(KIND=jplm),       INTENT(IN), OPTIONAL :: init
  LOGICAL(KIND=jplm),       INTENT(IN), OPTIONAL :: direct
!INTF_END
#include "rttov_errorreport.interface"
#include "rttov_init_auxrad_column.interface"
  LOGICAL(KIND=jplm) :: init1, direct1
  TRY
  init1 = .FALSE.
  IF (PRESENT(init)) init1 = init
  direct1 = .FALSE.
  IF (PRESENT(direct)) direct1 = direct

  IF (asw == 1) THEN
    CALL nullify_struct()

    ALLOCATE(auxrad_column%up(nlayers, 0:ncolumns, nchanprof),       &
             auxrad_column%down(nlayers, 0:ncolumns, nchanprof),     &
             auxrad_column%cloudy(0:ncolumns, nchanprof),            &
             auxrad_column%meanrad_up(0:ncolumns, nchanprof),        &
             auxrad_column%meanrad_down(0:ncolumns, nchanprof), STAT = err)
    THROWM(err .NE. 0 , "allocation of auxrad_column")

    IF (direct1) THEN
      ALLOCATE(auxrad_column%down_ref(nlayers, 0:ncolumns, nchanprof), STAT=err)
      THROWM(err .NE. 0 , "allocation of auxrad_column")
    ENDIF

    IF (opts%rt_all%do_lambertian .AND. &
        .NOT. ((opts%rt_ir%addaerosl .OR. opts%rt_ir%addclouds) .AND. &
               opts%rt_ir%ir_scatt_model == ir_scatt_dom)) THEN
      ALLOCATE(auxrad_column%down_p(nlayers, 0:ncolumns, nchanprof),     &
               auxrad_column%meanrad_down_p(0:ncolumns, nchanprof), STAT = err)
      THROWM(err .NE. 0 , "allocation of auxrad_column")

      IF (direct1) THEN
        ALLOCATE(auxrad_column%down_p_ref(nlayers, 0:ncolumns, nchanprof), STAT=err)
        THROWM(err .NE. 0 , "allocation of auxrad_column")
      ENDIF
    ENDIF

    IF (opts%rt_ir%addsolar .AND. .NOT. &
        ((opts%rt_ir%addaerosl .OR. opts%rt_ir%addclouds) .AND. &
         opts%rt_ir%vis_scatt_model == vis_scatt_mfasis)) THEN
      ALLOCATE(auxrad_column%up_solar(nlayers, 0:ncolumns, nchanprof),   &
               auxrad_column%meanrad_up_solar(0:ncolumns, nchanprof),    &
               auxrad_column%down_solar(nlayers, 0:ncolumns, nchanprof), &
               auxrad_column%meanrad_down_solar(0:ncolumns, nchanprof), STAT = err)
      THROWM(err .NE. 0 , "allocation of auxrad_column solar arrays")

      IF (direct1) THEN
        ALLOCATE(auxrad_column%down_ref_solar(nlayers, 0:ncolumns, nchanprof), STAT=err)
        THROWM(err .NE. 0 , "allocation of auxrad_column solar arrays")
      ENDIF

      IF (direct1 .AND. (opts%rt_ir%addaerosl .OR. opts%rt_ir%addclouds) .AND. &
          opts%rt_ir%vis_scatt_model == vis_scatt_single) THEN
        ALLOCATE(auxrad_column%fac1_2(0:1, nlayers, nchanprof), &
                 auxrad_column%fac3_2(nlayers, nchanprof),      &
                 auxrad_column%fac4_2(0:1, nlayers, nchanprof), &
                 auxrad_column%fac5_2(0:1, nlayers, nchanprof), &
                 auxrad_column%fac6_2(0:1, nlayers, nchanprof), &
                 auxrad_column%fac7_2(nlayers, nchanprof),      &
                 auxrad_column%fac4_3(0:ncolumns, nchanprof),   &
                 auxrad_column%fac5_3(0:ncolumns, nchanprof), STAT = err)
        THROWM(err .NE. 0 , "allocation of auxrad_column%fac")
      ENDIF
    ENDIF
    IF (init1) CALL rttov_init_auxrad_column(auxrad_column)
  ENDIF

  IF (asw == 0) THEN
    DEALLOCATE(auxrad_column%up,         &
               auxrad_column%down,       &
               auxrad_column%cloudy,     &
               auxrad_column%meanrad_up, &
               auxrad_column%meanrad_down, STAT = err)
    THROWM(err .NE. 0 , "deallocation of auxrad_column")

    IF (direct1) THEN
      DEALLOCATE(auxrad_column%down_ref, STAT=err)
      THROWM(err .NE. 0 , "deallocation of auxrad_column")
    ENDIF

    IF (opts%rt_all%do_lambertian .AND. &
        .NOT. ((opts%rt_ir%addaerosl .OR. opts%rt_ir%addclouds) .AND. &
               opts%rt_ir%ir_scatt_model == ir_scatt_dom)) THEN
      DEALLOCATE(auxrad_column%down_p,     &
                 auxrad_column%meanrad_down_p, STAT = err)
      THROWM(err .NE. 0 , "deallocation of auxrad_column")
      IF (direct1) THEN
        DEALLOCATE(auxrad_column%down_p_ref, STAT=err)
        THROWM(err .NE. 0 , "deallocation of auxrad_column")
      ENDIF
    ENDIF

    IF (opts%rt_ir%addsolar .AND. .NOT. &
        ((opts%rt_ir%addaerosl .OR. opts%rt_ir%addclouds) .AND. &
         opts%rt_ir%vis_scatt_model == vis_scatt_mfasis)) THEN
      DEALLOCATE(auxrad_column%up_solar,           &
                 auxrad_column%meanrad_up_solar,   &
                 auxrad_column%down_solar,         &
                 auxrad_column%meanrad_down_solar, STAT = err)
      THROWM(err .NE. 0 , "deallocation of auxrad_column solar arrays")

      IF (direct1) THEN
        DEALLOCATE(auxrad_column%down_ref_solar, STAT=err)
        THROWM(err .NE. 0 , "deallocation of auxrad_column solar arrays")
      ENDIF

      IF (direct1 .AND. (opts%rt_ir%addaerosl .OR. opts%rt_ir%addclouds) .AND. &
          opts%rt_ir%vis_scatt_model == vis_scatt_single) THEN
        DEALLOCATE (auxrad_column%fac1_2, &
                    auxrad_column%fac3_2, &
                    auxrad_column%fac4_2, &
                    auxrad_column%fac5_2, &
                    auxrad_column%fac6_2, &
                    auxrad_column%fac7_2, &
                    auxrad_column%fac4_3, &
                    auxrad_column%fac5_3, STAT = err)
        THROWM(err .NE. 0 , "deallocation of auxrad_column%fac")
      ENDIF
    ENDIF

    CALL nullify_struct()
  ENDIF
  CATCH
CONTAINS
  SUBROUTINE nullify_struct()
    NULLIFY (auxrad_column%up,                 &
             auxrad_column%down,               &
             auxrad_column%down_ref,           &
             auxrad_column%down_p,             &
             auxrad_column%down_p_ref,         &
             auxrad_column%cloudy,             &
             auxrad_column%meanrad_up,         &
             auxrad_column%meanrad_down,       &
             auxrad_column%meanrad_down_p,     &
             auxrad_column%up_solar,           &
             auxrad_column%down_solar,         &
             auxrad_column%meanrad_up_solar,   &
             auxrad_column%meanrad_down_solar, &
             auxrad_column%down_ref_solar,     &
             auxrad_column%fac1_2,             &
             auxrad_column%fac3_2,             &
             auxrad_column%fac4_2,             &
             auxrad_column%fac5_2,             &
             auxrad_column%fac6_2,             &
             auxrad_column%fac7_2,             &
             auxrad_column%fac4_3,             &
             auxrad_column%fac5_3)
  END SUBROUTINE nullify_struct
END SUBROUTINE 
