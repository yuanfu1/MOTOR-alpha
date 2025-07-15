! Description:
!> @file
!!   Reads BRDF atlas data. Data is loaded for the specified month.
!
!> @brief
!!   Reads BRDF atlas data. Data is loaded for the specified month.
!!
!! @details
!!   This subroutine initialises the BRDF atlas.
!!
!!   The coefs argument is used to pass the coefficients structure of
!!   an instrument that has been populated by calling rttov_read_coefs.
!!   This must be a v9 predictor coefficient file with visible/near-IR
!!   channels.
!!
!!   The coefs argument is optional: if it is omitted the initialised
!!   atlas can be used with any visible/IR sensor.
!!   If present the atlas data will be initialised for use with that
!!   specific instrument: the calls to rttov_get_brdf are then much
!!   faster. It is recommended to supply the coefs argument if you
!!   are only simulating one instrument and speed is important.
!!
!!   The atlas data for the specified month (and, where relevant,
!!   for the specified instrument) is read into the supplied
!!   atlas data structure argument. This allows you to initialise
!!   atlas data for multiple months (and instruments) simultaneously.
!!
!!   The atlas_id argument is used to choose between different
!!   atlases: currently there is only one BRDF atlas so this
!!   argument can be omitted.
!!
!! @param[out]    err        status on exit
!! @param[in]     opts       options to configure the simulations
!! @param[in]     imonth     month (1-12 => Jan-Dec) for which atlas should be initialised
!! @param[in,out] atlas      atlas data structure to initialise
!! @param[in]     atlas_id   ID of atlas to initialise, optional
!! @param[in]     path       path to directory containing BRDF atlas data, optional
!! @param[in]     coefs      coefficients structure for instrument, if present atlas is
!!                             initialiased specifically for this instrument, optional
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
SUBROUTINE rttov_setup_brdf_atlas( &
                err,      &! out
                opts,     &! in
                imonth,   &! in
                atlas,    &! inout
                atlas_id, &! in, optional
                path,     &! in, optional
                coefs)     ! in, optional
!INTF_OFF
#include "throw.h"
!INTF_ON
  USE parkind1, ONLY : jpim
  USE rttov_types, ONLY : rttov_coefs, rttov_options
  USE mod_rttov_brdf_atlas, ONLY : rttov_brdf_atlas_data
!INTF_OFF
  USE parkind1, ONLY : jplm
  USE rttov_const, ONLY : sensor_id_mw, sensor_id_po
  USE mod_rttov_brdf_atlas, ONLY : brdf_atlas_id
  USE mod_brdf_atlas, ONLY : rttov_visnirbrdf_init
!INTF_ON

  IMPLICIT NONE

  INTEGER(KIND=jpim),          INTENT(OUT)             :: err
  TYPE(rttov_options),         INTENT(IN)              :: opts
  INTEGER(KIND=jpim),          INTENT(IN)              :: imonth
  TYPE(rttov_brdf_atlas_data), INTENT(INOUT)           :: atlas
  INTEGER(KIND=jpim),          INTENT(IN),    OPTIONAL :: atlas_id
  CHARACTER(LEN=*),            INTENT(IN),    OPTIONAL :: path
  TYPE(rttov_coefs),           INTENT(IN),    OPTIONAL :: coefs
!INTF_END

#include "rttov_errorreport.interface"

  CHARACTER(LEN=300) :: fpath
  INTEGER(KIND=jpim) :: id
  LOGICAL(KIND=jplm) :: single_inst

!----------------------------------------------------------------------------
  TRY

  IF (imonth < 0 .OR. imonth > 12) THEN
    err = errorstatus_fatal
    THROWM(err .NE. 0, 'Month argument must lie between 1 and 12')
  ENDIF

  IF (atlas%init) THEN
    err = errorstatus_fatal
    THROWM(err .NE. 0, 'Atlas data already initialised: deallocate first')
  ENDIF

  fpath = './'
  IF (PRESENT(path)) fpath = TRIM(path)//'/'

  single_inst = PRESENT(coefs)

  id = -1
  IF (PRESENT(atlas_id)) id = atlas_id

  IF (id <= 0) id = brdf_atlas_id

  IF (id == brdf_atlas_id) THEN
    IF (single_inst) THEN
      IF (coefs%coef%id_sensor == sensor_id_mw .OR. &
          coefs%coef%id_sensor == sensor_id_po .OR. &
          ALL(coefs%coef%ss_val_chn(:) == 0_jpim)) THEN
        err = errorstatus_fatal
        THROWM(err .NE. 0, 'BRDF atlas only valid for sensors with VIS/NIR channels')
      ENDIF

      CALL rttov_visnirbrdf_init(TRIM(fpath),            &
                                 imonth,                 &
                                 opts%config%verbose,    &
                                 atlas%brdf_atlas,       &
                                 err,                    &
                                 coefs%coef%ff_cwn(:),   &
                                 coefs%coef%id_platform, &
                                 coefs%coef%id_sat,      &
                                 coefs%coef%id_inst)
    ELSE
      CALL rttov_visnirbrdf_init(TRIM(fpath),            &
                                 imonth,                 &
                                 opts%config%verbose,    &
                                 atlas%brdf_atlas,       &
                                 err)
    ENDIF
    THROWM(err .NE. 0, 'Error initialising BRDF atlas')
  ELSE
    err = errorstatus_fatal
    THROWM(err .NE. 0, 'Unknown BRDF atlas ID')
  ENDIF

  atlas%atlas_id = id
  atlas%init = .TRUE.

  CATCH

END SUBROUTINE rttov_setup_brdf_atlas
