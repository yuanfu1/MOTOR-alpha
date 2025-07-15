! Description:
!> @file
!!   Reads emissivity atlas data appropriate to the sensor
!!   type. Data is loaded for the specified month.
!
!> @brief
!!   Reads emissivity atlas data appropriate to the sensor
!!   type. Data is loaded for the specified month.
!!
!! @details
!!   This subroutine is called to initialise an IR or
!!   MW emissivity atlas: the type of atlas is specified
!!   in atlas_type (1 => MW, 2 => IR).
!!
!!   There are three IR atlases: the UW IR atlas and the CAMEL
!!   2007 and climatology atlases.
!!
!!   There are two MW atlases: TELSEM2 and the CNRM MW atlas.
!!
!!   The atlas to initialise in each case is specified by the
!!   atlas_id argument. If this is omitted the default atlas
!!   (UW IR atlas for IR sensors, TELSEM2 for MW sensors) is
!!   loaded.
!!
!!   The coefs argument is used to pass the coefficients structure of
!!   an instrument that has been populated by calling rttov_read_coefs.
!!
!!   This argument is mandatory for the CNRM MW atlas as this is always
!!   initialised for a specific instrument. The coefs argument is ignored
!!   for the TELSEM2 MW atlas as this can be used for any MW instrument.
!!
!!   For the UW and CAMEL IR atlases the coefs argument is optional.
!!   If omitted the initialised atlas may be used with any IR sensor.
!!   If present the atlas data will be initialised for use with that
!!   specific instrument: the calls to rttov_get_emis are then much
!!   faster. It is recommended to supply the coefs argument if you
!!   are only simulating one instrument and speed is important.
!!
!!   The atlas data for the specified month (and, where relevant,
!!   for the specified instrument) is read into the supplied
!!   atlas data structure argument. This allows you to initialise
!!   atlas data for multiple months (and instruments) simultaneously.
!!
!!   Data from different years are available for the CNRM MW atlas:
!!   the year of the data is indicated in the filenames and this is
!!   selected via the year argument.
!!
!!   There are additional options for the IR emissivity atlases:
!!   - the standard deviation information can optionally be
!!     returned by the atlas: do not select this if it is not
!!     required as the memory requirements are relatively large
!!   - the IR emissivites can optionally include a correction
!!     for the zenith angle. If this option is selected you
!!     must download the additional angular correction atlas
!!     files from the RTTOV website.
!!
!! @param[out]    err                         status on exit
!! @param[in]     opts                        options to configure the simulations
!! @param[in]     imonth                      month (1-12 => Jan-Dec) for which atlas should be initialised
!! @param[in]     atlas_type                  1 => MW atlas, 2 => IR atlas
!! @param[in,out] atlas                       atlas data structure to initialise
!! @param[in]     atlas_id                    ID of atlas to initialise, optional
!! @param[in]     path                        path to directory containing emissivity atlas data, optional
!! @param[in]     coefs                       instrument coefficients structure: required for CNRM atlas, ignored
!!                                              by TELSEM2, optional for IR atlases, if present initialises IR atlas
!!                                              for this specific instrument
!! @param[in]     ir_atlas_read_std           return standard deviations for IR emissivity atlas, optional
!! @param[in]     ir_atlas_ang_corr           apply zenith angle correction to IR emissivities, optional
!! @param[in]     year                        use atlas data for this year, CNRM atlas only, optional
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
SUBROUTINE rttov_setup_emis_atlas(          &
                err,                        &! out
                opts,                       &! in
                imonth,                     &! in
                atlas_type,                 &! in
                atlas,                      &! inout
                atlas_id,                   &! in, optional
                path,                       &! in, optional
                coefs,                      &! in, optional
                ir_atlas_read_std,          &! in, optional
                ir_atlas_ang_corr,          &! in, optional
                year)                        ! in, optional
!INTF_OFF
#include "throw.h"
!INTF_ON
  USE parkind1, ONLY : jpim, jplm
  USE rttov_types, ONLY : rttov_coefs, rttov_options
  USE mod_rttov_emis_atlas, ONLY : rttov_emis_atlas_data
!INTF_OFF
  USE parkind1, ONLY : jprb

  USE rttov_const, ONLY :      &
          sensor_id_mw,        &
          sensor_id_po,        &
          errorstatus_success, &
          errorstatus_fatal

  USE mod_rttov_emis_atlas, ONLY : &
    atlas_type_mw, atlas_type_ir,      &
    uwiremis_atlas_id, camel_atlas_id, camel_clim_atlas_id, &
    telsem2_atlas_id, cnrm_mw_atlas_id

  USE mod_uwiremis_atlas,   ONLY : rttov_uwiremis_init
  USE mod_camel_atlas,      ONLY : rttov_camel_init
  USE mod_camel_clim_atlas, ONLY : rttov_camel_clim_init
  USE mod_mwatlas_m2,       ONLY : rttov_read_telsem2_atlas => rttov_readmw_atlas
  USE mod_cnrm_mw_atlas,    ONLY : rttov_cnrmmwemis_init
!INTF_ON
  IMPLICIT NONE

  INTEGER(KIND=jpim),          INTENT(OUT)             :: err
  TYPE(rttov_options),         INTENT(IN)              :: opts
  INTEGER(KIND=jpim),          INTENT(IN)              :: imonth
  INTEGER(KIND=jpim),          INTENT(IN)              :: atlas_type
  TYPE(rttov_emis_atlas_data), INTENT(INOUT)           :: atlas
  INTEGER(KIND=jpim),          INTENT(IN),    OPTIONAL :: atlas_id
  CHARACTER(LEN=*),            INTENT(IN),    OPTIONAL :: path
  TYPE(rttov_coefs),           INTENT(IN),    OPTIONAL :: coefs
  LOGICAL(KIND=jplm),          INTENT(IN),    OPTIONAL :: ir_atlas_read_std
  LOGICAL(KIND=jplm),          INTENT(IN),    OPTIONAL :: ir_atlas_ang_corr
  INTEGER(KIND=jpim),          INTENT(IN),    OPTIONAL :: year
!INTF_END

#include "rttov_errorreport.interface"

  CHARACTER(LEN=300) :: fpath
  INTEGER(KIND=jpim) :: id, iyear
  LOGICAL(KIND=jplm) :: single_inst, ir_atlas_std_init, ir_atlas_do_ang_corr
  LOGICAL(KIND=jplm) :: htfrtc, sensor_mw
  INTEGER(KIND=jpim) :: id_platform, id_sat, id_inst
  REAL(KIND=jprb), POINTER :: chan_wvn(:)

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

  IF (atlas_type < 1 .OR. atlas_type > 2) THEN
    err = errorstatus_fatal
    THROWM(err .NE. 0, 'Invalid atlas_type: must be either 1 => MW or 2 => IR')
  ENDIF

  fpath = './'
  IF (PRESENT(path)) fpath = TRIM(path)//'/'

  single_inst = PRESENT(coefs)

  IF (single_inst) THEN
    ! Determine if user is calling HTFRTC
    htfrtc = ASSOCIATED(coefs%coef_htfrtc%freq)

    IF (htfrtc) THEN
      id_platform = 0  ! HTFRTC coefs don't contain these IDs (and in any case
      id_sat      = 0  !   they are irrelevant since the centroid frequencies
      id_inst     = 0  !   are the same regardless of the sensor).
      chan_wvn => coefs%coef_htfrtc%freq
      sensor_mw   = .FALSE.
    ELSE
      id_platform = coefs%coef%id_platform
      id_sat      = coefs%coef%id_sat
      id_inst     = coefs%coef%id_inst
      chan_wvn => coefs%coef%ff_cwn
      sensor_mw   = coefs%coef%id_sensor == sensor_id_mw .OR. &
                    coefs%coef%id_sensor == sensor_id_po
    ENDIF

    IF ((atlas_type == atlas_type_mw .AND. .NOT. sensor_mw) .OR. &
        (atlas_type == atlas_type_ir .AND. sensor_mw)) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0, 'Mismatch between atlas_type and coefficients')
    ENDIF
  ENDIF

  ir_atlas_std_init = .FALSE.
  IF (PRESENT(ir_atlas_read_std)) ir_atlas_std_init = ir_atlas_read_std

  ir_atlas_do_ang_corr = .FALSE.
  IF (PRESENT(ir_atlas_ang_corr)) ir_atlas_do_ang_corr = ir_atlas_ang_corr

  id = -1
  IF (PRESENT(atlas_id)) id = atlas_id

  iyear = 0
  IF (PRESENT(year)) iyear = year

  IF (atlas_type == atlas_type_mw) THEN

    ! MW atlas

    IF (id <= 0) id = telsem2_atlas_id

    IF (id == telsem2_atlas_id) THEN
      CALL rttov_read_telsem2_atlas(TRIM(fpath), imonth, atlas%telsem2_atlas, opts%config%verbose, err)
      THROWM(err .NE. 0, 'Error initialising TELSEM2 emissivity atlas')
    ELSEIF (id == cnrm_mw_atlas_id) THEN
      IF (.NOT. single_inst) THEN
        err = errorstatus_fatal
        THROWM(err.NE.0, 'The coefs argument is mandatory for the CNRM MW atlas')
      ENDIF
      CALL rttov_cnrmmwemis_init(TRIM(fpath), imonth, iyear, atlas%cnrm_mw_atlas, coefs%coef, opts%config%verbose, err)
      THROWM(err .NE. 0, 'Error initialising CNRM MW emissivity atlas')
    ELSE
      err = errorstatus_fatal
      THROWM(err .NE. 0, 'Unknown MW atlas ID')
    ENDIF

    atlas%is_mw = .TRUE.

  ELSE

    ! IR atlas

    IF (id <= 0) id = uwiremis_atlas_id

    IF (id == uwiremis_atlas_id) THEN
      IF (single_inst) THEN
        CALL rttov_uwiremis_init(TRIM(fpath),            &
                                 imonth,                 &
                                 opts%config%verbose,    &
                                 ir_atlas_std_init,      &
                                 ir_atlas_do_ang_corr,   &
                                 atlas%uwiremis_atlas,   &
                                 err,                    &
                                 chan_wvn(:),            &
                                 id_platform,            &
                                 id_sat,                 &
                                 id_inst)
      ELSE
        CALL rttov_uwiremis_init(TRIM(fpath),            &
                                 imonth,                 &
                                 opts%config%verbose,    &
                                 ir_atlas_std_init,      &
                                 ir_atlas_do_ang_corr,   &
                                 atlas%uwiremis_atlas,   &
                                 err)
      ENDIF
      THROWM(err .NE. 0, 'Error initialising UW IR emissivity atlas')
    ELSEIF (id == camel_atlas_id) THEN
      IF (single_inst) THEN
        CALL rttov_camel_init(TRIM(fpath),            &
                              imonth,                 &
                              opts%config%verbose,    &
                              ir_atlas_std_init,      &
                              ir_atlas_do_ang_corr,   &
                              atlas%camel_atlas,      &
                              err,                    &
                              chan_wvn(:),            &
                              id_platform,            &
                              id_sat,                 &
                              id_inst)
      ELSE
        CALL rttov_camel_init(TRIM(fpath),            &
                              imonth,                 &
                              opts%config%verbose,    &
                              ir_atlas_std_init,      &
                              ir_atlas_do_ang_corr,   &
                              atlas%camel_atlas,      &
                              err)
      ENDIF
      THROWM(err .NE. 0, 'Error initialising CAMEL IR emissivity atlas')
    ELSEIF (id == camel_clim_atlas_id) THEN
      IF (single_inst) THEN
        CALL rttov_camel_clim_init(TRIM(fpath),       &
                              imonth,                 &
                              opts%config%verbose,    &
                              ir_atlas_std_init,      &
                              ir_atlas_do_ang_corr,   &
                              atlas%camel_clim_atlas, &
                              err,                    &
                              chan_wvn(:),            &
                              id_platform,            &
                              id_sat,                 &
                              id_inst)
      ELSE
        CALL rttov_camel_clim_init(TRIM(fpath),       &
                              imonth,                 &
                              opts%config%verbose,    &
                              ir_atlas_std_init,      &
                              ir_atlas_do_ang_corr,   &
                              atlas%camel_clim_atlas, &
                              err)
      ENDIF
      THROWM(err .NE. 0, 'Error initialising CAMEL climatology IR emissivity atlas')
    ELSE
      err = errorstatus_fatal
      THROWM(err .NE. 0, 'Unknown IR atlas ID')
    ENDIF

    atlas%is_mw = .FALSE.

  ENDIF

  atlas%atlas_id = id
  atlas%init = .TRUE.

  CATCH

END SUBROUTINE rttov_setup_emis_atlas
