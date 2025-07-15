! Description:
!> @file
!!   Defines the subroutines which form the wrapper interface.
!
!> @brief
!!   Defines the subroutines which form the wrapper interface.
!!
!! @details
!!   Defines the subroutines which form the wrapper interface:
!!
!!   rttov_load_inst               - call to initialise an instrument (set options, read coefs)
!!   rttov_set_options             - call to change options for an instrument (after reading coefs)
!!   rttov_print_options           - call to print current options for an instrument
!!   rttov_call_direct             - call RTTOV direct
!!   rttov_visir_scatt_call_direct - call RTTOV direct with aer/cld optical properties (visible/IR)
!!   rttov_scatt_call_direct       - call RTTOV-SCATT direct (MW)
!!   rttov_call_k                  - call RTTOV K
!!   rttov_visir_scatt_call_k      - call RTTOV K with aer/cld optical properties (visible/IR)
!!   rttov_scatt_call_k            - call RTTOV-SCATT K (MW)
!!   rttov_drop_inst               - deallocation of a single instrument
!!   rttov_drop_all                - deallocation of all instruments and atlases
!!
!!   rttov_load_ir_emis_atlas      - initialise an IR emissivity atlas
!!   rttov_load_mw_emis_atlas      - initialise an MW emissivity atlas
!!   rttov_load_brdf_atlas         - initialise a BRDF atlas
!!   rttov_get_emisbrdf            - return emissivity/BRDF data from an atlas
!!   rttov_drop_atlas              - deallocation of an atlas
!!
!!   rttov_bpr                     - calculate the bpr parameter given a phase function
!!   rttov_legcoef                 - calculate Legendre coefficients given a phase function
!!
!!   rttov_get_*                   - call to return members of radiance, radiance2, transmission,
!!                                   and RTTOV-SCATT emissivity retrieval terms structures
!!
!!   rttov_get_coef_val_*          - call to return variables from RTTOV coefficients
!!                                   structure; not intended for general use
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

!> This subroutine is only here to ensure the Makefile.PL will create dependencies
SUBROUTINE rttov_wrapper_f2py
!INTF_END

  USE parkind1
  USE rttov_types
  USE rttov_const
  USE rttov_wrapper_handle
  USE rttov_wrapper_transfer

END SUBROUTINE rttov_wrapper_f2py

!> Initialises an instrument by setting the options and reading coefficients
!! @param[out]    inst_id     returned ID for instrument, if <=0 then initialisation failed
!! @param[in]     opts_str    string specifying options (see wrapper user guide)
!! @param[in]     nchannels   number of channels to load (not required in Python)
!! @param[in]     channels    list of channels to load, can be (/0/) in which case all channels are read
SUBROUTINE rttov_load_inst(inst_id, opts_str, nchannels, channels)
!f2py intent(hide):: nchannels=len(channels)

#include "throw.h"

  USE parkind1
  USE rttov_wrapper_handle

  IMPLICIT NONE

  INTEGER(jpim),    INTENT(OUT) :: inst_id
  CHARACTER(LEN=*), INTENT(IN)  :: opts_str
  INTEGER(jpim),    INTENT(IN)  :: nchannels
  INTEGER(jpim),    INTENT(IN)  :: channels(nchannels)

  INTEGER(jpim) :: err
  TYPE(rttovwrapperhandle_type), POINTER :: rth

#include "rttov_user_options_checkinput.interface"
#include "rttov_errorreport.interface"

!------------------------------------------------------------------------------

  TRY

  CALL rttov_wrapper_handle_new_inst(err, inst_id, rth)
  THROWM(err.NE.0, 'Error allocating new instrument')

  CALL rttov_set_options(err, inst_id, opts_str(1:LEN(opts_str)))
  THROWM(err.NE.0, 'Error setting options')

  CALL rttov_wrapper_handle_load(err, rth, channels)
  THROW(err.NE.0)

  CALL rttov_wrapper_nullify_structs(rth)

  rth%init = .TRUE.

  IF (rth%check_opts) THEN
    CALL rttov_user_options_checkinput(err, rth%opts, rth%coefs)
    THROW(err.NE.0)
  ENDIF

  IF (rth%verbose) THEN
    INFO("Load coefficients: ")
    INFO(rth%file_coef)
    IF (TRIM(rth%file_scaer) .NE. '') &
      INFO(rth%file_scaer)
    IF (TRIM(rth%file_sccld) .NE. '') &
      INFO(rth%file_sccld)
    IF (TRIM(rth%file_mfasis_cld) .NE. '') &
      INFO(rth%file_mfasis_cld)
    IF (TRIM(rth%file_mfasis_aer) .NE. '') &
      INFO(rth%file_mfasis_aer)
    IF (TRIM(rth%file_pccoef) .NE. '') &
      INFO(rth%file_pccoef)
    IF (TRIM(rth%file_hydrotable) .NE. '') &
      INFO(rth%file_hydrotable)
  ENDIF

  CATCH
  rth%init = .FALSE.
  inst_id = -1
END SUBROUTINE rttov_load_inst


!> Deallocates memory for an instrument
!! @param[out]    err         return status
!! @param[in]     inst_id     ID for instrument to deallocate
SUBROUTINE rttov_drop_inst(err, inst_id)

#include "throw.h"

  USE parkind1
  USE rttov_wrapper_handle

  IMPLICIT NONE
  INTEGER(jpim), INTENT(OUT) :: err
  INTEGER(jpim), INTENT(IN)  :: inst_id

  TYPE(rttovwrapperhandle_type), POINTER :: rth

#include "rttov_errorreport.interface"

!------------------------------------------------------------------------------

  TRY

  CALL rttov_wrapper_handle_get_inst(err, inst_id, rth)
  IF (err .NE. 0) THEN
    IF (rth%verbose) INFO('No coefficients loaded')
    RETURN
  ENDIF

  IF (rth%verbose) THEN
    INFO('Dropping coefficients: ')
    INFO(rth%file_coef)
    IF (TRIM(rth%file_scaer) .NE. '') &
      INFO(rth%file_scaer)
    IF (TRIM(rth%file_sccld) .NE. '') &
      INFO(rth%file_sccld)
    IF (TRIM(rth%file_mfasis_cld) .NE. '') &
      INFO(rth%file_mfasis_cld)
    IF (TRIM(rth%file_mfasis_aer) .NE. '') &
      INFO(rth%file_mfasis_aer)
    IF (TRIM(rth%file_pccoef) .NE. '') &
      INFO(rth%file_pccoef)
    IF (TRIM(rth%file_hydrotable) .NE. '') &
      INFO(rth%file_hydrotable)
  ENDIF

  CALL rttov_wrapper_handle_del_inst(err, rth)
  THROWM(err.NE.0, 'An error occured while dropping coefficients')

  CATCH
END SUBROUTINE rttov_drop_inst


!> Deallocates memory for all initialised instruments and atlases
!! @param[out]    err         return status
SUBROUTINE rttov_drop_all(err)

#include "throw.h"

  USE parkind1
  USE rttov_wrapper_handle

  IMPLICIT NONE
  INTEGER(jpim), INTENT(OUT) :: err

#include "rttov_errorreport.interface"

!------------------------------------------------------------------------------

  TRY

  CALL rttov_wrapper_handle_del_all_atlas(err)
  THROW(err.NE.0)

  CALL rttov_wrapper_handle_del_all_inst(err)
  THROW(err.NE.0)

  CATCH
END SUBROUTINE rttov_drop_all


!> Initialises an IR atlas
!! @param[out]    atlas_wrap_id     returned wrapper ID for atlas, if <=0 then initialisation failed
!! @param[in]     path              string specifying location of atlas data
!! @param[in]     month             month for which to load data (1-12)
!! @param[in]     atlas_id          selects the IR atlas to load (see user guide for the atlas IDs)
!! @param[in]     inst_id           if <=0 the atlas data are initialised for use with any instrument,
!!                                    otherwise they are initialised for use with given instrument
!! @param[in]     ang_corr          if non-zero the zenith angle correction will be included
SUBROUTINE rttov_load_ir_emis_atlas(atlas_wrap_id, path, month, atlas_id, inst_id, ang_corr)

#include "throw.h"

  USE parkind1
  USE rttov_types, ONLY : rttov_options
  USE rttov_wrapper_handle
  USE rttov_const, ONLY : sensor_id_mw, sensor_id_po
  USE mod_rttov_emis_atlas, ONLY : atlas_type_ir

  IMPLICIT NONE

  INTEGER(jpim),    INTENT(OUT) :: atlas_wrap_id
  CHARACTER(LEN=*), INTENT(IN)  :: path
  INTEGER(jpim),    INTENT(IN)  :: month
  INTEGER(jpim),    INTENT(IN)  :: atlas_id
  INTEGER(jpim),    INTENT(IN)  :: inst_id
  INTEGER(jpim),    INTENT(IN)  :: ang_corr

  INTEGER(jpim)       :: err
  LOGICAL(jplm)       :: do_ang_corr
  TYPE(rttov_options) :: opts

  TYPE(rttovwrapperhandle_type), POINTER :: rth
  TYPE(rttovwrapperatlas_type),  POINTER :: ath

#include "rttov_setup_emis_atlas.interface"
#include "rttov_errorreport.interface"

!------------------------------------------------------------------------------

  TRY

  CALL rttov_wrapper_handle_new_atlas(err, atlas_wrap_id, ath)
  THROW(err.NE.0)

  do_ang_corr = (ang_corr .NE. 0)

  IF (inst_id > 0) THEN
    CALL rttov_wrapper_handle_get_inst(err, inst_id, rth)
    THROWM(err.NE.0, 'inst_id is not valid: must be zero or the ID for an initialised instrument')

    IF (rth%coefs%coef%id_sensor == sensor_id_mw .OR. &
        rth%coefs%coef%id_sensor == sensor_id_po) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0, 'Specified instrument is not an IR sensor')
    ENDIF

    ! Initialise atlas for the specific instrument inst_id
    CALL rttov_setup_emis_atlas(          &
           err,                           &
           rth%opts,                      &
           month,                         &
           atlas_type_ir,                 &
           ath%emis_atlas,                &
           atlas_id = atlas_id,           &
           path = path(1:LEN(path)),      &
           coefs = rth%coefs,             &
           ir_atlas_read_std = .FALSE.,   &
           ir_atlas_ang_corr = do_ang_corr)
  ELSE
    ! Initialise atlas for multiple instruments; use dummy opts with verbose = false
    opts%config%verbose = .FALSE.
    CALL rttov_setup_emis_atlas(          &
           err,                           &
           opts,                          &
           month,                         &
           atlas_type_ir,                 &
           ath%emis_atlas,                &
           atlas_id = atlas_id,           &
           path = path(1:LEN(path)),      &
           ir_atlas_read_std = .FALSE.,   &
           ir_atlas_ang_corr = do_ang_corr)
  ENDIF
  THROW(err.NE.0)

  ath%init    = .TRUE.
  ath%is_emis = .TRUE.

  CATCH
  ath%init = .FALSE.
  atlas_wrap_id = -1
END SUBROUTINE rttov_load_ir_emis_atlas


!> Initialises a MW atlas
!! @param[out]    atlas_wrap_id     returned wrapper ID for atlas, if <=0 then initialisation failed
!! @param[in]     path              string specifying location of atlas data
!! @param[in]     month             month for which to load data (1-12)
!! @param[in]     atlas_id          selects the MW atlas to load (see user guide for the atlas IDs)
!! @param[in]     inst_id           ignored for TELSEM2, required for the CNRM MW atlas, specifies the
!!                                    instrument for which atlas data are required
!! @param[in]     year              specifies the year for which to load CNRM MW atlas data, ignored for TELSEM2
SUBROUTINE rttov_load_mw_emis_atlas(atlas_wrap_id, path, month, atlas_id, inst_id, year)

#include "throw.h"

  USE parkind1
  USE rttov_types, ONLY : rttov_options
  USE rttov_wrapper_handle
  USE rttov_const, ONLY : sensor_id_ir, sensor_id_hi
  USE mod_rttov_emis_atlas, ONLY : cnrm_mw_atlas_id, atlas_type_mw

  IMPLICIT NONE

  INTEGER(jpim),    INTENT(OUT) :: atlas_wrap_id
  CHARACTER(LEN=*), INTENT(IN)  :: path
  INTEGER(jpim),    INTENT(IN)  :: month
  INTEGER(jpim),    INTENT(IN)  :: atlas_id
  INTEGER(jpim),    INTENT(IN)  :: inst_id
  INTEGER(jpim),    INTENT(IN)  :: year

  INTEGER(jpim)       :: err
  TYPE(rttov_options) :: opts

  TYPE(rttovwrapperhandle_type), POINTER :: rth
  TYPE(rttovwrapperatlas_type),  POINTER :: ath

#include "rttov_setup_emis_atlas.interface"
#include "rttov_errorreport.interface"

!------------------------------------------------------------------------------

  TRY

  CALL rttov_wrapper_handle_new_atlas(err, atlas_wrap_id, ath)
  THROW(err.NE.0)

  IF (inst_id > 0) THEN
    CALL rttov_wrapper_handle_get_inst(err, inst_id, rth)
    THROWM(err.NE.0, 'inst_id is not valid: must be zero or the ID for an initialised instrument')

    IF (rth%coefs%coef%id_sensor == sensor_id_ir .OR. &
        rth%coefs%coef%id_sensor == sensor_id_hi) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0, 'Specified instrument is not a MW sensor')
    ENDIF

    ! Initialise atlas using the specific instrument inst_id
    ! This is required for CNRM atlas; makes no difference for TELSEM
    CALL rttov_setup_emis_atlas(     &
           err,                      &
           rth%opts,                 &
           month,                    &
           atlas_type_mw,            &
           ath%emis_atlas,           &
           atlas_id = atlas_id,      &
           path = path(1:LEN(path)), &
           coefs = rth%coefs,        &
           year = year)
  ELSE
    IF (atlas_id == cnrm_mw_atlas_id) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0, 'inst_id for a loaded instrument is required to initialise CNRM MW atlas')
    ENDIF

    ! Initialise TELSEM2 atlas using dummy opts with verbose = false (TELSEM2 doesn't use year argument)
    opts%config%verbose = .FALSE.
    CALL rttov_setup_emis_atlas( &
           err,                  &
           opts,                 &
           month,                &
           atlas_type_mw,        &
           ath%emis_atlas,       &
           atlas_id = atlas_id,  &
           path = path(1:LEN(path)))
  ENDIF
  THROW(err.NE.0)

  ath%init    = .TRUE.
  ath%is_emis = .TRUE.

  CATCH
  ath%init = .FALSE.
  atlas_wrap_id = -1
END SUBROUTINE rttov_load_mw_emis_atlas


!> Initialises a BRDF atlas
!! @param[out]    atlas_wrap_id     returned wrapper ID for atlas, if <=0 then initialisation failed
!! @param[in]     path              string specifying location of atlas data
!! @param[in]     month             month for which to load data (1-12)
!! @param[in]     atlas_id          selects the BRDF atlas to load (see user guide for the atlas IDs)
!! @param[in]     inst_id           if <=0 the atlas data are initialised for use with any instrument,
!!                                    otherwise they are initialised for use with given instrument
SUBROUTINE rttov_load_brdf_atlas(atlas_wrap_id, path, month, atlas_id, inst_id)

#include "throw.h"

  USE parkind1
  USE rttov_types, ONLY : rttov_options
  USE rttov_wrapper_handle
  USE rttov_const, ONLY : sensor_id_mw, sensor_id_po

  IMPLICIT NONE

  INTEGER(jpim),    INTENT(OUT) :: atlas_wrap_id
  CHARACTER(LEN=*), INTENT(IN)  :: path
  INTEGER(jpim),    INTENT(IN)  :: month
  INTEGER(jpim),    INTENT(IN)  :: atlas_id
  INTEGER(jpim),    INTENT(IN)  :: inst_id

  INTEGER(jpim)       :: err
  TYPE(rttov_options) :: opts

  TYPE(rttovwrapperhandle_type), POINTER :: rth
  TYPE(rttovwrapperatlas_type),  POINTER :: ath

#include "rttov_setup_brdf_atlas.interface"
#include "rttov_errorreport.interface"

!------------------------------------------------------------------------------

  TRY

  CALL rttov_wrapper_handle_new_atlas(err, atlas_wrap_id, ath)
  THROW(err.NE.0)

  IF (inst_id > 0) THEN
    CALL rttov_wrapper_handle_get_inst(err, inst_id, rth)
    THROWM(err.NE.0, 'inst_id is not valid: must be zero or the ID for an initialised instrument')

    IF (rth%coefs%coef%id_sensor == sensor_id_mw .OR. &
        rth%coefs%coef%id_sensor == sensor_id_po) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0, 'Specified instrument is not a VIS/IR sensor')
    ENDIF

    ! Initialise atlas for the specific instrument inst_id
    CALL rttov_setup_brdf_atlas(     &
           err,                      &
           rth%opts,                 &
           month,                    &
           ath%brdf_atlas,           &
           atlas_id = atlas_id,      &
           path = path(1:LEN(path)), &
           coefs = rth%coefs)
  ELSE
    ! Initialise atlas for multiple instruments using dummy opts with verbose = false
    opts%config%verbose = .FALSE.
    CALL rttov_setup_brdf_atlas(     &
           err,                      &
           opts,                     &
           month,                    &
           ath%brdf_atlas,           &
           atlas_id = atlas_id,      &
           path = path(1:LEN(path)))
  ENDIF
  THROW(err.NE.0)

  ath%init    = .TRUE.
  ath%is_emis = .FALSE.

  CATCH
  ath%init = .FALSE.
  atlas_wrap_id = -1
END SUBROUTINE rttov_load_brdf_atlas


!> Return emissivity or BRDF values from a loaded atlas for multiple profiles.
!! See user guide for definitions of profile variables and for information on
!! which variables are used by which atlases.
!! @param[out]    err               return status
!! @param[in]     atlas_wrap_id     wrapper ID for atlas
!! @param[in]     latitude          array of latitudes
!! @param[in]     longitude         array of longitudes
!! @param[in]     surftype          array of surftype values
!! @param[in]     watertype         array of watertype values
!! @param[in]     zenangle          array of zenith angles
!! @param[in]     azangle           array of azimuth angles
!! @param[in]     sunzenangle       array of solar zenith angles
!! @param[in]     sunazangle        array of solar azimuth angles
!! @param[in]     snow_fraction     array of snow fractions
!! @param[in]     inst_id           if the atlas data were initialised for use with a specific
!!                                    sensor this must be the ID for an appropriate instrument
!! @param[in]     channel_list      list of channels for which to return emissivities/BRDFs for each profile
!! @param[in,out] emisbrdf          returned emissivities/BRDFs
!! @param[in]     nchannels         number of channels (size of channel_list, not required in Python)
!! @param[in]     nprofiles         number of profiles (size of latitude/etc, not required in Python)
SUBROUTINE rttov_get_emisbrdf( &
              err,            &
              atlas_wrap_id,  &
              latitude,       &
              longitude,      &
              surftype,       &
              watertype,      &
              zenangle,       &
              azangle,        &
              sunzenangle,    &
              sunazangle,     &
              snow_fraction,  &
              inst_id,        &
              channel_list,   &
              emisbrdf,       &
              nchannels, nprofiles)
!f2py threadsafe
!f2py intent(hide):: nchannels=len(channel_list)
!f2py intent(hide):: nprofiles=len(latitude)

#include "throw.h"

  USE parkind1, ONLY : jpim, jprb
  USE rttov_wrapper_handle
  USE rttov_types, ONLY : rttov_chanprof, rttov_profile

  IMPLICIT NONE

  INTEGER(jpim), INTENT(OUT)   :: err
  INTEGER(jpim), INTENT(IN)    :: nchannels, nprofiles

  INTEGER(jpim), INTENT(IN)    :: atlas_wrap_id
  REAL(jprb),    INTENT(IN)    :: latitude(nprofiles)
  REAL(jprb),    INTENT(IN)    :: longitude(nprofiles)
  INTEGER(jpim), INTENT(IN)    :: surftype(nprofiles)
  INTEGER(jpim), INTENT(IN)    :: watertype(nprofiles)
  REAL(jprb),    INTENT(IN)    :: zenangle(nprofiles)
  REAL(jprb),    INTENT(IN)    :: azangle(nprofiles)
  REAL(jprb),    INTENT(IN)    :: sunzenangle(nprofiles)
  REAL(jprb),    INTENT(IN)    :: sunazangle(nprofiles)
  REAL(jprb),    INTENT(IN)    :: snow_fraction(nprofiles)
  INTEGER(jpim), INTENT(IN)    :: inst_id
  INTEGER(jpim), INTENT(IN)    :: channel_list(nchannels)
  REAL(jprb),    INTENT(INOUT) :: emisbrdf(nchannels,nprofiles)

  TYPE(rttov_chanprof) :: chanprof(nprofiles*nchannels)
  TYPE(rttov_profile)  :: profiles(nprofiles)
  INTEGER(jpim)        :: i, lo, hi

  TYPE(rttovwrapperhandle_type), POINTER :: rth
  TYPE(rttovwrapperatlas_type),  POINTER :: ath

#include "rttov_get_emis.interface"
#include "rttov_get_brdf.interface"
#include "rttov_errorreport.interface"

!------------------------------------------------------------------------------

  TRY

  CALL rttov_wrapper_handle_get_atlas(err, atlas_wrap_id, ath)
  THROW(err.NE.0)

  CALL rttov_wrapper_handle_get_inst(err, inst_id, rth)
  THROW(err.NE.0)

  DO i = 1, nprofiles
    lo = (i - 1) * nchannels + 1
    hi = lo + nchannels - 1

    chanprof(lo:hi)%prof = i
    chanprof(lo:hi)%chan = channel_list

    profiles(i)%latitude           = latitude(i)
    profiles(i)%longitude          = longitude(i)
    profiles(i)%zenangle           = zenangle(i)
    profiles(i)%azangle            = azangle(i)
    profiles(i)%sunzenangle        = sunzenangle(i)
    profiles(i)%sunazangle         = sunazangle(i)
    profiles(i)%skin%surftype      = surftype(i)
    profiles(i)%skin%watertype     = watertype(i)
    profiles(i)%skin%snow_fraction = snow_fraction(i)
  ENDDO

  IF (ath%is_emis) THEN
    CALL rttov_get_emis(    &
            err,            &
            rth%opts,       &
            chanprof,       &
            profiles,       &
            rth%coefs,      &
            ath%emis_atlas, &
            emisbrdf)
  ELSE
    CALL rttov_get_brdf(    &
            err,            &
            rth%opts,       &
            chanprof,       &
            profiles,       &
            rth%coefs,      &
            ath%brdf_atlas, &
            emisbrdf)
  ENDIF
  THROW(err.NE.0)

  CATCH
END SUBROUTINE rttov_get_emisbrdf


!> Deallocate data for an atlas
!! @param[out]    err               return status
!! @param[in]     atlas_wrap_id     wrapper ID for atlas to deallocate
SUBROUTINE rttov_drop_atlas(err, atlas_wrap_id)

#include "throw.h"

  USE parkind1, ONLY : jpim
  USE rttov_wrapper_handle

  IMPLICIT NONE

  INTEGER(jpim), INTENT(OUT) :: err
  INTEGER(jpim), INTENT(IN)  :: atlas_wrap_id

  TYPE(rttovwrapperatlas_type), POINTER :: ath

#include "rttov_errorreport.interface"

!------------------------------------------------------------------------------

  TRY

  CALL rttov_wrapper_handle_get_atlas(err, atlas_wrap_id, ath)
  THROW(err.NE.0)

  CALL rttov_wrapper_handle_del_atlas(err, ath)
  THROWM(err.NE.0, 'An error occured while dropping the atlas')

  CATCH
END SUBROUTINE rttov_drop_atlas


!> Call rttov_direct for clear-sky or visible/IR scattering simulations using
!! pre-defined particle types.
!! @param[out]    err               return status
!! @param[in]     inst_id           ID for instrument to simulate
!! @param[in]     channel_list      list of channels to simulate for for each profile
!! @param[in]     datetimes         profile dates/times: year, month, day, hour, min, sec (6,nprofiles)
!! @param[in]     angles            profile angles: zenang, aziang, sunzenang, sunaziang (4,nprofiles)
!! @param[in]     surfgeom          profile surface geometry: lat, lon, elevation (3,nprofiles)
!! @param[in]     surftype          profile surface type: surftype, watertype (2,nprofiles)
!! @param[in]     skin              profile skin data: skin T, salinity, snow_frac, foam_frac, fastem_coefsx5 (9,nprofiles)
!! @param[in]     s2m               profile s2m data: 2m p, 2m t, 2m q, 10m wind u, v, wind-fetch (6,nprofiles)
!! @param[in]     simplecloud       profile simple cloud data: ctp, cfraction (2,nprofiles)
!! @param[in]     clwscheme         profile clw scheme: clw_scheme, clwde_param (2,nprofiles)
!! @param[in]     icecloud          profile ice scheme: ice_scheme, icede_param (2,nprofiles)
!! @param[in]     zeeman            profile B-field data: Be, cosbk (2,nprofiles)
!! @param[in]     p                 pressure profiles (nlevels,nprofiles)
!! @param[in]     t                 temperature profiles (nlevels,nprofiles)
!! @param[in]     gas_units         profile gas units
!! @param[in]     mmr_cldaer        profile mmr_cldaer flag
!! @param[in]     gas_id            list of gas IDs specified in gases argument (ngases)
!! @param[in]     gases             gas, aerosol and cloud profiles (nlevels,nprofiles,ngases)
!! @param[in]     surfemisrefl      input/output surface emissivities, reflectances, specularities (nchannels,nprofiles,4)
!! @param[in,out] btrefl            returned BTs/reflectances (nchannels,nprofiles)
!! @param[in,out] rads              returned radiances (nchannels,nprofiles)
!! @param[in]     nchannels         number of channels (size of channel_list, not required in Python)
!! @param[in]     ngases            number of gases (size of gas_id, not required in Python)
!! @param[in]     nlevels           number of levels (size of p(:,1), not required in Python)
!! @param[in]     nprofiles         number of profiles (size of p(1,:), not required in Python)
SUBROUTINE rttov_call_direct( &
    err,               &
    inst_id,           &
    channel_list,      &
    datetimes,         &    ! profile dates/times                                     (6,nprofiles)
    angles,            &    ! satzen, satazi, sunzen, sunazi angles                   (4,nprofiles)
    surfgeom,          &    ! lat, lon, elevation                                     (3,nprofiles)
    surftype,          &    ! surftype, watertype                                     (2,nprofiles)
    skin,              &    ! skin T, salinity, snow_frac, foam_frac, fastem_coefsx5  (9,nprofiles)
    s2m,               &    ! 2m p, 2m t, 2m q, 10m wind u, v, wind-fetch             (6,nprofiles)
    simplecloud,       &    ! ctp, cfraction                                          (2,nprofiles)
    clwscheme,         &    ! clw_scheme, clwde_param                                 (2,nprofiles)
    icecloud,          &    ! ice_scheme, icede_param                                 (2,nprofiles)
    zeeman,            &    ! Be, cosbk                                               (2,nprofiles)
    p,                 &    ! pressure                                                (nlevels,nprofiles)
    t,                 &    ! temperature                                             (nlevels,nprofiles)
    gas_units,         &    ! units for gas profiles
    mmr_cldaer,        &    ! set mmr_cldaer flag
    gas_id,            &    ! gas ID list                                             (ngases)
    gases,             &    ! gas profiles                                            (nlevels,nprofiles,ngases)
    surfemisrefl,      &    ! surface emissivities, reflectances, specularities       (nchannels,nprofiles,4)
    btrefl,            &    ! output BTs/refls (for thermal/solar chans)              (nchannels,nprofiles)
    rads,              &    ! output radiances                                        (nchannels,nprofiles)
    nchannels, ngases, nlevels, nprofiles)
!f2py threadsafe
!f2py intent(hide):: nchannels=len(channel_list)
!f2py intent(hide):: ngases=len(gas_id)
!f2py intent(hide):: nlevels=shape(p, 1)
!f2py intent(hide):: nprofiles=shape(p, 2)

!
! Prepares input profiles, calls RTTOV and returns radiances
!

#include "throw.h"

  USE rttov_wrapper_handle

  USE rttov_wrapper_transfer

  USE parkind1, ONLY : jpim, jprb, jplm

  USE rttov_types, ONLY : &
      rttov_chanprof,     &
      rttov_emissivity,   &
      rttov_reflectance,  &
      rttov_profile,      &
      rttov_radiance,     &
      rttov_radiance2,    &
      rttov_transmission

  IMPLICIT NONE

  INTEGER(jpim),    INTENT(OUT)   :: err
  INTEGER(jpim),    INTENT(IN)    :: nchannels, ngases, nlevels, nprofiles

  INTEGER(jpim),    INTENT(IN)    :: inst_id
  INTEGER(jpim),    INTENT(IN)    :: channel_list(nchannels)
  INTEGER(jpim),    INTENT(IN)    :: datetimes(6,nprofiles)
  REAL(jprb),       INTENT(IN)    :: angles(4,nprofiles)
  REAL(jprb),       INTENT(IN)    :: surfgeom(3,nprofiles)
  INTEGER(jpim),    INTENT(IN)    :: surftype(2,nprofiles)
  REAL(jprb),       INTENT(IN)    :: skin(9,nprofiles)
  REAL(jprb),       INTENT(IN)    :: s2m(6,nprofiles)
  REAL(jprb),       INTENT(IN)    :: simplecloud(2,nprofiles)
  INTEGER(jpim),    INTENT(IN)    :: clwscheme(2,nprofiles)
  INTEGER(jpim),    INTENT(IN)    :: icecloud(2,nprofiles)
  REAL(jprb),       INTENT(IN)    :: zeeman(2,nprofiles)
  REAL(jprb),       INTENT(IN)    :: p(nlevels,nprofiles)
  REAL(jprb),       INTENT(IN)    :: t(nlevels,nprofiles)
  INTEGER(jpim),    INTENT(IN)    :: gas_units
  INTEGER(jpim),    INTENT(IN)    :: mmr_cldaer
  INTEGER(jpim),    INTENT(IN)    :: gas_id(ngases)
  REAL(jprb),       INTENT(IN)    :: gases(nlevels,nprofiles,ngases)

  REAL(jprb),       INTENT(INOUT) :: surfemisrefl(nchannels,nprofiles,4)

  REAL(jprb),       INTENT(INOUT) :: btrefl(nchannels,nprofiles)
  REAL(jprb),       INTENT(INOUT) :: rads(nchannels,nprofiles)

  INTEGER(jpim)      :: asw
  INTEGER(jpim)      :: i, k, iprof, iprof1, lochan, hichan
  INTEGER(jpim)      :: nchanprof, thisnchanprof, nprof, thisnprof, ncalls
  LOGICAL(jplm)      :: lmmr_cldaer
  CHARACTER(LEN=256) :: msg

  TYPE(rttov_chanprof),    POINTER :: chanprof(:)
  LOGICAL(jplm),           POINTER :: calcemis(:)
  TYPE(rttov_emissivity),  POINTER :: emissivity(:)
  LOGICAL(jplm),           POINTER :: calcrefl(:)
  TYPE(rttov_reflectance), POINTER :: reflectance(:)
  TYPE(rttov_profile),     POINTER :: profiles(:)
  TYPE(rttov_transmission)         :: transmission
  TYPE(rttov_radiance)             :: radiance
  TYPE(rttov_radiance2)            :: radiance2

  TYPE(rttovwrapperhandle_type), POINTER :: rth

#include "rttov_alloc_direct.interface"
#include "rttov_init_emis_refl.interface"
#include "rttov_direct.interface"
#include "rttov_parallel_direct.interface"
#include "rttov_errorreport.interface"

!------------------------------------------------------------------------------

  TRY

  NULLIFY(chanprof,    &
          calcemis,    &
          emissivity,  &
          calcrefl,    &
          reflectance, &
          profiles)

  CALL rttov_wrapper_handle_get_inst(err, inst_id, rth)
  THROW(err.NE.0)

  IF ((rth%opts%rt_ir%addclouds .AND. rth%opts%rt_ir%user_cld_opt_param) .OR. &
      (rth%opts%rt_ir%addaerosl .AND. rth%opts%rt_ir%user_aer_opt_param)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Use rttov_visir_scatt_call_direct when providing cloud/aerosol optical properties')
  ENDIF

  ! Deallocate (if necessary) and re-allocate the output structures
  CALL rttov_wrapper_handle_alloc(err, rth, nprofiles, nchannels, nlevels, .FALSE._jplm, 0_jpim)
  THROW(err.NE.0)
  CALL rttov_wrapper_handle_alloc(err, rth, nprofiles, nchannels, nlevels, .FALSE._jplm, 1_jpim)
  THROW(err.NE.0)

  ! nchannels     - input argument, number of channels simulated for every profile
  ! nprofiles     - input argument, total number of profiles passed to this subroutine
  ! ncalls        - number of calls to RTTOV
  ! nprof         - number of profiles per call to RTTOV
  ! thisnprof     - number of profiles in current call (may differ to nprof on last call)
  ! nchanprof     - size of chanprof per call to RTTOV
  ! thisnchanprof - size of chanprof for current call (may differ to nchanprof on last call)

  nprof = MIN(nprofiles, rth%nprofs_per_call)
  nchanprof = nprof * nchannels
  lmmr_cldaer = (mmr_cldaer .NE. 0)

  ! Allocate RTTOV data
  asw = 1_jpim
  CALL rttov_alloc_direct( &
              err,                       &
              asw,                       &
              nprof,                     &
              nchanprof,                 &
              nlevels,                   &
              chanprof,                  &
              rth%opts,                  &
              profiles,                  &
              rth%coefs,                 &
              transmission,              &
              radiance,                  &
              radiance2   = radiance2,   &
              calcemis    = calcemis,    &
              emissivity  = emissivity,  &
              calcrefl    = calcrefl,    &
              reflectance = reflectance, &
              init        = .TRUE._jplm)
  THROWM(err.NE.0, 'Error allocating RTTOV structures')

  ! Populate input arrays which don't change
  DO i = 1, nprof
    lochan = (i - 1) * nchannels + 1
    hichan = i * nchannels
    chanprof(lochan:hichan)%prof = i
    chanprof(lochan:hichan)%chan = channel_list(:)
  ENDDO

  IF (rth%verbose) THEN
    WRITE(msg, '(a,i4,a,i8)') &
      'Running RTTOV using nthreads = ', rth%nthreads, &
      ' and nprofs_per_call = ', rth%nprofs_per_call
    INFO(msg)
  ENDIF

  ! Main loop over batches of profiles

  ncalls = nprofiles / rth%nprofs_per_call
  IF (MOD(nprofiles, rth%nprofs_per_call) /= 0) ncalls = ncalls + 1
  thisnchanprof = nchanprof ! This may be overwritten in the final loop below
  thisnprof = rth%nprofs_per_call

  DO i = 1, ncalls

    ! Calculate low/high profile and chanprof indices
    iprof1 = (i - 1) * rth%nprofs_per_call + 1
    IF (iprof1 + rth%nprofs_per_call > nprofiles) THEN
      thisnprof = nprofiles - iprof1 + 1
      thisnchanprof = thisnprof * nchannels
    ENDIF

    ! Populate the profile structure
    CALL rttov_copy_to_profiles(  &
        rth,                      &
        profiles(1:thisnprof),    &
        iprof1,                   &
        datetimes,                &
        angles, surfgeom,         &
        surftype, skin, s2m,      &
        p, t, gas_units,          &
        gas_id, gases,            &
        lmmr_cldaer,              &
        simplecloud,              &
        clwscheme, icecloud, zeeman)

    ! Set surface emissivities and BRDFs: if user passed any non-negative
    ! emissivities/BRDFs then use these and set calcemis/calcrefl FALSE
    calcemis(:) = .TRUE.
    calcrefl(:) = .TRUE.
    CALL rttov_init_emis_refl(emissivity, reflectance)
    ! Simple cloud scheme cloud top albedo cannot currently be set by wrapper
    ! reflectance%refl_cloud_top = 0._jprb

    DO iprof = 1, thisnprof
      k = iprof1 + iprof - 1
      lochan = (iprof - 1) * nchannels + 1
      hichan = iprof * nchannels
      WHERE (surfemisrefl(:,k,1) >= 0._jprb)
        emissivity(lochan:hichan)%emis_in = surfemisrefl(:,k,1)
        calcemis(lochan:hichan) = .FALSE.
      ENDWHERE
      WHERE (surfemisrefl(:,k,2) >= 0._jprb)
        reflectance(lochan:hichan)%refl_in = surfemisrefl(:,k,2)
        calcrefl(lochan:hichan) = .FALSE.
      ENDWHERE
      WHERE (surfemisrefl(:,k,3) >= 0._jprb)
        reflectance(lochan:hichan)%diffuse_refl_in = surfemisrefl(:,k,3)
      ENDWHERE
      WHERE (surfemisrefl(:,k,4) >= 0._jprb)
        emissivity(lochan:hichan)%specularity = surfemisrefl(:,k,4)
      ENDWHERE
    ENDDO

    ! Run RTTOV
    IF (rth%nthreads <= 1) THEN
      CALL rttov_direct( &
          err,                                        &
          chanprof(1:thisnchanprof),                  &
          rth%opts,                                   &
          profiles(1:thisnprof),                      &
          rth%coefs,                                  &
          transmission,                               &
          radiance,                                   &
          radiance2,                                  &
          calcemis    = calcemis(1:thisnchanprof),    &
          emissivity  = emissivity(1:thisnchanprof),  &
          calcrefl    = calcrefl(1:thisnchanprof),    &
          reflectance = reflectance(1:thisnchanprof))
    ELSE
      CALL rttov_parallel_direct( &
          err,                                        &
          chanprof(1:thisnchanprof),                  &
          rth%opts,                                   &
          profiles(1:thisnprof),                      &
          rth%coefs,                                  &
          transmission,                               &
          radiance,                                   &
          radiance2,                                  &
          calcemis    = calcemis(1:thisnchanprof),    &
          emissivity  = emissivity(1:thisnchanprof),  &
          calcrefl    = calcrefl(1:thisnchanprof),    &
          reflectance = reflectance(1:thisnchanprof), &
          nthreads    = rth%nthreads)
    ENDIF
    THROWM(err.NE.0, 'Error running RTTOV')


    ! Store results
    DO iprof = 1, thisnprof
      k = iprof1 + iprof - 1
      lochan = (iprof - 1) * nchannels + 1
      hichan = iprof * nchannels
      WHERE (rth%coefs%coef%ss_val_chn(channel_list(:)) < 2)
        btrefl(:,k) = radiance%bt(lochan:hichan)
      ELSEWHERE
        btrefl(:,k) = radiance%refl(lochan:hichan)
      ENDWHERE
      rads(:,k) = radiance%total(lochan:hichan)
      surfemisrefl(:,k,1) = emissivity(lochan:hichan)%emis_out
      surfemisrefl(:,k,2) = reflectance(lochan:hichan)%refl_out
      surfemisrefl(:,k,3) = reflectance(lochan:hichan)%diffuse_refl_out
    ENDDO

    lochan = (iprof1 - 1) * nchannels + 1
    hichan = (iprof1 + thisnprof - 1) * nchannels

    IF (rth%store_trans) THEN
      rth%transmission%tau_total(lochan:hichan)             = transmission%tau_total(1:thisnchanprof)
      rth%transmission%tau_levels(:,lochan:hichan)          = transmission%tau_levels(:,1:thisnchanprof)
      rth%transmission%tausun_total_path2(lochan:hichan)    = transmission%tausun_total_path2(1:thisnchanprof)
      rth%transmission%tausun_levels_path2(:,lochan:hichan) = transmission%tausun_levels_path2(:,1:thisnchanprof)
      rth%transmission%tausun_total_path1(lochan:hichan)    = transmission%tausun_total_path1(1:thisnchanprof)
      rth%transmission%tausun_levels_path1(:,lochan:hichan) = transmission%tausun_levels_path1(:,1:thisnchanprof)
      rth%transmission%tau_total_cld(lochan:hichan)         = transmission%tau_total_cld(1:thisnchanprof)
      rth%transmission%tau_levels_cld(:,lochan:hichan)      = transmission%tau_levels_cld(:,1:thisnchanprof)
    ENDIF

    IF (rth%store_rad) THEN
      rth%radiance%clear(lochan:hichan)              = radiance%clear(1:thisnchanprof)
      rth%radiance%total(lochan:hichan)              = radiance%total(1:thisnchanprof)
      rth%radiance%bt_clear(lochan:hichan)           = radiance%bt_clear(1:thisnchanprof)
      rth%radiance%bt(lochan:hichan)                 = radiance%bt(1:thisnchanprof)
      rth%radiance%refl_clear(lochan:hichan)         = radiance%refl_clear(1:thisnchanprof)
      rth%radiance%refl(lochan:hichan)               = radiance%refl(1:thisnchanprof)
      rth%radiance%cloudy(lochan:hichan)             = radiance%cloudy(1:thisnchanprof)
      rth%radiance%overcast(:,lochan:hichan)         = radiance%overcast(:,1:thisnchanprof)
      rth%radiance%quality(lochan:hichan)            = radiance%quality(1:thisnchanprof)
      rth%radiance%plane_parallel                    = radiance%plane_parallel
      rth%radiance%geometric_height(:,lochan:hichan) = radiance%geometric_height(:,1:thisnchanprof)
    ENDIF

    IF (rth%store_rad2) THEN
      rth%radiance2%upclear(lochan:hichan)     = radiance2%upclear(1:thisnchanprof)
      rth%radiance2%dnclear(lochan:hichan)     = radiance2%dnclear(1:thisnchanprof)
      rth%radiance2%refldnclear(lochan:hichan) = radiance2%refldnclear(1:thisnchanprof)
      rth%radiance2%up(:,lochan:hichan)        = radiance2%up(:,1:thisnchanprof)
      rth%radiance2%down(:,lochan:hichan)      = radiance2%down(:,1:thisnchanprof)
      rth%radiance2%surf(:,lochan:hichan)      = radiance2%surf(:,1:thisnchanprof)
    ENDIF

  ENDDO ! loop over calls to RTTOV

  ! Deallocate RTTOV data
  asw = 0_jpim
  CALL rttov_alloc_direct( &
              err,                      &
              asw,                      &
              nprof,                    &
              nchanprof,                &
              nlevels,                  &
              chanprof,                 &
              rth%opts,                 &
              profiles,                 &
              rth%coefs,                &
              transmission,             &
              radiance,                 &
              radiance2   = radiance2,  &
              calcemis    = calcemis,   &
              emissivity  = emissivity, &
              calcrefl    = calcrefl,   &
              reflectance = reflectance)
  THROWM(err.NE.0, 'Error deallocating RTTOV structures')

  CATCH
END SUBROUTINE rttov_call_direct


!> Call rttov_direct for visible/IR scattering simulations using explicit optical properties
!! @param[out]    err               return status
!! @param[in]     inst_id           ID for instrument to simulate
!! @param[in]     channel_list      list of channels to simulate for for each profile
!! @param[in]     datetimes         profile dates/times: year, month, day, hour, min, sec (6,nprofiles)
!! @param[in]     angles            profile angles: zenang, aziang, sunzenang, sunaziang (4,nprofiles)
!! @param[in]     surfgeom          profile surface geometry: lat, lon, elevation (3,nprofiles)
!! @param[in]     surftype          profile surface type: surftype, watertype (2,nprofiles)
!! @param[in]     skin              profile skin data: skin T, salinity, snow_frac, foam_frac, fastem_coefsx5 (9,nprofiles)
!! @param[in]     s2m               profile s2m data: 2m p, 2m t, 2m q, 10m wind u, v, wind-fetch (6,nprofiles)
!! @param[in]     clwscheme         profile clw scheme: clw_scheme, clwde_param (2,nprofiles)
!! @param[in]     icecloud          profile ice scheme: ice_scheme, icede_param (2,nprofiles)
!! @param[in]     p                 pressure profiles (nlevels,nprofiles)
!! @param[in]     t                 temperature profiles (nlevels,nprofiles)
!! @param[in]     gas_units         profile gas units
!! @param[in]     mmr_cldaer        profile mmr_cldaer flag
!! @param[in]     gas_id            list of gas IDs specified in gases argument (ngases)
!! @param[in]     gases             gas, aerosol and cloud profiles (nlevels,nprofiles,ngases)
!! @param[in]     aer_phangle       aerosol phase fn angles (aer_nphangle)
!! @param[in]     aer_asb           aerosol abs, sca, bpr parameters (nlayers,nchannels,nprofiles,3)
!! @param[in]     aer_legcoef       aerosol phase function Legendre coefficients (aer_nmom+1,nlayers,nchannels,nprofiles)
!! @param[in]     aer_pha           aerosol phase fns (aer_nphangle,nlayers,nchannels,nprofiles)
!! @param[in]     cld_phangle       cloud phase fn angles (cld_nphangle)
!! @param[in]     cld_asb           cloud abs, sca, bpr parameters (nlayers,nchannels,nprofiles,3)
!! @param[in]     cld_legcoef       cloud phase function Legendre coefficients (cld_nmom+1,nlayers,nchannels,nprofiles)
!! @param[in]     cld_pha           cloud phase fns (cld_nphangle,nlayers,nchannels,nprofiles)
!! @param[in]     surfemisrefl      input/output surface emissivities, reflectances, specularities (nchannels,nprofiles,4)
!! @param[in,out] btrefl            returned BTs/reflectances (nchannels,nprofiles)
!! @param[in,out] rads              returned radiances (nchannels,nprofiles)
!! @param[in]     nchannels         number of channels (size of channel_list, not required in Python)
!! @param[in]     ngases            number of gases (size of gas_id, not required in Python)
!! @param[in]     nlevels           number of levels (size of p(:,1), not required in Python)
!! @param[in]     nprofiles         number of profiles (size of p(1,:), not required in Python)
!! @param[in]     aer_nphangle      number of aerosol phase angles (size of aer_phangle, not required in Python)
!! @param[in]     aer_nmom          number of aerosol Legendre moments excluding zeroth moment (size of 
!!                                      aer_legcoef(:,1,1,1)-1, not required in Python)
!! @param[in]     cld_nphangle      number of cloud phase angles (size of cld_phangle, not required in Python)
!! @param[in]     cld_nmom          number of cloud Legendre moments excluding zeroth moment (size of
!!                                      cld_legcoef(:,1,1,1)-1, not required in Python)
SUBROUTINE rttov_visir_scatt_call_direct( &
    err,               &
    inst_id,           &
    channel_list,      &
    datetimes,         &    ! profile dates/times                                     (6,nprofiles)
    angles,            &    ! satzen, satazi, sunzen, sunazi angles                   (4,nprofiles)
    surfgeom,          &    ! lat, lon, elevation                                     (3,nprofiles)
    surftype,          &    ! surftype, watertype                                     (2,nprofiles)
    skin,              &    ! skin T, salinity, snow_frac, foam_frac, fastem_coefsx5  (9,nprofiles)
    s2m,               &    ! 2m p, 2m t, 2m q, 10m wind u, v, wind-fetch             (6,nprofiles)
    clwscheme,         &    ! clw_scheme, clwde_param                                 (2,nprofiles)
    icecloud,          &    ! ice_scheme, icede_param                                 (2,nprofiles)
    p,                 &    ! pressure                                                (nlevels,nprofiles)
    t,                 &    ! temperature                                             (nlevels,nprofiles)
    gas_units,         &    ! units for gas profiles
    mmr_cldaer,        &    ! set mmr_cldaer flag
    gas_id,            &    ! gas ID list                                             (ngases)
    gases,             &    ! gas profiles                                            (nlevels,nprofiles,ngases)
    aer_phangle,       &    ! aerosol phase fn angles                                 (aer_nphangle)
    aer_asb,           &    ! aerosol abs, sca, bpr parameters                        (nlayers,nchannels,nprofiles,3)
    aer_legcoef,       &    ! aerosol phase fn Legendre coefficients                  (aer_nmom+1,nlayers,nchannels,nprofiles)
    aer_pha,           &    ! aerosol phase fns                                       (aer_nphangle,nlayers,nchannels,nprofiles)
    cld_phangle,       &    ! cloud phase fn angles                                   (cld_nphangle)
    cld_asb,           &    ! cloud abs, sca, bpr parameters                          (nlayers,nchannels,nprofiles,3)
    cld_legcoef,       &    ! cloud phase fn Legendre coefficients                    (cld_nmom+1,nlayers,nchannels,nprofiles)
    cld_pha,           &    ! cloud phase fns                                         (cld_nphangle,nlayers,nchannels,nprofiles)
    surfemisrefl,      &    ! surface emissivities, reflectances, specularities       (nchannels,nprofiles,4)
    btrefl,            &    ! output BTs/refls (for thermal/solar chans)              (nchannels,nprofiles)
    rads,              &    ! output radiances                                        (nchannels,nprofiles)
    nchannels, ngases, nlevels, nprofiles, aer_nphangle, aer_nmom, cld_nphangle, cld_nmom)
!f2py threadsafe
!f2py intent(hide):: nchannels=len(channel_list)
!f2py intent(hide):: ngases=len(gas_id)
!f2py intent(hide):: nlevels=shape(p, 1)
!f2py intent(hide):: nprofiles=shape(p, 2)
!f2py intent(hide):: aer_nphangle=len(aer_phangle)
!f2py intent(hide):: aer_nmom=shape(aer_legcoef, 1)-1
!f2py intent(hide):: cld_nphangle=len(cld_phangle)
!f2py intent(hide):: cld_nmom=shape(cld_legcoef, 1)-1

!
! Prepares input profiles, calls RTTOV and returns radiances
!

#include "throw.h"

  USE rttov_wrapper_handle

  USE rttov_wrapper_transfer

  USE parkind1, ONLY : jpim, jprb, jplm

  USE rttov_types, ONLY : &
      rttov_chanprof,     &
      rttov_emissivity,   &
      rttov_reflectance,  &
      rttov_profile,      &
      rttov_radiance,     &
      rttov_transmission, &
      rttov_opt_param

  USE rttov_const, ONLY : sensor_id_ir, sensor_id_hi

  IMPLICIT NONE

  INTEGER(jpim),    INTENT(OUT)   :: err
  INTEGER(jpim),    INTENT(IN)    :: nchannels, ngases, nlevels, nprofiles
  INTEGER(jpim),    INTENT(IN)    :: aer_nphangle, aer_nmom, cld_nphangle, cld_nmom

  INTEGER(jpim),    INTENT(IN)    :: inst_id
  INTEGER(jpim),    INTENT(IN)    :: channel_list(nchannels)
  INTEGER(jpim),    INTENT(IN)    :: datetimes(6,nprofiles)
  REAL(jprb),       INTENT(IN)    :: angles(4,nprofiles)
  REAL(jprb),       INTENT(IN)    :: surfgeom(3,nprofiles)
  INTEGER(jpim),    INTENT(IN)    :: surftype(2,nprofiles)
  REAL(jprb),       INTENT(IN)    :: skin(9,nprofiles)
  REAL(jprb),       INTENT(IN)    :: s2m(6,nprofiles)
  INTEGER(jpim),    INTENT(IN)    :: clwscheme(2,nprofiles)
  INTEGER(jpim),    INTENT(IN)    :: icecloud(2,nprofiles)
  REAL(jprb),       INTENT(IN)    :: p(nlevels,nprofiles)
  REAL(jprb),       INTENT(IN)    :: t(nlevels,nprofiles)
  INTEGER(jpim),    INTENT(IN)    :: gas_units
  INTEGER(jpim),    INTENT(IN)    :: mmr_cldaer
  INTEGER(jpim),    INTENT(IN)    :: gas_id(ngases)
  REAL(jprb),       INTENT(IN)    :: gases(nlevels,nprofiles,ngases)
  REAL(jprb),       INTENT(IN)    :: aer_phangle(aer_nphangle)
  REAL(jprb),       INTENT(IN)    :: aer_asb(nlevels-1,nchannels,nprofiles,3)
  REAL(jprb),       INTENT(IN)    :: aer_legcoef(aer_nmom+1,nlevels-1,nchannels,nprofiles)
  REAL(jprb),       INTENT(IN)    :: aer_pha(aer_nphangle,nlevels-1,nchannels,nprofiles)
  REAL(jprb),       INTENT(IN)    :: cld_phangle(cld_nphangle)
  REAL(jprb),       INTENT(IN)    :: cld_asb(nlevels-1,nchannels,nprofiles,3)
  REAL(jprb),       INTENT(IN)    :: cld_legcoef(cld_nmom+1,nlevels-1,nchannels,nprofiles)
  REAL(jprb),       INTENT(IN)    :: cld_pha(cld_nphangle,nlevels-1,nchannels,nprofiles)
  REAL(jprb),       INTENT(INOUT) :: surfemisrefl(nchannels,nprofiles,4)

  REAL(jprb),       INTENT(INOUT) :: btrefl(nchannels,nprofiles)
  REAL(jprb),       INTENT(INOUT) :: rads(nchannels,nprofiles)


  INTEGER(jpim)      :: asw
  INTEGER(jpim)      :: i, k, iprof, iprof1, lochan, hichan
  INTEGER(jpim)      :: nchanprof, thisnchanprof, nprof, thisnprof, ncalls
  LOGICAL(jplm)      :: lmmr_cldaer
  CHARACTER(LEN=256) :: msg

  TYPE(rttov_chanprof),    POINTER :: chanprof(:)
  LOGICAL(jplm),           POINTER :: calcemis(:)
  TYPE(rttov_emissivity),  POINTER :: emissivity(:)
  LOGICAL(jplm),           POINTER :: calcrefl(:)
  TYPE(rttov_reflectance), POINTER :: reflectance(:)
  TYPE(rttov_profile),     POINTER :: profiles(:)
  TYPE(rttov_transmission)         :: transmission
  TYPE(rttov_radiance)             :: radiance
  TYPE(rttov_opt_param)            :: aer_opt_param
  TYPE(rttov_opt_param)            :: cld_opt_param

  TYPE(rttovwrapperhandle_type), POINTER :: rth

#include "rttov_alloc_direct.interface"
#include "rttov_init_emis_refl.interface"
#include "rttov_init_opt_param.interface"
#include "rttov_direct.interface"
#include "rttov_parallel_direct.interface"
#include "rttov_errorreport.interface"

!------------------------------------------------------------------------------

  TRY

  NULLIFY(chanprof,    &
          calcemis,    &
          emissivity,  &
          calcrefl,    &
          reflectance, &
          profiles)

  CALL rttov_wrapper_handle_get_inst(err, inst_id, rth)
  THROW(err.NE.0)

  IF (.NOT. ((rth%opts%rt_ir%addclouds .AND. rth%opts%rt_ir%user_cld_opt_param) .OR. &
             (rth%opts%rt_ir%addaerosl .AND. rth%opts%rt_ir%user_aer_opt_param))) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'This subroutine is for user specified cloud and/or aerosol properties only')
  ENDIF

  IF (.NOT. (rth%coefs%coef%id_sensor == sensor_id_ir .OR. &
             rth%coefs%coef%id_sensor == sensor_id_hi)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'This subroutine is for visible/IR sensors only')
  ENDIF

  ! Deallocate (if necessary) and re-allocate the output structures
  CALL rttov_wrapper_handle_alloc(err, rth, nprofiles, nchannels, nlevels, .FALSE._jplm, 0_jpim)
  THROW(err.NE.0)
  CALL rttov_wrapper_handle_alloc(err, rth, nprofiles, nchannels, nlevels, .FALSE._jplm, 1_jpim)
  THROW(err.NE.0)

  ! nchannels     - input argument, number of channels simulated for every profile
  ! nprofiles     - input argument, total number of profiles passed to this subroutine
  ! ncalls        - number of calls to RTTOV
  ! nprof         - number of profiles per call to RTTOV
  ! thisnprof     - number of profiles in current call (may differ to nprof on last call)
  ! nchanprof     - size of chanprof per call to RTTOV
  ! thisnchanprof - size of chanprof for current call (may differ to nchanprof on last call)

  nprof = MIN(nprofiles, rth%nprofs_per_call)
  nchanprof = nprof * nchannels
  lmmr_cldaer = (mmr_cldaer .NE. 0)

  ! Allocate RTTOV data
  asw = 1_jpim
  CALL rttov_alloc_direct( &
              err,                           &
              asw,                           &
              nprof,                         &
              nchanprof,                     &
              nlevels,                       &
              chanprof,                      &
              rth%opts,                      &
              profiles,                      &
              rth%coefs,                     &
              transmission,                  &
              radiance,                      &
              calcemis      = calcemis,      &
              emissivity    = emissivity,    &
              calcrefl      = calcrefl,      &
              reflectance   = reflectance,   &
              aer_maxnmom   = aer_nmom,      &
              aer_nphangle  = aer_nphangle,  &
              aer_opt_param = aer_opt_param, &
              cld_maxnmom   = cld_nmom,      &
              cld_nphangle  = cld_nphangle,  &
              cld_opt_param = cld_opt_param, &
              init          = .TRUE._jplm)
  THROWM(err.NE.0, 'Error allocating RTTOV structures')

  ! Populate input arrays which don't change
  DO i = 1, nprof
    lochan = (i - 1) * nchannels + 1
    hichan = i * nchannels
    chanprof(lochan:hichan)%prof = i
    chanprof(lochan:hichan)%chan = channel_list(:)
  ENDDO

  IF (rth%verbose) THEN
    WRITE(msg, '(a,i4,a,i8)') &
      'Running RTTOV using nthreads = ', rth%nthreads, &
      ' and nprofs_per_call = ', rth%nprofs_per_call
    INFO(msg)
  ENDIF

  ! Main loop over batches of profiles

  ncalls = nprofiles / rth%nprofs_per_call
  IF (MOD(nprofiles, rth%nprofs_per_call) /= 0) ncalls = ncalls + 1
  thisnchanprof = nchanprof ! This may be overwritten in the final loop below
  thisnprof = rth%nprofs_per_call

  DO i = 1, ncalls

    ! Calculate low/high profile and chanprof indices
    iprof1 = (i - 1) * rth%nprofs_per_call + 1
    IF (iprof1 + rth%nprofs_per_call > nprofiles) THEN
      thisnprof = nprofiles - iprof1 + 1
      thisnchanprof = thisnprof * nchannels
    ENDIF

    ! Populate the profile structure
    CALL rttov_copy_to_profiles(  &
        rth,                      &
        profiles(1:thisnprof),    &
        iprof1,                   &
        datetimes,                &
        angles, surfgeom,         &
        surftype, skin, s2m,      &
        p, t, gas_units,          &
        gas_id, gases,            &
        lmmr_cldaer,              &
        clwscheme = clwscheme,    &
        icecloud = icecloud)

    CALL rttov_copy_to_opt_param( &
        rth,                                  &
        aer_opt_param, cld_opt_param,         &
        iprof1, nprofiles,                    &
        thisnprof, nchannels, nlevels-1_jpim, &
        aer_nmom, aer_nphangle,               &
        aer_phangle,                          &
        aer_asb, aer_legcoef, aer_pha,        &
        cld_nmom, cld_nphangle,               &
        cld_phangle,                          &
        cld_asb, cld_legcoef, cld_pha)

    IF (i == 1 .AND. rth%opts%rt_ir%addsolar) THEN
      IF (rth%opts%rt_ir%addaerosl) THEN
        CALL rttov_init_opt_param(err, rth%opts, aer_opt_param)
        THROWM(err.NE.0, 'Error in rttov_init_opt_param for aerosols')
      ENDIF
      IF (rth%opts%rt_ir%addclouds) THEN
        CALL rttov_init_opt_param(err, rth%opts, cld_opt_param)
        THROWM(err.NE.0, 'Error in rttov_init_opt_param for clouds')
      ENDIF
    ENDIF

    ! Set surface emissivities and BRDFs: if user passed any non-negative
    ! emissivities/BRDFs then use these and set calcemis/calcrefl FALSE
    calcemis(:) = .TRUE.
    calcrefl(:) = .TRUE.
    CALL rttov_init_emis_refl(emissivity, reflectance)

    DO iprof = 1, thisnprof
      k = iprof1 + iprof - 1
      lochan = (iprof - 1) * nchannels + 1
      hichan = iprof * nchannels
      WHERE (surfemisrefl(:,k,1) >= 0._jprb)
        emissivity(lochan:hichan)%emis_in = surfemisrefl(:,k,1)
        calcemis(lochan:hichan) = .FALSE.
      ENDWHERE
      WHERE (surfemisrefl(:,k,2) >= 0._jprb)
        reflectance(lochan:hichan)%refl_in = surfemisrefl(:,k,2)
        calcrefl(lochan:hichan) = .FALSE.
      ENDWHERE
      WHERE (surfemisrefl(:,k,3) >= 0._jprb)
        reflectance(lochan:hichan)%diffuse_refl_in = surfemisrefl(:,k,3)
      ENDWHERE
      WHERE (surfemisrefl(:,k,4) >= 0._jprb)
        emissivity(lochan:hichan)%specularity = surfemisrefl(:,k,4)
      ENDWHERE
    ENDDO

    ! Run RTTOV
    IF (rth%nthreads <= 1) THEN
      CALL rttov_direct( &
          err,                                          &
          chanprof(1:thisnchanprof),                    &
          rth%opts,                                     &
          profiles(1:thisnprof),                        &
          rth%coefs,                                    &
          transmission,                                 &
          radiance,                                     &
          calcemis      = calcemis(1:thisnchanprof),    &
          emissivity    = emissivity(1:thisnchanprof),  &
          calcrefl      = calcrefl(1:thisnchanprof),    &
          reflectance   = reflectance(1:thisnchanprof), &
          aer_opt_param = aer_opt_param,                &
          cld_opt_param = cld_opt_param)
    ELSE
      CALL rttov_parallel_direct( &
          err,                                          &
          chanprof(1:thisnchanprof),                    &
          rth%opts,                                     &
          profiles(1:thisnprof),                        &
          rth%coefs,                                    &
          transmission,                                 &
          radiance,                                     &
          calcemis      = calcemis(1:thisnchanprof),    &
          emissivity    = emissivity(1:thisnchanprof),  &
          calcrefl      = calcrefl(1:thisnchanprof),    &
          reflectance   = reflectance(1:thisnchanprof), &
          aer_opt_param = aer_opt_param,                &
          cld_opt_param = cld_opt_param,                &
          nthreads      = rth%nthreads)
    ENDIF
    THROWM(err.NE.0, 'Error running RTTOV')


    ! Store results
    DO iprof = 1, thisnprof
      k = iprof1 + iprof - 1
      lochan = (iprof - 1) * nchannels + 1
      hichan = iprof * nchannels
      WHERE (rth%coefs%coef%ss_val_chn(channel_list(:)) < 2)
        btrefl(:,k) = radiance%bt(lochan:hichan)
      ELSEWHERE
        btrefl(:,k) = radiance%refl(lochan:hichan)
      ENDWHERE
      rads(:,k) = radiance%total(lochan:hichan)
      surfemisrefl(:,k,1) = emissivity(lochan:hichan)%emis_out
      surfemisrefl(:,k,2) = reflectance(lochan:hichan)%refl_out
      surfemisrefl(:,k,3) = reflectance(lochan:hichan)%diffuse_refl_out
    ENDDO

    lochan = (iprof1 - 1) * nchannels + 1
    hichan = (iprof1 + thisnprof - 1) * nchannels

    IF (rth%store_trans) THEN
      rth%transmission%tau_total(lochan:hichan)             = transmission%tau_total(1:thisnchanprof)
      rth%transmission%tau_levels(:,lochan:hichan)          = transmission%tau_levels(:,1:thisnchanprof)
      rth%transmission%tausun_total_path2(lochan:hichan)    = transmission%tausun_total_path2(1:thisnchanprof)
      rth%transmission%tausun_levels_path2(:,lochan:hichan) = transmission%tausun_levels_path2(:,1:thisnchanprof)
      rth%transmission%tausun_total_path1(lochan:hichan)    = transmission%tausun_total_path1(1:thisnchanprof)
      rth%transmission%tausun_levels_path1(:,lochan:hichan) = transmission%tausun_levels_path1(:,1:thisnchanprof)
      rth%transmission%tau_total_cld(lochan:hichan)         = transmission%tau_total_cld(1:thisnchanprof)
      rth%transmission%tau_levels_cld(:,lochan:hichan)      = transmission%tau_levels_cld(:,1:thisnchanprof)
    ENDIF

    IF (rth%store_rad) THEN
      rth%radiance%clear(lochan:hichan)              = radiance%clear(1:thisnchanprof)
      rth%radiance%total(lochan:hichan)              = radiance%total(1:thisnchanprof)
      rth%radiance%bt_clear(lochan:hichan)           = radiance%bt_clear(1:thisnchanprof)
      rth%radiance%bt(lochan:hichan)                 = radiance%bt(1:thisnchanprof)
      rth%radiance%refl_clear(lochan:hichan)         = radiance%refl_clear(1:thisnchanprof)
      rth%radiance%refl(lochan:hichan)               = radiance%refl(1:thisnchanprof)
      rth%radiance%cloudy(lochan:hichan)             = radiance%cloudy(1:thisnchanprof)
      rth%radiance%overcast(:,lochan:hichan)         = radiance%overcast(:,1:thisnchanprof)
      rth%radiance%quality(lochan:hichan)            = radiance%quality(1:thisnchanprof)
      rth%radiance%plane_parallel                    = radiance%plane_parallel
      rth%radiance%geometric_height(:,lochan:hichan) = radiance%geometric_height(:,1:thisnchanprof)
    ENDIF

  ENDDO ! loop over calls to RTTOV

  ! Deallocate RTTOV data
  asw = 0_jpim
  CALL rttov_alloc_direct( &
              err,                           &
              asw,                           &
              nprof,                         &
              nchanprof,                     &
              nlevels,                       &
              chanprof,                      &
              rth%opts,                      &
              profiles,                      &
              rth%coefs,                     &
              transmission,                  &
              radiance,                      &
              calcemis      = calcemis,      &
              emissivity    = emissivity,    &
              calcrefl      = calcrefl,      &
              reflectance   = reflectance,   &
              aer_maxnmom   = aer_nmom,      &
              aer_nphangle  = aer_nphangle,  &
              aer_opt_param = aer_opt_param, &
              cld_maxnmom   = cld_nmom,      &
              cld_nphangle  = cld_nphangle,  &
              cld_opt_param = cld_opt_param, &
              init          = .TRUE._jplm)
  THROWM(err.NE.0, 'Error deallocating RTTOV structures')

  CATCH
END SUBROUTINE rttov_visir_scatt_call_direct


!> Call RTTOV-SCATT direct model.
!! @param[out]    err               return status
!! @param[in]     inst_id           ID for instrument to simulate
!! @param[in]     channel_list      list of channels to simulate for for each profile
!! @param[in]     datetimes         profile dates/times: year, month, day, hour, min, sec (6,nprofiles)
!! @param[in]     angles            profile angles: zenang, aziang (2,nprofiles)
!! @param[in]     surfgeom          profile surface geometry: lat, lon, elevation (3,nprofiles)
!! @param[in]     surftype          profile surface type (nprofiles)
!! @param[in]     skin              profile skin data: skin T, salinity, foam_frac, fastem_coefsx5 (8,nprofiles)
!! @param[in]     s2m               profile s2m data: 2m p, 2m t, 2m q, 10m wind u, v (5,nprofiles)
!! @param[in]     zeeman            profile B-field data: Be, cosbk (2,nprofiles)
!! @param[in]     p                 pressure profiles (nlevels,nprofiles)
!! @param[in]     t                 temperature profiles (nlevels,nprofiles)
!! @param[in]     gas_units         profile gas units
!! @param[in]     gas_id            list of gas IDs specified in gases argument (ngases)
!! @param[in]     gases             gas, cloud and hydrometeor profiles (nlevels,nprofiles,ngases)
!! @param[in]     ph                pressure half-level profiles (nlevels+1,nprofiles)
!! @param[in]     cfrac             user cloud fraction (nprofiles)
!! @param[in]     multi_hydro_frac  false => single hydro_frac profile, true => one hydro_frac profile per hydrometeor
!! @param[in]     calc_zef          enable/disable the radar reflectivity calculations
!! @param[in]     surfemis          input/output surface emissivities (nchannels,nprofiles)
!! @param[in,out] bt                returned BTs (nchannels,nprofiles)
!! @param[in]     nchannels         number of channels (size of channel_list, not required in Python)
!! @param[in]     ngases            number of gases (size of gas_id, not required in Python)
!! @param[in]     nlevels           number of levels (size of p(:,1), not required in Python)
!! @param[in]     nprofiles         number of profiles (size of p(1,:), not required in Python)
SUBROUTINE rttov_scatt_call_direct( &
    err,               &
    inst_id,           &
    channel_list,      &
    datetimes,         &    ! profile dates/times                                     (6,nprofiles)
    angles,            &    ! satzen, satazi angles                                   (2,nprofiles)
    surfgeom,          &    ! lat, lon, elevation                                     (3,nprofiles)
    surftype,          &    ! surftype                                                (nprofiles)
    skin,              &    ! skin T, salinity, foam_frac, fastem_coefsx5             (8,nprofiles)
    s2m,               &    ! 2m p, 2m t, 2m q, 10m wind u, v                         (5,nprofiles)
    zeeman,            &    ! Be, cosbk                                               (2,nprofiles)
    p,                 &    ! pressure                                                (nlevels,nprofiles)
    t,                 &    ! temperature                                             (nlevels,nprofiles)
    gas_units,         &    ! units for gas profiles
    gas_id,            &    ! gas ID list                                             (ngases)
    gases,             &    ! gas profiles                                            (nlevels,nprofiles,ngases)
    ph,                &    ! pressure half-levels                                    (nlevels+1,nprofiles)
    cfrac,             &    ! user cloud fraction                                     (nprofiles)
    multi_hydro_frac,  &    ! false => single hydro_frac profile, true => one hydro_frac profile per hydrometeor
    calc_zef,          &    ! enable/disable the radar reflectivity calculations
    surfemis,          &    ! surface emissivities                                    (nchannels,nprofiles)
    bt,                &    ! output BTs                                              (nchannels,nprofiles)
    nchannels, ngases, nlevels, nprofiles)
!f2py threadsafe
!f2py intent(hide):: nchannels=len(channel_list)
!f2py intent(hide):: ngases=len(gas_id)
!f2py intent(hide):: nlevels=shape(p, 1)
!f2py intent(hide):: nprofiles=shape(p, 2)

!
! Prepares input profiles, calls RTTOV-SCATT and returns radiances
!

#include "throw.h"

  USE rttov_wrapper_handle

  USE rttov_wrapper_transfer

  USE parkind1, ONLY : jpim, jprb, jplm

  USE rttov_types, ONLY :  &
      rttov_chanprof,      &
      rttov_emissivity,    &
      rttov_profile,       &
      rttov_profile_cloud, &
      rttov_radiance,      &
      rttov_reflectivity,  &
      rttov_scatt_emis_retrieval_type

  USE rttov_const, ONLY : sensor_id_mw, sensor_id_po

  IMPLICIT NONE

  INTEGER(jpim),    INTENT(OUT)   :: err
  INTEGER(jpim),    INTENT(IN)    :: nchannels, ngases, nlevels, nprofiles

  INTEGER(jpim),    INTENT(IN)    :: inst_id
  INTEGER(jpim),    INTENT(IN)    :: channel_list(nchannels)
  INTEGER(jpim),    INTENT(IN)    :: datetimes(6,nprofiles)
  REAL(jprb),       INTENT(IN)    :: angles(2,nprofiles)
  REAL(jprb),       INTENT(IN)    :: surfgeom(3,nprofiles)
  INTEGER(jpim),    INTENT(IN)    :: surftype(nprofiles)
  REAL(jprb),       INTENT(IN)    :: skin(8,nprofiles)
  REAL(jprb),       INTENT(IN)    :: s2m(5,nprofiles)
  REAL(jprb),       INTENT(IN)    :: zeeman(2,nprofiles)
  REAL(jprb),       INTENT(IN)    :: p(nlevels,nprofiles)
  REAL(jprb),       INTENT(IN)    :: t(nlevels,nprofiles)
  INTEGER(jpim),    INTENT(IN)    :: gas_units
  INTEGER(jpim),    INTENT(IN)    :: gas_id(ngases)
  REAL(jprb),       INTENT(IN)    :: gases(nlevels,nprofiles,ngases)
  REAL(jprb),       INTENT(IN)    :: ph(nlevels+1,nprofiles)
  REAL(jprb),       INTENT(INOUT) :: cfrac(nprofiles)
  LOGICAL(jpim),    INTENT(IN)    :: multi_hydro_frac
  LOGICAL(jpim),    INTENT(IN)    :: calc_zef

  REAL(jprb),       INTENT(INOUT) :: surfemis(nchannels,nprofiles)

  REAL(jprb),       INTENT(INOUT) :: bt(nchannels,nprofiles)


  INTEGER(jpim)      :: asw
  INTEGER(jpim)      :: i, k, iprof, iprof1, lochan, hichan
  INTEGER(jpim)      :: nchanprof, thisnchanprof, nprof, thisnprof, ncalls
  INTEGER(jpim)      :: nhydro_frac
  CHARACTER(LEN=256) :: msg

  TYPE(rttov_chanprof),         POINTER :: chanprof(:)
  INTEGER(KIND=jpim),           POINTER :: frequencies(:)
  LOGICAL(KIND=jplm),           POINTER :: use_chan(:,:)
  LOGICAL(jplm),                POINTER :: calcemis(:)
  TYPE(rttov_emissivity),       POINTER :: emissivity(:)
  TYPE(rttov_profile),          POINTER :: profiles(:)
  TYPE(rttov_profile_cloud),    POINTER :: cld_profiles(:)
  TYPE(rttov_radiance)                  :: radiance
  TYPE(rttov_reflectivity) ,    POINTER :: reflectivity
  REAL(jprb)                            :: cfrac_out(nprofiles)
  TYPE(rttov_scatt_emis_retrieval_type) :: emis_terms

  TYPE(rttovwrapperhandle_type), POINTER :: rth

#include "rttov_alloc_direct.interface"
#include "rttov_scatt_setupindex.interface"
#include "rttov_init_emis_refl.interface"
#include "rttov_scatt.interface"
#include "rttov_parallel_scatt.interface"
#include "rttov_errorreport.interface"

!------------------------------------------------------------------------------

  TRY

  NULLIFY(chanprof,     &
          frequencies,  &
          use_chan,     &
          calcemis,     &
          emissivity,   &
          profiles,     &
          cld_profiles, &
          reflectivity)

  ! If reflectivity is nullified, it tells the routines below that radar
  ! calculations are turned off
  IF (calc_zef) ALLOCATE(reflectivity)

  CALL rttov_wrapper_handle_get_inst(err, inst_id, rth)
  THROW(err.NE.0)

  IF (.NOT. (rth%coefs%coef%id_sensor == sensor_id_mw .OR. &
             rth%coefs%coef%id_sensor == sensor_id_po)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'RTTOV-SCATT is for MW sensors only')
  ENDIF

  ! Deallocate (if necessary) and re-allocate the output structures
  CALL rttov_wrapper_handle_alloc(err, rth, nprofiles, nchannels, nlevels, calc_zef, 0_jpim)
  THROW(err.NE.0)
  CALL rttov_wrapper_handle_alloc(err, rth, nprofiles, nchannels, nlevels, calc_zef, 1_jpim)
  THROW(err.NE.0)

  ! nchannels     - input argument, number of channels simulated for every profile
  ! nprofiles     - input argument, total number of profiles passed to this subroutine
  ! ncalls        - number of calls to RTTOV
  ! nprof         - number of profiles per call to RTTOV
  ! thisnprof     - number of profiles in current call (may differ to nprof on last call)
  ! nchanprof     - size of chanprof per call to RTTOV
  ! thisnchanprof - size of chanprof for current call (may differ to nchanprof on last call)

  nprof = MIN(nprofiles, rth%nprofs_per_call)
  nchanprof = nprof * nchannels

  nhydro_frac = 1
  IF (multi_hydro_frac) nhydro_frac = rth%coefs_scatt%nhydro

  ! Allocate RTTOV data
  asw = 1_jpim
  CALL rttov_alloc_direct( &
              err,                                    &
              asw,                                    &
              nprof,                                  &
              nchanprof,                              &
              nlevels,                                &
              chanprof,                               &
              rth%opts,                               &
              profiles,                               &
              rth%coefs,                              &
              radiance             = radiance,        &
              calcemis             = calcemis,        &
              emissivity           = emissivity,      &
              frequencies          = frequencies,     &
              coef_scatt           = rth%coefs_scatt, &
              nhydro_frac          = nhydro_frac,     &
              cld_profiles         = cld_profiles,    &
              emis_retrieval_terms = emis_terms,      &
              reflectivity         = reflectivity,    &
              init                 = .TRUE._jplm)
  THROWM(err.NE.0, 'Error allocating RTTOV structures')

  ! use_chan array is dimensioned by the total number of instrument channels
  ALLOCATE(use_chan(nprof,rth%coefs%coef%fmv_chn))

  ! Set use_chan to .TRUE. only for required channels
  use_chan(:,:) = .FALSE._jplm
  DO i = 1, nprof
    use_chan(i,channel_list(1:nchannels)) = .TRUE._jplm
  ENDDO

  ! Populate chanprof and frequencies arrays
  CALL rttov_scatt_setupindex ( &
      err,                    &
      nprof,                  &
      rth%coefs%coef%fmv_chn, &
      rth%coefs,              &
      rth%coefs_scatt,        &
      nchanprof,              &
      chanprof,               &
      frequencies,            &
      use_chan)
  THROWM(err.NE.0, 'Error determining channels, frequencies and polarisations from coefficient files')

  IF (rth%verbose) THEN
    WRITE(msg, '(a,i4,a,i8)') &
      'Running RTTOV-SCATT using nthreads = ', rth%nthreads, &
      ' and nprofs_per_call = ', rth%nprofs_per_call
    INFO(msg)
  ENDIF

  ! Main loop over batches of profiles

  ncalls = nprofiles / rth%nprofs_per_call
  IF (MOD(nprofiles, rth%nprofs_per_call) /= 0) ncalls = ncalls + 1
  thisnchanprof = nchanprof ! This may be overwritten in the final loop below
  thisnprof = rth%nprofs_per_call

  DO i = 1, ncalls

    ! Calculate low/high profile and chanprof indices
    iprof1 = (i - 1) * rth%nprofs_per_call + 1
    IF (iprof1 + rth%nprofs_per_call > nprofiles) THEN
      thisnprof = nprofiles - iprof1 + 1
      thisnchanprof = thisnprof * nchannels
    ENDIF

    ! Populate the profile structure
    CALL rttov_scatt_copy_to_profiles(  &
        rth,                       &
        profiles(1:thisnprof),     &
        cld_profiles(1:thisnprof), &
        iprof1,                    &
        datetimes,                 &
        angles, surfgeom,          &
        surftype, skin, s2m,       &
        zeeman,                    &
        p, t, gas_units,           &
        gas_id, gases, ph, cfrac)

    ! Set surface emissivities: if user passed any non-negative
    ! emissivities then use these and set calcemis FALSE
    calcemis(:) = .TRUE.
    CALL rttov_init_emis_refl(emissivity)

    DO iprof = 1, thisnprof
      k = iprof1 + iprof - 1
      lochan = (iprof - 1) * nchannels + 1
      hichan = iprof * nchannels
      WHERE (surfemis(:,k) >= 0._jprb)
        emissivity(lochan:hichan)%emis_in = surfemis(:,k)
        calcemis(lochan:hichan) = .FALSE.
      ENDWHERE
    ENDDO

    ! Run RTTOV-SCATT
    IF (rth%nthreads <= 1) THEN
      CALL rttov_scatt ( &
          err,                            &
          rth%opts_scatt,                 &
          nlevels,                        &
          chanprof(1:thisnchanprof),      &
          frequencies(1:thisnchanprof),   &
          profiles(1:thisnprof),          &
          cld_profiles(1:thisnprof),      &
          rth%coefs,                      &
          rth%coefs_scatt,                &
          calcemis(1:thisnchanprof),      &
          emissivity(1:thisnchanprof),    &
          radiance,                       &
          cfrac_out(1:thisnprof),         &
          emis_terms,                     &
          reflectivity)
    ELSE
      IF (calc_zef) THEN
        CALL rttov_parallel_scatt ( &
            err,                            &
            rth%opts_scatt,                 &
            nlevels,                        &
            chanprof(1:thisnchanprof),      &
            frequencies(1:thisnchanprof),   &
            profiles(1:thisnprof),          &
            cld_profiles(1:thisnprof),      &
            rth%coefs,                      &
            rth%coefs_scatt,                &
            calcemis(1:thisnchanprof),      &
            emissivity(1:thisnchanprof),    &
            radiance,                       &
            cfrac_out(1:thisnprof),         &
            reflectivity = reflectivity,    &
            nthreads = rth%nthreads)
      ELSE
        CALL rttov_parallel_scatt ( &
            err,                            &
            rth%opts_scatt,                 &
            nlevels,                        &
            chanprof(1:thisnchanprof),      &
            frequencies(1:thisnchanprof),   &
            profiles(1:thisnprof),          &
            cld_profiles(1:thisnprof),      &
            rth%coefs,                      &
            rth%coefs_scatt,                &
            calcemis(1:thisnchanprof),      &
            emissivity(1:thisnchanprof),    &
            radiance,                       &
            cfrac_out(1:thisnprof),         &
            emis_terms,                     &
            nthreads = rth%nthreads)
      ENDIF
    ENDIF
    THROWM(err.NE.0, 'Error running RTTOV-SCATT')

    ! Store results
    DO iprof = 1, thisnprof
      k = iprof1 + iprof - 1
      lochan = (iprof - 1) * nchannels + 1
      hichan = iprof * nchannels
      bt(:,k) = radiance%bt(lochan:hichan)
      surfemis(:,k) = emissivity(lochan:hichan)%emis_out
      IF (.NOT. rth%opts_scatt%lusercfrac) cfrac(k) = cfrac_out(iprof)
    ENDDO

    lochan = (iprof1 - 1) * nchannels + 1
    hichan = (iprof1 + thisnprof - 1) * nchannels

    IF (rth%store_rad) THEN
      rth%radiance%clear(lochan:hichan)              = radiance%clear(1:thisnchanprof)
      rth%radiance%total(lochan:hichan)              = radiance%total(1:thisnchanprof)
      rth%radiance%bt_clear(lochan:hichan)           = radiance%bt_clear(1:thisnchanprof)
      rth%radiance%bt(lochan:hichan)                 = radiance%bt(1:thisnchanprof)
      rth%radiance%quality(lochan:hichan)            = radiance%quality(1:thisnchanprof)
      rth%radiance%plane_parallel                    = radiance%plane_parallel
      rth%radiance%geometric_height(:,lochan:hichan) = radiance%geometric_height(:,1:thisnchanprof)
      IF (ASSOCIATED(reflectivity)) THEN
        rth%reflectivity%zef(:,lochan:hichan)        = reflectivity%zef(:,1:thisnchanprof)
        rth%reflectivity%azef(:,lochan:hichan)       = reflectivity%azef(:,1:thisnchanprof)
      ENDIF
    ENDIF

    IF (rth%store_emis_terms) THEN
      rth%emis_terms%cfrac(lochan:hichan)    = emis_terms%cfrac(1:thisnchanprof)
      rth%emis_terms%bsfc(lochan:hichan)     = emis_terms%bsfc(1:thisnchanprof)
      rth%emis_terms%tau_cld(lochan:hichan)  = emis_terms%tau_cld(1:thisnchanprof)
      rth%emis_terms%up_cld(lochan:hichan)   = emis_terms%up_cld(1:thisnchanprof)
      rth%emis_terms%down_cld(lochan:hichan) = emis_terms%down_cld(1:thisnchanprof)
      rth%emis_terms%tau_clr(lochan:hichan)  = emis_terms%tau_clr(1:thisnchanprof)
      rth%emis_terms%up_clr(lochan:hichan)   = emis_terms%up_clr(1:thisnchanprof)
      rth%emis_terms%down_clr(lochan:hichan) = emis_terms%down_clr(1:thisnchanprof)
    ENDIF

  ENDDO ! loop over calls to RTTOV-SCATT

  ! Deallocate RTTOV data
  asw = 0_jpim
  CALL rttov_alloc_direct( &
              err,                                    &
              asw,                                    &
              nprof,                                  &
              nchanprof,                              &
              nlevels,                                &
              chanprof,                               &
              rth%opts,                               &
              profiles,                               &
              rth%coefs,                              &
              radiance             = radiance,        &
              calcemis             = calcemis,        &
              emissivity           = emissivity,      &
              frequencies          = frequencies,     &
              coef_scatt           = rth%coefs_scatt, &
              nhydro_frac          = nhydro_frac,     &
              cld_profiles         = cld_profiles,    &
              emis_retrieval_terms = emis_terms,      &
              reflectivity         = reflectivity)
  THROWM(err.NE.0, 'Error deallocating RTTOV structures')

  IF (calc_zef) DEALLOCATE(reflectivity)

  DEALLOCATE(use_chan, stat=err)
  THROWM(err.NE.0, 'Error dellocating RTTOV-SCATT arrays')

  CATCH
END SUBROUTINE rttov_scatt_call_direct


!> Call rttov_k for clear-sky or visible/IR scattering simulations using
!! pre-defined particle types.
!! @param[out]    err               return status
!! @param[in]     inst_id           ID for instrument to simulate
!! @param[in]     channel_list      list of channels to simulate for for each profile
!! @param[in]     datetimes         profile dates/times: year, month, day, hour, min, sec (6,nprofiles)
!! @param[in]     angles            profile angles: zenang, aziang, sunzenang, sunaziang (4,nprofiles)
!! @param[in]     surfgeom          profile surface geometry: lat, lon, elevation (3,nprofiles)
!! @param[in]     surftype          profile surface type: surftype, watertype (2,nprofiles)
!! @param[in]     skin              profile skin data: skin T, salinity, snow_frac, foam_frac, fastem_coefsx5 (9,nprofiles)
!! @param[in,out] skin_k            Jacobians for skin data (9,nchannels,nprofiles)
!! @param[in]     s2m               profile s2m data: 2m p, 2m t, 2m q, 10m wind u, v, wind-fetch (6,nprofiles)
!! @param[in,out] s2m_k             Jacobians for s2m data (6,nchannels,nprofiles)
!! @param[in]     simplecloud       profile simple cloud data: ctp, cfraction (2,nprofiles)
!! @param[in,out] simplecloud_k     Jacobians for simple cloud data (2,nchannels,nprofiles)
!! @param[in]     clwscheme         profile clw scheme: clw_scheme, clwde_param (2,nprofiles)
!! @param[in]     icecloud          profile ice scheme: ice_scheme, icede_param (2,nprofiles)
!! @param[in]     zeeman            profile B-field data: Be, cosbk (2,nprofiles)
!! @param[in]     p                 pressure profiles (nlevels,nprofiles)
!! @param[in,out] p_k               pressure Jacobians (nlevels,nchannels,nprofiles)
!! @param[in]     t                 temperature profiles (nlevels,nprofiles)
!! @param[in,out] t_k               temperature Jacobians (nlevels,nchannels,nprofiles)
!! @param[in]     gas_units         profile gas units
!! @param[in]     mmr_cldaer        profile mmr_cldaer flag
!! @param[in]     gas_id            list of gas IDs specified in gases argument (ngases)
!! @param[in]     gases             gas, aerosol and cloud profiles (nlevels,nprofiles,ngases)
!! @param[in,out] gases_k           Jacobians for contents of gases array (nlevels,nchannels,nprofiles,ngases)
!! @param[in]     surfemisrefl      input/output surface emissivities, reflectances, specularities (nchannels,nprofiles,4)
!! @param[in,out] surfemisrefl_k    emissivity/reflectance/specularity Jacobians (nchannels,nprofiles,4)
!! @param[in,out] btrefl            returned BTs/reflectances (nchannels,nprofiles)
!! @param[in,out] rads              returned radiances (nchannels,nprofiles)
!! @param[in]     bt_k              input BT perturbations (nchannels,nprofiles)
!! @param[in]     rads_k            input radiance perturbations (nchannels,nprofiles)
!! @param[in]     nchannels         number of channels (size of channel_list, not required in Python)
!! @param[in]     ngases            number of gases (size of gas_id, not required in Python)
!! @param[in]     nlevels           number of levels (size of p(:,1), not required in Python)
!! @param[in]     nprofiles         number of profiles (size of p(1,:), not required in Python)
SUBROUTINE rttov_call_k( &
    err,               &
    inst_id,           &
    channel_list,      &
    datetimes,         &    ! profile dates/times                                     (6,nprofiles)
    angles,            &    ! satzen, satazi, sunzen, sunazi angles                   (4,nprofiles)
    surfgeom,          &    ! lat, lon, elevation                                     (3,nprofiles)
    surftype,          &    ! surftype, watertype                                     (2,nprofiles)
    skin,              &    ! skin T, salinity, snow_frac, foam_frac, fastem_coefsx5  (9,nprofiles)
    skin_k,            &    ! skin K                                                  (9,nchannels,nprofiles)
    s2m,               &    ! 2m p, 2m t, 2m q, 10m wind u, v, wind-fetch             (6,nprofiles)
    s2m_k,             &    ! 2m K                                                    (6,nchannels,nprofiles)
    simplecloud,       &    ! ctp, cfraction                                          (2,nprofiles)
    simplecloud_k,     &    ! ctp, cfraction K                                        (2,nchannels,nprofiles)
    clwscheme,         &    ! clw_scheme, clwde_param                                 (2,nprofiles)
    icecloud,          &    ! ice_scheme, icede_param                                 (2,nprofiles)
    zeeman,            &    ! Be, cosbk                                               (2,nprofiles)
    p,                 &    ! pressure                                                (nlevels,nprofiles)
    p_k,               &    ! output pressure K                                       (nlevels,nchannels,nprofiles)
    t,                 &    ! temperature                                             (nlevels,nprofiles)
    t_k,               &    ! output temperature K                                    (nlevels,nchannels,nprofiles)
    gas_units,         &    ! units for gas profiles
    mmr_cldaer,        &    ! set mmr_cldaer flag
    gas_id,            &    ! gas ID list                                             (ngases)
    gases,             &    ! gas profiles                                            (nlevels,nprofiles,ngases)
    gases_k,           &    ! output gas profiles K                                   (nlevels,nchannels,nprofiles,ngases)
    surfemisrefl,      &    ! surface emissivities, reflectances, specularities       (nchannels,nprofiles,4)
    surfemisrefl_k,    &    ! output surface emis, refl, spec K                       (nchannels,nprofiles,4)
    btrefl,            &    ! output BTs/refls (for thermal/solar chans)              (nchannels,nprofiles)
    rads,              &    ! output radiances                                        (nchannels,nprofiles)
    bt_k,              &    ! input BT perturbations (thermal chans only)             (nchannels,nprofiles)
    rads_k,            &    ! input radiance perturbations                            (nchannels,nprofiles)
    nchannels, ngases, nlevels, nprofiles)
!f2py threadsafe
!f2py intent(hide):: nchannels=len(channel_list)
!f2py intent(hide):: ngases=len(gas_id)
!f2py intent(hide):: nlevels=shape(p, 1)
!f2py intent(hide):: nprofiles=shape(p, 2)

!
! Prepares input profiles, calls RTTOV K and returns radiances and profiles_k
!

#include "throw.h"

  USE rttov_wrapper_handle

  USE rttov_wrapper_transfer

  USE parkind1, ONLY : jpim, jprb, jplm

  USE rttov_types, ONLY : &
      rttov_chanprof,     &
      rttov_emissivity,   &
      rttov_reflectance,  &
      rttov_profile,      &
      rttov_radiance,     &
      rttov_radiance2,    &
      rttov_transmission

  IMPLICIT NONE

  INTEGER(jpim),    INTENT(OUT)   :: err

  INTEGER(jpim),    INTENT(IN)    :: nchannels, ngases, nlevels, nprofiles

  INTEGER(jpim),    INTENT(IN)    :: inst_id
  INTEGER(jpim),    INTENT(IN)    :: channel_list(nchannels)
  INTEGER(jpim),    INTENT(IN)    :: datetimes(6,nprofiles)
  REAL(jprb),       INTENT(IN)    :: angles(4,nprofiles)
  REAL(jprb),       INTENT(IN)    :: surfgeom(3,nprofiles)
  INTEGER(jpim),    INTENT(IN)    :: surftype(2,nprofiles)
  REAL(jprb),       INTENT(IN)    :: skin(9,nprofiles)
  REAL(jprb),       INTENT(INOUT) :: skin_k(9,nchannels,nprofiles)
  REAL(jprb),       INTENT(IN)    :: s2m(6,nprofiles)
  REAL(jprb),       INTENT(INOUT) :: s2m_k(6,nchannels,nprofiles)
  REAL(jprb),       INTENT(IN)    :: simplecloud(2,nprofiles)
  REAL(jprb),       INTENT(INOUT) :: simplecloud_k(2,nchannels,nprofiles)
  INTEGER(jpim),    INTENT(IN)    :: clwscheme(2,nprofiles)
  INTEGER(jpim),    INTENT(IN)    :: icecloud(2,nprofiles)
  REAL(jprb),       INTENT(IN)    :: zeeman(2,nprofiles)
  REAL(jprb),       INTENT(IN)    :: p(nlevels,nprofiles)
  REAL(jprb),       INTENT(INOUT) :: p_k(nlevels,nchannels,nprofiles)
  REAL(jprb),       INTENT(IN)    :: t(nlevels,nprofiles)
  REAL(jprb),       INTENT(INOUT) :: t_k(nlevels,nchannels,nprofiles)
  INTEGER(jpim),    INTENT(IN)    :: gas_units
  INTEGER(jpim),    INTENT(IN)    :: mmr_cldaer
  INTEGER(jpim),    INTENT(IN)    :: gas_id(ngases)
  REAL(jprb),       INTENT(IN)    :: gases(nlevels,nprofiles,ngases)
  REAL(jprb),       INTENT(INOUT) :: gases_k(nlevels,nchannels,nprofiles,ngases)

  REAL(jprb),       INTENT(INOUT) :: surfemisrefl(nchannels,nprofiles,4)
  REAL(jprb),       INTENT(INOUT) :: surfemisrefl_k(nchannels,nprofiles,4)

  REAL(jprb),       INTENT(INOUT) :: btrefl(nchannels,nprofiles)
  REAL(jprb),       INTENT(INOUT) :: rads(nchannels,nprofiles)
  REAL(jprb),       INTENT(IN)    :: bt_k(nchannels,nprofiles)
  REAL(jprb),       INTENT(IN)    :: rads_k(nchannels,nprofiles)


  INTEGER(jpim)      :: asw
  INTEGER(jpim)      :: i, j, k, iprof, iprof1, lochan, hichan
  INTEGER(jpim)      :: nchanprof, thisnchanprof, nprof, thisnprof, ncalls
  LOGICAL(jplm)      :: lmmr_cldaer
  CHARACTER(LEN=128) :: msg

  TYPE(rttov_chanprof),    POINTER :: chanprof(:)
  LOGICAL(jplm),           POINTER :: calcemis(:)
  TYPE(rttov_emissivity),  POINTER :: emissivity(:)
  TYPE(rttov_emissivity),  POINTER :: emissivity_k(:)
  LOGICAL(jplm),           POINTER :: calcrefl(:)
  TYPE(rttov_reflectance), POINTER :: reflectance(:)
  TYPE(rttov_reflectance), POINTER :: reflectance_k(:)
  TYPE(rttov_profile),     POINTER :: profiles(:)
  TYPE(rttov_profile),     POINTER :: profiles_k(:)
  TYPE(rttov_transmission)         :: transmission
  TYPE(rttov_transmission)         :: transmission_k
  TYPE(rttov_radiance)             :: radiance
  TYPE(rttov_radiance)             :: radiance_k
  TYPE(rttov_radiance2)            :: radiance2

  TYPE(rttovwrapperhandle_type), POINTER :: rth

#include "rttov_alloc_k.interface"
#include "rttov_init_emis_refl.interface"
#include "rttov_init_rad.interface"
#include "rttov_init_transmission.interface"
#include "rttov_init_prof.interface"
#include "rttov_k.interface"
#include "rttov_parallel_k.interface"
#include "rttov_errorreport.interface"

!------------------------------------------------------------------------------

  TRY

  NULLIFY(chanprof,      &
          calcemis,      &
          emissivity,    &
          emissivity_k,  &
          calcrefl,      &
          reflectance,   &
          reflectance_k, &
          profiles,      &
          profiles_k)

  CALL rttov_wrapper_handle_get_inst(err, inst_id, rth)
  THROW(err.NE.0)

  IF ((rth%opts%rt_ir%addclouds .AND. rth%opts%rt_ir%user_cld_opt_param) .OR. &
      (rth%opts%rt_ir%addaerosl .AND. rth%opts%rt_ir%user_aer_opt_param)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Use rttov_visir_scatt_call_k when providing cloud/aerosol optical properties')
  ENDIF

  ! Deallocate (if necessary) and re-allocate the output structures
  CALL rttov_wrapper_handle_alloc(err, rth, nprofiles, nchannels, nlevels, .FALSE._jplm, 0_jpim)
  THROW(err.NE.0)
  CALL rttov_wrapper_handle_alloc(err, rth, nprofiles, nchannels, nlevels, .FALSE._jplm, 1_jpim)
  THROW(err.NE.0)

  ! nchannels     - input argument, number of channels simulated for every profile
  ! nprofiles     - input argument, total number of profiles passed to this subroutine
  ! ncalls        - number of calls to RTTOV
  ! nprof         - number of profiles per call to RTTOV
  ! thisnprof     - number of profiles in current call (may differ to nprof on last call)
  ! iprof1        - index of first profile in current batch
  ! nchanprof     - size of chanprof per call to RTTOV
  ! thisnchanprof - size of chanprof for current call (may differ to nchanprof on last call)

  nprof = MIN(nprofiles, rth%nprofs_per_call)
  nchanprof = nprof * nchannels
  lmmr_cldaer = (mmr_cldaer .NE. 0)

  ! Allocate RTTOV data
  asw = 1_jpim
  CALL rttov_alloc_k( &
              err,            &
              asw,            &
              nprof,          &
              nchanprof,      &
              nlevels,        &
              chanprof,       &
              rth%opts,       &
              profiles,       &
              profiles_k,     &
              rth%coefs,      &
              transmission,   &
              transmission_k, &
              radiance,       &
              radiance_k,     &
              radiance2,      &
              calcemis,       &
              emissivity,     &
              emissivity_k,   &
              calcrefl,       &
              reflectance,    &
              reflectance_k,  &
              init = .TRUE._jplm)
  THROWM(err.NE.0, 'Error allocating RTTOV structures')

  ! Populate input arrays which don't change
  DO j = 1, nprof
    lochan = (j - 1) * nchannels + 1
    hichan = j * nchannels
    chanprof(lochan:hichan)%prof = j
    chanprof(lochan:hichan)%chan = channel_list(:)
  ENDDO

  IF (rth%verbose) THEN
    WRITE(msg, '(a,i4,a,i8)') &
      'Running RTTOV K using nthreads = ', rth%nthreads, &
      ' and nprofs_per_call = ', rth%nprofs_per_call
    INFO(msg)
  ENDIF

  ! Main loop over batches of profiles

  ncalls = nprofiles / rth%nprofs_per_call
  IF (MOD(nprofiles, rth%nprofs_per_call) /= 0) ncalls = ncalls + 1
  thisnchanprof = nchanprof
  thisnprof = rth%nprofs_per_call

  DO i = 1, ncalls

    ! Calculate low/high profile and chanprof indices
    iprof1 = (i - 1) * rth%nprofs_per_call + 1
    IF (iprof1 + rth%nprofs_per_call > nprofiles) THEN
      thisnprof = nprofiles - iprof1 + 1
      thisnchanprof = thisnprof * nchannels
    ENDIF

    ! Populate the profile structure
    CALL rttov_copy_to_profiles(  &
        rth,                      &
        profiles(1:thisnprof),    &
        iprof1,                   &
        datetimes,                &
        angles, surfgeom,         &
        surftype, skin, s2m,      &
        p, t, gas_units,          &
        gas_id, gases,            &
        lmmr_cldaer,              &
        simplecloud,              &
        clwscheme, icecloud, zeeman)

    ! Set surface emissivities and BRDFs: if user passed any non-negative
    ! emissivities/BRDFs then use these and set calcemis/calcrefl FALSE
    calcemis(:) = .TRUE.
    calcrefl(:) = .TRUE.
    CALL rttov_init_emis_refl(emissivity, reflectance)
    ! Simple cloud scheme cloud top albedo cannot currently be set by wrapper
    ! reflectance%refl_cloud_top = 0._jprb

    DO iprof = 1, thisnprof
      k = iprof1 + iprof - 1
      lochan = (iprof - 1) * nchannels + 1
      hichan = iprof * nchannels
      WHERE (surfemisrefl(:,k,1) >= 0._jprb)
        emissivity(lochan:hichan)%emis_in = surfemisrefl(:,k,1)
        calcemis(lochan:hichan) = .FALSE.
      ENDWHERE
      WHERE (surfemisrefl(:,k,2) >= 0._jprb)
        reflectance(lochan:hichan)%refl_in = surfemisrefl(:,k,2)
        calcrefl(lochan:hichan) = .FALSE.
      ENDWHERE
      WHERE (surfemisrefl(:,k,3) >= 0._jprb)
        reflectance(lochan:hichan)%diffuse_refl_in = surfemisrefl(:,k,3)
      ENDWHERE
      WHERE (surfemisrefl(:,k,4) >= 0._jprb)
        emissivity(lochan:hichan)%specularity = surfemisrefl(:,k,4)
      ENDWHERE
    ENDDO

    ! Populate the radiance_k structure
    CALL rttov_init_rad(radiance_k)
    CALL rttov_copy_to_radiance_k( &
        radiance_k,                &
        iprof1,                    &
        thisnprof,                 &
        bt_k, rads_k)

    ! Initialise other K arrays
    CALL rttov_init_prof(profiles_k)
    CALL rttov_init_transmission(transmission_k)
    CALL rttov_init_emis_refl(emissivity_k, reflectance_k)

    ! Run RTTOV K
    IF (rth%nthreads <= 1) THEN
      CALL rttov_k( &
          err,                                            &
          chanprof(1:thisnchanprof),                      &
          rth%opts,                                       &
          profiles(1:thisnprof),                          &
          profiles_k(1:thisnchanprof),                    &
          rth%coefs,                                      &
          transmission,                                   &
          transmission_k,                                 &
          radiance,                                       &
          radiance_k,                                     &
          radiance2     = radiance2,                      &
          calcemis      = calcemis(1:thisnchanprof),      &
          emissivity    = emissivity(1:thisnchanprof),    &
          emissivity_k  = emissivity_k(1:thisnchanprof),  &
          calcrefl      = calcrefl(1:thisnchanprof),      &
          reflectance   = reflectance(1:thisnchanprof),   &
          reflectance_k = reflectance_k(1:thisnchanprof))
    ELSE
      CALL rttov_parallel_k( &
          err,                                            &
          chanprof(1:thisnchanprof),                      &
          rth%opts,                                       &
          profiles(1:thisnprof),                          &
          profiles_k(1:thisnchanprof),                    &
          rth%coefs,                                      &
          transmission,                                   &
          transmission_k,                                 &
          radiance,                                       &
          radiance_k,                                     &
          radiance2     = radiance2,                      &
          calcemis      = calcemis(1:thisnchanprof),      &
          emissivity    = emissivity(1:thisnchanprof),    &
          emissivity_k  = emissivity_k(1:thisnchanprof),  &
          calcrefl      = calcrefl(1:thisnchanprof),      &
          reflectance   = reflectance(1:thisnchanprof),   &
          reflectance_k = reflectance_k(1:thisnchanprof), &
          nthreads      = rth%nthreads)
    ENDIF
    THROWM(err.NE.0, 'Error running RTTOV')

    ! Store results
    DO iprof = 1, thisnprof
      k = iprof1 + iprof - 1
      lochan = (iprof - 1) * nchannels + 1
      hichan = iprof * nchannels
      WHERE (rth%coefs%coef%ss_val_chn(channel_list(:)) < 2)
        btrefl(:,k) = radiance%bt(lochan:hichan)
      ELSEWHERE
        btrefl(:,k) = radiance%refl(lochan:hichan)
      ENDWHERE
      rads(:,k) = radiance%total(lochan:hichan)
      surfemisrefl(:,k,1) = emissivity(lochan:hichan)%emis_out
      surfemisrefl(:,k,2) = reflectance(lochan:hichan)%refl_out
      surfemisrefl(:,k,3) = reflectance(lochan:hichan)%diffuse_refl_out
      surfemisrefl_k(:,k,1) = emissivity_k(lochan:hichan)%emis_in
      surfemisrefl_k(:,k,2) = reflectance_k(lochan:hichan)%refl_in
      surfemisrefl_k(:,k,3) = reflectance_k(lochan:hichan)%diffuse_refl_in
      surfemisrefl_k(:,k,4) = emissivity_k(lochan:hichan)%specularity
    ENDDO
    CALL rttov_copy_from_profiles_k( &
        rth,                           &
        profiles_k(1:thisnchanprof),   &
        iprof1,                        &
        skin_k, s2m_k,                 &
        t_k, gas_id, gases_k,          &
        simplecloud_k = simplecloud_k, &
        p_k = p_k)

    lochan = (iprof1 - 1) * nchannels + 1
    hichan = (iprof1 + thisnprof - 1) * nchannels

    IF (rth%store_trans) THEN
      rth%transmission%tau_total(lochan:hichan)             = transmission%tau_total(1:thisnchanprof)
      rth%transmission%tau_levels(:,lochan:hichan)          = transmission%tau_levels(:,1:thisnchanprof)
      rth%transmission%tausun_total_path2(lochan:hichan)    = transmission%tausun_total_path2(1:thisnchanprof)
      rth%transmission%tausun_levels_path2(:,lochan:hichan) = transmission%tausun_levels_path2(:,1:thisnchanprof)
      rth%transmission%tausun_total_path1(lochan:hichan)    = transmission%tausun_total_path1(1:thisnchanprof)
      rth%transmission%tausun_levels_path1(:,lochan:hichan) = transmission%tausun_levels_path1(:,1:thisnchanprof)
      rth%transmission%tau_total_cld(lochan:hichan)         = transmission%tau_total_cld(1:thisnchanprof)
      rth%transmission%tau_levels_cld(:,lochan:hichan)      = transmission%tau_levels_cld(:,1:thisnchanprof)
    ENDIF

    IF (rth%store_rad) THEN
      rth%radiance%clear(lochan:hichan)              = radiance%clear(1:thisnchanprof)
      rth%radiance%total(lochan:hichan)              = radiance%total(1:thisnchanprof)
      rth%radiance%bt_clear(lochan:hichan)           = radiance%bt_clear(1:thisnchanprof)
      rth%radiance%bt(lochan:hichan)                 = radiance%bt(1:thisnchanprof)
      rth%radiance%refl_clear(lochan:hichan)         = radiance%refl_clear(1:thisnchanprof)
      rth%radiance%refl(lochan:hichan)               = radiance%refl(1:thisnchanprof)
      rth%radiance%cloudy(lochan:hichan)             = radiance%cloudy(1:thisnchanprof)
      rth%radiance%overcast(:,lochan:hichan)         = radiance%overcast(:,1:thisnchanprof)
      rth%radiance%quality(lochan:hichan)            = radiance%quality(1:thisnchanprof)
      rth%radiance%plane_parallel                    = radiance%plane_parallel
      rth%radiance%geometric_height(:,lochan:hichan) = radiance%geometric_height(:,1:thisnchanprof)
    ENDIF

    IF (rth%store_rad2) THEN
      rth%radiance2%upclear(lochan:hichan)     = radiance2%upclear(1:thisnchanprof)
      rth%radiance2%dnclear(lochan:hichan)     = radiance2%dnclear(1:thisnchanprof)
      rth%radiance2%refldnclear(lochan:hichan) = radiance2%refldnclear(1:thisnchanprof)
      rth%radiance2%up(:,lochan:hichan)        = radiance2%up(:,1:thisnchanprof)
      rth%radiance2%down(:,lochan:hichan)      = radiance2%down(:,1:thisnchanprof)
      rth%radiance2%surf(:,lochan:hichan)      = radiance2%surf(:,1:thisnchanprof)
    ENDIF

  ENDDO ! loop over calls to RTTOV

  ! Deallocate RTTOV data
  asw = 0_jpim
  CALL rttov_alloc_k( &
              err,            &
              asw,            &
              nprof,          &
              nchanprof,      &
              nlevels,        &
              chanprof,       &
              rth%opts,       &
              profiles,       &
              profiles_k,     &
              rth%coefs,      &
              transmission,   &
              transmission_k, &
              radiance,       &
              radiance_k,     &
              radiance2,      &
              calcemis,       &
              emissivity,     &
              emissivity_k,   &
              calcrefl,       &
              reflectance,    &
              reflectance_k)
  THROWM(err.NE.0, 'Error deallocating RTTOV structures')

  CATCH
END SUBROUTINE rttov_call_k


!> Call rttov_k for visible/IR scattering simulations using explicit optical properties
!! @param[out]    err               return status
!! @param[in]     inst_id           ID for instrument to simulate
!! @param[in]     channel_list      list of channels to simulate for for each profile
!! @param[in]     datetimes         profile dates/times: year, month, day, hour, min, sec (6,nprofiles)
!! @param[in]     angles            profile angles: zenang, aziang, sunzenang, sunaziang (4,nprofiles)
!! @param[in]     surfgeom          profile surface geometry: lat, lon, elevation (3,nprofiles)
!! @param[in]     surftype          profile surface type: surftype, watertype (2,nprofiles)
!! @param[in]     skin              profile skin data: skin T, salinity, snow_frac, foam_frac, fastem_coefsx5 (9,nprofiles)
!! @param[in,out] skin_k            Jacobians for skin data (9,nchannels,nprofiles)
!! @param[in]     s2m               profile s2m data: 2m p, 2m t, 2m q, 10m wind u, v, wind-fetch (6,nprofiles)
!! @param[in,out] s2m_k             Jacobians for s2m data (6,nchannels,nprofiles)
!! @param[in]     clwscheme         profile clw scheme: clw_scheme, clwde_param (2,nprofiles)
!! @param[in]     icecloud          profile ice scheme: ice_scheme, icede_param (2,nprofiles)
!! @param[in]     p                 pressure profiles (nlevels,nprofiles)
!! @param[in,out] p_k               pressure Jacobians (nlevels,nchannels,nprofiles)
!! @param[in]     t                 temperature profiles (nlevels,nprofiles)
!! @param[in,out] t_k               temperature Jacobians (nlevels,nchannels,nprofiles)
!! @param[in]     gas_units         profile gas units
!! @param[in]     mmr_cldaer        profile mmr_cldaer flag
!! @param[in]     gas_id            list of gas IDs specified in gases argument (ngases)
!! @param[in]     gases             gas, aerosol and cloud profiles (nlevels,nprofiles,ngases)
!! @param[in,out] gases_k           Jacobians for contents of gases array (nlevels,nchannels,nprofiles,ngases)
!! @param[in]     aer_phangle       aerosol phase fn angles (aer_nphangle)
!! @param[in]     aer_asb           aerosol abs, sca, bpr parameters (nlayers,nchannels,nprofiles,3)
!! @param[in]     aer_legcoef       aerosol phase function Legendre coefficients (aer_nmom+1,nlayers,nchannels,nprofiles)
!! @param[in]     aer_pha           aerosol phase fns (aer_nphangle,nlayers,nchannels,nprofiles)
!! @param[in]     cld_phangle       cloud phase fn angles (cld_nphangle)
!! @param[in]     cld_asb           cloud abs, sca, bpr parameters (nlayers,nchannels,nprofiles,3)
!! @param[in]     cld_legcoef       cloud phase function Legendre coefficients (cld_nmom+1,nlayers,nchannels,nprofiles)
!! @param[in]     cld_pha           cloud phase fns (cld_nphangle,nlayers,nchannels,nprofiles)
!! @param[in]     surfemisrefl      input/output surface emissivities, reflectances, specularities (nchannels,nprofiles,4)
!! @param[in,out] surfemisrefl_k    emissivity/reflectance/specularity Jacobians (nchannels,nprofiles,4)
!! @param[in,out] btrefl            returned BTs/reflectances (nchannels,nprofiles)
!! @param[in,out] rads              returned radiances (nchannels,nprofiles)
!! @param[in]     bt_k              input BT perturbations (nchannels,nprofiles)
!! @param[in]     rads_k            input radiance perturbations (nchannels,nprofiles)
!! @param[in]     nchannels         number of channels (size of channel_list, not required in Python)
!! @param[in]     ngases            number of gases (size of gas_id, not required in Python)
!! @param[in]     nlevels           number of levels (size of p(:,1), not required in Python)
!! @param[in]     nprofiles         number of profiles (size of p(1,:), not required in Python)
!! @param[in]     aer_nphangle      number of aerosol phase angles (size of aer_phangle, not required in Python)
!! @param[in]     aer_nmom          number of aerosol Legendre moments excluding zeroth moment (size of
!!                                      aer_legcoef(:,1,1,1)-1, not required in Python)
!! @param[in]     cld_nphangle      number of cloud phase angles (size of cld_phangle, not required in Python)
!! @param[in]     cld_nmom          number of cloud Legendre moments excluding zeroth moment (size of
!!                                      cld_legcoef(:,1,1,1)-1, not required in Python)
SUBROUTINE rttov_visir_scatt_call_k( &
    err,               &
    inst_id,           &
    channel_list,      &
    datetimes,         &    ! profile dates/times                                     (6,nprofiles)
    angles,            &    ! satzen, satazi, sunzen, sunazi angles                   (4,nprofiles)
    surfgeom,          &    ! lat, lon, elevation                                     (3,nprofiles)
    surftype,          &    ! surftype, watertype                                     (2,nprofiles)
    skin,              &    ! skin T, salinity, snow_frac, foam_frac, fastem_coefsx5  (9,nprofiles)
    skin_k,            &    ! skin K                                                  (9,nchannels,nprofiles)
    s2m,               &    ! 2m p, 2m t, 2m q, 10m wind u, v, wind-fetch             (6,nprofiles)
    s2m_k,             &    ! 2m K                                                    (6,nchannels,nprofiles)
    clwscheme,         &    ! clw_scheme, clwde_param                                 (2,nprofiles)
    icecloud,          &    ! ice_scheme, icede_param                                 (2,nprofiles)
    p,                 &    ! pressure                                                (nlevels,nprofiles)
    p_k,               &    ! output pressure K                                       (nlevels,nchannels,nprofiles)
    t,                 &    ! temperature                                             (nlevels,nprofiles)
    t_k,               &    ! output temperature K                                    (nlevels,nchannels,nprofiles)
    gas_units,         &    ! units for gas profiles
    mmr_cldaer,        &    ! set mmr_cldaer flag
    gas_id,            &    ! gas ID list                                             (ngases)
    gases,             &    ! gas profiles                                            (nlevels,nprofiles,ngases)
    gases_k,           &    ! output gas profiles K                                   (nlevels,nchannels,nprofiles,ngases)
    aer_phangle,       &    ! aerosol phase fn angles                                 (aer_nphangle)
    aer_asb,           &    ! aerosol abs, sca, bpr parameters                        (nlayers,nchannels,nprofiles,3)
    aer_legcoef,       &    ! aerosol phase fn Legendre coefficients                  (aer_nmom+1,nlayers,nchannels,nprofiles)
    aer_pha,           &    ! aerosol phase fns                                       (aer_nphangle,nlayers,nchannels,nprofiles)
    cld_phangle,       &    ! cloud phase fn angles                                   (cld_nphangle)
    cld_asb,           &    ! cloud abs, sca, bpr parameters                          (nlayers,nchannels,nprofiles,3)
    cld_legcoef,       &    ! cloud phase fn Legendre coefficients                    (cld_nmom+1,nlayers,nchannels,nprofiles)
    cld_pha,           &    ! cloud phase fns                                         (cld_nphangle,nlayers,nchannels,nprofiles)
    surfemisrefl,      &    ! surface emissivities, reflectances, specularities       (nchannels,nprofiles,4)
    surfemisrefl_k,    &    ! output surface emis, refl, spec K                       (nchannels,nprofiles,4)
    btrefl,            &    ! output BTs/refls (for thermal/solar chans)              (nchannels,nprofiles)
    rads,              &    ! output radiances                                        (nchannels,nprofiles)
    bt_k,              &    ! input BT perturbations (thermal chans only)             (nchannels,nprofiles)
    rads_k,            &    ! input radiance perturbations                            (nchannels,nprofiles)
    nchannels, ngases, nlevels, nprofiles, aer_nphangle, aer_nmom, cld_nphangle, cld_nmom)
!f2py threadsafe
!f2py intent(hide):: nchannels=len(channel_list)
!f2py intent(hide):: ngases=len(gas_id)
!f2py intent(hide):: nlevels=shape(p, 1)
!f2py intent(hide):: nprofiles=shape(p, 2)
!f2py intent(hide):: aer_nphangle=len(aer_phangle)
!f2py intent(hide):: aer_nmom=shape(aer_legcoef, 1)-1
!f2py intent(hide):: cld_nphangle=len(cld_phangle)
!f2py intent(hide):: cld_nmom=shape(cld_legcoef, 1)-1

!
! Prepares input profiles, calls RTTOV K and returns radiances and profiles_k
!

#include "throw.h"

  USE rttov_wrapper_handle

  USE rttov_wrapper_transfer

  USE parkind1, ONLY : jpim, jprb, jplm

  USE rttov_types, ONLY : &
      rttov_chanprof,     &
      rttov_emissivity,   &
      rttov_reflectance,  &
      rttov_profile,      &
      rttov_radiance,     &
      rttov_transmission, &
      rttov_opt_param

  USE rttov_const, ONLY : sensor_id_ir, sensor_id_hi

  IMPLICIT NONE

  INTEGER(jpim),    INTENT(OUT)   :: err

  INTEGER(jpim),    INTENT(IN)    :: nchannels, ngases, nlevels, nprofiles

  INTEGER(jpim),    INTENT(IN)    :: inst_id
  INTEGER(jpim),    INTENT(IN)    :: channel_list(nchannels)
  INTEGER(jpim),    INTENT(IN)    :: datetimes(6,nprofiles)
  REAL(jprb),       INTENT(IN)    :: angles(4,nprofiles)
  REAL(jprb),       INTENT(IN)    :: surfgeom(3,nprofiles)
  INTEGER(jpim),    INTENT(IN)    :: surftype(2,nprofiles)
  REAL(jprb),       INTENT(IN)    :: skin(9,nprofiles)
  REAL(jprb),       INTENT(INOUT) :: skin_k(9,nchannels,nprofiles)
  REAL(jprb),       INTENT(IN)    :: s2m(6,nprofiles)
  REAL(jprb),       INTENT(INOUT) :: s2m_k(6,nchannels,nprofiles)
  INTEGER(jpim),    INTENT(IN)    :: clwscheme(2,nprofiles)
  INTEGER(jpim),    INTENT(IN)    :: icecloud(2,nprofiles)
  REAL(jprb),       INTENT(IN)    :: p(nlevels,nprofiles)
  REAL(jprb),       INTENT(INOUT) :: p_k(nlevels,nchannels,nprofiles)
  REAL(jprb),       INTENT(IN)    :: t(nlevels,nprofiles)
  REAL(jprb),       INTENT(INOUT) :: t_k(nlevels,nchannels,nprofiles)
  INTEGER(jpim),    INTENT(IN)    :: gas_units
  INTEGER(jpim),    INTENT(IN)    :: mmr_cldaer
  INTEGER(jpim),    INTENT(IN)    :: gas_id(ngases)
  REAL(jprb),       INTENT(IN)    :: gases(nlevels,nprofiles,ngases)
  REAL(jprb),       INTENT(INOUT) :: gases_k(nlevels,nchannels,nprofiles,ngases)
  INTEGER(jpim),    INTENT(IN)    :: aer_nphangle, aer_nmom, cld_nphangle, cld_nmom
  REAL(jprb),       INTENT(IN)    :: aer_phangle(aer_nphangle)
  REAL(jprb),       INTENT(IN)    :: aer_asb(nlevels-1,nchannels,nprofiles,3)
  REAL(jprb),       INTENT(IN)    :: aer_legcoef(aer_nmom+1,nlevels-1,nchannels,nprofiles)
  REAL(jprb),       INTENT(IN)    :: aer_pha(aer_nphangle,nlevels-1,nchannels,nprofiles)
  REAL(jprb),       INTENT(IN)    :: cld_phangle(cld_nphangle)
  REAL(jprb),       INTENT(IN)    :: cld_asb(nlevels-1,nchannels,nprofiles,3)
  REAL(jprb),       INTENT(IN)    :: cld_legcoef(cld_nmom+1,nlevels-1,nchannels,nprofiles)
  REAL(jprb),       INTENT(IN)    :: cld_pha(cld_nphangle,nlevels-1,nchannels,nprofiles)

  REAL(jprb),       INTENT(INOUT) :: surfemisrefl(nchannels,nprofiles,4)
  REAL(jprb),       INTENT(INOUT) :: surfemisrefl_k(nchannels,nprofiles,4)

  REAL(jprb),       INTENT(INOUT) :: btrefl(nchannels,nprofiles)
  REAL(jprb),       INTENT(INOUT) :: rads(nchannels,nprofiles)
  REAL(jprb),       INTENT(IN)    :: bt_k(nchannels,nprofiles)
  REAL(jprb),       INTENT(IN)    :: rads_k(nchannels,nprofiles)


  INTEGER(jpim)      :: asw
  INTEGER(jpim)      :: i, j, k, iprof, iprof1, lochan, hichan
  INTEGER(jpim)      :: nchanprof, thisnchanprof, nprof, thisnprof, ncalls
  LOGICAL(jplm)      :: lmmr_cldaer
  CHARACTER(LEN=128) :: msg

  TYPE(rttov_chanprof),    POINTER :: chanprof(:)
  LOGICAL(jplm),           POINTER :: calcemis(:)
  TYPE(rttov_emissivity),  POINTER :: emissivity(:)
  TYPE(rttov_emissivity),  POINTER :: emissivity_k(:)
  LOGICAL(jplm),           POINTER :: calcrefl(:)
  TYPE(rttov_reflectance), POINTER :: reflectance(:)
  TYPE(rttov_reflectance), POINTER :: reflectance_k(:)
  TYPE(rttov_profile),     POINTER :: profiles(:)
  TYPE(rttov_profile),     POINTER :: profiles_k(:)
  TYPE(rttov_transmission)         :: transmission
  TYPE(rttov_transmission)         :: transmission_k
  TYPE(rttov_radiance)             :: radiance
  TYPE(rttov_radiance)             :: radiance_k
  TYPE(rttov_opt_param)            :: aer_opt_param
  TYPE(rttov_opt_param)            :: cld_opt_param

  TYPE(rttovwrapperhandle_type), POINTER :: rth

#include "rttov_alloc_k.interface"
#include "rttov_init_emis_refl.interface"
#include "rttov_init_rad.interface"
#include "rttov_init_transmission.interface"
#include "rttov_init_prof.interface"
#include "rttov_init_opt_param.interface"
#include "rttov_k.interface"
#include "rttov_parallel_k.interface"
#include "rttov_errorreport.interface"

!------------------------------------------------------------------------------

  TRY

  NULLIFY(chanprof,      &
          calcemis,      &
          emissivity,    &
          emissivity_k,  &
          calcrefl,      &
          reflectance,   &
          reflectance_k, &
          profiles,      &
          profiles_k)

  CALL rttov_wrapper_handle_get_inst(err, inst_id, rth)
  THROW(err.NE.0)

  IF (.NOT. ((rth%opts%rt_ir%addclouds .AND. rth%opts%rt_ir%user_cld_opt_param) .OR. &
             (rth%opts%rt_ir%addaerosl .AND. rth%opts%rt_ir%user_aer_opt_param))) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'This subroutine is for user specified cloud and/or aerosol properties only')
  ENDIF

  IF (.NOT. (rth%coefs%coef%id_sensor == sensor_id_ir .OR. &
             rth%coefs%coef%id_sensor == sensor_id_hi)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'This subroutine is for visible/IR sensors only')
  ENDIF

  ! Deallocate (if necessary) and re-allocate the output structures
  CALL rttov_wrapper_handle_alloc(err, rth, nprofiles, nchannels, nlevels, .FALSE._jplm, 0_jpim)
  THROW(err.NE.0)
  CALL rttov_wrapper_handle_alloc(err, rth, nprofiles, nchannels, nlevels, .FALSE._jplm, 1_jpim)
  THROW(err.NE.0)

  ! nchannels     - input argument, number of channels simulated for every profile
  ! nprofiles     - input argument, total number of profiles passed to this subroutine
  ! ncalls        - number of calls to RTTOV
  ! nprof         - number of profiles per call to RTTOV
  ! thisnprof     - number of profiles in current call (may differ to nprof on last call)
  ! iprof1        - index of first profile in current batch
  ! nchanprof     - size of chanprof per call to RTTOV
  ! thisnchanprof - size of chanprof for current call (may differ to nchanprof on last call)

  nprof = MIN(nprofiles, rth%nprofs_per_call)
  nchanprof = nprof * nchannels
  lmmr_cldaer = (mmr_cldaer .NE. 0)

  ! Allocate RTTOV data
  asw = 1_jpim
  CALL rttov_alloc_k( &
              err,                           &
              asw,                           &
              nprof,                         &
              nchanprof,                     &
              nlevels,                       &
              chanprof,                      &
              rth%opts,                      &
              profiles,                      &
              profiles_k,                    &
              rth%coefs,                     &
              transmission,                  &
              transmission_k,                &
              radiance,                      &
              radiance_k,                    &
              calcemis      = calcemis,      &
              emissivity    = emissivity,    &
              emissivity_k  = emissivity_k,  &
              calcrefl      = calcrefl,      &
              reflectance   = reflectance,   &
              reflectance_k = reflectance_k, &
              aer_maxnmom   = aer_nmom,      &
              aer_nphangle  = aer_nphangle,  &
              aer_opt_param = aer_opt_param, &
              cld_maxnmom   = cld_nmom,      &
              cld_nphangle  = cld_nphangle,  &
              cld_opt_param = cld_opt_param, &
              init          = .TRUE._jplm)
  THROWM(err.NE.0, 'Error allocating RTTOV structures')

  ! Populate input arrays which don't change
  DO j = 1, nprof
    lochan = (j - 1) * nchannels + 1
    hichan = j * nchannels
    chanprof(lochan:hichan)%prof = j
    chanprof(lochan:hichan)%chan = channel_list(:)
  ENDDO

  IF (rth%verbose) THEN
    WRITE(msg, '(a,i4,a,i8)') &
      'Running RTTOV K using nthreads = ', rth%nthreads, &
      ' and nprofs_per_call = ', rth%nprofs_per_call
    INFO(msg)
  ENDIF

  ! Main loop over batches of profiles

  ncalls = nprofiles / rth%nprofs_per_call
  IF (MOD(nprofiles, rth%nprofs_per_call) /= 0) ncalls = ncalls + 1
  thisnchanprof = nchanprof
  thisnprof = rth%nprofs_per_call

  DO i = 1, ncalls

    ! Calculate low/high profile and chanprof indices
    iprof1 = (i - 1) * rth%nprofs_per_call + 1
    IF (iprof1 + rth%nprofs_per_call > nprofiles) THEN
      thisnprof = nprofiles - iprof1 + 1
      thisnchanprof = thisnprof * nchannels
    ENDIF

    ! Populate the profile structure
    CALL rttov_copy_to_profiles(  &
        rth,                      &
        profiles(1:thisnprof),    &
        iprof1,                   &
        datetimes,                &
        angles, surfgeom,         &
        surftype, skin, s2m,      &
        p, t, gas_units,          &
        gas_id, gases,            &
        lmmr_cldaer,              &
        clwscheme = clwscheme,    &
        icecloud = icecloud)

    CALL rttov_copy_to_opt_param( &
        rth,                                  &
        aer_opt_param, cld_opt_param,         &
        iprof1, nprofiles,                    &
        thisnprof, nchannels, nlevels-1_jpim, &
        aer_nmom, aer_nphangle,               &
        aer_phangle,                          &
        aer_asb, aer_legcoef, aer_pha,        &
        cld_nmom, cld_nphangle,               &
        cld_phangle,                          &
        cld_asb, cld_legcoef, cld_pha)

    IF (i == 1 .AND. rth%opts%rt_ir%addsolar) THEN
      IF (rth%opts%rt_ir%addaerosl) THEN
        CALL rttov_init_opt_param(err, rth%opts, aer_opt_param)
        THROWM(err.NE.0, 'Error in rttov_init_opt_param for aerosols')
      ENDIF
      IF (rth%opts%rt_ir%addclouds) THEN
        CALL rttov_init_opt_param(err, rth%opts, cld_opt_param)
        THROWM(err.NE.0, 'Error in rttov_init_opt_param for clouds')
      ENDIF
    ENDIF

    ! Set surface emissivities and BRDFs: if user passed any non-negative
    ! emissivities/BRDFs then use these and set calcemis/calcrefl FALSE
    calcemis(:) = .TRUE.
    calcrefl(:) = .TRUE.
    CALL rttov_init_emis_refl(emissivity, reflectance)

    DO iprof = 1, thisnprof
      k = iprof1 + iprof - 1
      lochan = (iprof - 1) * nchannels + 1
      hichan = iprof * nchannels
      WHERE (surfemisrefl(:,k,1) >= 0._jprb)
        emissivity(lochan:hichan)%emis_in = surfemisrefl(:,k,1)
        calcemis(lochan:hichan) = .FALSE.
      ENDWHERE
      WHERE (surfemisrefl(:,k,2) >= 0._jprb)
        reflectance(lochan:hichan)%refl_in = surfemisrefl(:,k,2)
        calcrefl(lochan:hichan) = .FALSE.
      ENDWHERE
      WHERE (surfemisrefl(:,k,3) >= 0._jprb)
        reflectance(lochan:hichan)%diffuse_refl_in = surfemisrefl(:,k,3)
      ENDWHERE
      WHERE (surfemisrefl(:,k,4) >= 0._jprb)
        emissivity(lochan:hichan)%specularity = surfemisrefl(:,k,4)
      ENDWHERE
    ENDDO

    ! Populate the radiance_k structure
    CALL rttov_init_rad(radiance_k)
    CALL rttov_copy_to_radiance_k( &
        radiance_k,                &
        iprof1,                    &
        thisnprof,                 &
        bt_k, rads_k)

    ! Initialise other K arrays
    CALL rttov_init_prof(profiles_k)
    CALL rttov_init_transmission(transmission_k)
    CALL rttov_init_emis_refl(emissivity_k, reflectance_k)

    ! Run RTTOV K
    IF (rth%nthreads <= 1) THEN
      CALL rttov_k( &
          err,                                            &
          chanprof(1:thisnchanprof),                      &
          rth%opts,                                       &
          profiles(1:thisnprof),                          &
          profiles_k(1:thisnchanprof),                    &
          rth%coefs,                                      &
          transmission,                                   &
          transmission_k,                                 &
          radiance,                                       &
          radiance_k,                                     &
          calcemis      = calcemis(1:thisnchanprof),      &
          emissivity    = emissivity(1:thisnchanprof),    &
          emissivity_k  = emissivity_k(1:thisnchanprof),  &
          calcrefl      = calcrefl(1:thisnchanprof),      &
          reflectance   = reflectance(1:thisnchanprof),   &
          reflectance_k = reflectance_k(1:thisnchanprof), &
          aer_opt_param = aer_opt_param,                  &
          cld_opt_param = cld_opt_param)

    ELSE
      CALL rttov_parallel_k( &
          err,                                            &
          chanprof(1:thisnchanprof),                      &
          rth%opts,                                       &
          profiles(1:thisnprof),                          &
          profiles_k(1:thisnchanprof),                    &
          rth%coefs,                                      &
          transmission,                                   &
          transmission_k,                                 &
          radiance,                                       &
          radiance_k,                                     &
          calcemis      = calcemis(1:thisnchanprof),      &
          emissivity    = emissivity(1:thisnchanprof),    &
          emissivity_k  = emissivity_k(1:thisnchanprof),  &
          calcrefl      = calcrefl(1:thisnchanprof),      &
          reflectance   = reflectance(1:thisnchanprof),   &
          reflectance_k = reflectance_k(1:thisnchanprof), &
          aer_opt_param = aer_opt_param,                  &
          cld_opt_param = cld_opt_param,                  &
          nthreads      = rth%nthreads)
    ENDIF
    THROWM(err.NE.0, 'Error running RTTOV')

    ! Store results
    DO iprof = 1, thisnprof
      k = iprof1 + iprof - 1
      lochan = (iprof - 1) * nchannels + 1
      hichan = iprof * nchannels
      WHERE (rth%coefs%coef%ss_val_chn(channel_list(:)) < 2)
        btrefl(:,k) = radiance%bt(lochan:hichan)
      ELSEWHERE
        btrefl(:,k) = radiance%refl(lochan:hichan)
      ENDWHERE
      rads(:,k) = radiance%total(lochan:hichan)
      surfemisrefl(:,k,1) = emissivity(lochan:hichan)%emis_out
      surfemisrefl(:,k,2) = reflectance(lochan:hichan)%refl_out
      surfemisrefl(:,k,3) = reflectance(lochan:hichan)%diffuse_refl_out
      surfemisrefl_k(:,k,1) = emissivity_k(lochan:hichan)%emis_in
      surfemisrefl_k(:,k,2) = reflectance_k(lochan:hichan)%refl_in
      surfemisrefl_k(:,k,3) = reflectance_k(lochan:hichan)%diffuse_refl_in
      surfemisrefl_k(:,k,4) = emissivity_k(lochan:hichan)%specularity
    ENDDO
    CALL rttov_copy_from_profiles_k( &
        rth,                         &
        profiles_k(1:thisnchanprof), &
        iprof1,                      &
        skin_k, s2m_k,               &
        t_k, gas_id, gases_k,        &
        p_k = p_k)

    lochan = (iprof1 - 1) * nchannels + 1
    hichan = (iprof1 + thisnprof - 1) * nchannels

    IF (rth%store_trans) THEN
      rth%transmission%tau_total(lochan:hichan)             = transmission%tau_total(1:thisnchanprof)
      rth%transmission%tau_levels(:,lochan:hichan)          = transmission%tau_levels(:,1:thisnchanprof)
      rth%transmission%tausun_total_path2(lochan:hichan)    = transmission%tausun_total_path2(1:thisnchanprof)
      rth%transmission%tausun_levels_path2(:,lochan:hichan) = transmission%tausun_levels_path2(:,1:thisnchanprof)
      rth%transmission%tausun_total_path1(lochan:hichan)    = transmission%tausun_total_path1(1:thisnchanprof)
      rth%transmission%tausun_levels_path1(:,lochan:hichan) = transmission%tausun_levels_path1(:,1:thisnchanprof)
      rth%transmission%tau_total_cld(lochan:hichan)         = transmission%tau_total_cld(1:thisnchanprof)
      rth%transmission%tau_levels_cld(:,lochan:hichan)      = transmission%tau_levels_cld(:,1:thisnchanprof)
    ENDIF

    IF (rth%store_rad) THEN
      rth%radiance%clear(lochan:hichan)              = radiance%clear(1:thisnchanprof)
      rth%radiance%total(lochan:hichan)              = radiance%total(1:thisnchanprof)
      rth%radiance%bt_clear(lochan:hichan)           = radiance%bt_clear(1:thisnchanprof)
      rth%radiance%bt(lochan:hichan)                 = radiance%bt(1:thisnchanprof)
      rth%radiance%refl_clear(lochan:hichan)         = radiance%refl_clear(1:thisnchanprof)
      rth%radiance%refl(lochan:hichan)               = radiance%refl(1:thisnchanprof)
      rth%radiance%cloudy(lochan:hichan)             = radiance%cloudy(1:thisnchanprof)
      rth%radiance%overcast(:,lochan:hichan)         = radiance%overcast(:,1:thisnchanprof)
      rth%radiance%quality(lochan:hichan)            = radiance%quality(1:thisnchanprof)
      rth%radiance%plane_parallel                    = radiance%plane_parallel
      rth%radiance%geometric_height(:,lochan:hichan) = radiance%geometric_height(:,1:thisnchanprof)
    ENDIF
  ENDDO ! loop over calls to RTTOV

  ! Deallocate RTTOV data
  asw = 0_jpim
  CALL rttov_alloc_k( &
              err,                           &
              asw,                           &
              nprof,                         &
              nchanprof,                     &
              nlevels,                       &
              chanprof,                      &
              rth%opts,                      &
              profiles,                      &
              profiles_k,                    &
              rth%coefs,                     &
              transmission,                  &
              transmission_k,                &
              radiance,                      &
              radiance_k,                    &
              calcemis      = calcemis,      &
              emissivity    = emissivity,    &
              emissivity_k  = emissivity_k,  &
              calcrefl      = calcrefl,      &
              reflectance   = reflectance,   &
              reflectance_k = reflectance_k, &
              aer_maxnmom   = aer_nmom,      &
              aer_nphangle  = aer_nphangle,  &
              aer_opt_param = aer_opt_param, &
              cld_maxnmom   = cld_nmom,      &
              cld_nphangle  = cld_nphangle,  &
              cld_opt_param = cld_opt_param)
  THROWM(err.NE.0, 'Error deallocating RTTOV structures')

  CATCH
END SUBROUTINE rttov_visir_scatt_call_k


!> Call RTTOV-SCATT K model
!! @param[out]    err               return status
!! @param[in]     inst_id           ID for instrument to simulate
!! @param[in]     channel_list      list of channels to simulate for for each profile
!! @param[in]     datetimes         profile dates/times: year, month, day, hour, min, sec (6,nprofiles)
!! @param[in]     angles            profile angles: zenang, aziang (2,nprofiles)
!! @param[in]     surfgeom          profile surface geometry: lat, lon, elevation (3,nprofiles)
!! @param[in]     surftype          profile surface type (nprofiles)
!! @param[in]     skin              profile skin data: skin T, salinity, foam_frac, fastem_coefsx5 (8,nprofiles)
!! @param[in,out] skin_k            Jacobians for skin data (8,nchannels,nprofiles)
!! @param[in]     s2m               profile s2m data: 2m p, 2m t, 2m q, 10m wind u, v (5,nprofiles)
!! @param[in,out] s2m_k             Jacobians for s2m data (5,nchannels,nprofiles)
!! @param[in]     zeeman            profile B-field data: Be, cosbk (2,nprofiles)
!! @param[in]     p                 pressure profiles (nlevels,nprofiles)
!! @param[in,out] p_k               pressure Jacobians (nlevels,nchannels,nprofiles)
!! @param[in]     t                 temperature profiles (nlevels,nprofiles)
!! @param[in,out] t_k               temperature Jacobians (nlevels,nchannels,nprofiles)
!! @param[in]     gas_units         profile gas units
!! @param[in]     gas_id            list of gas IDs specified in gases argument (ngases)
!! @param[in]     gases             gas, cloud and hydrometeor profiles (nlevels,nprofiles,ngases)
!! @param[in,out] gases_k           Jacobians for contents of gases array (nlevels,nchannels,nprofiles,ngases)
!! @param[in]     ph                pressure half-level profiles (nlevels+1,nprofiles)
!! @param[in,out] ph_k              pressure half-level Jacobians (nlevels+1,nchannels,nprofiles)
!! @param[in]     cfrac             user cloud fraction (nprofiles)
!! @param[in,out] cfrac_k           user cloud fraction Jacobians (nchannels,nprofiles)
!! @param[in]     multi_hydro_frac  false => single hydro_frac profile, true => one hydro_frac profile per hydrometeor
!! @param[in]     calc_zef          enable/disable the radar reflectivity calculations
!! @param[in]     surfemis          input/output surface emissivities (nchannels,nprofiles)
!! @param[in,out] surfemis_k        emissivity Jacobians (nchannels,nprofiles)
!! @param[in,out] bt                returned BTs (nchannels,nprofiles)
!! @param[in]     bt_k              input BT perturbations (nchannels,nprofiles)
!! @param[in]     zef_k             input radar reflectivity perturbations (nlevels,nchannels,nprofiles)
!! @param[in]     nchannels         number of channels (size of channel_list, not required in Python)
!! @param[in]     ngases            number of gases (size of gas_id, not required in Python)
!! @param[in]     nlevels           number of levels (size of p(:,1), not required in Python)
!! @param[in]     nprofiles         number of profiles (size of p(1,:), not required in Python)
SUBROUTINE rttov_scatt_call_k( &
    err,               &
    inst_id,           &
    channel_list,      &
    datetimes,         &    ! profile dates/times                                     (6,nprofiles)
    angles,            &    ! satzen, satazi angles                                   (2,nprofiles)
    surfgeom,          &    ! lat, lon, elevation                                     (3,nprofiles)
    surftype,          &    ! surftype                                                (nprofiles)
    skin,              &    ! skin T, salinity, foam_frac, fastem_coefsx5             (8,nprofiles)
    skin_k,            &    ! skin K                                                  (8,nchannels,nprofiles)
    s2m,               &    ! 2m p, 2m t, 2m q, 10m wind u, v                         (5,nprofiles)
    s2m_k,             &    ! 2m K                                                    (5,nchannels,nprofiles)
    zeeman,            &    ! Be, cosbk                                               (2,nprofiles)
    p,                 &    ! pressure                                                (nlevels,nprofiles)
    p_k,               &    ! output pressure K                                       (nlevels,nchannels,nprofiles)
    t,                 &    ! temperature                                             (nlevels,nprofiles)
    t_k,               &    ! output temperature K                                    (nlevels,nchannels,nprofiles)
    gas_units,         &    ! units for gas profiles
    gas_id,            &    ! gas ID list                                             (ngases)
    gases,             &    ! gas profiles                                            (nlevels,nprofiles,ngases)
    gases_k,           &    ! output gas profiles K                                   (nlevels,nchannels,nprofiles,ngases)
    ph,                &    ! pressure half-levels                                    (nlevels+1,nprofiles)
    ph_k,              &    ! pressure half-levels K                                  (nlevels+1,nchannels,nprofiles)
    cfrac,             &    ! cfrac                                                   (nprofiles)
    cfrac_k,           &    ! cfrac K                                                 (nchannels,nprofiles)
    multi_hydro_frac,  &    ! false => single hydro_frac profile, true => one hydro_frac profile per hydrometeor
    calc_zef,          &    ! enable/disable the radar reflectivity calculations
    surfemis,          &    ! surface emissivities                                    (nchannels,nprofiles)
    surfemis_k,        &    ! output surface emissivities K                           (nchannels,nprofiles)
    bt,                &    ! output BTs                                              (nchannels,nprofiles)
    bt_k,              &    ! input BT perturbations                                  (nchannels,nprofiles)
    zef_k,             &    ! input Z reflectivity perturbations                      (nlevels,nchannels,nprofiles)
    nchannels, ngases, nlevels, nprofiles)
!f2py threadsafe
!f2py intent(hide):: nchannels=len(channel_list)
!f2py intent(hide):: ngases=len(gas_id)
!f2py intent(hide):: nlevels=shape(p, 1)
!f2py intent(hide):: nprofiles=shape(p, 2)

!
! Prepares input profiles, calls RTTOV-SCATT K and returns radiances, profiles_k and cld_profiles_k
!

#include "throw.h"

  USE rttov_wrapper_handle

  USE rttov_wrapper_transfer

  USE parkind1, ONLY : jpim, jprb, jplm

  USE rttov_types, ONLY :  &
      rttov_chanprof,      &
      rttov_emissivity,    &
      rttov_profile,       &
      rttov_profile_cloud, &
      rttov_radiance

  USE rttov_const, ONLY : sensor_id_mw, sensor_id_po

  IMPLICIT NONE

  INTEGER(jpim),    INTENT(OUT)   :: err

  INTEGER(jpim),    INTENT(IN)    :: nchannels, ngases, nlevels, nprofiles

  INTEGER(jpim),    INTENT(IN)    :: inst_id
  INTEGER(jpim),    INTENT(IN)    :: channel_list(nchannels)
  INTEGER(jpim),    INTENT(IN)    :: datetimes(6,nprofiles)
  REAL(jprb),       INTENT(IN)    :: angles(2,nprofiles)
  REAL(jprb),       INTENT(IN)    :: surfgeom(3,nprofiles)
  INTEGER(jpim),    INTENT(IN)    :: surftype(nprofiles)
  REAL(jprb),       INTENT(IN)    :: skin(8,nprofiles)
  REAL(jprb),       INTENT(INOUT) :: skin_k(8,nchannels,nprofiles)
  REAL(jprb),       INTENT(IN)    :: s2m(5,nprofiles)
  REAL(jprb),       INTENT(INOUT) :: s2m_k(5,nchannels,nprofiles)
  REAL(jprb),       INTENT(IN)    :: zeeman(2,nprofiles)
  REAL(jprb),       INTENT(IN)    :: p(nlevels,nprofiles)
  REAL(jprb),       INTENT(INOUT) :: p_k(nlevels,nchannels,nprofiles)
  REAL(jprb),       INTENT(IN)    :: t(nlevels,nprofiles)
  REAL(jprb),       INTENT(INOUT) :: t_k(nlevels,nchannels,nprofiles)
  INTEGER(jpim),    INTENT(IN)    :: gas_units
  INTEGER(jpim),    INTENT(IN)    :: gas_id(ngases)
  REAL(jprb),       INTENT(IN)    :: gases(nlevels,nprofiles,ngases)
  REAL(jprb),       INTENT(INOUT) :: gases_k(nlevels,nchannels,nprofiles,ngases)
  REAL(jprb),       INTENT(IN)    :: ph(nlevels+1,nprofiles)
  REAL(jprb),       INTENT(INOUT) :: ph_k(nlevels+1,nchannels,nprofiles)
  REAL(jprb),       INTENT(IN)    :: cfrac(nprofiles)
  REAL(jprb),       INTENT(INOUT) :: cfrac_k(nchannels,nprofiles)
  LOGICAL(jpim),    INTENT(IN)    :: multi_hydro_frac
  LOGICAL(jpim),    INTENT(IN)    :: calc_zef

  REAL(jprb),       INTENT(INOUT) :: surfemis(nchannels,nprofiles)
  REAL(jprb),       INTENT(INOUT) :: surfemis_k(nchannels,nprofiles)

  REAL(jprb),       INTENT(INOUT) :: bt(nchannels,nprofiles)
  REAL(jprb),       INTENT(IN)    :: bt_k(nchannels,nprofiles)
  REAL(jprb),       INTENT(IN)    :: zef_k(nlevels,nchannels,nprofiles)

  INTEGER(jpim)      :: asw
  INTEGER(jpim)      :: i, k, iprof, iprof1, lochan, hichan
  INTEGER(jpim)      :: nchanprof, thisnchanprof, nprof, thisnprof, ncalls
  INTEGER(jpim)      :: nhydro_frac
  CHARACTER(LEN=128) :: msg

  TYPE(rttov_chanprof),      POINTER :: chanprof(:)
  INTEGER(KIND=jpim),        POINTER :: frequencies(:)
  LOGICAL(KIND=jplm),        POINTER :: use_chan(:,:)
  LOGICAL(jplm),             POINTER :: calcemis(:)
  TYPE(rttov_emissivity),    POINTER :: emissivity(:)
  TYPE(rttov_emissivity),    POINTER :: emissivity_k(:)
  TYPE(rttov_profile),       POINTER :: profiles(:)
  TYPE(rttov_profile),       POINTER :: profiles_k(:)
  TYPE(rttov_profile_cloud), POINTER :: cld_profiles(:)
  TYPE(rttov_profile_cloud), POINTER :: cld_profiles_k(:)
  TYPE(rttov_radiance)               :: radiance
  TYPE(rttov_radiance)               :: radiance_k
  TYPE(rttov_reflectivity),  POINTER :: reflectivity
  TYPE(rttov_reflectivity),  POINTER :: reflectivity_k

  TYPE(rttovwrapperhandle_type), POINTER :: rth

#include "rttov_alloc_k.interface"
#include "rttov_scatt_setupindex.interface"
#include "rttov_init_emis_refl.interface"
#include "rttov_init_rad.interface"
#include "rttov_init_prof.interface"
#include "rttov_init_scatt_prof.interface"
#include "rttov_scatt_ad.interface"
#include "rttov_parallel_scatt_ad.interface"
#include "rttov_errorreport.interface"

!------------------------------------------------------------------------------

  TRY

  NULLIFY(chanprof,       &
          frequencies,    &
          use_chan,       &
          calcemis,       &
          emissivity,     &
          emissivity_k,   &
          profiles,       &
          profiles_k,     &
          cld_profiles,   &
          cld_profiles_k, &
          reflectivity,   &
          reflectivity_k)

  ! If reflectivity is nullified, it tells the routines below that radar
  ! calculations are turned off
  IF (calc_zef) ALLOCATE(reflectivity, reflectivity_k)

  CALL rttov_wrapper_handle_get_inst(err, inst_id, rth)
  THROW(err.NE.0)

  IF (.NOT. (rth%coefs%coef%id_sensor == sensor_id_mw .OR. &
             rth%coefs%coef%id_sensor == sensor_id_po)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'RTTOV-SCATT is for MW sensors only')
  ENDIF

  ! Deallocate (if necessary) and re-allocate the output structures
  CALL rttov_wrapper_handle_alloc(err, rth, nprofiles, nchannels, nlevels, calc_zef, 0_jpim)
  THROW(err.NE.0)
  CALL rttov_wrapper_handle_alloc(err, rth, nprofiles, nchannels, nlevels, calc_zef, 1_jpim)
  THROW(err.NE.0)

  ! nchannels     - input argument, number of channels simulated for every profile
  ! nprofiles     - input argument, total number of profiles passed to this subroutine
  ! ncalls        - number of calls to RTTOV
  ! nprof         - number of profiles per call to RTTOV
  ! thisnprof     - number of profiles in current call (may differ to nprof on last call)
  ! iprof1        - index of first profile in current batch
  ! nchanprof     - size of chanprof per call to RTTOV
  ! thisnchanprof - size of chanprof for current call (may differ to nchanprof on last call)

  nprof = MIN(nprofiles, rth%nprofs_per_call)
  nchanprof = nprof * nchannels

  nhydro_frac = 1
  IF (multi_hydro_frac) nhydro_frac = rth%coefs_scatt%nhydro

  ! Allocate RTTOV data
  asw = 1_jpim
  CALL rttov_alloc_k( &
              err,                              &
              asw,                              &
              nprof,                            &
              nchanprof,                        &
              nlevels,                          &
              chanprof,                         &
              rth%opts,                         &
              profiles,                         &
              profiles_k,                       &
              rth%coefs,                        &
              radiance       = radiance,        &
              radiance_k     = radiance_k,      &
              calcemis       = calcemis,        &
              emissivity     = emissivity,      &
              emissivity_k   = emissivity_k,    &
              frequencies    = frequencies,     &
              coef_scatt     = rth%coefs_scatt, &
              nhydro_frac    = nhydro_frac,     &
              cld_profiles   = cld_profiles,    &
              cld_profiles_k = cld_profiles_k,  &
              reflectivity   = reflectivity,    &
              reflectivity_k = reflectivity_k,  &
              init           = .TRUE._jplm)
  THROWM(err.NE.0, 'Error allocating RTTOV structures')

  ! use_chan array is dimensioned by the total number of instrument channels
  ALLOCATE(use_chan(nprof,rth%coefs%coef%fmv_chn))

  ! Set use_chan to .TRUE. only for required channels
  use_chan(:,:) = .FALSE._jplm
  DO i = 1, nprof
    use_chan(i,channel_list(1:nchannels)) = .TRUE._jplm
  ENDDO

  ! Populate chanprof and frequencies arrays
  CALL rttov_scatt_setupindex ( &
      err,                    &
      nprof,                  &
      rth%coefs%coef%fmv_chn, &
      rth%coefs,              &
      rth%coefs_scatt,        &
      nchanprof,              &
      chanprof,               &
      frequencies,            &
      use_chan)
  THROWM(err.NE.0, 'Error determining channels, frequencies and polarisations from coefficient files')

  IF (rth%verbose) THEN
    WRITE(msg, '(a,i4,a,i8)') &
      'Running RTTOV-SCATT K using nthreads = ', rth%nthreads, &
      ' and nprofs_per_call = ', rth%nprofs_per_call
    INFO(msg)
  ENDIF

  ! Main loop over batches of profiles

  ncalls = nprofiles / rth%nprofs_per_call
  IF (MOD(nprofiles, rth%nprofs_per_call) /= 0) ncalls = ncalls + 1
  thisnchanprof = nchanprof
  thisnprof = rth%nprofs_per_call

  DO i = 1, ncalls

    ! Calculate low/high profile and chanprof indices
    iprof1 = (i - 1) * rth%nprofs_per_call + 1
    IF (iprof1 + rth%nprofs_per_call > nprofiles) THEN
      thisnprof = nprofiles - iprof1 + 1
      thisnchanprof = thisnprof * nchannels
    ENDIF

    ! Populate the profile structure
    CALL rttov_scatt_copy_to_profiles(  &
        rth,                       &
        profiles(1:thisnprof),     &
        cld_profiles(1:thisnprof), &
        iprof1,                    &
        datetimes,                 &
        angles, surfgeom,          &
        surftype, skin, s2m,       &
        zeeman,                    &
        p, t, gas_units,           &
        gas_id, gases, ph, cfrac)

    ! Set surface emissivities: if user passed any non-negative
    ! emissivities then use these and set calcemis FALSE
    calcemis(:) = .TRUE.
    CALL rttov_init_emis_refl(emissivity)

    DO iprof = 1, thisnprof
      k = iprof1 + iprof - 1
      lochan = (iprof - 1) * nchannels + 1
      hichan = iprof * nchannels
      WHERE (surfemis(:,k) >= 0._jprb)
        emissivity(lochan:hichan)%emis_in = surfemis(:,k)
        calcemis(lochan:hichan) = .FALSE.
      ENDWHERE
    ENDDO

    ! Populate the radiance_k structure
    CALL rttov_init_rad(radiance_k)
    CALL rttov_copy_to_radiance_k( &
        radiance_k,                    &
        iprof1,                        &
        thisnprof,                     &
        bt_k,                          &
        reflectivity_k=reflectivity_k, &
        zef_k=zef_k)

    ! Initialise other K arrays
    CALL rttov_init_prof(profiles_k)
    CALL rttov_init_scatt_prof(cld_profiles_k)
    CALL rttov_init_emis_refl(emissivity_k)

    ! Run RTTOV-SCATT K
    IF (rth%nthreads <= 1) THEN
      CALL rttov_scatt_ad ( &
          err,                             &
          rth%opts_scatt,                  &
          nlevels,                         &
          chanprof(1:thisnchanprof),       &
          frequencies(1:thisnchanprof),    &
          profiles(1:thisnprof),           &
          cld_profiles(1:thisnprof),       &
          rth%coefs,                       &
          rth%coefs_scatt,                 &
          calcemis(1:thisnchanprof),       &
          emissivity(1:thisnchanprof),     &
          profiles_k(1:thisnchanprof),     &
          cld_profiles_k(1:thisnchanprof), &
          emissivity_k(1:thisnchanprof),   &
          radiance,                        &
          radiance_k,                      &
          reflectivity=reflectivity,       &
          reflectivity_ad=reflectivity_k)
    ELSE
      CALL rttov_parallel_scatt_ad ( &
          err,                             &
          rth%opts_scatt,                  &
          nlevels,                         &
          chanprof(1:thisnchanprof),       &
          frequencies(1:thisnchanprof),    &
          profiles(1:thisnprof),           &
          cld_profiles(1:thisnprof),       &
          rth%coefs,                       &
          rth%coefs_scatt,                 &
          calcemis(1:thisnchanprof),       &
          emissivity(1:thisnchanprof),     &
          profiles_k(1:thisnchanprof),     &
          cld_profiles_k(1:thisnchanprof), &
          emissivity_k(1:thisnchanprof),   &
          radiance,                        &
          radiance_k,                      &
          reflectivity=reflectivity,       &
          reflectivity_ad=reflectivity_k,  &
          nthreads = rth%nthreads)
    ENDIF
    THROWM(err.NE.0, 'Error running RTTOV-SCATT K')

    ! Store results
    DO iprof = 1, thisnprof
      k = iprof1 + iprof - 1
      lochan = (iprof - 1) * nchannels + 1
      hichan = iprof * nchannels
      bt(:,k) = radiance%bt(lochan:hichan)
      surfemis(:,k) = emissivity(lochan:hichan)%emis_out
      surfemis_k(:,k) = emissivity_k(lochan:hichan)%emis_in
    ENDDO

    CALL rttov_copy_from_profiles_k( &
        rth,                         &
        profiles_k(1:thisnchanprof), &
        iprof1,                      &
        skin_k, s2m_k,               &
        t_k, gas_id, gases_k,        &
        p_k = p_k)

    CALL rttov_copy_from_cld_profiles_k( &
        cld_profiles_k(1:thisnchanprof), &
        iprof1,                          &
        gas_id, gases_k,                 &
        ph_k, cfrac_k)

    lochan = (iprof1 - 1) * nchannels + 1
    hichan = (iprof1 + thisnprof - 1) * nchannels

    IF (rth%store_rad) THEN
      rth%radiance%clear(lochan:hichan)              = radiance%clear(1:thisnchanprof)
      rth%radiance%total(lochan:hichan)              = radiance%total(1:thisnchanprof)
      rth%radiance%bt_clear(lochan:hichan)           = radiance%bt_clear(1:thisnchanprof)
      rth%radiance%bt(lochan:hichan)                 = radiance%bt(1:thisnchanprof)
      rth%radiance%quality(lochan:hichan)            = radiance%quality(1:thisnchanprof)
      rth%radiance%plane_parallel                    = radiance%plane_parallel
      rth%radiance%geometric_height(:,lochan:hichan) = radiance%geometric_height(:,1:thisnchanprof)
      IF (ASSOCIATED(reflectivity)) THEN
        rth%reflectivity%zef(:,lochan:hichan)        = reflectivity%zef(:,1:thisnchanprof)
        rth%reflectivity%azef(:,lochan:hichan)       = reflectivity%azef(:,1:thisnchanprof)
      ENDIF
    ENDIF

  ENDDO ! loop over calls to RTTOV

  ! Deallocate RTTOV data
  asw = 0_jpim
  CALL rttov_alloc_k( &
              err,                              &
              asw,                              &
              nprof,                            &
              nchanprof,                        &
              nlevels,                          &
              chanprof,                         &
              rth%opts,                         &
              profiles,                         &
              profiles_k,                       &
              rth%coefs,                        &
              radiance       = radiance,        &
              radiance_k     = radiance_k,      &
              calcemis       = calcemis,        &
              emissivity     = emissivity,      &
              emissivity_k   = emissivity_k,    &
              frequencies    = frequencies,     &
              coef_scatt     = rth%coefs_scatt, &
              nhydro_frac    = nhydro_frac,     &
              cld_profiles   = cld_profiles,    &
              cld_profiles_k = cld_profiles_k,  &
              reflectivity   = reflectivity,    &
              reflectivity_k = reflectivity_k)
  THROWM(err.NE.0, 'Error deallocating RTTOV structures')

  IF (calc_zef) DEALLOCATE(reflectivity, reflectivity_k)

  DEALLOCATE(use_chan, stat=err)
  THROWM(err.NE.0, 'Error dellocating RTTOV-SCATT arrays')

  CATCH
END SUBROUTINE rttov_scatt_call_k


!> Sets/updates the options associated with an instrument
!! @param[out]    err         return status
!! @param[in]     inst_id     ID for instrument
!! @param[in]     opts_str    string specifying options (see wrapper user guide)
SUBROUTINE rttov_set_options(err, inst_id, opts_str)
!
! Set RTTOV options (once instrument has been initialised
!   with a call to rttov_load_inst). This also sets values
!   in the RTTOV-SCATT options structure.
!
#include "throw.h"

  USE rttov_wrapper_handle
  USE rttov_getoptions
  USE parkind1, ONLY : jprb, jpim, jplm

  IMPLICIT NONE

  INTEGER(jpim),    INTENT(OUT) :: err
  INTEGER(jpim),    INTENT(IN)  :: inst_id
  CHARACTER(LEN=*), INTENT(IN)  :: opts_str

  INTEGER(jpim)       :: ival
  REAL(jprb)          :: rval
  CHARACTER(LEN=256)  :: strval
  LOGICAL(jplm)       :: exists
  TYPE(rttovwrapperhandle_type), POINTER :: rth

#include "rttov_user_options_checkinput.interface"
#include "rttov_errorreport.interface"

!------------------------------------------------------------------------------

  TRY

  CALL rttov_wrapper_handle_get_inst(err, inst_id, rth)
  THROW(err.NE.0)

  CALL initoptionsstr(opts_str(1:LEN(opts_str)))


  ! General config options

  CALL getoption('opts%config%apply_reg_limits', ival, exists = exists)
  IF (exists) THEN
    rth%opts%config%apply_reg_limits = (ival .NE. 0)
    rth%opts_scatt%config%apply_reg_limits = (ival .NE. 0)
  ENDIF

  CALL getoption('opts_scatt%config%apply_reg_limits', ival, exists = exists)
  IF (exists) THEN
    rth%opts%config%apply_reg_limits = (ival .NE. 0)
    rth%opts_scatt%config%apply_reg_limits = (ival .NE. 0)
  ENDIF

  CALL getoption('opts%config%verbose', ival, exists = exists)
  IF (exists) THEN
    rth%opts%config%verbose = (ival .NE. 0)
    rth%opts_scatt%config%verbose = (ival .NE. 0)
  ENDIF

  CALL getoption('opts_scatt%config%verbose', ival, exists = exists)
  IF (exists) THEN
    rth%opts%config%verbose = (ival .NE. 0)
    rth%opts_scatt%config%verbose = (ival .NE. 0)
  ENDIF

  CALL getoption('opts%config%do_checkinput', ival, exists = exists)
  IF (exists) THEN
    rth%opts%config%do_checkinput = (ival .NE. 0)
    rth%opts_scatt%config%do_checkinput = (ival .NE. 0)
  ENDIF

  CALL getoption('opts_scatt%config%do_checkinput', ival, exists = exists)
  IF (exists) THEN
    rth%opts%config%do_checkinput = (ival .NE. 0)
    rth%opts_scatt%config%do_checkinput = (ival .NE. 0)
  ENDIF

  CALL getoption('opts%config%fix_hgpl', ival, exists = exists)
  IF (exists) THEN
    rth%opts%config%fix_hgpl = (ival .NE. 0)
    rth%opts_scatt%config%fix_hgpl = (ival .NE. 0)
  ENDIF

  CALL getoption('opts_scatt%config%fix_hgpl', ival, exists = exists)
  IF (exists) THEN
    rth%opts%config%fix_hgpl = (ival .NE. 0)
    rth%opts_scatt%config%fix_hgpl = (ival .NE. 0)
  ENDIF


  ! Interpolation options

  CALL getoption('opts%interpolation%addinterp', ival, exists = exists)
  IF (exists) rth%opts%interpolation%addinterp = (ival .NE. 0)

  CALL getoption('opts%interpolation%interp_mode', ival, exists = exists)
  IF (exists) THEN
    rth%opts%interpolation%interp_mode = ival
    rth%opts_scatt%interp_mode = ival
  ENDIF

  CALL getoption('opts_scatt%interp_mode', ival, exists = exists)
  IF (exists) THEN
    rth%opts%interpolation%interp_mode = ival
    rth%opts_scatt%interp_mode = ival
  ENDIF

  CALL getoption('opts%interpolation%lgradp', ival, exists = exists)
  IF (exists) THEN
    rth%opts%interpolation%lgradp = (ival .NE. 0)
    rth%opts_scatt%lgradp = (ival .NE. 0)
  ENDIF

  CALL getoption('opts_scatt%lgradp', ival, exists = exists)
  IF (exists) THEN
    rth%opts%interpolation%lgradp = (ival .NE. 0)
    rth%opts_scatt%lgradp = (ival .NE. 0)
  ENDIF

  CALL getoption('opts%interpolation%spacetop', ival, exists = exists)
  IF (exists) rth%opts%interpolation%spacetop = (ival .NE. 0)

  CALL getoption('opts%interpolation%reg_limit_extrap', ival, exists = exists)
  IF (exists) THEN
    rth%opts%interpolation%reg_limit_extrap = (ival .NE. 0)
    rth%opts_scatt%reg_limit_extrap = (ival .NE. 0)
  ENDIF

  CALL getoption('opts_scatt%reg_limit_extrap', ival, exists = exists)
  IF (exists) THEN
    rth%opts%interpolation%reg_limit_extrap = (ival .NE. 0)
    rth%opts_scatt%reg_limit_extrap = (ival .NE. 0)
  ENDIF


  ! General RT options

  CALL getoption('opts%rt_all%addrefrac', ival, exists = exists)
  IF (exists) THEN
    rth%opts%rt_all%addrefrac = (ival .NE. 0)
    rth%opts_scatt%addrefrac = (ival .NE. 0)
  ENDIF

  CALL getoption('opts_scatt%addrefrac', ival, exists = exists)
  IF (exists) THEN
    rth%opts%rt_all%addrefrac = (ival .NE. 0)
    rth%opts_scatt%addrefrac = (ival .NE. 0)
  ENDIF

  CALL getoption('opts%rt_all%plane_parallel', ival, exists = exists)
  IF (exists) rth%opts%rt_all%plane_parallel = (ival .NE. 0)

  CALL getoption('opts%rt_all%switchrad', ival, exists = exists)
  IF (exists) rth%opts%rt_all%switchrad = (ival .NE. 0)

  CALL getoption('opts%rt_all%use_t2m_opdep', ival, exists = exists)
  IF (exists) THEN
    rth%opts%rt_all%use_t2m_opdep = (ival .NE. 0)
    rth%opts_scatt%use_t2m_opdep = (ival .NE. 0)
  ENDIF

  CALL getoption('opts_scatt%use_t2m_opdep', ival, exists = exists)
  IF (exists) THEN
    rth%opts%rt_all%use_t2m_opdep = (ival .NE. 0)
    rth%opts_scatt%use_t2m_opdep = (ival .NE. 0)
  ENDIF

  CALL getoption('opts%rt_all%use_q2m', ival, exists = exists)
  IF (exists) THEN
    rth%opts%rt_all%use_q2m = (ival .NE. 0)
    rth%opts_scatt%use_q2m = (ival .NE. 0)
  ENDIF

  CALL getoption('opts_scatt%use_q2m', ival, exists = exists)
  IF (exists) THEN
    rth%opts%rt_all%use_q2m = (ival .NE. 0)
    rth%opts_scatt%use_q2m = (ival .NE. 0)
  ENDIF

  CALL getoption('opts%rt_all%do_lambertian', ival, exists = exists)
  IF (exists) rth%opts%rt_all%do_lambertian = (ival .NE. 0)

  CALL getoption('opts%rt_all%lambertian_fixed_angle', ival, exists = exists)
  IF (exists) rth%opts%rt_all%lambertian_fixed_angle = (ival .NE. 0)

  CALL getoption('opts%rt_all%rad_down_lin_tau', ival, exists = exists)
  IF (exists) THEN
    rth%opts%rt_all%rad_down_lin_tau = (ival .NE. 0)
    rth%opts_scatt%rad_down_lin_tau = (ival .NE. 0)
  ENDIF

  CALL getoption('opts_scatt%rad_down_lin_tau', ival, exists = exists)
  IF (exists) THEN
    rth%opts%rt_all%rad_down_lin_tau = (ival .NE. 0)
    rth%opts_scatt%rad_down_lin_tau = (ival .NE. 0)
  ENDIF

  CALL getoption('opts%rt_all%dtau_test', ival, exists = exists)
  IF (exists) THEN
    rth%opts%rt_all%dtau_test = (ival .NE. 0)
    rth%opts_scatt%dtau_test = (ival .NE. 0)
  ENDIF

  CALL getoption('opts_scatt%dtau_test', ival, exists = exists)
  IF (exists) THEN
    rth%opts%rt_all%dtau_test = (ival .NE. 0)
    rth%opts_scatt%dtau_test = (ival .NE. 0)
  ENDIF


  CALL getoption('opts%rt_all%ozone_data', ival, exists = exists)
  IF (exists) THEN
    rth%opts%rt_all%ozone_data = (ival .NE. 0)
    rth%opts_scatt%ozone_data = (ival .NE. 0)
  ENDIF

  CALL getoption('opts_scatt%ozone_data', ival, exists = exists)
  IF (exists) THEN
    rth%opts%rt_all%ozone_data = (ival .NE. 0)
    rth%opts_scatt%ozone_data = (ival .NE. 0)
  ENDIF

  CALL getoption('opts%rt_all%co2_data', ival, exists = exists)
  IF (exists) rth%opts%rt_all%co2_data = (ival .NE. 0)

  CALL getoption('opts%rt_all%n2o_data', ival, exists = exists)
  IF (exists) rth%opts%rt_all%n2o_data = (ival .NE. 0)

  CALL getoption('opts%rt_all%co_data', ival, exists = exists)
  IF (exists) rth%opts%rt_all%co_data = (ival .NE. 0)

  CALL getoption('opts%rt_all%ch4_data', ival, exists = exists)
  IF (exists) rth%opts%rt_all%ch4_data = (ival .NE. 0)

  CALL getoption('opts%rt_all%so2_data', ival, exists = exists)
  IF (exists) rth%opts%rt_all%so2_data = (ival .NE. 0)



  ! MW RT options

  CALL getoption('opts%rt_mw%fastem_version', ival, exists = exists)
  IF (exists) THEN
    rth%opts%rt_mw%fastem_version = ival
    rth%opts_scatt%fastem_version = ival
  ENDIF

  CALL getoption('opts_scatt%fastem_version', ival, exists = exists)
  IF (exists) THEN
    rth%opts%rt_mw%fastem_version = ival
    rth%opts_scatt%fastem_version = ival
  ENDIF

  CALL getoption('opts%rt_mw%supply_foam_fraction', ival, exists = exists)
  IF (exists) THEN
    rth%opts%rt_mw%supply_foam_fraction = (ival .NE. 0)
    rth%opts_scatt%supply_foam_fraction = (ival .NE. 0)
  ENDIF

  CALL getoption('opts_scatt%supply_foam_fraction', ival, exists = exists)
  IF (exists) THEN
    rth%opts%rt_mw%supply_foam_fraction = (ival .NE. 0)
    rth%opts_scatt%supply_foam_fraction = (ival .NE. 0)
  ENDIF

  CALL getoption('opts%rt_mw%clw_data', ival, exists = exists)
  IF (exists) rth%opts%rt_mw%clw_data = (ival .NE. 0)

  CALL getoption('opts%rt_mw%clw_scheme', ival, exists = exists)
  IF (exists) rth%opts%rt_mw%clw_scheme = ival

  CALL getoption('opts%rt_mw%clw_cloud_top', rval, exists = exists)
  IF (exists) rth%opts%rt_mw%clw_cloud_top = rval


  ! RTTOV-SCATT-only options

  CALL getoption('opts_scatt%ice_polarisation', rval, exists = exists)
  IF (exists) rth%opts_scatt%ice_polarisation = rval

  CALL getoption('opts_scatt%lusercfrac', ival, exists = exists)
  IF (exists) rth%opts_scatt%lusercfrac = (ival .NE. 0)

  CALL getoption('opts_scatt%cc_threshold', rval, exists = exists)
  IF (exists) rth%opts_scatt%cc_threshold = rval

  CALL getoption('opts_scatt%hydro_cfrac_tlad', ival, exists = exists)
  IF (exists) rth%opts_scatt%hydro_cfrac_tlad = (ival .NE. 0)

  CALL getoption('opts_scatt%zero_hydro_tlad', ival, exists = exists)
  IF (exists) rth%opts_scatt%zero_hydro_tlad = (ival .NE. 0)


  ! VIS/IR RT options

  CALL getoption('opts%rt_ir%solar_sea_brdf_model', ival, exists = exists)
  IF (exists) rth%opts%rt_ir%solar_sea_brdf_model = ival

  CALL getoption('opts%rt_ir%ir_sea_emis_model', ival, exists = exists)
  IF (exists) rth%opts%rt_ir%ir_sea_emis_model = ival

  CALL getoption('opts%rt_ir%addsolar', ival, exists = exists)
  IF (exists) rth%opts%rt_ir%addsolar = (ival .NE. 0)

  CALL getoption('opts%rt_ir%rayleigh_max_wavelength', rval, exists = exists)
  IF (exists) rth%opts%rt_ir%rayleigh_max_wavelength = rval

  CALL getoption('opts%rt_ir%rayleigh_min_pressure', rval, exists = exists)
  IF (exists) rth%opts%rt_ir%rayleigh_min_pressure = rval

  CALL getoption('opts%rt_ir%rayleigh_single_scatt', ival, exists = exists)
  IF (exists) rth%opts%rt_ir%rayleigh_single_scatt = (ival .NE. 0)

  CALL getoption('opts%rt_ir%do_nlte_correction', ival, exists = exists)
  IF (exists) rth%opts%rt_ir%do_nlte_correction = (ival .NE. 0)

  CALL getoption('opts%rt_ir%addaerosl', ival, exists = exists)
  IF (exists) rth%opts%rt_ir%addaerosl = (ival .NE. 0)

  CALL getoption('opts%rt_ir%user_aer_opt_param', ival, exists = exists)
  IF (exists) rth%opts%rt_ir%user_aer_opt_param = (ival .NE. 0)

  CALL getoption('opts%rt_ir%addclouds', ival, exists = exists)
  IF (exists) rth%opts%rt_ir%addclouds = (ival .NE. 0)

  CALL getoption('opts%rt_ir%user_cld_opt_param', ival, exists = exists)
  IF (exists) rth%opts%rt_ir%user_cld_opt_param = (ival .NE. 0)

  CALL getoption('opts%rt_ir%grid_box_avg_cloud', ival, exists = exists)
  IF (exists) rth%opts%rt_ir%grid_box_avg_cloud = (ival .NE. 0)

  CALL getoption('opts%rt_ir%cldcol_threshold', rval, exists = exists)
  IF (exists) rth%opts%rt_ir%cldcol_threshold = rval

  CALL getoption('opts%rt_ir%cloud_overlap', ival, exists = exists)
  IF (exists) rth%opts%rt_ir%cloud_overlap = ival

  CALL getoption('opts%rt_ir%cc_low_cloud_top', rval, exists = exists)
  IF (exists) rth%opts%rt_ir%cc_low_cloud_top = rval

  CALL getoption('opts%rt_ir%ir_scatt_model', ival, exists = exists)
  IF (exists) rth%opts%rt_ir%ir_scatt_model = ival

  CALL getoption('opts%rt_ir%vis_scatt_model', ival, exists = exists)
  IF (exists) rth%opts%rt_ir%vis_scatt_model = ival

  CALL getoption('opts%rt_ir%dom_nstreams', ival, exists = exists)
  IF (exists) rth%opts%rt_ir%dom_nstreams = ival

  CALL getoption('opts%rt_ir%dom_accuracy', rval, exists = exists)
  IF (exists) rth%opts%rt_ir%dom_accuracy = rval

  CALL getoption('opts%rt_ir%dom_opdep_threshold', rval, exists = exists)
  IF (exists) rth%opts%rt_ir%dom_opdep_threshold = rval

  CALL getoption('opts%rt_ir%dom_rayleigh', ival, exists = exists)
  IF (exists) rth%opts%rt_ir%dom_rayleigh = (ival .NE. 0)


  ! PC options

  CALL getoption('opts%rt_ir%pc%addpc', ival, exists = exists)
  IF (exists) rth%opts%rt_ir%pc%addpc = (ival .NE. 0)

  CALL getoption('opts%rt_ir%pc%addradrec', ival, exists = exists)
  IF (exists) rth%opts%rt_ir%pc%addradrec = (ival .NE. 0)

  CALL getoption('opts%rt_ir%pc%npcscores', ival, exists = exists)
  IF (exists) rth%opts%rt_ir%pc%npcscores = ival

  CALL getoption('opts%rt_ir%pc%ipcbnd', ival, exists = exists)
  IF (exists) rth%opts%rt_ir%pc%ipcbnd = ival

  CALL getoption('opts%rt_ir%pc%ipcreg', ival, exists = exists)
  IF (exists) rth%opts%rt_ir%pc%ipcreg = ival


  ! HTFRTC options

  CALL getoption('opts%htfrtc_opts%htfrtc', ival, exists = exists)
  IF (exists) rth%opts%htfrtc_opts%htfrtc = (ival .NE. 0)

  CALL getoption('opts%htfrtc_opts%reconstruct', ival, exists = exists)
  IF (exists) rth%opts%htfrtc_opts%reconstruct = (ival .NE. 0)

  CALL getoption('opts%htfrtc_opts%n_pc_in', ival, exists = exists)
  IF (exists) rth%opts%htfrtc_opts%n_pc_in = ival

  CALL getoption('opts%htfrtc_opts%simple_cloud', ival, exists = exists)
  IF (exists) rth%opts%htfrtc_opts%simple_cloud = (ival .NE. 0)

  CALL getoption('opts%htfrtc_opts%overcast', ival, exists = exists)
  IF (exists) rth%opts%htfrtc_opts%overcast = (ival .NE. 0)


  ! Developer options

  CALL getoption('opts%dev%do_opdep_calc', ival, exists = exists)
  IF (exists) rth%opts%dev%do_opdep_calc = (ival .NE. 0)


  ! Wrapper options

  CALL getoption('nthreads', ival, exists = exists)
  IF (exists) rth%nthreads = ival

  CALL getoption('nprofs_per_call', ival, exists = exists)
  IF (exists) rth%nprofs_per_call = ival
  IF (rth%nprofs_per_call < 1) rth%nprofs_per_call = 1

  CALL getoption('verbose_wrapper', ival, exists = exists)
  IF (exists) rth%verbose = (ival .NE. 0)

  CALL getoption('check_opts', ival, exists = exists)
  IF (exists) rth%check_opts = (ival .NE. 0)

  CALL getoption('store_trans', ival, exists = exists)
  IF (exists) rth%store_trans = (ival .NE. 0)

  CALL getoption('store_rad', ival, exists = exists)
  IF (exists) rth%store_rad = (ival .NE. 0)

  CALL getoption('store_rad2', ival, exists = exists)
  IF (exists) rth%store_rad2 = (ival .NE. 0)

  ! Radiance structure must be allocated if radiance2 is allocated
  IF (rth%store_rad2) rth%store_rad = .TRUE.

  CALL getoption('store_emis_terms', ival, exists = exists)
  IF (exists) rth%store_emis_terms = (ival .NE. 0)


  CALL getoption('file_coef', strval, exists = exists)
  IF (exists) rth%file_coef = TRIM(strval)
  CALL getoption('file_scaer', strval, exists = exists)
  IF (exists) rth%file_scaer = TRIM(strval)
  CALL getoption('file_sccld', strval, exists = exists)
  IF (exists) rth%file_sccld = TRIM(strval)
  CALL getoption('file_mfasis_cld', strval, exists = exists)
  IF (exists) rth%file_mfasis_cld = TRIM(strval)
!   CALL getoption('file_mfasis_aer', strval, exists = exists)
!   IF (exists) rth%file_mfasis_aer = TRIM(strval)
!   CALL getoption('file_pccoef', strval, exists = exists)
!   IF (exists) rth%file_pccoef = TRIM(strval)
  CALL getoption('file_hydrotable', strval, exists = exists)
  IF (exists) rth%file_hydrotable = TRIM(strval)


  CALL checkoptions()

  ! When this is called from rttov_load_inst the coefs haven't been loaded (init = FALSE)
  ! so we shouldn't call rttov_user_options_checkinput at this point: only call it when
  ! user is calling rttov_set_options after initialising an intrument.
  IF (rth%init .AND. rth%check_opts) THEN
    CALL rttov_user_options_checkinput(err, rth%opts, rth%coefs)
    THROW(err.NE.0)
  ENDIF

  CATCH
END SUBROUTINE rttov_set_options

!> Prints the options associated with an instrument
!! @param[out]    err         return status
!! @param[in]     inst_id     ID for instrument
SUBROUTINE rttov_print_options(err, inst_id)
!
! Print current RTTOV/RTTOV-SCATT options for an instrument
!
#include "throw.h"

  USE rttov_wrapper_handle
  USE parkind1, ONLY : jpim
  USE rttov_global, ONLY : error_unit

  IMPLICIT NONE

  INTEGER(jpim),    INTENT(OUT) :: err
  INTEGER(jpim),    INTENT(IN)  :: inst_id

  CHARACTER(LEN=20)   :: tmp_text ! temporary string for formatting output
  TYPE(rttovwrapperhandle_type), POINTER :: rth

#include "rttov_print_opts.interface"
#include "rttov_print_opts_scatt.interface"
#include "rttov_errorreport.interface"

!------------------------------------------------------------------------------

  TRY

  CALL rttov_wrapper_handle_get_inst(err, inst_id, rth)
  THROW(err.NE.0)

  IF (rth%file_hydrotable .NE. '') THEN
    CALL rttov_print_opts_scatt(rth%opts_scatt)
  ELSE
    CALL rttov_print_opts(rth%opts)
  ENDIF

  WRITE(error_unit,'(a)')           "Wrapper options"
  WRITE(error_unit,'(2x,a,a)')      "Coef file            ", TRIM(rth%file_coef)
  IF (rth%file_scaer .NE. '') &
    WRITE(error_unit,'(2x,a,a)')      "Aerosol coef file    ", TRIM(rth%file_scaer)
  IF (rth%file_sccld .NE. '') &
    WRITE(error_unit,'(2x,a,a)')      "Cloud coef file      ", TRIM(rth%file_sccld)
  IF (rth%file_mfasis_cld .NE. '') &
    WRITE(error_unit,'(2x,a,a)')      "MFASIS cloud file    ", TRIM(rth%file_mfasis_cld)
  IF (rth%file_mfasis_aer .NE. '') &
    WRITE(error_unit,'(2x,a,a)')      "MFASIS aerosol file  ", TRIM(rth%file_mfasis_aer)
  IF (rth%file_pccoef .NE. '') &
    WRITE(error_unit,'(2x,a,a)')      "PC-RTTOV coef file   ", rth%file_pccoef
  IF (rth%file_hydrotable .NE. '') &
    WRITE(error_unit,'(2x,a,a)')      "RTTOV-SCATT hydrotable ", rth%file_hydrotable
  WRITE(tmp_text,'(i6)') rth%nprofs_per_call
  WRITE(error_unit,'(2x,a,a)')      "nprofs_per_call      ", TRIM(ADJUSTL(tmp_text))
  WRITE(tmp_text,'(i5)') rth%nthreads
  WRITE(error_unit,'(2x,a,a)')      "nthreads             ", TRIM(ADJUSTL(tmp_text))
  WRITE(error_unit,'(2x,a,l1)')     "verbose_wrapper      ", rth%verbose
  WRITE(error_unit,'(2x,a,l1)')     "check_opts           ", rth%check_opts
  WRITE(error_unit,'(2x,a,l1)')     "store_rad            ", rth%store_rad
  IF (rth%file_hydrotable .NE. '') THEN
    WRITE(error_unit,'(2x,a,l1)')     "store_emis_terms     ", rth%store_emis_terms
  ELSE
    WRITE(error_unit,'(2x,a,l1)')     "store_rad2           ", rth%store_rad2
    WRITE(error_unit,'(2x,a,l1)')     "store_trans          ", rth%store_trans
  ENDIF
  WRITE(error_unit,'(a)',advance='yes')

  CATCH
END SUBROUTINE rttov_print_options


! The following subroutines wrap the RTTOV bpr and Legendre coef calculation
! subroutines

!> Calculate bpr values for Chou-scaling IR scattering parameterisation
!! @param[out]    err         return status
!! @param[in]     phangle     phase function angle grid (nphangle)
!! @param[in]     pha         phase function (nphangle)
!! @param[out]    bpr         computed bpr value
!! @param[in]     nthreads    number of threads to use (RTTOV must have been compiled with OpenMP)
!! @param[in]     nphangle    size of phangle (not required in Python)
SUBROUTINE rttov_bpr(err, phangle, pha, bpr, nthreads, nphangle)
!f2py intent(hide):: nphangle=len(phangle)
#include "throw.h"
  USE rttov_wrapper_handle
  USE parkind1, ONLY : jpim, jprb

  IMPLICIT NONE

  INTEGER(jpim), INTENT(OUT)   :: err
  INTEGER(jpim), INTENT(IN)    :: nthreads
  INTEGER(jpim), INTENT(IN)    :: nphangle
  REAL(jprb),    INTENT(IN)    :: phangle(nphangle)
  REAL(jprb),    INTENT(IN)    :: pha(nphangle)
  REAL(jprb),    INTENT(OUT)   :: bpr

#include "rttov_errorreport.interface"
#include "rttov_bpr_init.interface"
#include "rttov_bpr_calc.interface"
#include "rttov_bpr_dealloc.interface"

  TRY

  CALL rttov_bpr_init(err, phangle)
  THROWM(err.NE.0, 'Error initialising bpr LUT')

  CALL rttov_bpr_calc(err, pha, phangle, bpr, nthreads)
  THROWM(err.NE.0, 'Error calculating bpr')

  CALL rttov_bpr_dealloc(err)
  THROWM(err.NE.0, 'Error deallocating bpr LUT')

  CATCH
END SUBROUTINE rttov_bpr

!> Calculate Legendre coefficients for DOM visible/IR scattering solver
!! @param[out]    err         return status
!! @param[in]     phangle     phase function angle grid (nphangle)
!! @param[in]     pha         phase function (nphangle)
!! @param[in,out] legcoef     computed Legendre coefficients (nmom+1)
!! @param[in]     ngauss      size of Gaussian quadrature to use (uses default if <=nmom)
!! @param[in]     nphangle    size of phangle (not required in Python)
!! @param[in]     nmom        size of legcoef minus one (not required in Python)
SUBROUTINE rttov_legcoef(err, phangle, pha, legcoef, ngauss, nphangle, nmom)
!f2py intent(hide):: nphangle=len(phangle)
!f2py intent(hide):: nmom=len(legcoef)-1
#include "throw.h"
  USE rttov_wrapper_handle
  USE parkind1, ONLY : jpim, jprb

  IMPLICIT NONE

  INTEGER(jpim), INTENT(OUT)   :: err
  INTEGER(jpim), INTENT(IN)    :: nphangle
  INTEGER(jpim), INTENT(IN)    :: nmom
  REAL(jprb),    INTENT(IN)    :: phangle(nphangle)
  REAL(jprb),    INTENT(IN)    :: pha(nphangle)
  REAL(jprb),    INTENT(INOUT) :: legcoef(nmom+1)
  INTEGER(jpim), INTENT(IN)    :: ngauss

#include "rttov_errorreport.interface"
#include "rttov_legcoef_calc.interface"

  TRY

  IF (ngauss >= nmom) THEN
    CALL rttov_legcoef_calc(err, pha, phangle, nmom, legcoef, ngauss)
  ELSE
    CALL rttov_legcoef_calc(err, pha, phangle, nmom, legcoef)
  ENDIF
  THROWM(err.NE.0, 'Error calculating Legendre coefficients')

  CATCH
END SUBROUTINE rttov_legcoef


! The following subroutines return RTTOV output arrays: they can be called
! after RTTOV has been called.

!> Return clear-sky radiances
!! @param[out]    err         return status
!! @param[in]     inst_id     ID for instrument
!! @param[in,out] rad_clear   clear-sky radiances (nchanprof)
!! @param[in]     nchanprof   size of rad_clear (not required in Python)
SUBROUTINE rttov_get_rad_clear(err, inst_id, rad_clear, nchanprof)
!f2py intent(hide):: nchanprof=len(rad_clear)
#include "throw.h"
  USE rttov_wrapper_handle
  USE parkind1, ONLY : jpim, jprb

  IMPLICIT NONE

  INTEGER(jpim), INTENT(OUT)   :: err
  INTEGER(jpim), INTENT(IN)    :: inst_id
  INTEGER(jpim), INTENT(IN)    :: nchanprof
  REAL(jprb),    INTENT(INOUT) :: rad_clear(nchanprof)

  TYPE(rttovwrapperhandle_type), POINTER :: rth

#include "rttov_errorreport.interface"

  TRY
  CALL rttov_wrapper_handle_get_inst(err, inst_id, rth)
  THROW(err.NE.0)

  IF (.NOT. ASSOCIATED(rth%radiance%clear)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Wrapper radiance not allocated: has RTTOV been called?')
  ENDIF

  rad_clear(:) = rth%radiance%clear(:)
  CATCH
END SUBROUTINE rttov_get_rad_clear

!> Return total radiances
!! @param[out]    err         return status
!! @param[in]     inst_id     ID for instrument
!! @param[in,out] rad_total   total radiances (nchanprof)
!! @param[in]     nchanprof   size of rad_total (not required in Python)
SUBROUTINE rttov_get_rad_total(err, inst_id, rad_total, nchanprof)
!f2py intent(hide):: nchanprof=len(rad_total)
#include "throw.h"
  USE rttov_wrapper_handle
  USE parkind1, ONLY : jpim, jprb

  IMPLICIT NONE

  INTEGER(jpim), INTENT(OUT)   :: err
  INTEGER(jpim), INTENT(IN)    :: inst_id
  INTEGER(jpim), INTENT(IN)    :: nchanprof
  REAL(jprb),    INTENT(INOUT) :: rad_total(nchanprof)

  TYPE(rttovwrapperhandle_type), POINTER :: rth

#include "rttov_errorreport.interface"

  TRY
  CALL rttov_wrapper_handle_get_inst(err, inst_id, rth)
  THROW(err.NE.0)

  IF (.NOT. ASSOCIATED(rth%radiance%total)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Wrapper radiance not allocated: has RTTOV been called?')
  ENDIF

  rad_total(:) = rth%radiance%total(:)
  CATCH
END SUBROUTINE rttov_get_rad_total

!> Return 100% cloudy radiances
!! @param[out]    err         return status
!! @param[in]     inst_id     ID for instrument
!! @param[in,out] rad_cloudy  100% cloudy radiances (nchanprof)
!! @param[in]     nchanprof   size of rad_cloudy (not required in Python)
SUBROUTINE rttov_get_rad_cloudy(err, inst_id, rad_cloudy, nchanprof)
!f2py intent(hide):: nchanprof=len(rad_cloudy)
#include "throw.h"
  USE rttov_wrapper_handle
  USE parkind1, ONLY : jpim, jprb

  IMPLICIT NONE

  INTEGER(jpim), INTENT(OUT)   :: err
  INTEGER(jpim), INTENT(IN)    :: inst_id
  INTEGER(jpim), INTENT(IN)    :: nchanprof
  REAL(jprb),    INTENT(INOUT) :: rad_cloudy(nchanprof)

  TYPE(rttovwrapperhandle_type), POINTER :: rth

#include "rttov_errorreport.interface"

  TRY
  CALL rttov_wrapper_handle_get_inst(err, inst_id, rth)
  THROW(err.NE.0)

  IF (.NOT. ASSOCIATED(rth%radiance%cloudy)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Wrapper radiance not allocated: has RTTOV been called?')
  ENDIF

  rad_cloudy(:) = rth%radiance%cloudy(:)
  CATCH
END SUBROUTINE rttov_get_rad_cloudy

!> Return clear-sky BTs
!! @param[out]    err         return status
!! @param[in]     inst_id     ID for instrument
!! @param[in,out] bt_clear    clear-sky BTs (nchanprof)
!! @param[in]     nchanprof   size of bt_clear (not required in Python)
SUBROUTINE rttov_get_bt_clear(err, inst_id, bt_clear, nchanprof)
!f2py intent(hide):: nchanprof=len(bt_clear)
#include "throw.h"
  USE rttov_wrapper_handle
  USE parkind1, ONLY : jpim, jprb

  IMPLICIT NONE

  INTEGER(jpim), INTENT(OUT)   :: err
  INTEGER(jpim), INTENT(IN)    :: inst_id
  INTEGER(jpim), INTENT(IN)    :: nchanprof
  REAL(jprb),    INTENT(INOUT) :: bt_clear(nchanprof)

  TYPE(rttovwrapperhandle_type), POINTER :: rth

#include "rttov_errorreport.interface"

  TRY
  CALL rttov_wrapper_handle_get_inst(err, inst_id, rth)
  THROW(err.NE.0)

  IF (.NOT. ASSOCIATED(rth%radiance%bt_clear)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Wrapper radiance not allocated: has RTTOV been called?')
  ENDIF

  bt_clear(:) = rth%radiance%bt_clear(:)
  CATCH
END SUBROUTINE rttov_get_bt_clear

!> Return total BTs
!! @param[out]    err         return status
!! @param[in]     inst_id     ID for instrument
!! @param[in,out] bt          total BTs (nchanprof)
!! @param[in]     nchanprof   size of bt (not required in Python)
SUBROUTINE rttov_get_bt(err, inst_id, bt, nchanprof)
!f2py intent(hide):: nchanprof=len(bt)
#include "throw.h"
  USE rttov_wrapper_handle
  USE parkind1, ONLY : jpim, jprb

  IMPLICIT NONE

  INTEGER(jpim), INTENT(OUT)   :: err
  INTEGER(jpim), INTENT(IN)    :: inst_id
  INTEGER(jpim), INTENT(IN)    :: nchanprof
  REAL(jprb),    INTENT(INOUT) :: bt(nchanprof)

  TYPE(rttovwrapperhandle_type), POINTER :: rth

#include "rttov_errorreport.interface"

  TRY
  CALL rttov_wrapper_handle_get_inst(err, inst_id, rth)
  THROW(err.NE.0)

  IF (.NOT. ASSOCIATED(rth%radiance%bt)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Wrapper radiance not allocated: has RTTOV been called?')
  ENDIF

  bt(:) = rth%radiance%bt(:)
  CATCH
END SUBROUTINE rttov_get_bt

!> Return clear-sky reflectances
!! @param[out]    err         return status
!! @param[in]     inst_id     ID for instrument
!! @param[in,out] refl_clear  clear-sky reflectances (nchanprof)
!! @param[in]     nchanprof   size of refl_clear (not required in Python)
SUBROUTINE rttov_get_refl_clear(err, inst_id, refl_clear, nchanprof)
!f2py intent(hide):: nchanprof=len(refl_clear)
#include "throw.h"
  USE rttov_wrapper_handle
  USE parkind1, ONLY : jpim, jprb

  IMPLICIT NONE

  INTEGER(jpim), INTENT(OUT)   :: err
  INTEGER(jpim), INTENT(IN)    :: inst_id
  INTEGER(jpim), INTENT(IN)    :: nchanprof
  REAL(jprb),    INTENT(INOUT) :: refl_clear(nchanprof)

  TYPE(rttovwrapperhandle_type), POINTER :: rth

#include "rttov_errorreport.interface"

  TRY
  CALL rttov_wrapper_handle_get_inst(err, inst_id, rth)
  THROW(err.NE.0)

  IF (.NOT. ASSOCIATED(rth%radiance%refl_clear)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Wrapper radiance not allocated: has RTTOV been called?')
  ENDIF

  refl_clear(:) = rth%radiance%refl_clear(:)
  CATCH
END SUBROUTINE rttov_get_refl_clear

!> Return total reflectances
!! @param[out]    err         return status
!! @param[in]     inst_id     ID for instrument
!! @param[in,out] refl        total reflectances (nchanprof)
!! @param[in]     nchanprof   size of refl (not required in Python)
SUBROUTINE rttov_get_refl(err, inst_id, refl, nchanprof)
!f2py intent(hide):: nchanprof=len(refl)
#include "throw.h"
  USE rttov_wrapper_handle
  USE parkind1, ONLY : jpim, jprb

  IMPLICIT NONE

  INTEGER(jpim), INTENT(OUT)   :: err
  INTEGER(jpim), INTENT(IN)    :: inst_id
  INTEGER(jpim), INTENT(IN)    :: nchanprof
  REAL(jprb),    INTENT(INOUT) :: refl(nchanprof)

  TYPE(rttovwrapperhandle_type), POINTER :: rth

#include "rttov_errorreport.interface"

  TRY
  CALL rttov_wrapper_handle_get_inst(err, inst_id, rth)
  THROW(err.NE.0)

  IF (.NOT. ASSOCIATED(rth%radiance%refl)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Wrapper radiance not allocated: has RTTOV been called?')
  ENDIF

  refl(:) = rth%radiance%refl(:)
  CATCH
END SUBROUTINE rttov_get_refl

!> Return overcast radiances
!! @param[out]    err         return status
!! @param[in]     inst_id     ID for instrument
!! @param[in,out] overcast    total radiances (nlayers,nchanprof)
!! @param[in]     nchanprof   size of overcast(1,:) (not required in Python)
!! @param[in]     nlayers     size of overcast(:,1) (not required in Python)
SUBROUTINE rttov_get_overcast(err, inst_id, overcast, nchanprof, nlayers)
!f2py intent(hide):: nchanprof=len(overcast)
!f2py intent(hide):: nlayers=len(overcast[0])
#include "throw.h"
  USE rttov_wrapper_handle
  USE parkind1, ONLY : jpim, jprb

  IMPLICIT NONE

  INTEGER(jpim), INTENT(OUT)   :: err
  INTEGER(jpim), INTENT(IN)    :: inst_id
  INTEGER(jpim), INTENT(IN)    :: nchanprof
  INTEGER(jpim), INTENT(IN)    :: nlayers
  REAL(jprb),    INTENT(INOUT) :: overcast(nlayers,nchanprof)

  TYPE(rttovwrapperhandle_type), POINTER :: rth

#include "rttov_errorreport.interface"

  TRY
  CALL rttov_wrapper_handle_get_inst(err, inst_id, rth)
  THROW(err.NE.0)

  IF (.NOT. ASSOCIATED(rth%radiance%overcast)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Wrapper radiance not allocated: has RTTOV been called?')
  ENDIF

  overcast(:,:) = rth%radiance%overcast(:,:)
  CATCH
END SUBROUTINE rttov_get_overcast

!> Return radiance quality flags
!! @param[out]    err         return status
!! @param[in]     inst_id     ID for instrument
!! @param[in,out] rad_quality quality flag (nchanprof)
!! @param[in]     nchanprof   size of rad_quality (not required in Python)
SUBROUTINE rttov_get_rad_quality(err, inst_id, rad_quality, nchanprof)
!f2py intent(hide):: nchanprof=len(rad_quality)
#include "throw.h"
  USE rttov_wrapper_handle
  USE parkind1, ONLY : jpim

  IMPLICIT NONE

  INTEGER(jpim), INTENT(OUT)   :: err
  INTEGER(jpim), INTENT(IN)    :: inst_id
  INTEGER(jpim), INTENT(IN)    :: nchanprof
  INTEGER(jpim), INTENT(INOUT) :: rad_quality(nchanprof)

  TYPE(rttovwrapperhandle_type), POINTER :: rth

#include "rttov_errorreport.interface"

  TRY
  CALL rttov_wrapper_handle_get_inst(err, inst_id, rth)
  THROW(err.NE.0)

  IF (.NOT. ASSOCIATED(rth%radiance%quality)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Wrapper radiance not allocated: has RTTOV been called?')
  ENDIF

  rad_quality(:) = rth%radiance%quality(:)
  CATCH
END SUBROUTINE rttov_get_rad_quality

!> Return plane_parallel flag
!! @param[out]    err             return status
!! @param[in]     inst_id         ID for instrument
!! @param[out]    plane_parallel  value of flag (non-zero => strict plane_parallel)
SUBROUTINE rttov_get_plane_parallel(err, inst_id, plane_parallel)
#include "throw.h"
  USE rttov_wrapper_handle
  USE parkind1, ONLY : jpim

  IMPLICIT NONE

  INTEGER(jpim), INTENT(OUT)   :: err
  INTEGER(jpim), INTENT(IN)    :: inst_id
  INTEGER(jpim), INTENT(OUT)   :: plane_parallel

  TYPE(rttovwrapperhandle_type), POINTER :: rth

#include "rttov_errorreport.interface"

  TRY
  CALL rttov_wrapper_handle_get_inst(err, inst_id, rth)
  THROW(err.NE.0)

  IF (.NOT. ASSOCIATED(rth%radiance%total)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Wrapper radiance not allocated: has RTTOV been called?')
  ENDIF

  IF (rth%radiance%plane_parallel) THEN
    plane_parallel = 1
  ELSE
    plane_parallel = 0
  ENDIF
  CATCH
END SUBROUTINE rttov_get_plane_parallel

!> Return level geometric heights
!! @param[out]    err               return status
!! @param[in]     inst_id           ID for instrument
!! @param[in,out] geometric_height  geometric heights (nlevels,nchanprof)
!! @param[in]     nchanprof         size of geometric_height(1,:) (not required in Python)
!! @param[in]     nlevels           size of geometric_height(:,1) (not required in Python)
SUBROUTINE rttov_get_geometric_height(err, inst_id, geometric_height, nchanprof, nlevels)
!f2py intent(hide):: nchanprof=len(geometric_height)
!f2py intent(hide):: nlevels=len(geometric_height[0])
#include "throw.h"
  USE rttov_wrapper_handle
  USE parkind1, ONLY : jpim, jprb

  IMPLICIT NONE

  INTEGER(jpim), INTENT(OUT)   :: err
  INTEGER(jpim), INTENT(IN)    :: inst_id
  INTEGER(jpim), INTENT(IN)    :: nchanprof
  INTEGER(jpim), INTENT(IN)    :: nlevels
  REAL(jprb),    INTENT(INOUT) :: geometric_height(nlevels,nchanprof)

  TYPE(rttovwrapperhandle_type), POINTER :: rth

#include "rttov_errorreport.interface"

  TRY
  CALL rttov_wrapper_handle_get_inst(err, inst_id, rth)
  THROW(err.NE.0)

  IF (.NOT. ASSOCIATED(rth%radiance%geometric_height)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Wrapper radiance not allocated: has RTTOV been called?')
  ENDIF

  geometric_height(:,:) = rth%radiance%geometric_height(:,:)

  CATCH
END SUBROUTINE rttov_get_geometric_height


!> Return upclear radiances
!! @param[out]    err           return status
!! @param[in]     inst_id       ID for instrument
!! @param[in,out] rad2_upclear  upclear radiances (nchanprof)
!! @param[in]     nchanprof     size of rad2_upclear (not required in Python)
SUBROUTINE rttov_get_rad2_upclear(err, inst_id, rad2_upclear, nchanprof)
!f2py intent(hide):: nchanprof=len(rad2_upclear)
#include "throw.h"
  USE rttov_wrapper_handle
  USE parkind1, ONLY : jpim, jprb

  IMPLICIT NONE

  INTEGER(jpim), INTENT(OUT)   :: err
  INTEGER(jpim), INTENT(IN)    :: inst_id
  INTEGER(jpim), INTENT(IN)    :: nchanprof
  REAL(jprb),    INTENT(INOUT) :: rad2_upclear(nchanprof)

  TYPE(rttovwrapperhandle_type), POINTER :: rth

#include "rttov_errorreport.interface"

  TRY
  CALL rttov_wrapper_handle_get_inst(err, inst_id, rth)
  THROW(err.NE.0)

  IF (.NOT. ASSOCIATED(rth%radiance2%upclear)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Wrapper radiance2 not allocated: has RTTOV been called?')
  ENDIF

  rad2_upclear(:) = rth%radiance2%upclear(:)
  CATCH
END SUBROUTINE rttov_get_rad2_upclear

!> Return dnclear radiances
!! @param[out]    err           return status
!! @param[in]     inst_id       ID for instrument
!! @param[in,out] rad2_dnclear  dnclear radiances (nchanprof)
!! @param[in]     nchanprof     size of rad2_dnclear (not required in Python)
SUBROUTINE rttov_get_rad2_dnclear(err, inst_id, rad2_dnclear, nchanprof)
!f2py intent(hide):: nchanprof=len(rad2_dnclear)
#include "throw.h"
  USE rttov_wrapper_handle
  USE parkind1, ONLY : jpim, jprb

  IMPLICIT NONE

  INTEGER(jpim), INTENT(OUT)   :: err
  INTEGER(jpim), INTENT(IN)    :: inst_id
  INTEGER(jpim), INTENT(IN)    :: nchanprof
  REAL(jprb),    INTENT(INOUT) :: rad2_dnclear(nchanprof)

  TYPE(rttovwrapperhandle_type), POINTER :: rth

#include "rttov_errorreport.interface"

  TRY
  CALL rttov_wrapper_handle_get_inst(err, inst_id, rth)
  THROW(err.NE.0)

  IF (.NOT. ASSOCIATED(rth%radiance2%dnclear)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Wrapper radiance2 not allocated: has RTTOV been called?')
  ENDIF

  rad2_dnclear(:) = rth%radiance2%dnclear(:)
  CATCH
END SUBROUTINE rttov_get_rad2_dnclear

!> Return refldnclear radiances
!! @param[out]    err               return status
!! @param[in]     inst_id           ID for instrument
!! @param[in,out] rad2_refldnclear  refldnclear radiances (nchanprof)
!! @param[in]     nchanprof         size of rad2_refldnclear (not required in Python)
SUBROUTINE rttov_get_rad2_refldnclear(err, inst_id, rad2_refldnclear, nchanprof)
!f2py intent(hide):: nchanprof=len(rad2_refldnclear)
#include "throw.h"
  USE rttov_wrapper_handle
  USE parkind1, ONLY : jpim, jprb

  IMPLICIT NONE

  INTEGER(jpim), INTENT(OUT)   :: err
  INTEGER(jpim), INTENT(IN)    :: inst_id
  INTEGER(jpim), INTENT(IN)    :: nchanprof
  REAL(jprb),    INTENT(INOUT) :: rad2_refldnclear(nchanprof)

  TYPE(rttovwrapperhandle_type), POINTER :: rth

#include "rttov_errorreport.interface"

  TRY
  CALL rttov_wrapper_handle_get_inst(err, inst_id, rth)
  THROW(err.NE.0)

  IF (.NOT. ASSOCIATED(rth%radiance2%refldnclear)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Wrapper radiance2 not allocated: has RTTOV been called?')
  ENDIF

  rad2_refldnclear(:) = rth%radiance2%refldnclear(:)
  CATCH
END SUBROUTINE rttov_get_rad2_refldnclear

!> Return upwelling radiances
!! @param[out]    err         return status
!! @param[in]     inst_id     ID for instrument
!! @param[in,out] rad2_up     upwelling radiances (nlayers,nchanprof)
!! @param[in]     nchanprof   size of rad2_up(1,:) (not required in Python)
!! @param[in]     nlayers     size of rad2_up(:,1) (not required in Python)
SUBROUTINE rttov_get_rad2_up(err, inst_id, rad2_up, nchanprof, nlayers)
!f2py intent(hide):: nchanprof=len(rad2_up)
!f2py intent(hide):: nlayers=len(rad2_up[0])
#include "throw.h"
  USE rttov_wrapper_handle
  USE parkind1, ONLY : jpim, jprb

  IMPLICIT NONE

  INTEGER(jpim), INTENT(OUT)   :: err
  INTEGER(jpim), INTENT(IN)    :: inst_id
  INTEGER(jpim), INTENT(IN)    :: nchanprof
  INTEGER(jpim), INTENT(IN)    :: nlayers
  REAL(jprb),    INTENT(INOUT) :: rad2_up(nlayers,nchanprof)

  TYPE(rttovwrapperhandle_type), POINTER :: rth

#include "rttov_errorreport.interface"

  TRY
  CALL rttov_wrapper_handle_get_inst(err, inst_id, rth)
  THROW(err.NE.0)

  IF (.NOT. ASSOCIATED(rth%radiance2%up)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Wrapper radiance2 not allocated: has RTTOV been called?')
  ENDIF

  rad2_up(:,:) = rth%radiance2%up(:,:)
  CATCH
END SUBROUTINE rttov_get_rad2_up

!> Return downwelling radiances
!! @param[out]    err         return status
!! @param[in]     inst_id     ID for instrument
!! @param[in,out] rad2_down   downwelling radiances (nlayers,nchanprof)
!! @param[in]     nchanprof   size of rad2_down(1,:) (not required in Python)
!! @param[in]     nlayers     size of rad2_down(:,1) (not required in Python)
SUBROUTINE rttov_get_rad2_down(err, inst_id, rad2_down, nchanprof, nlayers)
!f2py intent(hide):: nchanprof=len(rad2_down)
!f2py intent(hide):: nlayers=len(rad2_down[0])
#include "throw.h"
  USE rttov_wrapper_handle
  USE parkind1, ONLY : jpim, jprb

  IMPLICIT NONE

  INTEGER(jpim), INTENT(OUT)   :: err
  INTEGER(jpim), INTENT(IN)    :: inst_id
  INTEGER(jpim), INTENT(IN)    :: nchanprof
  INTEGER(jpim), INTENT(IN)    :: nlayers
  REAL(jprb),    INTENT(INOUT) :: rad2_down(nlayers,nchanprof)

  TYPE(rttovwrapperhandle_type), POINTER :: rth

#include "rttov_errorreport.interface"

  TRY
  CALL rttov_wrapper_handle_get_inst(err, inst_id, rth)
  THROW(err.NE.0)

  IF (.NOT. ASSOCIATED(rth%radiance2%down)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Wrapper radiance2 not allocated: has RTTOV been called?')
  ENDIF

  rad2_down(:,:) = rth%radiance2%down(:,:)
  CATCH
END SUBROUTINE rttov_get_rad2_down

!> Return Planck radiance at each level
!! @param[out]    err         return status
!! @param[in]     inst_id     ID for instrument
!! @param[in,out] rad2_surf   Planck radiances (nlayers,nchanprof)
!! @param[in]     nchanprof   size of rad2_surf(1,:) (not required in Python)
!! @param[in]     nlayers     size of rad2_surf(:,1) (not required in Python)
SUBROUTINE rttov_get_rad2_surf(err, inst_id, rad2_surf, nchanprof, nlayers)
!f2py intent(hide):: nchanprof=len(rad2_surf)
!f2py intent(hide):: nlayers=len(rad2_surf[0])
#include "throw.h"
  USE rttov_wrapper_handle
  USE parkind1, ONLY : jpim, jprb

  IMPLICIT NONE

  INTEGER(jpim), INTENT(OUT)   :: err
  INTEGER(jpim), INTENT(IN)    :: inst_id
  INTEGER(jpim), INTENT(IN)    :: nchanprof
  INTEGER(jpim), INTENT(IN)    :: nlayers
  REAL(jprb),    INTENT(INOUT) :: rad2_surf(nlayers,nchanprof)

  TYPE(rttovwrapperhandle_type), POINTER :: rth

#include "rttov_errorreport.interface"

  TRY
  CALL rttov_wrapper_handle_get_inst(err, inst_id, rth)
  THROW(err.NE.0)

  IF (.NOT. ASSOCIATED(rth%radiance2%surf)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Wrapper radiance2 not allocated: has RTTOV been called?')
  ENDIF

  rad2_surf(:,:) = rth%radiance2%surf(:,:)
  CATCH
END SUBROUTINE rttov_get_rad2_surf


!> Return total transmittances for surface-satellite path for emitted radiation
!! @param[out]    err               return status
!! @param[in]     inst_id           ID for instrument
!! @param[in,out] tau_total         total transmittances (nchanprof)
!! @param[in]     nchanprof         size of tau_total (not required in Python)
SUBROUTINE rttov_get_tau_total(err, inst_id, tau_total, nchanprof)
!f2py intent(hide):: nchanprof=len(tau_total)
#include "throw.h"
  USE rttov_wrapper_handle
  USE parkind1, ONLY : jpim, jprb

  IMPLICIT NONE

  INTEGER(jpim), INTENT(OUT)   :: err
  INTEGER(jpim), INTENT(IN)    :: inst_id
  INTEGER(jpim), INTENT(IN)    :: nchanprof
  REAL(jprb),    INTENT(INOUT) :: tau_total(nchanprof)

  TYPE(rttovwrapperhandle_type), POINTER :: rth

#include "rttov_errorreport.interface"

  TRY
  CALL rttov_wrapper_handle_get_inst(err, inst_id, rth)
  THROW(err.NE.0)

  IF (.NOT. ASSOCIATED(rth%transmission%tau_total)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Wrapper transmission not allocated: has RTTOV been called?')
  ENDIF

  tau_total(:) = rth%transmission%tau_total(:)
  CATCH
END SUBROUTINE rttov_get_tau_total

!> Return level-space transmittances for surface-satellite path for emitted radiation
!! @param[out]    err               return status
!! @param[in]     inst_id           ID for instrument
!! @param[in,out] tau_levels        transmittances (nlevels,nchanprof)
!! @param[in]     nchanprof         size of tau_levels(1,:) (not required in Python)
!! @param[in]     nlevels           size of tau_levels(:,1) (not required in Python)
SUBROUTINE rttov_get_tau_levels(err, inst_id, tau_levels, nchanprof, nlevels)
!f2py intent(hide):: nchanprof=len(tau_levels)
!f2py intent(hide):: nlevels=len(tau_levels[0])
#include "throw.h"
  USE rttov_wrapper_handle
  USE parkind1, ONLY : jpim, jprb

  IMPLICIT NONE

  INTEGER(jpim), INTENT(OUT)   :: err
  INTEGER(jpim), INTENT(IN)    :: inst_id
  INTEGER(jpim), INTENT(IN)    :: nchanprof
  INTEGER(jpim), INTENT(IN)    :: nlevels
  REAL(jprb),    INTENT(INOUT) :: tau_levels(nlevels,nchanprof)

  TYPE(rttovwrapperhandle_type), POINTER :: rth

#include "rttov_errorreport.interface"

  TRY
  CALL rttov_wrapper_handle_get_inst(err, inst_id, rth)
  THROW(err.NE.0)

  IF (.NOT. ASSOCIATED(rth%transmission%tau_levels)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Wrapper transmission not allocated: has RTTOV been called?')
  ENDIF

  tau_levels(:,:) = rth%transmission%tau_levels(:,:)
  CATCH
END SUBROUTINE rttov_get_tau_levels

!> Return total transmittances for sun-surface-satellite path for solar radiation
!! @param[out]    err                 return status
!! @param[in]     inst_id             ID for instrument
!! @param[in,out] tausun_total_path2  total transmittances (nchanprof)
!! @param[in]     nchanprof           size of tausun_total_path2 (not required in Python)
SUBROUTINE rttov_get_tausun_total_path2(err, inst_id, tausun_total_path2, nchanprof)
!f2py intent(hide):: nchanprof=len(tausun_total_path2)
#include "throw.h"
  USE rttov_wrapper_handle
  USE parkind1, ONLY : jpim, jprb

  IMPLICIT NONE

  INTEGER(jpim), INTENT(OUT)   :: err
  INTEGER(jpim), INTENT(IN)    :: inst_id
  INTEGER(jpim), INTENT(IN)    :: nchanprof
  REAL(jprb),    INTENT(INOUT) :: tausun_total_path2(nchanprof)

  TYPE(rttovwrapperhandle_type), POINTER :: rth

#include "rttov_errorreport.interface"

  TRY
  CALL rttov_wrapper_handle_get_inst(err, inst_id, rth)
  THROW(err.NE.0)

  IF (.NOT. ASSOCIATED(rth%transmission%tausun_total_path2)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Wrapper transmission not allocated: has RTTOV been called?')
  ENDIF

  tausun_total_path2(:) = rth%transmission%tausun_total_path2(:)
  CATCH
END SUBROUTINE rttov_get_tausun_total_path2

!> Return level-space transmittances for sun-surface-satellite path for solar radiation
!! @param[out]    err                   return status
!! @param[in]     inst_id               ID for instrument
!! @param[in,out] tausun_levels_path2   transmittances (nlevels,nchanprof)
!! @param[in]     nchanprof             size of tausun_levels_path2(1,:) (not required in Python)
!! @param[in]     nlevels               size of tausun_levels_path2(:,1) (not required in Python)
SUBROUTINE rttov_get_tausun_levels_path2(err, inst_id, tausun_levels_path2, nchanprof, nlevels)
!f2py intent(hide):: nchanprof=len(tausun_levels_path2)
!f2py intent(hide):: nlevels=len(tausun_levels_path2[0])
#include "throw.h"
  USE rttov_wrapper_handle
  USE parkind1, ONLY : jpim, jprb

  IMPLICIT NONE

  INTEGER(jpim), INTENT(OUT)   :: err
  INTEGER(jpim), INTENT(IN)    :: inst_id
  INTEGER(jpim), INTENT(IN)    :: nchanprof
  INTEGER(jpim), INTENT(IN)    :: nlevels
  REAL(jprb),    INTENT(INOUT) :: tausun_levels_path2(nlevels,nchanprof)

  TYPE(rttovwrapperhandle_type), POINTER :: rth

#include "rttov_errorreport.interface"

  TRY
  CALL rttov_wrapper_handle_get_inst(err, inst_id, rth)
  THROW(err.NE.0)

  IF (.NOT. ASSOCIATED(rth%transmission%tausun_levels_path2)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Wrapper transmission not allocated: has RTTOV been called?')
  ENDIF

  tausun_levels_path2(:,:) = rth%transmission%tausun_levels_path2(:,:)
  CATCH
END SUBROUTINE rttov_get_tausun_levels_path2

!> Return total transmittances for surface-satellite path for solar radiation
!! @param[out]    err                 return status
!! @param[in]     inst_id             ID for instrument
!! @param[in,out] tausun_total_path1  total transmittances (nchanprof)
!! @param[in]     nchanprof           size of tausun_total_path1 (not required in Python)
SUBROUTINE rttov_get_tausun_total_path1(err, inst_id, tausun_total_path1, nchanprof)
!f2py intent(hide):: nchanprof=len(tausun_total_path1)
#include "throw.h"
  USE rttov_wrapper_handle
  USE parkind1, ONLY : jpim, jprb

  IMPLICIT NONE

  INTEGER(jpim), INTENT(OUT)   :: err
  INTEGER(jpim), INTENT(IN)    :: inst_id
  INTEGER(jpim), INTENT(IN)    :: nchanprof
  REAL(jprb),    INTENT(INOUT) :: tausun_total_path1(nchanprof)

  TYPE(rttovwrapperhandle_type), POINTER :: rth

#include "rttov_errorreport.interface"

  TRY
  CALL rttov_wrapper_handle_get_inst(err, inst_id, rth)
  THROW(err.NE.0)

  IF (.NOT. ASSOCIATED(rth%transmission%tausun_total_path1)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Wrapper transmission not allocated: has RTTOV been called?')
  ENDIF

  tausun_total_path1(:) = rth%transmission%tausun_total_path1(:)
  CATCH
END SUBROUTINE rttov_get_tausun_total_path1

!> Return level-space transmittances for surface-satellite path for solar radiation
!! @param[out]    err                   return status
!! @param[in]     inst_id               ID for instrument
!! @param[in,out] tausun_levels_path1   transmittances (nlevels,nchanprof)
!! @param[in]     nchanprof             size of tausun_levels_path1(1,:) (not required in Python)
!! @param[in]     nlevels               size of tausun_levels_path1(:,1) (not required in Python)
SUBROUTINE rttov_get_tausun_levels_path1(err, inst_id, tausun_levels_path1, nchanprof, nlevels)
!f2py intent(hide):: nchanprof=len(tausun_levels_path1)
!f2py intent(hide):: nlevels=len(tausun_levels_path1[0])
#include "throw.h"
  USE rttov_wrapper_handle
  USE parkind1, ONLY : jpim, jprb

  IMPLICIT NONE

  INTEGER(jpim), INTENT(OUT)   :: err
  INTEGER(jpim), INTENT(IN)    :: inst_id
  INTEGER(jpim), INTENT(IN)    :: nchanprof
  INTEGER(jpim), INTENT(IN)    :: nlevels
  REAL(jprb),    INTENT(INOUT) :: tausun_levels_path1(nlevels,nchanprof)

  TYPE(rttovwrapperhandle_type), POINTER :: rth

#include "rttov_errorreport.interface"

  TRY
  CALL rttov_wrapper_handle_get_inst(err, inst_id, rth)
  THROW(err.NE.0)

  IF (.NOT. ASSOCIATED(rth%transmission%tausun_levels_path1)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Wrapper transmission not allocated: has RTTOV been called?')
  ENDIF

  tausun_levels_path1(:,:) = rth%transmission%tausun_levels_path1(:,:)
  CATCH
END SUBROUTINE rttov_get_tausun_levels_path1

!> Return total cloud transmittances for surface-satellite path
!! @param[out]    err               return status
!! @param[in]     inst_id           ID for instrument
!! @param[in,out] tau_total_cld     total cloud transmittances (nchanprof)
!! @param[in]     nchanprof         size of tau_total_cld (not required in Python)
SUBROUTINE rttov_get_tau_total_cld(err, inst_id, tau_total_cld, nchanprof)
!f2py intent(hide):: nchanprof=len(tau_total_cld)
#include "throw.h"
  USE rttov_wrapper_handle
  USE parkind1, ONLY : jpim, jprb

  IMPLICIT NONE

  INTEGER(jpim), INTENT(OUT)   :: err
  INTEGER(jpim), INTENT(IN)    :: inst_id
  INTEGER(jpim), INTENT(IN)    :: nchanprof
  REAL(jprb),    INTENT(INOUT) :: tau_total_cld(nchanprof)

  TYPE(rttovwrapperhandle_type), POINTER :: rth

#include "rttov_errorreport.interface"

  TRY
  CALL rttov_wrapper_handle_get_inst(err, inst_id, rth)
  THROW(err.NE.0)

  IF (.NOT. ASSOCIATED(rth%transmission%tau_total_cld)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Wrapper transmission not allocated: has RTTOV been called?')
  ENDIF

  tau_total_cld(:) = rth%transmission%tau_total_cld(:)
  CATCH
END SUBROUTINE rttov_get_tau_total_cld

!> Return level-space cloud transmittances for surface-satellite path
!! @param[out]    err               return status
!! @param[in]     inst_id           ID for instrument
!! @param[in,out] tau_levels_cld    cloud transmittances (nlevels,nchanprof)
!! @param[in]     nchanprof         size of tau_levels_cld(1,:) (not required in Python)
!! @param[in]     nlevels           size of tau_levels_cld(:,1) (not required in Python)
SUBROUTINE rttov_get_tau_levels_cld(err, inst_id, tau_levels_cld, nchanprof, nlevels)
!f2py intent(hide):: nchanprof=len(tau_levels_cld)
!f2py intent(hide):: nlevels=len(tau_levels_cld[0])
#include "throw.h"
  USE rttov_wrapper_handle
  USE parkind1, ONLY : jpim, jprb

  IMPLICIT NONE

  INTEGER(jpim), INTENT(OUT)   :: err
  INTEGER(jpim), INTENT(IN)    :: inst_id
  INTEGER(jpim), INTENT(IN)    :: nchanprof
  INTEGER(jpim), INTENT(IN)    :: nlevels
  REAL(jprb),    INTENT(INOUT) :: tau_levels_cld(nlevels,nchanprof)

  TYPE(rttovwrapperhandle_type), POINTER :: rth

#include "rttov_errorreport.interface"

  TRY
  CALL rttov_wrapper_handle_get_inst(err, inst_id, rth)
  THROW(err.NE.0)

  IF (.NOT. ASSOCIATED(rth%transmission%tau_levels_cld)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Wrapper transmission not allocated: has RTTOV been called?')
  ENDIF

  tau_levels_cld(:,:) = rth%transmission%tau_levels_cld(:,:)
  CATCH
END SUBROUTINE rttov_get_tau_levels_cld

!> Return RTTOV-SCATT emissivity retrieval cfrac output
!! @param[out]    err         return status
!! @param[in]     inst_id     ID for instrument
!! @param[in,out] cfrac       cloud fraction (nchanprof)
!! @param[in]     nchanprof   size of cfrac (not required in Python)
SUBROUTINE rttov_get_emis_terms_cfrac(err, inst_id, cfrac, nchanprof)
!f2py intent(hide):: nchanprof=len(cfrac)
#include "throw.h"
  USE rttov_wrapper_handle
  USE parkind1, ONLY : jpim, jprb

  IMPLICIT NONE

  INTEGER(jpim), INTENT(OUT)   :: err
  INTEGER(jpim), INTENT(IN)    :: inst_id
  INTEGER(jpim), INTENT(IN)    :: nchanprof
  REAL(jprb),    INTENT(INOUT) :: cfrac(nchanprof)

  TYPE(rttovwrapperhandle_type), POINTER :: rth

#include "rttov_errorreport.interface"

  TRY
  CALL rttov_wrapper_handle_get_inst(err, inst_id, rth)
  THROW(err.NE.0)

  IF (.NOT. ASSOCIATED(rth%emis_terms%cfrac)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Wrapper emis_terms not allocated: has RTTOV been called?')
  ENDIF

  cfrac(:) = rth%emis_terms%cfrac(:)
  CATCH
END SUBROUTINE rttov_get_emis_terms_cfrac

!> Return RTTOV-SCATT emissivity retrieval bsfc output
!! @param[out]    err         return status
!! @param[in]     inst_id     ID for instrument
!! @param[in,out] bsfc        surface black-body radiance (nchanprof)
!! @param[in]     nchanprof   size of bsfc (not required in Python)
SUBROUTINE rttov_get_emis_terms_bsfc(err, inst_id, bsfc, nchanprof)
!f2py intent(hide):: nchanprof=len(bsfc)
#include "throw.h"
  USE rttov_wrapper_handle
  USE parkind1, ONLY : jpim, jprb

  IMPLICIT NONE

  INTEGER(jpim), INTENT(OUT)   :: err
  INTEGER(jpim), INTENT(IN)    :: inst_id
  INTEGER(jpim), INTENT(IN)    :: nchanprof
  REAL(jprb),    INTENT(INOUT) :: bsfc(nchanprof)

  TYPE(rttovwrapperhandle_type), POINTER :: rth

#include "rttov_errorreport.interface"

  TRY
  CALL rttov_wrapper_handle_get_inst(err, inst_id, rth)
  THROW(err.NE.0)

  IF (.NOT. ASSOCIATED(rth%emis_terms%bsfc)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Wrapper emis_terms not allocated: has RTTOV been called?')
  ENDIF

  bsfc(:) = rth%emis_terms%bsfc(:)
  CATCH
END SUBROUTINE rttov_get_emis_terms_bsfc

!> Return RTTOV-SCATT emissivity retrieval tau_cld output
!! @param[out]    err         return status
!! @param[in]     inst_id     ID for instrument
!! @param[in,out] tau_cld     cloudy transmittance (nchanprof)
!! @param[in]     nchanprof   size of tau_cld (not required in Python)
SUBROUTINE rttov_get_emis_terms_tau_cld(err, inst_id, tau_cld, nchanprof)
!f2py intent(hide):: nchanprof=len(tau_cld)
#include "throw.h"
  USE rttov_wrapper_handle
  USE parkind1, ONLY : jpim, jprb

  IMPLICIT NONE

  INTEGER(jpim), INTENT(OUT)   :: err
  INTEGER(jpim), INTENT(IN)    :: inst_id
  INTEGER(jpim), INTENT(IN)    :: nchanprof
  REAL(jprb),    INTENT(INOUT) :: tau_cld(nchanprof)

  TYPE(rttovwrapperhandle_type), POINTER :: rth

#include "rttov_errorreport.interface"

  TRY
  CALL rttov_wrapper_handle_get_inst(err, inst_id, rth)
  THROW(err.NE.0)

  IF (.NOT. ASSOCIATED(rth%emis_terms%tau_cld)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Wrapper emis_terms not allocated: has RTTOV been called?')
  ENDIF

  tau_cld(:) = rth%emis_terms%tau_cld(:)
  CATCH
END SUBROUTINE rttov_get_emis_terms_tau_cld

!> Return RTTOV-SCATT emissivity retrieval up_cld output
!! @param[out]    err         return status
!! @param[in]     inst_id     ID for instrument
!! @param[in,out] up_cld      cloudy upwelling radiance (nchanprof)
!! @param[in]     nchanprof   size of up_cld (not required in Python)
SUBROUTINE rttov_get_emis_terms_up_cld(err, inst_id, up_cld, nchanprof)
!f2py intent(hide):: nchanprof=len(up_cld)
#include "throw.h"
  USE rttov_wrapper_handle
  USE parkind1, ONLY : jpim, jprb

  IMPLICIT NONE

  INTEGER(jpim), INTENT(OUT)   :: err
  INTEGER(jpim), INTENT(IN)    :: inst_id
  INTEGER(jpim), INTENT(IN)    :: nchanprof
  REAL(jprb),    INTENT(INOUT) :: up_cld(nchanprof)

  TYPE(rttovwrapperhandle_type), POINTER :: rth

#include "rttov_errorreport.interface"

  TRY
  CALL rttov_wrapper_handle_get_inst(err, inst_id, rth)
  THROW(err.NE.0)

  IF (.NOT. ASSOCIATED(rth%emis_terms%up_cld)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Wrapper emis_terms not allocated: has RTTOV been called?')
  ENDIF

  up_cld(:) = rth%emis_terms%up_cld(:)
  CATCH
END SUBROUTINE rttov_get_emis_terms_up_cld

!> Return RTTOV-SCATT emissivity retrieval down_cld output
!! @param[out]    err         return status
!! @param[in]     inst_id     ID for instrument
!! @param[in,out] down_cld    cloudy downwelling radiance (nchanprof)
!! @param[in]     nchanprof   size of down_cld (not required in Python)
SUBROUTINE rttov_get_emis_terms_down_cld(err, inst_id, down_cld, nchanprof)
!f2py intent(hide):: nchanprof=len(down_cld)
#include "throw.h"
  USE rttov_wrapper_handle
  USE parkind1, ONLY : jpim, jprb

  IMPLICIT NONE

  INTEGER(jpim), INTENT(OUT)   :: err
  INTEGER(jpim), INTENT(IN)    :: inst_id
  INTEGER(jpim), INTENT(IN)    :: nchanprof
  REAL(jprb),    INTENT(INOUT) :: down_cld(nchanprof)

  TYPE(rttovwrapperhandle_type), POINTER :: rth

#include "rttov_errorreport.interface"

  TRY
  CALL rttov_wrapper_handle_get_inst(err, inst_id, rth)
  THROW(err.NE.0)

  IF (.NOT. ASSOCIATED(rth%emis_terms%down_cld)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Wrapper emis_terms not allocated: has RTTOV been called?')
  ENDIF

  down_cld(:) = rth%emis_terms%down_cld(:)
  CATCH
END SUBROUTINE rttov_get_emis_terms_down_cld

!> Return RTTOV-SCATT emissivity retrieval tau_clr output
!! @param[out]    err         return status
!! @param[in]     inst_id     ID for instrument
!! @param[in,out] tau_clr     clear transmittance (nchanprof)
!! @param[in]     nchanprof   size of tau_clr (not required in Python)
SUBROUTINE rttov_get_emis_terms_tau_clr(err, inst_id, tau_clr, nchanprof)
!f2py intent(hide):: nchanprof=len(tau_clr)
#include "throw.h"
  USE rttov_wrapper_handle
  USE parkind1, ONLY : jpim, jprb

  IMPLICIT NONE

  INTEGER(jpim), INTENT(OUT)   :: err
  INTEGER(jpim), INTENT(IN)    :: inst_id
  INTEGER(jpim), INTENT(IN)    :: nchanprof
  REAL(jprb),    INTENT(INOUT) :: tau_clr(nchanprof)

  TYPE(rttovwrapperhandle_type), POINTER :: rth

#include "rttov_errorreport.interface"

  TRY
  CALL rttov_wrapper_handle_get_inst(err, inst_id, rth)
  THROW(err.NE.0)

  IF (.NOT. ASSOCIATED(rth%emis_terms%tau_clr)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Wrapper emis_terms not allocated: has RTTOV been called?')
  ENDIF

  tau_clr(:) = rth%emis_terms%tau_clr(:)
  CATCH
END SUBROUTINE rttov_get_emis_terms_tau_clr

!> Return RTTOV-SCATT emissivity retrieval up_clr output
!! @param[out]    err         return status
!! @param[in]     inst_id     ID for instrument
!! @param[in,out] up_clr      clear upwelling radiance (nchanprof)
!! @param[in]     nchanprof   size of up_clr (not required in Python)
SUBROUTINE rttov_get_emis_terms_up_clr(err, inst_id, up_clr, nchanprof)
!f2py intent(hide):: nchanprof=len(up_clr)
#include "throw.h"
  USE rttov_wrapper_handle
  USE parkind1, ONLY : jpim, jprb

  IMPLICIT NONE

  INTEGER(jpim), INTENT(OUT)   :: err
  INTEGER(jpim), INTENT(IN)    :: inst_id
  INTEGER(jpim), INTENT(IN)    :: nchanprof
  REAL(jprb),    INTENT(INOUT) :: up_clr(nchanprof)

  TYPE(rttovwrapperhandle_type), POINTER :: rth

#include "rttov_errorreport.interface"

  TRY
  CALL rttov_wrapper_handle_get_inst(err, inst_id, rth)
  THROW(err.NE.0)

  IF (.NOT. ASSOCIATED(rth%emis_terms%up_clr)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Wrapper emis_terms not allocated: has RTTOV been called?')
  ENDIF

  up_clr(:) = rth%emis_terms%up_clr(:)
  CATCH
END SUBROUTINE rttov_get_emis_terms_up_clr

!> Return RTTOV-SCATT emissivity retrieval down_clr output
!! @param[out]    err         return status
!! @param[in]     inst_id     ID for instrument
!! @param[in,out] down_clr    clear downwelling radiance (nchanprof)
!! @param[in]     nchanprof   size of down_clr (not required in Python)
SUBROUTINE rttov_get_emis_terms_down_clr(err, inst_id, down_clr, nchanprof)
!f2py intent(hide):: nchanprof=len(down_clr)
#include "throw.h"
  USE rttov_wrapper_handle
  USE parkind1, ONLY : jpim, jprb

  IMPLICIT NONE

  INTEGER(jpim), INTENT(OUT)   :: err
  INTEGER(jpim), INTENT(IN)    :: inst_id
  INTEGER(jpim), INTENT(IN)    :: nchanprof
  REAL(jprb),    INTENT(INOUT) :: down_clr(nchanprof)

  TYPE(rttovwrapperhandle_type), POINTER :: rth

#include "rttov_errorreport.interface"

  TRY
  CALL rttov_wrapper_handle_get_inst(err, inst_id, rth)
  THROW(err.NE.0)

  IF (.NOT. ASSOCIATED(rth%emis_terms%down_clr)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Wrapper emis_terms not allocated: has RTTOV been called?')
  ENDIF

  down_clr(:) = rth%emis_terms%down_clr(:)
  CATCH
END SUBROUTINE rttov_get_emis_terms_down_clr


!> Return RTTOV-SCATT reflectivities output
!! @param[out]    err         return status
!! @param[in]     inst_id     ID for instrument
!! @param[in,out] zef         reflectivities (nlevels, nchanprof)
!! @param[in]     nchanprof   size of zef (not required in Python)
!! @param[in]     nlevels     size of zef (not required in Python)
SUBROUTINE rttov_get_zef(err, inst_id, zef, nchanprof, nlevels)
!f2py intent(hide):: nchanprof=nchanprof
!f2py intent(hide):: nlevels=nlevels
#include "throw.h"
  USE rttov_wrapper_handle
  USE parkind1, ONLY : jpim, jprb

  IMPLICIT NONE

  INTEGER(jpim), INTENT(OUT)   :: err
  INTEGER(jpim), INTENT(IN)    :: inst_id
  INTEGER(jpim), INTENT(IN)    :: nchanprof
  INTEGER(jpim), INTENT(IN)    :: nlevels
  REAL(jprb),    INTENT(INOUT) :: zef(nlevels, nchanprof)

  TYPE(rttovwrapperhandle_type), POINTER :: rth

#include "rttov_errorreport.interface"

  TRY
  CALL rttov_wrapper_handle_get_inst(err, inst_id, rth)
  THROW(err.NE.0)

  IF (.NOT. ASSOCIATED(rth%reflectivity%zef)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Wrapper reflectivity not allocated: has RTTOV been called?')
  ENDIF

  zef(:,:) = rth%reflectivity%zef(:,:)
  CATCH
END SUBROUTINE rttov_get_zef

!> Return RTTOV-SCATT attenuated reflectivities output
!! @param[out]    err         return status
!! @param[in]     inst_id     ID for instrument
!! @param[in,out] azef        attenuated reflectivities (nlevels, nchanprof)
!! @param[in]     nchanprof   size of azef (not required in Python)
!! @param[in]     nlevels     size of azef (not required in Python)
SUBROUTINE rttov_get_azef(err, inst_id, azef, nchanprof, nlevels)
!f2py intent(hide):: nchanprof=nchanprof
!f2py intent(hide):: nlevels=nlevels
#include "throw.h"
  USE rttov_wrapper_handle
  USE parkind1, ONLY : jpim, jprb

  IMPLICIT NONE

  INTEGER(jpim), INTENT(OUT)   :: err
  INTEGER(jpim), INTENT(IN)    :: inst_id
  INTEGER(jpim), INTENT(IN)    :: nchanprof
  INTEGER(jpim), INTENT(IN)    :: nlevels
  REAL(jprb),    INTENT(INOUT) :: azef(nlevels, nchanprof)

  TYPE(rttovwrapperhandle_type), POINTER :: rth

#include "rttov_errorreport.interface"

  TRY
  CALL rttov_wrapper_handle_get_inst(err, inst_id, rth)
  THROW(err.NE.0)

  IF (.NOT. ASSOCIATED(rth%reflectivity%azef)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Wrapper azef not allocated: has RTTOV been called?')
  ENDIF

  azef(:,:) = rth%reflectivity%azef(:,:)
  CATCH
END SUBROUTINE rttov_get_azef

! The following subroutines return information from an RTTOV coefficients
! structure once an instrument has been initialised. These are not intended
! for general use.

!> Return integer variable from optical depth coefficient structure
!! @param[out]    err                   return status
!! @param[in]     inst_id               ID for instrument
!! @param[in]     varch                 name of variable to return
!! @param[out]    i0                    returned variable value
SUBROUTINE RTTOV_GET_COEF_VAL_I0(err, inst_id, varch, I0)
!
#include "throw.h"

USE parkind1, ONLY : jpim
!INTF_OFF
USE rttov_wrapper_handle
USE rttov_const, ONLY : ngases_max
!INTF_ON
  IMPLICIT NONE
  INTEGER(jpim),      INTENT(OUT) :: err
  INTEGER(jpim),      INTENT(IN)  :: inst_id
  CHARACTER(LEN=*),   INTENT(IN)  :: varch
  INTEGER(jpim),      INTENT(OUT) :: I0
!INTF_END

  TYPE(rttovwrapperhandle_type), POINTER :: rth

#include "rttov_errorreport.interface"
!
TRY

  CALL rttov_wrapper_handle_get_inst(err, inst_id, rth)
  THROW(err.NE.0)

  SELECT CASE (TRIM(VARCH))
    CASE ("ID_PLATFORM")
      I0 = rth%COEFS%COEF%ID_PLATFORM
    CASE ("ID_SAT")
      I0 = rth%COEFS%COEF%ID_SAT
    CASE ("ID_INST")
      I0 = rth%COEFS%COEF%ID_INST
    CASE ("ID_SENSOR")
      I0 = rth%COEFS%COEF%ID_SENSOR
    CASE ("ID_COMP_LVL")
      I0 = rth%COEFS%COEF%ID_COMP_LVL
    CASE ("ID_COMP_PC")
      I0 = rth%COEFS%COEF%ID_COMP_PC
    CASE ("FMV_MODEL_VER")
      I0 = rth%COEFS%COEF%FMV_MODEL_VER
    CASE ("FMV_CHN")
      I0 = rth%COEFS%COEF%FMV_CHN
    CASE ("FMV_GAS")
      I0 = rth%COEFS%COEF%FMV_GAS
    CASE ("NLEVELS")
      I0 = rth%COEFS%COEF%NLEVELS
    CASE ("NLAYERS")
      I0 = rth%COEFS%COEF%NLAYERS
    CASE ("NMIXED")
      I0 = rth%COEFS%COEF%NMIXED
    CASE ("NWATER")
      I0 = rth%COEFS%COEF%NWATER
    CASE ("NOZONE")
      I0 = rth%COEFS%COEF%NOZONE
    CASE ("NWVCONT")
      I0 = rth%COEFS%COEF%NWVCONT
    CASE ("NCO2")
      I0 = rth%COEFS%COEF%NCO2
    CASE ("NN2O")
      I0 = rth%COEFS%COEF%NN2O
    CASE ("NCO")
      I0 = rth%COEFS%COEF%NCO
    CASE ("NCH4")
      I0 = rth%COEFS%COEF%NCH4
    CASE ("NSO2")
      I0 = rth%COEFS%COEF%NSO2
    CASE ("INCZEEMAN")
      IF( rth%COEFS%COEF%INCZEEMAN ) THEN
        I0 = 1
      ELSE
        I0 = 0
      ENDIF
    CASE ("FMV_PC_BANDS")
      I0 = rth%COEFS%COEF_PCCOMP%FMV_PC_BANDS
    CASE ("FMV_PC_CLD")
      I0 = rth%COEFS%COEF_PCCOMP%FMV_PC_CLD
    CASE ("FMV_PC_MNUM")
      I0 = rth%COEFS%COEF_PCCOMP%FMV_PC_MNUM
    CASE ("SOLARCOEF")
      IF( rth%COEFS%COEF%SOLARCOEF ) THEN
        I0 = 1
      ELSE
        I0 = 0
      ENDIF
    CASE ("NLTECOEF")
      IF( rth%COEFS%COEF%NLTECOEF) THEN
        I0 = 1
      ELSE
        I0 = 0
      ENDIF
    CASE ("PMC_SHIFT")
      IF( rth%COEFS%COEF%PMC_SHIFT ) THEN
        I0 = 1
      ELSE
        I0 = 0
      ENDIF

! water optical constant
! wave spectrum
! fastem
! ssirem                 all omited

! aliases
    CASE ("NCHANNELS")
      I0 = rth%COEFS%COEF%FMV_CHN

! array sizes
    CASE ("SIZE_FMV_GAS_ID" ,&
        & "SIZE_FMV_VAR"    ,&
        & "SIZE_FMV_COE"    ,&
        & "SIZE_FMV_LVL"    )
      I0 = rth%COEFS%COEF%FMV_GAS

    CASE ("SIZE_FMV_GAS_POS")
      I0 = ngases_max

    CASE ("SIZE_FF_ORI_CHN" ,&
        & "SIZE_FF_VAL_CHN" ,&
        & "SIZE_FF_CWN" ,&
        & "SIZE_FF_BCO" ,&
        & "SIZE_FF_BCS" ,&
        & "SIZE_FF_GAM" ,&
        & "SIZE_TT_VAL_CHN" ,&
        & "SIZE_TT_A0"  ,&
        & "SIZE_TT_A1"  ,&
        & "SIZE_PW_VAL_CHN" ,&
        & "SIZE_SS_VAL_CHN" ,&
        & "SIZE_SS_SOLAR_SPECTRUM" )
      I0 = rth%COEFS%COEF%FMV_CHN

    CASE ("SIZE_REF_PRFL_P" ,&
        & "SIZE_LIM_PRFL_P" ,&
        & "SIZE_LIM_PRFL_TMAX" ,&
        & "SIZE_LIM_PRFL_TMIN" )
      I0 = rth%COEFS%COEF%NLEVELS

    CASE ("SIZE_NOISE_IN")
      I0=rth%COEFS%COEF_PCCOMP%FMV_PC_NCHN
    CASE ("SIZE_FMV_PC_SETS")
      I0 = rth%COEFS%COEF_PCCOMP%FMV_PC_BANDS

!    CASE ("SIZE_LINE_BY_LINE")
!      I0 = 20

! 2 dimensions arrays, should be specified in I0 and I1
    CASE ("SIZE_FMV_PC_NPRED")
      I0 = 2

! array size for aliases
    CASE ("SIZE_WAVENUMBERS")
      I0 = rth%COEFS%COEF%FMV_CHN

    CASE ("SIZE_REF_PRESSURE" ,&
        & "SIZE_REF_TEMPERATURE" ,&
        & "SIZE_REF_WATERVAPOR" ,&
        & "SIZE_REF_OZONE" ,&
        & "SIZE_REF_CO2" ,&
        & "SIZE_REF_N2O" ,&
        & "SIZE_REF_CO" ,&
        & "SIZE_REF_CH4" ,&
        & "SIZE_REF_SO2" ,&
        & "SIZE_MIN_TEMPERATURE" ,&
        & "SIZE_MIN_WATERVAPOR" ,&
        & "SIZE_MIN_OZONE" ,&
        & "SIZE_MIN_CO2" ,&
        & "SIZE_MIN_N2O" ,&
        & "SIZE_MIN_CO" ,&
        & "SIZE_MIN_CH4" ,&
        & "SIZE_MIN_SO2" ,&
        & "SIZE_MAX_TEMPERATURE" ,&
        & "SIZE_MAX_WATERVAPOR" ,&
        & "SIZE_MAX_OZONE" ,&
        & "SIZE_MAX_CO2" ,&
        & "SIZE_MAX_N2O" ,&
        & "SIZE_MAX_CO" ,&
        & "SIZE_MAX_CH4" ,&
        & "SIZE_MAX_SO2")
      I0 = rth%COEFS%COEF%NLEVELS

    CASE DEFAULT
      err = errorstatus_fatal
      THROWM(ERR.NE.0,"INTERGER SCALAR VARIABLE "//TRIM(VARCH)//" NOT FOUND")

  END SELECT

CATCH

END SUBROUTINE

!> Return real variable from optical depth coefficient structure
!! @param[out]    err                   return status
!! @param[in]     inst_id               ID for instrument
!! @param[in]     varch                 name of variable to return
!! @param[out]    r0                    returned variable value
SUBROUTINE RTTOV_GET_COEF_VAL_R0(err, inst_id, varch, R0)
!
#include "throw.h"

USE parkind1, ONLY : jpim, jprb
!INTF_OFF
USE rttov_wrapper_handle
USE rttov_const, ONLY : speedl
!INTF_ON
  IMPLICIT NONE
  INTEGER(jpim),      INTENT(OUT) :: err
  INTEGER(jpim),      INTENT(IN)  :: inst_id
  CHARACTER(LEN=*),   INTENT(IN)  :: varch
  REAL(jprb),         INTENT(OUT) :: R0
!INTF_END

  TYPE(rttovwrapperhandle_type), POINTER :: rth

#include "rttov_errorreport.interface"
!
TRY

  CALL rttov_wrapper_handle_get_inst(err, inst_id, rth)
  THROW(err.NE.0)

  SELECT CASE (TRIM(VARCH))

    CASE ("FC_SPEEDL")
      R0 = speedl
    CASE ("FC_PLANCK_C1")
      R0 = rth%COEFS%COEF%FC_PLANCK_C1
    CASE ("FC_PLANCK_C2")
      R0 = rth%COEFS%COEF%FC_PLANCK_C2
    CASE ("FC_SAT_HEIGHT")
      R0 = rth%COEFS%COEF%FC_SAT_HEIGHT

    CASE DEFAULT
      err = errorstatus_fatal
      THROWM(ERR.NE.0,"REAL SCALAR VARIABLE "//TRIM(VARCH)//" NOT FOUND")

  END SELECT

CATCH

END SUBROUTINE

!> Return character string variable from optical depth coefficient structure
!! @param[out]    err                   return status
!! @param[in]     inst_id               ID for instrument
!! @param[in]     varch                 name of variable to return
!! @param[out]    c0                    returned variable value
SUBROUTINE RTTOV_GET_COEF_VAL_C0(err, inst_id, varch, C0)
!
#include "throw.h"

USE parkind1, ONLY : jpim
!INTF_OFF
USE rttov_wrapper_handle
!INTF_ON
  IMPLICIT NONE
  INTEGER(jpim),      INTENT(OUT) :: err
  INTEGER(jpim),      INTENT(IN)  :: inst_id
  CHARACTER(LEN=*),   INTENT(IN)  :: varch
  CHARACTER(LEN=80),  INTENT(OUT) :: C0
!INTF_END

  TYPE(rttovwrapperhandle_type), POINTER :: rth

#include "rttov_errorreport.interface"
!
TRY

  CALL rttov_wrapper_handle_get_inst(err, inst_id, rth)
  THROW(err.NE.0)

  SELECT CASE (TRIM(VARCH))

    CASE ("ID_CREATION")
      C0 = rth%COEFS%COEF%ID_CREATION
    CASE ("ID_COMMON_NAME")
      C0 = rth%COEFS%COEF%ID_COMMON_NAME
    CASE ("FMV_MODEL_DEF")
      C0 = rth%COEFS%COEF%FMV_MODEL_DEF

    CASE DEFAULT
      err = errorstatus_fatal
      THROWM(ERR.NE.0,"CHARACTER STRING VARIABLE "//TRIM(VARCH)//" NOT FOUND")

  END SELECT

CATCH

END SUBROUTINE

!> Return 1D integer array variable from optical depth coefficient structure
!! @param[out]    err                   return status
!! @param[in]     inst_id               ID for instrument
!! @param[in]     varch                 name of variable to return
!! @param[in]     m                     size of i1
!! @param[out]    i1                    returned variable value
SUBROUTINE RTTOV_GET_COEF_VAL_I1(err, inst_id, varch, M, I1)
!
#include "throw.h"

USE parkind1, ONLY : jpim
!INTF_OFF
USE rttov_wrapper_handle
!INTF_ON
  IMPLICIT NONE
  INTEGER(jpim),      INTENT(OUT) :: err
  INTEGER(jpim),      INTENT(IN)  :: inst_id
  CHARACTER(LEN=*),   INTENT(IN)  :: varch
  INTEGER(jpim),      INTENT(IN)  :: M
  INTEGER(jpim),      INTENT(OUT) :: I1(M)
!INTF_END

  TYPE(rttovwrapperhandle_type), POINTER :: rth

#include "rttov_errorreport.interface"
!
TRY

  CALL rttov_wrapper_handle_get_inst(err, inst_id, rth)
  THROW(err.NE.0)

  SELECT CASE (TRIM(VARCH))
    CASE ("FMV_GAS_ID")
      I1 = rth%COEFS%COEF%FMV_GAS_ID
    CASE ("FMV_GAS_POS")
      I1 = rth%COEFS%COEF%FMV_GAS_POS
    CASE ("FMV_VAR")
      I1 = rth%COEFS%COEF%FMV_VAR
    CASE ("FMV_COE")
      I1 = rth%COEFS%COEF%FMV_COE
    CASE ("FMV_LVL")
      I1 = rth%COEFS%COEF%FMV_LVL
    CASE ("FF_ORI_CHN")
      I1 = rth%COEFS%COEF%FF_ORI_CHN
    CASE ("FF_VAL_CHN")
      I1 = rth%COEFS%COEF%FF_VAL_CHN
    CASE ("TT_VAL_CHN")
      I1 = rth%COEFS%COEF%TT_VAL_CHN
    CASE ("PW_VAL_CHN")
      I1 = rth%COEFS%COEF%PW_VAL_CHN
    CASE ("SS_VAL_CHN")
      I1 = rth%COEFS%COEF%SS_VAL_CHN
    CASE ("FMV_PC_SETS")
      I1 = rth%COEFS%COEF_PCCOMP%FMV_PC_SETS

    CASE ("SIZE_FMV_PC_NPRED")
      I1 = (/rth%COEFS%COEF_PCCOMP%FMV_PC_BANDS ,rth%COEFS%COEF_PCCOMP%FMV_PC_SETS(1)/)

    CASE DEFAULT
      err = errorstatus_fatal
      THROWM(ERR.NE.0,"INTEGER ARRAY VARIABLE "//TRIM(VARCH)//" NOT FOUND")

  END SELECT

CATCH

END SUBROUTINE

!> Return 2D integer array variable from optical depth coefficient structure
!! @param[out]    err                   return status
!! @param[in]     inst_id               ID for instrument
!! @param[in]     varch                 name of variable to return
!! @param[in]     m                     size of i2(:,1)
!! @param[in]     n                     size of i2(1,:)
!! @param[out]    i2                    returned variable value
SUBROUTINE RTTOV_GET_COEF_VAL_I2(err, inst_id, varch, M, N, I2)
!
#include "throw.h"

USE parkind1, ONLY : jpim
!INTF_OFF
USE rttov_wrapper_handle
!INTF_ON
  IMPLICIT NONE
  INTEGER(jpim),      INTENT(OUT) :: err
  INTEGER(jpim),      INTENT(IN)  :: inst_id
  CHARACTER(LEN=*),   INTENT(IN)  :: varch
  INTEGER(jpim),      INTENT(IN)  :: M
  INTEGER(jpim),      INTENT(IN)  :: N
  INTEGER(jpim),      INTENT(OUT) :: I2(M, N)
!INTF_END

  TYPE(rttovwrapperhandle_type), POINTER :: rth

#include "rttov_errorreport.interface"
!
TRY

  CALL rttov_wrapper_handle_get_inst(err, inst_id, rth)
  THROW(err.NE.0)

  SELECT CASE (TRIM(VARCH))

    CASE ("FMV_PC_NPRED")
      I2(:,:) = rth%COEFS%COEF_PCCOMP%PCREG(:,:)%FMV_PC_NPRED

    CASE DEFAULT
      err = errorstatus_fatal
      THROWM(ERR.NE.0,"INTEGER ARRAY VARIABLE "//TRIM(VARCH)//" NOT FOUND")

  END SELECT

CATCH

END SUBROUTINE

!> Return 1D real array variable from optical depth coefficient structure
!! @param[out]    err                   return status
!! @param[in]     inst_id               ID for instrument
!! @param[in]     varch                 name of variable to return
!! @param[in]     m                     size of r1
!! @param[out]    r1                    returned variable value
SUBROUTINE RTTOV_GET_COEF_VAL_R1(err, inst_id, varch, M, R1)

#include "throw.h"

USE parkind1, ONLY : jpim, jprb
!INTF_OFF
USE rttov_wrapper_handle
USE rttov_const, ONLY : GID_MIXED => GAS_ID_MIXED, &
                     &  GID_OZONE => GAS_ID_OZONE, &
                     &  GID_WATERVAPOUR => GAS_ID_WATERVAPOUR, &
                     &  GID_CO2 => GAS_ID_CO2,     &
                     &  GID_N2O => GAS_ID_N2O,     &
                     &  GID_CO  => GAS_ID_CO,      &
                     &  GID_CH4 => GAS_ID_CH4,     &
                     &  GID_SO2 => GAS_ID_SO2
!INTF_ON
  IMPLICIT NONE
  INTEGER(jpim),      INTENT(OUT) :: err
  INTEGER(jpim),      INTENT(IN)  :: inst_id
  CHARACTER(LEN=*),   INTENT(IN)  :: varch
  INTEGER(jpim),      INTENT(IN)  :: M
  REAL(jprb),         INTENT(OUT) :: R1(M)
!INTF_END

  TYPE(rttovwrapperhandle_type), POINTER :: rth
  INTEGER(jpim) :: N

#include "rttov_errorreport.interface"
!
TRY

  CALL rttov_wrapper_handle_get_inst(err, inst_id, rth)
  THROW(err.NE.0)

  SELECT CASE (TRIM(VARCH))

    CASE ("FF_CWN")
      R1 = rth%COEFS%COEF%FF_CWN
    CASE ("FF_BCO")
      R1 = rth%COEFS%COEF%FF_BCO
    CASE ("FF_BCS")
      R1 = rth%COEFS%COEF%FF_BCS
    CASE ("FF_GAM")
      R1 = rth%COEFS%COEF%FF_GAM
    CASE ("TT_A0")
      R1 = rth%COEFS%COEF%TT_A0
    CASE ("TT_A1")
      R1 = rth%COEFS%COEF%TT_A1
    CASE ("SS_SOLAR_SPECTRUM")
      R1 = rth%COEFS%COEF%SS_SOLAR_SPECTRUM
    CASE ("REF_PRFL_P")
      R1 = rth%COEFS%COEF%REF_PRFL_P
    CASE ("LIM_PRFL_P")
      R1 = rth%COEFS%COEF%LIM_PRFL_P
    CASE ("LIM_PRFL_TMAX")
      R1 = rth%COEFS%COEF%LIM_PRFL_TMAX
    CASE ("LIM_PRFL_TMIN")
      R1 = rth%COEFS%COEF%LIM_PRFL_TMIN
    CASE ("NOISE_IN")
      R1 = rth%COEFS%COEF_PCCOMP%NOISE_IN

! Aliases
    CASE ("WAVENUMBERS")
      R1 = rth%COEFS%COEF%FF_CWN

    CASE ("REF_PRESSURE")
      R1 = rth%COEFS%COEF%REF_PRFL_P

    CASE ("REF_TEMPERATURE")
      N  = rth%COEFS%COEF%FMV_GAS_POS(GID_MIXED)
      R1 = rth%COEFS%COEF%REF_PRFL_T(:,N)

    CASE ("REF_WATERVAPOR")
      IF( SIZE(R1) .NE. rth%COEFS%COEF%NLEVELS) THEN
        err = errorstatus_fatal
      ENDIF
      THROWM(ERR.NE.0,"REAL ARRAY ARGUMENT DOES NOT HAVE THE CORRECT DIMENSION")
      N  = rth%COEFS%COEF%FMV_GAS_POS(GID_WATERVAPOUR)
      R1 = rth%COEFS%COEF%REF_PRFL_MR(:,N)

     CASE ("REF_OZONE")
      IF( rth%COEFS%COEF%NOZONE .LE. 0 ) THEN
        err = errorstatus_fatal
      ENDIF
      THROWM(ERR.NE.0,"DOES NOT CONTAIN OZONE")
      N  = rth%COEFS%COEF%FMV_GAS_POS(GID_OZONE)
      R1 = rth%COEFS%COEF%REF_PRFL_MR(:,N)

     CASE ("REF_CO2")
      IF( rth%COEFS%COEF%NCO2 .LE. 0 ) THEN
        err = errorstatus_fatal
      ENDIF
      THROWM(ERR.NE.0,"DOES NOT CONTAIN CO2")
      N  = rth%COEFS%COEF%FMV_GAS_POS(GID_CO2)
      R1 = rth%COEFS%COEF%REF_PRFL_MR(:,N)

     CASE ("REF_N2O")
      IF( rth%COEFS%COEF%NN2O .LE. 0 ) THEN
        err = errorstatus_fatal
      ENDIF
      THROWM(ERR.NE.0,"DOES NOT CONTAIN N2O")
      N  = rth%COEFS%COEF%FMV_GAS_POS(GID_N2O)
      R1 = rth%COEFS%COEF%REF_PRFL_MR(:,N)

     CASE ("REF_CO")
      IF( rth%COEFS%COEF%NCO .LE. 0 ) THEN
        err = errorstatus_fatal
      ENDIF
      THROWM(ERR.NE.0,"DOES NOT CONTAIN CO")
      N  = rth%COEFS%COEF%FMV_GAS_POS(GID_CO)
      R1 = rth%COEFS%COEF%REF_PRFL_MR(:,N)

     CASE ("REF_CH4")
      IF( rth%COEFS%COEF%NCH4 .LE. 0 ) THEN
        err = errorstatus_fatal
      ENDIF
      THROWM(ERR.NE.0,"DOES NOT CONTAIN CH4")
      N  = rth%COEFS%COEF%FMV_GAS_POS(GID_CH4)
      R1 = rth%COEFS%COEF%REF_PRFL_MR(:,N)

     CASE ("REF_SO2")
      IF( rth%COEFS%COEF%NSO2 .LE. 0 ) THEN
        err = errorstatus_fatal
      ENDIF
      THROWM(ERR.NE.0,"DOES NOT CONTAIN SO2")
      N  = rth%COEFS%COEF%FMV_GAS_POS(GID_SO2)
      R1 = rth%COEFS%COEF%REF_PRFL_MR(:,N)


    CASE ("MIN_TEMPERATURE")
      R1 = rth%COEFS%COEF%LIM_PRFL_TMIN(:)

    CASE ("MIN_WATERVAPOR")
      IF( SIZE(R1) .NE. rth%COEFS%COEF%NLEVELS) THEN
        err = errorstatus_fatal
      ENDIF
      THROWM(ERR.NE.0,"REAL ARRAY ARGUMENT DOES NOT HAVE THE CORRECT DIMENSION")
      N  = rth%COEFS%COEF%FMV_GAS_POS(GID_WATERVAPOUR)
      R1 = rth%COEFS%COEF%LIM_PRFL_GMIN(:,N)

     CASE ("MIN_OZONE")
      IF( rth%COEFS%COEF%NOZONE .LE. 0 ) THEN
        err = errorstatus_fatal
      ENDIF
      THROWM(ERR.NE.0,"DOES NOT CONTAIN OZONE")
      N  = rth%COEFS%COEF%FMV_GAS_POS(GID_OZONE)
      R1 = rth%COEFS%COEF%LIM_PRFL_GMIN(:,N)

     CASE ("MIN_CO2")
      IF( rth%COEFS%COEF%NCO2 .LE. 0 ) THEN
        err = errorstatus_fatal
      ENDIF
      THROWM(ERR.NE.0,"DOES NOT CONTAIN CO2")
      N  = rth%COEFS%COEF%FMV_GAS_POS(GID_CO2)
      R1 = rth%COEFS%COEF%LIM_PRFL_GMIN(:,N)

     CASE ("MIN_N2O")
      IF( rth%COEFS%COEF%NN2O .LE. 0 ) THEN
        err = errorstatus_fatal
      ENDIF
      THROWM(ERR.NE.0,"DOES NOT CONTAIN N2O")
      N  = rth%COEFS%COEF%FMV_GAS_POS(GID_N2O)
      R1 = rth%COEFS%COEF%LIM_PRFL_GMIN(:,N)

     CASE ("MIN_CO")
      IF( rth%COEFS%COEF%NCO .LE. 0 ) THEN
        err = errorstatus_fatal
      ENDIF
      THROWM(ERR.NE.0,"DOES NOT CONTAIN CO")
      N  = rth%COEFS%COEF%FMV_GAS_POS(GID_CO)
      R1 = rth%COEFS%COEF%LIM_PRFL_GMIN(:,N)

     CASE ("MIN_CH4")
      IF( rth%COEFS%COEF%NCH4 .LE. 0 ) THEN
        err = errorstatus_fatal
      ENDIF
      THROWM(ERR.NE.0,"DOES NOT CONTAIN CH4")
      N  = rth%COEFS%COEF%FMV_GAS_POS(GID_CH4)
      R1 = rth%COEFS%COEF%LIM_PRFL_GMIN(:,N)

     CASE ("MIN_SO2")
      IF( rth%COEFS%COEF%NSO2 .LE. 0 ) THEN
        err = errorstatus_fatal
      ENDIF
      THROWM(ERR.NE.0,"DOES NOT CONTAIN SO2")
      N  = rth%COEFS%COEF%FMV_GAS_POS(GID_SO2)
      R1 = rth%COEFS%COEF%LIM_PRFL_GMIN(:,N)


    CASE ("MAX_TEMPERATURE")
      R1 = rth%COEFS%COEF%LIM_PRFL_TMAX(:)

    CASE ("MAX_WATERVAPOR")
      IF( SIZE(R1) .NE. rth%COEFS%COEF%NLEVELS) THEN
        err = errorstatus_fatal
      ENDIF
      THROWM(ERR.NE.0,"REAL ARRAY ARGUMENT DOES NOT HAVE THE CORRECT DIMENSION")
      N  = rth%COEFS%COEF%FMV_GAS_POS(GID_WATERVAPOUR)
      R1 = rth%COEFS%COEF%LIM_PRFL_GMAX(:,N)

     CASE ("MAX_OZONE")
      IF( rth%COEFS%COEF%NOZONE .LE. 0 ) THEN
        err = errorstatus_fatal
      ENDIF
      THROWM(ERR.NE.0,"DOES NOT CONTAIN OZONE")
      N  = rth%COEFS%COEF%FMV_GAS_POS(GID_OZONE)
      R1 = rth%COEFS%COEF%LIM_PRFL_GMAX(:,N)

     CASE ("MAX_CO2")
      IF( rth%COEFS%COEF%NCO2 .LE. 0 ) THEN
        err = errorstatus_fatal
      ENDIF
      THROWM(ERR.NE.0,"DOES NOT CONTAIN CO2")
      N  = rth%COEFS%COEF%FMV_GAS_POS(GID_CO2)
      R1 = rth%COEFS%COEF%LIM_PRFL_GMAX(:,N)

     CASE ("MAX_N2O")
      IF( rth%COEFS%COEF%NN2O .LE. 0 ) THEN
        err = errorstatus_fatal
      ENDIF
      THROWM(ERR.NE.0,"DOES NOT CONTAIN N2O")
      N  = rth%COEFS%COEF%FMV_GAS_POS(GID_N2O)
      R1 = rth%COEFS%COEF%LIM_PRFL_GMAX(:,N)

     CASE ("MAX_CO")
      IF( rth%COEFS%COEF%NCO .LE. 0 ) THEN
        err = errorstatus_fatal
      ENDIF
      THROWM(ERR.NE.0,"DOES NOT CONTAIN CO")
      N  = rth%COEFS%COEF%FMV_GAS_POS(GID_CO)
      R1 = rth%COEFS%COEF%LIM_PRFL_GMAX(:,N)

     CASE ("MAX_CH4")
      IF( rth%COEFS%COEF%NCH4 .LE. 0 ) THEN
        err = errorstatus_fatal
      ENDIF
      THROWM(ERR.NE.0,"DOES NOT CONTAIN CH4")
      N  = rth%COEFS%COEF%FMV_GAS_POS(GID_CH4)
      R1 = rth%COEFS%COEF%LIM_PRFL_GMAX(:,N)

     CASE ("MAX_SO2")
      IF( rth%COEFS%COEF%NSO2 .LE. 0 ) THEN
        err = errorstatus_fatal
      ENDIF
      THROWM(ERR.NE.0,"DOES NOT CONTAIN SO2")
      N  = rth%COEFS%COEF%FMV_GAS_POS(GID_SO2)
      R1 = rth%COEFS%COEF%LIM_PRFL_GMAX(:,N)


    CASE DEFAULT
      err = errorstatus_fatal
      THROWM(ERR.NE.0,"REAL ARRAY VARIABLE "//TRIM(VARCH)//" NOT FOUND")

  END SELECT

CATCH

END SUBROUTINE
