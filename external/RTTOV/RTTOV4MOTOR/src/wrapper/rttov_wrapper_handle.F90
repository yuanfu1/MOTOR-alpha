! Description:
!> @file
!!   Defines internal data structures, constants and subroutines 
!!   for the wrapper.
!
!> @brief
!!   Defines internal data structures, constants and subroutines 
!!   for the wrapper.
!!
!! @details
!!   Defines internal data structure for the information associated 
!!   with each loaded instrument.
!!
!!   Defines the "gas" constants which are used to pass gas, aerosol
!!   and hydrometeor arrays to the wrapper.
!!
!!   Also contains subroutines for memory allocation/deallocation
!!   and reading coefs.
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
MODULE rttov_wrapper_handle

#include "throw.h"

  USE parkind1
  USE rttov_types
  USE rttov_const, ONLY : nwcl_max, naer_opac, naer_cams, sensor_id_ir, sensor_id_hi
  USE mod_rttov_emis_atlas, ONLY : rttov_emis_atlas_data
  USE mod_rttov_brdf_atlas, ONLY : rttov_brdf_atlas_data

  IMPLICIT NONE

  PUBLIC
  
  INTEGER(jpim), PRIVATE :: iii

  ! --- Gas ID definitions ---
  ! Gases
  INTEGER(jpim), PARAMETER :: gas_id_q   = 1
  INTEGER(jpim), PARAMETER :: gas_id_o3  = 2
  INTEGER(jpim), PARAMETER :: gas_id_co2 = 3
  INTEGER(jpim), PARAMETER :: gas_id_n2o = 4
  INTEGER(jpim), PARAMETER :: gas_id_co  = 5
  INTEGER(jpim), PARAMETER :: gas_id_ch4 = 6
  INTEGER(jpim), PARAMETER :: gas_id_so2 = 7

  ! CLW as absorber (non-scattering MW calculations)
  INTEGER(jpim), PARAMETER :: gas_id_clw = 15

  ! VIS/IR cloud profile inputs
  INTEGER(jpim), PARAMETER :: gas_id_cfrac = 20
  INTEGER(jpim), PARAMETER :: gas_id_lwc(1:nwcl_max) = (/ (iii, iii = 21, 25) /)
  INTEGER(jpim), PARAMETER :: gas_id_iwc   = 30
  INTEGER(jpim), PARAMETER :: gas_id_icede = 31
  INTEGER(jpim), PARAMETER :: gas_id_clwde = 32

  ! VIS/IR aerosol profile inputs
  INTEGER(jpim), PARAMETER :: gas_id_aer_opac(1:naer_opac) = (/ (iii, iii = 41, 40+naer_opac) /)
  INTEGER(jpim), PARAMETER :: gas_id_aer_cams(1:naer_cams) = (/ (iii, iii = 81, 80+naer_cams) /)
  INTEGER(jpim), PARAMETER :: gas_id_aer_user_min = 101
  INTEGER(jpim), PARAMETER :: gas_id_aer_user_max = 130

  ! RTTOV-SCATT (MW) cloud profile inputs
  INTEGER(jpim), PARAMETER :: gas_id_scatt_hydro_frac     = 60
  INTEGER(jpim), PARAMETER :: gas_id_scatt_clw            = 61
  INTEGER(jpim), PARAMETER :: gas_id_scatt_ciw            = 62
  INTEGER(jpim), PARAMETER :: gas_id_scatt_rain           = 63
  INTEGER(jpim), PARAMETER :: gas_id_scatt_snow           = 64
  INTEGER(jpim), PARAMETER :: gas_id_scatt_graupel        = 65
  INTEGER(jpim), PARAMETER :: gas_id_scatt_hydro_min      = 201
  INTEGER(jpim), PARAMETER :: gas_id_scatt_hydro_max      = 230
  INTEGER(jpim), PARAMETER :: gas_id_scatt_hydro_frac_min = 301
  INTEGER(jpim), PARAMETER :: gas_id_scatt_hydro_frac_max = 330
  ! --------------------------

  !> Data structure for a loaded instrument
  TYPE rttovwrapperhandle_type
    LOGICAL(jplm)                          :: init             = .FALSE.
    CHARACTER(LEN=256)                     :: file_coef        = ''
    CHARACTER(LEN=256)                     :: file_scaer       = ''
    CHARACTER(LEN=256)                     :: file_sccld       = ''
    CHARACTER(LEN=256)                     :: file_mfasis_cld  = ''
    CHARACTER(LEN=256)                     :: file_mfasis_aer  = ''
    CHARACTER(LEN=256)                     :: file_pccoef      = ''
    CHARACTER(LEN=256)                     :: file_hydrotable  = ''

    LOGICAL(jplm)                          :: verbose          = .FALSE.
    LOGICAL(jplm)                          :: check_opts       = .TRUE.
    LOGICAL(jplm)                          :: store_rad        = .FALSE.
    LOGICAL(jplm)                          :: store_rad2       = .FALSE.
    LOGICAL(jplm)                          :: store_trans      = .FALSE.
    LOGICAL(jplm)                          :: store_emis_terms = .FALSE.
    INTEGER(jpim)                          :: nthreads         = 1
    INTEGER(jpim)                          :: nprofs_per_call  = 1

    TYPE(rttov_options)                    :: opts
    TYPE(rttov_options_scatt)              :: opts_scatt
    TYPE(rttov_coefs)                      :: coefs
    TYPE(rttov_scatt_coef)                 :: coefs_scatt

    TYPE(rttov_transmission)               :: transmission
    TYPE(rttov_radiance)                   :: radiance
    TYPE(rttov_radiance2)                  :: radiance2
    TYPE(rttov_reflectivity)               :: reflectivity
    TYPE(rttov_scatt_emis_retrieval_type)  :: emis_terms

    INTEGER(jpim)                          :: inst_id

    TYPE(rttovwrapperhandle_type), POINTER :: next => NULL()
  ENDTYPE rttovwrapperhandle_type

  TYPE(rttovwrapperhandle_type), POINTER, SAVE, PRIVATE :: rth_first => NULL()

  !> Data structure for a loaded atlas data
  TYPE rttovwrapperatlas_type
    LOGICAL(jplm)               :: init = .FALSE.
    LOGICAL(jplm)               :: is_emis
    TYPE(rttov_emis_atlas_data) :: emis_atlas
    TYPE(rttov_brdf_atlas_data) :: brdf_atlas

    INTEGER(jpim)               :: atlas_id

    TYPE(rttovwrapperatlas_type), POINTER :: next => NULL()
  ENDTYPE rttovwrapperatlas_type

  TYPE(rttovwrapperatlas_type), POINTER, SAVE, PRIVATE :: ath_first => NULL()

CONTAINS


!> Creates a new instrument
!! @param[out]  err      return status
!! @param[out]  inst_id  ID of new instrument
!! @param       rth      pointer to handle of inst_id
SUBROUTINE rttov_wrapper_handle_new_inst(err, inst_id, rth)

  INTEGER(jpim), INTENT(OUT) :: err
  INTEGER(jpim), INTENT(OUT) :: inst_id
  TYPE(rttovwrapperhandle_type), POINTER :: rth

  TYPE(rttovwrapperhandle_type), POINTER :: next

#include "rttov_errorreport.interface"

  TRY

  IF (.NOT. ASSOCIATED(rth_first)) THEN
    ALLOCATE(rth_first, STAT=err)
    THROWM(err.NE.0, 'Error allocating new instrument handle')
    rth => rth_first
    rth%inst_id = 1
    inst_id = 1
  ELSE
    ! Check if rth_first has inst_id 1: if not we'll insert the new one at the beginning of the list
    IF (rth_first%inst_id > 1) THEN
      ALLOCATE(rth, STAT=err)
      THROWM(err.NE.0, 'Error allocating new instrument handle')
      rth%inst_id = 1
      inst_id = 1
      rth%next => rth_first
      rth_first => rth
    ELSE
      ! Find first missing inst_id (they are strictly increasing) or stop at end of list
      rth => rth_first
      inst_id = 1
      DO WHILE (ASSOCIATED(rth%next))
        inst_id = inst_id + 1
        IF (inst_id .NE. rth%next%inst_id) EXIT
        rth => rth%next
      ENDDO

      IF (.NOT. ASSOCIATED(rth%next)) THEN ! Add a new instrument at the end of the list
        ALLOCATE(rth%next, STAT=err)                              ! Allocate the new inst in next
        THROWM(err.NE.0, 'Error allocating new instrument handle')
        inst_id = rth%inst_id + 1                                 ! Increment the inst_id
        rth => rth%next                                           ! Move to the new one
      ELSE ! Insert instrument in the middle of the list at the first free inst_id
        next => rth%next                                          ! Store pointer to next, we are inserting after rth
        ALLOCATE(rth%next, STAT=err)                              ! Allocate the new inst in next
        THROWM(err.NE.0, 'Error allocating new instrument handle')
        rth => rth%next                                           ! Move to the new one
        rth%next => next                                          ! Link next to the one that was stored above
      ENDIF
      rth%inst_id = inst_id
    ENDIF
  ENDIF

  CATCH
END SUBROUTINE rttov_wrapper_handle_new_inst


!> Returns the handle for the given instrument ID
!! @param[out]  err      return status
!! @param[in]   inst_id  ID of instrument to return
!! @param       rth      pointer to handle of inst_id
SUBROUTINE rttov_wrapper_handle_get_inst(err, inst_id, rth)

  INTEGER(jpim),                 INTENT(OUT) :: err
  INTEGER(jpim),                 INTENT(IN)  :: inst_id
  TYPE(rttovwrapperhandle_type), POINTER     :: rth

#include "rttov_errorreport.interface"

  TRY

  rth => rth_first
  IF (.NOT. ASSOCIATED(rth)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'No instruments allocated')
  ELSE
    DO WHILE (ASSOCIATED(rth))
      IF (rth%inst_id == inst_id) RETURN
      rth => rth%next
    ENDDO
  ENDIF
  NULLIFY(rth)
  err = errorstatus_fatal
  THROWM(err.NE.0, 'Specified inst ID not found')

  CATCH
END SUBROUTINE rttov_wrapper_handle_get_inst


!> Deletes the specified instrument ID
!! @param[out]  err      return status
!! @param       rth      pointer to handle of instrument to delete
SUBROUTINE rttov_wrapper_handle_del_inst(err, rth)

  INTEGER(jpim), INTENT(OUT) :: err
  TYPE(rttovwrapperhandle_type), POINTER :: rth

  TYPE(rttovwrapperhandle_type), POINTER :: prev

#include "rttov_errorreport.interface"
#include "rttov_dealloc_coefs.interface"
#include "rttov_dealloc_scattcoeffs.interface"

  TRY

  IF (rth%init) THEN
    CALL rttov_dealloc_coefs(err, rth%coefs)
    THROWM(err.NE.0, 'Error deallocating coefficients')

    CALL rttov_dealloc_scattcoeffs(rth%coefs_scatt) ! These were nullified on load so this is safe

    ! Deallocate RTTOV output structures - for deallocation the actual values
    !   of nprofiles, nchannels and nlevels arguments don't matter
    CALL rttov_wrapper_handle_alloc(err, rth, 1_jpim, 1_jpim, 1_jpim, .TRUE._jplm, 0_jpim)
    THROW(err.NE.0)

    rth%init = .FALSE.
  ENDIF

  IF (rth%inst_id == rth_first%inst_id) THEN ! First in the list
    IF (ASSOCIATED(rth%next)) THEN ! There are others in the list
      rth_first => rth%next
    ELSE ! There is only one in the list
      NULLIFY(rth_first)
    ENDIF
  ELSE ! Not the first in the list
    ! Find the handle of the previous instrument
    prev => rth_first
    DO WHILE (rth%inst_id .NE. prev%next%inst_id)
      prev => prev%next
    ENDDO

    IF (ASSOCIATED(rth%next)) THEN ! In the middle of the list
      prev%next => rth%next
    ELSE ! At the end of the list
      NULLIFY(prev%next)
    ENDIF
  ENDIF

  DEALLOCATE(rth, STAT=err)
  THROWM(err.NE.0, 'Deallocation of instrument handle')
  NULLIFY(rth)

  CATCH
END SUBROUTINE rttov_wrapper_handle_del_inst


!> Delete all loaded instrument handle structures
!! @param[out]    err        return status
SUBROUTINE rttov_wrapper_handle_del_all_inst(err)

  INTEGER(jpim), INTENT(OUT) :: err

  TYPE(rttovwrapperhandle_type), POINTER :: rth

#include "rttov_errorreport.interface"

  TRY

  rth => rth_first
  DO WHILE(ASSOCIATED(rth))
    CALL rttov_wrapper_handle_del_inst(err, rth)
    THROW(err.NE.0)
  ENDDO

  CATCH
END SUBROUTINE rttov_wrapper_handle_del_all_inst


!> Creates a new atlas
!! @param[out]  err       return status
!! @param[out]  atlas_id  ID of new atlas
!! @param       ath       pointer to handle of atlas_id
SUBROUTINE rttov_wrapper_handle_new_atlas(err, atlas_id, ath)

  INTEGER(jpim), INTENT(OUT) :: err
  INTEGER(jpim), INTENT(OUT) :: atlas_id
  TYPE(rttovwrapperatlas_type), POINTER :: ath

  TYPE(rttovwrapperatlas_type), POINTER :: next

#include "rttov_errorreport.interface"

  TRY

  IF (.NOT. ASSOCIATED(ath_first)) THEN
    ALLOCATE(ath_first, STAT=err)
    THROWM(err.NE.0, 'Error allocating new atlas handle')
    ath => ath_first
    ath%atlas_id = 1
    atlas_id = 1
  ELSE
    ! Check if ath_first has atlas_id 1: if not we'll insert the new one at the beginning of the list
    IF (ath_first%atlas_id > 1) THEN
      ALLOCATE(ath, STAT=err)
      THROWM(err.NE.0, 'Error allocating new atlas handle')
      ath%atlas_id = 1
      atlas_id = 1
      ath%next => ath_first
      ath_first => ath
    ELSE
      ! Find first missing atlas_id (they are strictly increasing) or stop at end of list
      ath => ath_first
      atlas_id = 1
      DO WHILE (ASSOCIATED(ath%next))
        atlas_id = atlas_id + 1
        IF (atlas_id .NE. ath%next%atlas_id) EXIT
        ath => ath%next
      ENDDO

      IF (.NOT. ASSOCIATED(ath%next)) THEN ! Add a new atlas at the end of the list
        ALLOCATE(ath%next, STAT=err)                              ! Allocate the new atlas in next
        THROWM(err.NE.0, 'Error allocating new atlas handle')
        atlas_id = ath%atlas_id + 1                               ! Increment the atlas_id
        ath => ath%next                                           ! Move to the new one
      ELSE ! Insert atlas in the middle of the list at the first free atlas_id
        next => ath%next                                          ! Store pointer to next, we are inserting after ath
        ALLOCATE(ath%next, STAT=err)                              ! Allocate the new inst in next
        THROWM(err.NE.0, 'Error allocating new atlas handle')
        ath => ath%next                                           ! Move to the new one
        ath%next => next                                          ! Link next to the one that was stored above
      ENDIF
      ath%atlas_id = atlas_id
    ENDIF
  ENDIF

  CATCH
END SUBROUTINE rttov_wrapper_handle_new_atlas


!> Returns the handle for the given atlas ID
!! @param[out]  err       return status
!! @param[in]   atlas_id  ID of instrument to return
!! @param       ath       pointer to handle of atlas_id
SUBROUTINE rttov_wrapper_handle_get_atlas(err, atlas_id, ath)

  INTEGER(jpim),                INTENT(OUT) :: err
  INTEGER(jpim),                INTENT(IN)  :: atlas_id
  TYPE(rttovwrapperatlas_type), POINTER     :: ath

#include "rttov_errorreport.interface"

  TRY

  ath => ath_first
  IF (.NOT. ASSOCIATED(ath)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'No atlases allocated')
  ELSE
    DO WHILE (ASSOCIATED(ath))
      IF (ath%atlas_id == atlas_id) RETURN
      ath => ath%next
    ENDDO
  ENDIF
  NULLIFY(ath)
  err = errorstatus_fatal
  THROWM(err.NE.0, 'Specified atlas ID not found')

  CATCH
END SUBROUTINE rttov_wrapper_handle_get_atlas


!> Deletes the specified atlas ID
!! @param[out]  err      return status
!! @param       ath      pointer to handle of atlas to delete
SUBROUTINE rttov_wrapper_handle_del_atlas(err, ath)

  INTEGER(jpim), INTENT(OUT) :: err
  TYPE(rttovwrapperatlas_type), POINTER :: ath

  TYPE(rttovwrapperatlas_type), POINTER :: prev

#include "rttov_errorreport.interface"
#include "rttov_deallocate_emis_atlas.interface"
#include "rttov_deallocate_brdf_atlas.interface"

  TRY

  IF (ath%init) THEN
    IF (ath%is_emis) THEN
      CALL rttov_deallocate_emis_atlas(ath%emis_atlas)
    ELSE
      CALL rttov_deallocate_brdf_atlas(ath%brdf_atlas)
    ENDIF
    ath%init = .FALSE.
  ENDIF

  IF (ath%atlas_id == ath_first%atlas_id) THEN ! First in the list
    IF (ASSOCIATED(ath%next)) THEN ! There are others in the list
      ath_first => ath%next
    ELSE ! There is only one in the list
      NULLIFY(ath_first)
    ENDIF
  ELSE ! Not the first in the list
    ! Find the handle of the previous atlas
    prev => ath_first
    DO WHILE (ath%atlas_id .NE. prev%next%atlas_id)
      prev => prev%next
    ENDDO

    IF (ASSOCIATED(ath%next)) THEN ! In the middle of the list
      prev%next => ath%next
    ELSE ! At the end of the list
      NULLIFY(prev%next)
    ENDIF
  ENDIF

  DEALLOCATE(ath, STAT=err)
  THROWM(err.NE.0, 'Deallocation of atlas handle')
  NULLIFY(ath)

  CATCH
END SUBROUTINE rttov_wrapper_handle_del_atlas


!> Delete all loaded atlas handle structures
!! @param[out]    err        return status
SUBROUTINE rttov_wrapper_handle_del_all_atlas(err)

  INTEGER(jpim), INTENT(OUT) :: err

  TYPE(rttovwrapperatlas_type), POINTER :: ath

#include "rttov_errorreport.interface"

  TRY

  ath => ath_first
  DO WHILE(ASSOCIATED(ath))
    CALL rttov_wrapper_handle_del_atlas(err, ath)
    THROW(err.NE.0)
  ENDDO

  CATCH
END SUBROUTINE rttov_wrapper_handle_del_all_atlas


!> Reads coefficient files
!! @param[out]    err       return status
!! @param         rth       instrument structure
!! @param[in]     channels  list of channels to read ([0] => read all channels)
SUBROUTINE rttov_wrapper_handle_load(err, rth, channels)

  USE rttov_const, ONLY : vis_scatt_mfasis

  INTEGER(jpim),                 INTENT(OUT)   :: err
  TYPE(rttovwrapperhandle_type), POINTER       :: rth
  INTEGER(jpim),                 INTENT(IN)    :: channels(:)

#include "rttov_errorreport.interface"
#include "rttov_read_coefs.interface"
#include "rttov_read_scattcoeffs.interface"
#include "rttov_nullify_scattcoeffs.interface"

  LOGICAL(jplm) :: allchannels, addclouds, addaerosl
  INTEGER(jpim) :: vis_scatt_model

  TRY

  allchannels = (SIZE(channels) == 1 .AND. channels(1) <= 0) .OR. (rth%file_hydrotable .NE. '')

  ! Warn if addclouds/addaerosl set but cld/aer optical property files are missing
  IF (rth%opts%rt_ir%addclouds .AND. &
      .NOT. rth%opts%rt_ir%user_cld_opt_param .AND. &
      .NOT. rth%file_sccld .NE. '') THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'sccldcoef file must be specified if addclouds is true')
  ENDIF
  IF (rth%opts%rt_ir%addaerosl .AND. &
      .NOT. rth%opts%rt_ir%user_aer_opt_param .AND. &
      .NOT. rth%file_scaer .NE. '') THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'scaercoef file must be specified if addaerosl is true')
  ENDIF

  ! Ensure cloud/aerosol coefs are read if filenames were specified
  addclouds = rth%opts%rt_ir%addclouds
  addaerosl = rth%opts%rt_ir%addaerosl
  rth%opts%rt_ir%addclouds = (rth%file_sccld .NE. '')
  rth%opts%rt_ir%addaerosl = (rth%file_scaer .NE. '')

  ! Ensure MFASIS LUT is read if filename was specified
  vis_scatt_model = rth%opts%rt_ir%vis_scatt_model
  IF (rth%file_mfasis_cld .NE. '' .OR. rth%file_mfasis_aer .NE. '') &
    rth%opts%rt_ir%vis_scatt_model = vis_scatt_mfasis

  IF (allchannels) THEN
    CALL rttov_read_coefs(err, rth%coefs, rth%opts, &
                          file_coef       = rth%file_coef,       &
                          file_scaer      = rth%file_scaer,      &
                          file_sccld      = rth%file_sccld,      &
!                           file_mfasis_aer = rth%file_mfasis_aer, &
                          file_mfasis_cld = rth%file_mfasis_cld, &
                          file_pccoef     = rth%file_pccoef)
  ELSE
    CALL rttov_read_coefs(err, rth%coefs, rth%opts, &
                          file_coef       = rth%file_coef,       &
                          file_scaer      = rth%file_scaer,      &
                          file_sccld      = rth%file_sccld,      &
!                           file_mfasis_aer = rth%file_mfasis_aer, &
                          file_mfasis_cld = rth%file_mfasis_cld, &
                          file_pccoef     = rth%file_pccoef,     &
                          channels        = channels)
  ENDIF

  ! Revert the options to original values
  rth%opts%rt_ir%addclouds = addclouds
  rth%opts%rt_ir%addaerosl = addaerosl
  rth%opts%rt_ir%vis_scatt_model = vis_scatt_model

  THROWM(err.NE.0, 'Error reading coefficient files')

  CALL rttov_nullify_scattcoeffs(rth%coefs_scatt)
  IF (rth%file_hydrotable .NE. '') THEN
    CALL rttov_read_scattcoeffs(err, rth%opts_scatt, rth%coefs, rth%coefs_scatt, &
                                file_coef = rth%file_hydrotable)
    THROWM(err.NE.0, 'Error reading hydrotable file')
  ENDIF

  CATCH
END SUBROUTINE rttov_wrapper_handle_load


!> Allocate/deallocate wrapper handle output structures
!! @param[out]    err        return status
!! @param         rth        instrument structure
!! @param[in]     nprofiles  number of profiles
!! @param[in]     nchannels  number of channels simulated PER PROFILE
!! @param[in]     nlevels    number of levels in input profiles
!! @param[in]     calc_zef   flag for radar calculations
!! @param[in]     asw        1 => allocate; 0 => deallocate
SUBROUTINE rttov_wrapper_handle_alloc(err, rth, nprofiles, nchannels, nlevels, calc_zef, asw)

  INTEGER(jpim),                 INTENT(OUT)   :: err
  TYPE(rttovwrapperhandle_type), POINTER       :: rth
  INTEGER(jpim),                 INTENT(IN)    :: nprofiles
  INTEGER(jpim),                 INTENT(IN)    :: nchannels
  INTEGER(jpim),                 INTENT(IN)    :: nlevels
  LOGICAL(jplm),                 INTENT(IN)    :: calc_zef
  INTEGER(jpim),                 INTENT(IN)    :: asw

  INTEGER(jpim) :: nchanprof

#include "rttov_errorreport.interface"
#include "rttov_alloc_rad.interface"
#include "rttov_alloc_transmission.interface"
#include "rttov_alloc_reflectivity.interface"
#include "rttov_alloc_emis_ret_terms.interface"

  TRY

  nchanprof = nchannels * nprofiles

  IF (asw == 1_jpim) THEN

    IF (rth%store_rad) THEN
      IF (rth%store_rad2) THEN
        CALL rttov_alloc_rad(err, nchanprof, rth%radiance, nlevels, asw, rth%radiance2)
      ELSE
        CALL rttov_alloc_rad(err, nchanprof, rth%radiance, nlevels, asw)
      ENDIF
      THROWM(err.NE.0, 'Error allocating wrapper radiance')
    ENDIF

    IF (rth%store_trans) THEN
      CALL rttov_alloc_transmission(err, rth%transmission, nlevels, nchanprof, asw)
      THROWM(err.NE.0, 'Error allocating wrapper transmission')
    ENDIF

    IF (rth%store_emis_terms) THEN
      CALL rttov_alloc_emis_ret_terms(err, nchanprof, rth%emis_terms, asw)
      THROWM(err.NE.0, 'Error allocating wrapper emis_terms')
    ENDIF

    IF (calc_zef) THEN
      CALL rttov_alloc_reflectivity(err, nchanprof, rth%reflectivity, nlevels, asw)
      THROWM(err.NE.0, 'Error allocating wrapper radar reflectivity')
    ENDIF

  ELSE IF (asw == 0_jpim) THEN

    IF (rth%store_rad) THEN
      IF (ASSOCIATED(rth%radiance%clear)) THEN
        IF (rth%store_rad2) THEN
          CALL rttov_alloc_rad(err, nchanprof, rth%radiance, nlevels, asw, rth%radiance2)
        ELSE
          CALL rttov_alloc_rad(err, nchanprof, rth%radiance, nlevels, asw)
        ENDIF
        THROWM(err.NE.0, 'Error deallocating wrapper radiance')
      ENDIF
    ENDIF

    IF (rth%store_trans) THEN
      IF (ASSOCIATED(rth%transmission%tau_total)) THEN
        CALL rttov_alloc_transmission(err, rth%transmission, nlevels, nchanprof, asw)
        THROWM(err.NE.0, 'Error deallocating wrapper transmission')
      ENDIF
    ENDIF

    IF (rth%store_emis_terms) THEN
      IF (ASSOCIATED(rth%emis_terms%cfrac)) THEN
        CALL rttov_alloc_emis_ret_terms(err, nchanprof, rth%emis_terms, asw)
        THROWM(err.NE.0, 'Error deallocating wrapper emis_terms')
      ENDIF
    ENDIF

    IF (ASSOCIATED(rth%reflectivity%zef)) THEN
      CALL rttov_alloc_reflectivity(err, nchanprof, rth%reflectivity, nlevels, asw)
      THROWM(err.NE.0, 'Error deallocating wrapper radar reflectivity')
    ENDIF

  ENDIF

  CATCH
END SUBROUTINE rttov_wrapper_handle_alloc


!> Nullify output structures
!! @param         rth        instrument structure
SUBROUTINE rttov_wrapper_nullify_structs(rth)

  TYPE(rttovwrapperhandle_type), POINTER :: rth

  NULLIFY(rth%radiance%clear,      &
          rth%radiance%total,      &
          rth%radiance%bt_clear,   &
          rth%radiance%bt,         &
          rth%radiance%refl,       &
          rth%radiance%refl_clear, &
          rth%radiance%cloudy,     &
          rth%radiance%overcast,   &
          rth%radiance%quality,    &
          rth%radiance%geometric_height)
  NULLIFY(rth%radiance2%up,      &
          rth%radiance2%down,    &
          rth%radiance2%surf,    &
          rth%radiance2%upclear, &
          rth%radiance2%dnclear, &
          rth%radiance2%refldnclear)
  NULLIFY(rth%reflectivity%zef, &
          rth%reflectivity%azef)
  NULLIFY(rth%transmission%tau_total,           &
          rth%transmission%tau_levels,          &
          rth%transmission%tausun_total_path2,  &
          rth%transmission%tausun_levels_path2, &
          rth%transmission%tausun_total_path1,  &
          rth%transmission%tausun_levels_path1, &
          rth%transmission%tau_total_cld,       &
          rth%transmission%tau_levels_cld)
  NULLIFY(rth%emis_terms%cfrac,    &
          rth%emis_terms%bsfc,     &
          rth%emis_terms%tau_cld,  &
          rth%emis_terms%up_cld,   &
          rth%emis_terms%down_cld, &
          rth%emis_terms%tau_clr,  &
          rth%emis_terms%up_clr,   &
          rth%emis_terms%down_clr)
END SUBROUTINE

END MODULE rttov_wrapper_handle
