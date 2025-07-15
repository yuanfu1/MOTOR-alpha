! Description:
!> @file
!!   Subroutines for METEO-FRANCE CNRM MW emissivity atlas
!
!> @brief
!!   Subroutines for METEO-FRANCE CNRM MW emissivity atlas
!!
!! @details
!!   It is intended that this atlas be used via the RTTOV interface
!!   rather than by calling these subroutines directly.
!!
!!   Atlas data are stored in HDF5 files, one file per month/year/instrument
!!
!!   References:
!!   Karbou, F., E.Gérard, and F. Rabier, 2006: Microwave land emissivity and
!!   skin temperature for AMSU-A and -B assimilation over land, Q.J.R.
!!   Meteorol. Soc. 132, No. 620, Part A, pp. 2333-2355(23), doi:10.1256/qj.05.216
!!
!!   Karbou, F., E. Gérard, and F. Rabier, 2010: Global 4DVAR assimilation and
!!   forecast experiments using AMSU observations over land. Part I: Impacts of
!!   various land surface emissivity parameterizations. Wea. Forecasting, 25, 5-19.
!!
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
MODULE mod_cnrm_mw_atlas

#include "throw.h"

  USE parkind1, ONLY : jpim, jpis, jprb, jplm
  USE rttov_types, ONLY : rttov_coef
  USE rttov_const, ONLY :  &
        errorstatus_fatal, &
        inst_id_amsua,     &
        inst_id_amsub,     &
        inst_id_mhs,       &
        inst_id_ssmis,     &
        inst_id_atms,      &
        inst_name

#ifdef _RTTOV_HDF
  USE rttov_hdf_mod, ONLY : &
    open_hdf,             &
    close_hdf,            &
    is_hdf_open,          &
    is_hdf_64bit_reals
#endif

  IMPLICIT NONE

#include "rttov_errorreport.interface"

  PRIVATE
  PUBLIC :: cnrm_mw_atlas_data,    &
            rttov_cnrmmwemis_init, &
            rttov_cnrmmwemis,      &
            rttov_cnrm_close_atlas

  !> Data type for CNRM MW atlas data
  TYPE cnrm_mw_atlas_data
    PRIVATE
    INTEGER(KIND=jpis), POINTER :: emissivity(:,:,:,:)
    REAL(KIND=jprb),    POINTER :: frequencies(:)
    INTEGER(KIND=jpim), POINTER :: polarisation(:)
    REAL(KIND=jprb),    POINTER :: freq_range_min(:)
    REAL(KIND=jprb),    POINTER :: freq_range_max(:)
    REAL(KIND=jprb),    POINTER :: zenith_angles(:)
    REAL(KIND=jprb)             :: emis_offset
    REAL(KIND=jprb)             :: emis_scale
    REAL(KIND=jprb)             :: lat_start
    REAL(KIND=jprb)             :: lat_step
    REAL(KIND=jprb)             :: lat_stop
    REAL(KIND=jprb)             :: lon_start
    REAL(KIND=jprb)             :: lon_step
    REAL(KIND=jprb)             :: lon_stop
    INTEGER(KIND=jpim)          :: month
    INTEGER(KIND=jpim)          :: version_id
    INTEGER(KIND=jpim)          :: inst_id
    INTEGER(KIND=jpim)          :: year
    CHARACTER(LEN=32)           :: instrument
  END TYPE cnrm_mw_atlas_data

CONTAINS

!------------------------------------------
! Routines for initialising database
!------------------------------------------

  !> Initialise a CNRM MW atlas data structure. The atlas data are initialised
  !! for the instrument specified via the coefficients structure. Atlas data
  !! are available for multiple years: the year is specified when reading the
  !! data.
  !! @param[in]       path            path to atlas data files
  !! @param[in]       imonth          month of data to read (1-12)
  !! @param[in]       year            year of data to read (see atlas data file names)
  !! @param[in,out]   atlas           CNRM MW atlas data structure to initialise
  !! @param[in]       coef            RTTOV optical depth coefficient structure
  !! @param[in]       verbose         flag to turn verbose output on/off
  !! @param[out]      err             status on exit
  SUBROUTINE rttov_cnrmmwemis_init( &
                      path,       &
                      imonth,     &
                      year,       &
                      atlas,      &
                      coef,       &
                      verbose,    &
                      err         )

   USE rttov_unix_env, ONLY : rttov_upper_case

    CHARACTER(LEN=*),         INTENT(IN)    :: path
    INTEGER(KIND=jpim),       INTENT(IN)    :: imonth
    INTEGER(KIND=jpim),       INTENT(IN)    :: year
    TYPE(cnrm_mw_atlas_data), INTENT(INOUT) :: atlas
    TYPE(rttov_coef),         INTENT(IN)    :: coef
    LOGICAL(KIND=jplm),       INTENT(IN)    :: verbose
    INTEGER(KIND=jpim),       INTENT(OUT)   :: err

    CHARACTER(LEN=300) :: fn
    CHARACTER(LEN=2)   :: cmonth
    CHARACTER(LEN=4)   :: cyear
    CHARACTER(LEN=9)   :: iname
    CHARACTER(LEN=256) :: msg

    LOGICAL(KIND=jplm) :: file_exists
#ifdef _RTTOV_HDF
    LOGICAL(KIND=jplm) :: hdf_was_open, hdf_was_64bit_reals
#endif
    TRY
    CALL rttov_cnrm_nullify_pointers(atlas)

#ifdef _RTTOV_HDF
    hdf_was_open = is_hdf_open
    hdf_was_64bit_reals = is_hdf_64bit_reals
    CALL open_hdf(.TRUE._jplm, err)
    THROW(err.NE.0)
#else
    err = errorstatus_fatal
    THROWM(err.NE.0,'RTTOV must be compiled with HDF5 to use the CNRM atlas')
#endif

    IF (coef%id_inst == inst_id_amsua .OR. &
        coef%id_inst == inst_id_amsub .OR. &
        coef%id_inst == inst_id_mhs   .OR. &
        coef%id_inst == inst_id_ssmis .OR. &
        coef%id_inst == inst_id_atms) THEN

      IF (year > 1970) THEN
        WRITE(cyear, '(i4)') year
      ELSE
        cyear = '2015'
      ENDIF
      WRITE(cmonth, '(i2.2)') imonth
      CALL rttov_upper_case(iname, inst_name(coef%id_inst))
      IF (TRIM(iname) == 'MHS') iname = 'AMSUB'

      !----------------------------------------------------------------------------
      ! reading the data of CNRM MW Land Surface Global Emissivity
      !----------------------------------------------------------------------------
      fn = TRIM(path)//'CNRM_MWEMIS_'//TRIM(iname)//'_'//cyear//cmonth//'.H5'
      INQUIRE(FILE=fn, EXIST=file_exists)
      IF (file_exists) THEN
#ifdef _RTTOV_HDF
        CALL rttov_cnrm_read_coefs(err, TRIM(fn), atlas)
#endif
        THROW(err.NE.0)
      ELSE
        err = errorstatus_fatal
        THROWM(err.NE.errorstatus_success,'CNRM atlas file not found: '//TRIM(fn))
      ENDIF
      IF (verbose) INFO('Using CNRM coefs: '//TRIM(fn))

      IF (coef%id_inst == inst_id_mhs .AND. atlas%inst_id == inst_id_amsub) THEN
        atlas%inst_id = inst_id_mhs
      ENDIF

      IF (atlas%inst_id .NE. coef%id_inst) THEN
        err = errorstatus_fatal
        WRITE(msg,'(a)') 'Atlas instrument does not match coefficients: '// &
                         TRIM(inst_name(coef%id_inst))//'  '//TRIM(inst_name(atlas%inst_id))
        THROWM(err.NE.0,msg)
      ENDIF
    ELSE
      err = errorstatus_fatal
      THROWM(err.NE.0,'Instrument not supported '//TRIM(inst_name(coef%id_inst)))
    END IF

#ifdef _RTTOV_HDF
    ! If HDF5 lib was open before rttov_read_coefs was called, make sure the real
    ! kind is the same as it was previously. Otherwise close the library.
    IF (hdf_was_open) THEN
      CALL open_hdf(hdf_was_64bit_reals, err)
      THROW(err.NE.0)
    ELSE
      CALL close_hdf(err)
      THROW(err.NE.0)
    ENDIF
#endif

    CATCH
  END SUBROUTINE rttov_cnrmmwemis_init


#ifdef _RTTOV_HDF
  !> Read CNRM MW atlas data file
  !! @param[out]     err     status on exit
  !! @param[in]      fn      full path of file to read
  !! @param[in,out]  atlas   CNRM MW atlas data structure
  SUBROUTINE rttov_cnrm_read_coefs(err, fn, atlas)

    INTEGER(KIND=jpim),       INTENT(OUT)   :: err
    CHARACTER(LEN=*),         INTENT(IN)    :: fn
    TYPE(cnrm_mw_atlas_data), INTENT(INOUT) :: atlas

#include "rttov_hdf_load.interface"
    TRY

    CALL rttov_hdf_load(err, fn, "/", sname="EMISSIVITY", PIS4=atlas%emissivity)
    THROWM(err.NE.0,"Cannot load cnrm emissivity from "//TRIM(fn))

    CALL rttov_hdf_load(err, fn, "/", sname="FREQUENCIES", PR1=atlas%frequencies)
    THROWM(err.NE.0,"Cannot load cnrm FREQUENCIES from "//TRIM(fn))
    CALL rttov_hdf_load(err, fn, "/", sname="FREQ_RANGE_MIN", PR1=atlas%freq_range_min)
    THROWM(err.NE.0,"Cannot load cnrm FREQUENCY RANGE MIN from "//TRIM(fn))
    CALL rttov_hdf_load(err, fn, "/", sname="FREQ_RANGE_MAX", PR1=atlas%freq_range_max)
    THROWM(err.NE.0,"Cannot load cnrm FREQUENCY RANGE MAX from "//TRIM(fn))

    CALL rttov_hdf_load(err, fn, "/", sname="POLARISATION", PI1=atlas%polarisation)
    THROWM(err.NE.0,"Cannot load cnrm POLARISATION from "//TRIM(fn))

    CALL rttov_hdf_load(err, fn, "/", sname="ZENITH_ANGLES", PR1=atlas%zenith_angles)
    THROWM(err.NE.0,"Cannot load cnrm ZENITH_ANGLES from "//TRIM(fn))

    CALL rttov_hdf_load(err, fn, "/", sname="EMIS_OFFSET", R0=atlas%emis_offset)
    THROWM(err.NE.0,"Cannot load cnrm EMIS_OFFSET from "//TRIM(fn))
    CALL rttov_hdf_load(err, fn, "/", sname="EMIS_SCALE", R0=atlas%emis_scale)
    THROWM(err.NE.0,"Cannot load cnrm EMIS_SCALE from "//TRIM(fn))

    CALL rttov_hdf_load(err, fn, "/", sname="INSTRUMENT", C0=atlas%instrument)
    THROWM(err.NE.0,"Cannot load cnrm INSTRUMENT from "//TRIM(fn))
    CALL rttov_hdf_load(err, fn, "/", sname="INST_ID", I0=atlas%inst_id)
    THROWM(err.NE.0,"Cannot load cnrm INST_ID from "//TRIM(fn))

    CALL rttov_hdf_load(err, fn, "/", sname="LAT_START", R0=atlas%lat_start)
    THROWM(err.NE.0,"Cannot load cnrm LAT_START from "//TRIM(fn))
    CALL rttov_hdf_load(err, fn, "/", sname="LAT_STEP", R0=atlas%lat_step)
    THROWM(err.NE.0,"Cannot load cnrm LAT_STEP from "//TRIM(fn))
    CALL rttov_hdf_load(err, fn, "/", sname="LAT_STOP", R0=atlas%lat_stop)
    THROWM(err.NE.0,"Cannot load cnrm LAT_STOP from "//TRIM(fn))

    CALL rttov_hdf_load(err, fn, "/", sname="LON_START", R0=atlas%lon_start)
    THROWM(err.NE.0,"Cannot load cnrm LON_START from "//TRIM(fn))
    CALL rttov_hdf_load(err, fn, "/", sname="LON_STEP", R0=atlas%lon_step)
    THROWM(err.NE.0,"Cannot load cnrm LON_STEP from "//TRIM(fn))
    CALL rttov_hdf_load(err, fn, "/", sname="LON_STOP", R0=atlas%lon_stop)
    THROWM(err.NE.0,"Cannot load cnrm LON_STOP from "//TRIM(fn))

    CALL rttov_hdf_load(err, fn, "/", sname="MONTH", I0=atlas%month)
    THROWM(err.NE.0,"Cannot load cnrm MONTH from "//TRIM(fn))
    CALL rttov_hdf_load(err, fn, "/", sname="VERSION_ID", I0=atlas%version_id)
    THROWM(err.NE.0,"Cannot load cnrm VERSION_ID from "//TRIM(fn))
    CALL rttov_hdf_load(err, fn, "/", sname="YEAR", I0=atlas%year)
    THROWM(err.NE.0,"Cannot load cnrm YEAR from "//TRIM(fn))

    CATCH
  END SUBROUTINE rttov_cnrm_read_coefs
#endif

  !> Main subroutine to return emissivities for the given location, satellite
  !! geometry, and frequencies. The id_inst argument is used as a safety check
  !! to ensure the atlas data were initialised for this instrument.
  !! @param[out]      err                 status on exit
  !! @param[in]       id_inst             instrument ID from associated rtcoef structure
  !! @param[in]       nchan               number of channels for which to return emissivities
  !! @param[in]       freq                channel frequency list (GHz) for required emissivities
  !! @param[in]       polar               channel polarisation list (FASTEM convention) for required emissivities
  !! @param[in]       atlas               CNRM MW atlas data structure
  !! @param[in]       latitude            latitude for which to return emissivities
  !! @param[in]       longitude_in        longitude for which to return emissivities
  !! @param[in]       zenangle            viewing zenith angle
  !! @param[out]      emis_out            calculated emissivities values
  SUBROUTINE rttov_cnrmmwemis(    &
                  err,            &
                  id_inst,        &
                  nchan,          &
                  freq,           &
                  polar,          &
                  atlas,          &
                  latitude,       &
                  longitude_in,   &
                  zenangle,       &
                  emis_out)

    INTEGER(KIND=jpim),       INTENT(OUT) :: err
    INTEGER(KIND=jpim),       INTENT(IN)  :: id_inst
    INTEGER(KIND=jpim),       INTENT(IN)  :: nchan
    REAL(KIND=jprb),          INTENT(IN)  :: freq(nchan)
    INTEGER(KIND=jpim),       INTENT(IN)  :: polar(nchan)
    TYPE(cnrm_mw_atlas_data), INTENT(IN)  :: atlas
    REAL(KIND=jprb),          INTENT(IN)  :: latitude
    REAL(KIND=jprb),          INTENT(IN)  :: longitude_in
    REAL(KIND=jprb),          INTENT(IN)  :: zenangle
    REAL(KIND=jprb),          INTENT(OUT) :: emis_out(nchan)

    REAL(KIND=jprb)    :: longitude ! in (-180,180]
    INTEGER(KIND=jpim) :: ilat, ilon, izen, ifreq, iz, ifr, ic

    TRY

    IF (id_inst .NE. atlas%inst_id) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0,"Atlas data initialised for different instrument to coefficients")
    ENDIF

    longitude = MODULO(longitude_in, 360._jprb)
    IF (longitude > 180._jprb) longitude = longitude - 360._jprb

    emis_out(:) = 0.0_jprb

    ilat = NINT((latitude  - atlas%lat_start) / atlas%lat_step + 1_jpim)
    ilon = NINT((longitude - atlas%lon_start) / atlas%lon_step + 1_jpim)

    izen = 0
    DO iz = 2, SIZE(atlas%zenith_angles)
      IF (zenangle >= atlas%zenith_angles(iz-1) .AND. zenangle < atlas%zenith_angles(iz)) THEN
        izen = iz - 1
      ENDIF
    ENDDO

    IF (izen > 0) THEN
      DO ic = 1, nchan
        ifreq = 0
        DO ifr = 1, SIZE(atlas%frequencies)
          IF (atlas%polarisation(ifr) < 0_jpim) THEN
            ! no test for polarisation
            IF (freq(ic) >= atlas%freq_range_min(ifr) .AND. &
                freq(ic) <= atlas%freq_range_max(ifr)) THEN
              ifreq = ifr
            ENDIF
          ELSE
            IF (freq(ic)  >= atlas%freq_range_min(ifr) .AND. &
                freq(ic)  <= atlas%freq_range_max(ifr) .AND. &
                polar(ic) == atlas%polarisation(ifr)) THEN
              ifreq = ifr
            ENDIF
          ENDIF
        ENDDO
        IF (ifreq > 0) THEN
          emis_out(ic) = atlas%emissivity(ifreq, izen, ilon, ilat) * 1._jprb / atlas%emis_scale - atlas%emis_offset
        ENDIF
      ENDDO
    ENDIF

    CATCH
  END SUBROUTINE rttov_cnrmmwemis

!------------------------------------
! Routine to deallocate atlas arrays
!------------------------------------
  !> Deallocate data in CNRM MW atlas data structure
  !! @param[in,out]   atlas   CNRM MW atlas data structure to deallocate
  SUBROUTINE rttov_cnrm_close_atlas(atlas)
    TYPE(cnrm_mw_atlas_data), INTENT(INOUT) :: atlas

    IF (ASSOCIATED(atlas%emissivity   )) DEALLOCATE(atlas%emissivity)
    IF (ASSOCIATED(atlas%frequencies  )) DEALLOCATE(atlas%frequencies)
    IF (ASSOCIATED(atlas%zenith_angles)) DEALLOCATE(atlas%zenith_angles)
    IF (ASSOCIATED(atlas%freq_range_min)) DEALLOCATE(atlas%freq_range_min)
    IF (ASSOCIATED(atlas%freq_range_max)) DEALLOCATE(atlas%freq_range_max)
    IF (ASSOCIATED(atlas%polarisation  )) DEALLOCATE(atlas%polarisation)

    CALL rttov_cnrm_nullify_pointers(atlas)

  END SUBROUTINE rttov_cnrm_close_atlas

  !> Nullify pointers in CNRM MW atlas data structure
  !! @param[in,out]   atlas   CNRM MW atlas data structure to nullify
  SUBROUTINE rttov_cnrm_nullify_pointers(atlas)
    TYPE(cnrm_mw_atlas_data), INTENT(INOUT) :: atlas

    NULLIFY(atlas%emissivity)
    NULLIFY(atlas%frequencies)
    NULLIFY(atlas%freq_range_min)
    NULLIFY(atlas%freq_range_max)
    NULLIFY(atlas%zenith_angles)
    NULLIFY(atlas%polarisation)

    atlas%emis_offset = 0.0_jprb
    atlas%emis_scale  = 0.0_jprb
    atlas%instrument  = "xxxx"
    atlas%lat_start   = 0.0_jprb
    atlas%lat_step    = 0.0_jprb
    atlas%lat_stop    = 0.0_jprb
    atlas%lon_start   = 0.0_jprb
    atlas%lon_step    = 0.0_jprb
    atlas%lon_stop    = 0.0_jprb
    atlas%month       = 0_jpim
    atlas%version_id  = 0_jpim
    atlas%inst_id     = 0_jpim
    atlas%year        = 0_jpim

  END SUBROUTINE rttov_cnrm_nullify_pointers

END MODULE mod_cnrm_mw_atlas
