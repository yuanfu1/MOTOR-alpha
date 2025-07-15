! Description:
!> @file
!!   Executable for converting RTTOV v12-compatible
!!   coefficient files to v10/v11 format.
!
!> @brief
!!   Executable for converting RTTOV v12-compatible
!!   coefficient files to v10/v11 format.
!!
!! @details
!!   Usage:
!!   $ rttov11_conv_coef_12to11.exe \-\-coef-in ... \-\-coef-out ...
!!
!!   where \-\-coef-in specifies the input v12-compatible
!!   file and \-\-coef-out specifies the output v10/v11-
!!   compatible file. Input file must be ASCII or HDF5 format.
!!
!!   Note that certain types of v12 file cannot be converted
!!   (in particular files with variable SO2).
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
PROGRAM rttov11_conv_coef_12to11

#include "throw.h"

  USE parkind1, ONLY : jpim, jplm
  USE rttov_types, ONLY : rttov_coef
  USE rttov_unix_env, ONLY : rttov_exit
  USE rttov_getoptions, ONLY : initoptions, getoption, checkoptions
  USE rttov_coef_io_mod, ONLY : getlun, closelun
#ifdef _RTTOV_HDF
  USE rttov_hdf_mod, ONLY : open_hdf, close_hdf
#endif

  IMPLICIT NONE

  CHARACTER(LEN=256) :: f_coef_in, f_coef_out, f1
  CHARACTER(LEN=32)  :: file_format
  INTEGER(jpim)      :: err
  INTEGER(jpim)      :: file_id
  TYPE(rttov_coef)   :: coef

#include "rttov_errorreport.interface"
#include "rttov_cmpuc.interface"
#include "rttov_read_ascii_coef.interface"
#include "rttov11_write_ascii_coef.interface"
#include "rttov11_write_hdf5_coef.interface"
#include "rttov_init_coef.interface"
#include "rttov_nullify_coef.interface"
#include "rttov_dealloc_coef.interface"
#ifdef _RTTOV_HDF
#include "rttov_hdf_load.interface"
#endif

  TRY

  CALL initoptions( "This program converts RTTOV v12 optical depth coefficient files to v10/v11 format")
#ifndef _RTTOV_HDF
  CALL getoption("--coef-in", f_coef_in, mnd=.TRUE._jplm, use="Input ASCII v12 optical depth coefficient file")
  CALL getoption("--coef-out", f_coef_out, mnd=.TRUE._jplm, use="Output ASCII v10/v11 optical depth coefficient file")
#else
  CALL getoption("--coef-in", f_coef_in, mnd=.TRUE._jplm, use="Input ASCII/HDF5 v12 optical depth coefficient file")
  CALL getoption("--coef-out", f_coef_out, mnd=.TRUE._jplm, use="Output ASCII/HDF5 v10/v11 optical depth coefficient file")
#endif
  CALL checkoptions()


  ! ---------------------------------------------------------------------------
  ! Read input coefficient file
  ! ---------------------------------------------------------------------------

  CALL rttov_nullify_coef(coef) 

  CALL getlun(err, file_id, f1, file_format, .FALSE._jplm, 'rtcoef', f=TRIM(f_coef_in))
  THROW(err.NE.0)

  IF (rttov_cmpuc(file_format, 'unformatted')) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Invalid file format: must be ASCII or HDF5')
  ELSEIF (rttov_cmpuc(file_format, 'formatted')) THEN
    CALL rttov_read_ascii_coef(err, coef, file_id)
    THROWM(err.NE.0, 'Cannot open ASCII coefficient file '//TRIM(f_coef_in))
  ELSEIF (rttov_cmpuc(file_format, 'hdf5')) THEN
#ifndef _RTTOV_HDF
    err = errorstatus_fatal
    THROWM(err.NE.0, 'This program is not compiled with HDF5 capability; use RTTOV_HDF=1 with Makefile.PL')
#else
    CALL open_hdf(.TRUE._jplm, err)
    THROWM(err.NE.0, 'Error opening HDF5 interface')

    CALL rttov_hdf_load(err, f1, "/COEF", coef=coef)
    THROWM(err.NE.0, 'Cannot open HDF5 coefficient file '//TRIM(f_coef_in))
#endif
  ELSE
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Unknown coefficient file format '//TRIM(file_format))
  ENDIF

  CALL closelun(err, file_id)
  THROW(err.NE.0)

  CALL rttov_init_coef(err, coef)
  THROW(err.NE.0)

  ! Check for v13 coefficient files: these cannot be converted
  IF (coef%id_comp_lvl > 12) THEN
    CALL rttov_dealloc_coef(err, coef)
    err = errorstatus_fatal
    THROWM(err.NE.0,'v13-format rtcoef files cannot be converted, only v12-compatible files')
  ENDIF

  ! ---------------------------------------------------------------------------
  ! Write output coefficient file
  ! ---------------------------------------------------------------------------

  IF (TRIM(file_format) == 'formatted') THEN
    OPEN(file_id, file=TRIM(f_coef_out), form='formatted', iostat=err)
    THROWM(err.NE.0, 'Error opening output ASCII file')

    CALL rttov11_write_ascii_coef(err, file_id, coef)
    THROWM(err.NE.0, 'Error writing ASCII coefficients')

    CLOSE(file_id, iostat=err)
    THROWM(err.NE.0, 'Error closing output ASCII file')
  ELSEIF (TRIM(file_format) == 'hdf5') THEN
#ifdef _RTTOV_HDF
!     CALL rttov11_write_hdf5_coef(err, f_coef_out, "/COEF", create=.TRUE._jplm, coef=coef, compress=.TRUE._jplm)
    CALL rttov11_write_hdf5_coef(err, TRIM(f_coef_out), coef)
    THROWM(err.NE.0, 'Error writing HDF5 coefficients')

    CALL close_hdf(err)
    THROWM(err.NE.0, 'Error closing HDF5 interface')
#endif
  ENDIF

  CALL rttov_dealloc_coef(err, coef)
  THROWM(err.NE.0, 'Error deallocating coef')

  PCATCH
END PROGRAM rttov11_conv_coef_12to11
