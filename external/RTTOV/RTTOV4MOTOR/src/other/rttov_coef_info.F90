! Description:
!> @file
!!   Executable for displaying information about an optical
!!   depth (rtcoef) coefficient file.
!
!> @brief
!!   Executable for displaying information about an optical
!!   depth (rtcoef) coefficient file.
!!
!! @details
!!   Usage:
!!   $ rttov_coef_info.exe \-\-coef ... [\-\-format ...] [\-\-verbose]
!!
!!   Only the \-\-coef argument is mandatory and is used to specify
!!   the rtcoef file to read.
!!
!!   RTTOV can usually detect the file format so the
!!   \-\-format argument is not normally required. If specified the
!!   format should be one of FORMATTED, UNFORMATTED or HDF5.
!!
!!   If supplied \-\-verbose will cause additional information
!!   to be printed about each channel.
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
PROGRAM rttov_coef_info

#include "throw.h"

  USE rttov_types, ONLY : rttov_coefs, rttov_options
  USE parkind1, ONLY : jpim, jplm
  USE rttov_getoptions, ONLY : initoptions, getoption, checkoptions
  USE rttov_unix_env, ONLY : rttov_exit
#ifdef _RTTOV_HDF
  USE rttov_hdf_mod, ONLY : open_hdf, close_hdf
#endif

  IMPLICIT NONE

#include "rttov_errorreport.interface"
#include "rttov_read_coefs.interface"
#include "rttov_dealloc_coefs.interface"
#include "rttov_print_info.interface"

  INTEGER(KIND=jpim)  :: err
  TYPE(rttov_coefs)   :: coefs
  TYPE(rttov_options) :: opts
  CHARACTER(LEN=256)  :: f_coef = ""
  CHARACTER(LEN=256)  :: coef_format = ""
  LOGICAL(KIND=jplm)  :: verbose = .FALSE.

  !- End of header --------------------------------------------------------
  TRY

  CALL initoptions("This program prints out information on RTTOV optical depth coefficient files")
  CALL getoption("--coef", f_coef, mnd=.TRUE._jplm, use="coefficient file name")
#ifdef _RTTOV_HDF
  CALL getoption("--format", coef_format, use="FORMATTED|UNFORMATTED|HDF5")
#else
  CALL getoption("--format", coef_format, use="FORMATTED|UNFORMATTED")
#endif
  CALL getoption("--verbose", verbose, use="more verbose output")
  CALL checkoptions()

#ifdef _RTTOV_HDF
  CALL open_hdf(.FALSE._jplm, err)
  THROW(err .NE. 0)
#endif

  CALL rttov_read_coefs(err, coefs, opts, form_coef=coef_format, file_coef=f_coef)
  THROWM(err .NE. 0, "Failure reading coefficient file")

#ifdef _RTTOV_HDF
  CALL close_hdf(err)
  THROW(err .NE. 0)
#endif

  CALL rttov_print_info(coefs, verbose=verbose)

  CALL rttov_dealloc_coefs(err, coefs)
  THROWM(err .NE. 0, "Failure deallocating coefficients")

  PCATCH
END PROGRAM
