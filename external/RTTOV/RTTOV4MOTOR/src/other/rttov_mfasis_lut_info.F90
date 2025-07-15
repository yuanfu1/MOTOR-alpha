! Description:
!> @file
!!   Executable for displaying information about an MFASIS LUT file.
!
!> @brief
!!   Executable for displaying information about an MFASIS LUT file.
!!
!! @details
!!   Usage:
!!   $ rttov_coef_info.exe \-\-mfasis_lut ... [\-\-verbose]
!!
!!   Only the \-\-mfasis_lut argument is mandatory and is used to specify
!!   the MFASIS LUT file to read. The file must be in HDF5 format.
!!
!!   If supplied \-\-verbose will cause the README_LUT section to be printed
!!   out which may contain additional information.
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
!    Copyright 2018, EUMETSAT, All Rights Reserved.
!
PROGRAM rttov_mfasis_lut_info

#include "throw.h"

  USE rttov_types, ONLY : rttov_coef_mfasis
  USE parkind1, ONLY : jpim, jplm
  USE rttov_getoptions, ONLY : initoptions, getoption, checkoptions
  USE rttov_unix_env, ONLY : rttov_exit
  USE rttov_const, ONLY : &
    mfasis_dim_albedo,    &
    mfasis_dim_opdp,      &
    mfasis_dim_effdia,    &
    mfasis_dim_scaangle
#ifdef _RTTOV_HDF
  USE rttov_hdf_mod, ONLY : open_hdf, close_hdf
#endif

  IMPLICIT NONE

#include "rttov_errorreport.interface"
#include "rttov_dealloc_coef_mfasis.interface"
#ifdef _RTTOV_HDF
#include "rttov_hdf_load.interface"
#endif

  INTEGER(KIND=jpim)      :: err
  TYPE(rttov_coef_mfasis) :: coef_mfasis
  CHARACTER(LEN=256)      :: f_mfasis_lut = ""
  LOGICAL(KIND=jplm)      :: verbose = .FALSE.

  INTEGER(jpim)     :: l, d
  LOGICAL(jplm)     :: clw_opdp, clw_deff
  CHARACTER(LEN=64) :: msg
  !- End of header --------------------------------------------------------
  TRY

#ifndef _RTTOV_HDF
  err = errorstatus_fatal
  THROWM(err.NE.0,"This program requires RTTOV to be compiled against the HDF5 library")
#endif

  CALL initoptions("This program prints out information from an RTTOV/MFASIS LUT file")
  CALL getoption("--mfasis_lut", f_mfasis_lut, mnd=.TRUE._jplm, use="MFASIS LUT file name")
  CALL getoption("--verbose", verbose, use="more verbose output")
  CALL checkoptions()

#ifdef _RTTOV_HDF
  CALL open_hdf(.FALSE._jplm, err)
  THROW(err.NE.0)

  CALL rttov_hdf_load(err, f_mfasis_lut, "/MFASIS_CLD", mfasislutcld=coef_mfasis)
  THROWM(err.NE.0, "Failure reading MFASIS LUT file")
#endif

  WRITE(*,'(2x,a)') 'RTTOV/MFASIS cloud LUT:'
  WRITE(*,'(a)',advance='yes')

  WRITE(*,'(4x,a,i6)') 'LUT format version = ', coef_mfasis%version
  WRITE(*,'(a)',advance='yes')

  WRITE(*,'(4x,a,i6)') 'Number of channels in corresponding rtcoef file = ', coef_mfasis%nchannels_coef
  WRITE(*,'(4x,a,i6)') 'Number of channels supported = ', coef_mfasis%nchannels
  WRITE(*,'(4x,a,10i6)') 'Channels supported: ', coef_mfasis%channel_list
  WRITE(*,'(4x,a)')      'Number of LUTs per channel (>1 => variable water vapour):'
  WRITE(*,'(4x,20x,10i6)') coef_mfasis%lut(:)%nluts
  WRITE(*,'(a)',advance='yes')

  WRITE(*,'(4x,a)') 'Cloud properties used in training:'
  WRITE(*,'(6x,a,i2)') 'Cloud liquid water scheme = ', coef_mfasis%clw_scheme
  WRITE(*,'(6x,a,i2)') 'Cloud ice water scheme    = ', coef_mfasis%ice_scheme
  WRITE(*,'(a)',advance='yes')

  WRITE(*,'(4x,a)') 'LUT dimension min/max limits:'
  clw_opdp = .FALSE.  ! Assume liquid water opdp and deff axes come before ice (see rttov_mfasis)
  clw_deff = .FALSE.
  DO d = 1, coef_mfasis%ndims
    msg = ''
    SELECT CASE (coef_mfasis%lut_axes(d)%dim_type)
    CASE (mfasis_dim_albedo)
      msg = 'Albedo min/max                      ='
    CASE (mfasis_dim_scaangle)
      msg = 'Scattering angle min/max (deg)      ='
    CASE (mfasis_dim_opdp)
      IF (.NOT. clw_opdp) THEN
        msg = 'CLW total optical depth min/max     ='
        clw_opdp = .TRUE.
      ELSE
        msg = 'Ice total optical depth min/max     ='
      ENDIF
    CASE (mfasis_dim_effdia)
      IF (.NOT. clw_deff) THEN
        msg = 'CLW effective diameter min/max (um) ='
        clw_deff = .TRUE.
      ELSE
        msg = 'Ice effective diameter min/max (um) ='
      ENDIF
    END SELECT
    IF (msg .NE. '') THEN
      WRITE(*,'(6x,a,2f10.1)') TRIM(msg), coef_mfasis%lut_axes(d)%values(1), &
                               coef_mfasis%lut_axes(d)%values(coef_mfasis%lut_axes(d)%nvalues)
    ENDIF
  ENDDO
  IF (coef_mfasis%version > 0) THEN
    WRITE(*,'(6x,a,f10.1)') 'Maximum zenith angle (deg)          =', coef_mfasis%maxzenangle
  ENDIF
  WRITE(*,'(a)',advance='yes')

  IF (verbose) THEN
    WRITE(*,'(2x,a)') "README_LUT"
    DO l = 1, SIZE(coef_mfasis%readme_lut)
      IF (coef_mfasis%readme_lut(l) .EQ. 'xxxx') EXIT
      IF (coef_mfasis%readme_lut(l) .EQ. '') EXIT
      WRITE(*,'(4x,a)') TRIM(coef_mfasis%readme_lut(l))
    ENDDO
    WRITE(*,'(a)',advance='yes')
  ENDIF

  CALL rttov_dealloc_coef_mfasis(err, coef_mfasis)
  THROW(err.NE.0)

#ifdef _RTTOV_HDF
  CALL close_hdf(err)
  THROW(err.NE.0)
#endif

  PCATCH
END PROGRAM rttov_mfasis_lut_info
