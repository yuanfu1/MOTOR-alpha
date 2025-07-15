! Description:
!> @file
!!   Executable for converting coefficient file formats or
!!   creating new coefficient files with a subset of channels.
!
!> @brief
!!   Executable for converting coefficient file formats or
!!   creating new coefficient files with a subset of channels.
!!
!! @details
!!   For usage details see user guide or run:
!!   $ rttov_conv_coef.exe \-\-help
!!
!!   The \-\-coef-in argument is mandatory and the \-\-format-out
!!   argument should normally be specified.
!!
!!   RTTOV can usually detect the input file format so the
!!   \-\-format-in argument is not normally required.
!!
!!   Note that after extracting a subset of n channels to a
!!   new file these will indexed as 1...n rather than by
!!   their original channel numbering.
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
PROGRAM rttov_conv_coef

#include "throw.h"

  USE parkind1, ONLY : jpim, jplm
  USE rttov_types, ONLY : rttov_coefs, rttov_options
  USE rttov_const, ONLY : vis_scatt_mfasis
  USE rttov_getoptions, ONLY : initoptions, getoption, checkoptions
  USE rttov_unix_env, ONLY : rttov_exit

#ifdef _RTTOV_HDF
  USE rttov_hdf_mod, ONLY : open_hdf, close_hdf
#endif

  IMPLICIT NONE

#include "rttov_read_coefs.interface"
#include "rttov_write_coefs.interface"
#include "rttov_dealloc_coefs.interface"
#include "rttov_errorreport.interface"

  TYPE(rttov_coefs)           :: coefs
  TYPE(rttov_options)         :: opts
  INTEGER(KIND=jpim)          :: err
  INTEGER(KIND=jpim), POINTER :: channels(:) => NULL()
  INTEGER(KIND=jpim), POINTER :: channels_rec(:) => NULL()

  CHARACTER(LEN=256) :: f_coef_in  = "", f_scaer_in  = "", f_sccld_in  = "", f_pccoef_in  = ""
!   CHARACTER(LEN=256) :: f_mfasis_aer_in  = "", 
  CHARACTER(LEN=256) :: f_mfasis_cld_in  = ""
  CHARACTER(LEN=256) :: f_coef_out = "", f_scaer_out = "", f_sccld_out = "", f_pccoef_out = ""
!   CHARACTER(LEN=256) :: f_mfasis_aer_out = ""
  CHARACTER(LEN=256) :: f_mfasis_cld_out = ""
  CHARACTER(LEN=256) :: coef_format_in = "", coef_format_out = ""
  LOGICAL(KIND=jplm)  :: AllInOne = .FALSE.
  LOGICAL(KIND=jplm)  :: HDF5REALS32 = .FALSE.
  LOGICAL(KIND=jplm)  :: HDF5REALS64 = .FALSE.
  LOGICAL(KIND=jplm)  :: compress = .FALSE.
  LOGICAL(KIND=jplm)  :: force_single = .FALSE.
  LOGICAL(KIND=jplm)  :: force_double
  !- End of header --------------------------------------------------------

  TRY

  CALL initoptions( "This program converts RTTOV coefficient files")

#ifdef _RTTOV_HDF
  CALL getoption( "--format-in",  coef_format_in,  USE="FORMATTED|UNFORMATTED|HDF5")
  CALL getoption( "--format-out", coef_format_out, mnd=.TRUE._jplm, USE="FORMATTED|UNFORMATTED|HDF5")
  CALL getoption( "--hdf5-reals32", HDF5REALS32, USE="Store reals in HDF5 32 bits ; (default is KIND jprb)")
  CALL getoption( "--all-in-one", AllInOne, USE="output all coefs in one file (HDF5 only)")
  CALL getoption( "--compress", compress, USE="use HDF5 internal GZIP compression (useless for low spectral res. instruments)")
  CALL getoption( "--force-single", force_single, USE="use 32bits reals for storage of 2D arrays and more")
#else
  CALL getoption( "--format-in",  coef_format_in,  USE="FORMATTED|UNFORMATTED")
  CALL getoption( "--format-out", coef_format_out, mnd=.TRUE._jplm, USE="FORMATTED|UNFORMATTED")
#endif

  CALL getoption( "--channels",  channels ,  USE="list of channels to extract")
  !CALL getoption( "--channels-rec", channels_rec, USE="list of channels_rec for PC")

  CALL getoption( "--coef-in",       f_coef_in, mnd=.TRUE._jplm, USE="input coefficient file name")
  CALL getoption( "--scaer-in",      f_scaer_in,      USE="input aerosol coefficient file name")
  CALL getoption( "--sccld-in",      f_sccld_in,      USE="input cloud coefficient file name")
!   CALL getoption( "--mfasis_aer-in", f_mfasis_aer_in, USE="input MFASIS aerosol file name")
  CALL getoption( "--mfasis_cld-in", f_mfasis_cld_in, USE="input MFASIS cloud file name")
  CALL getoption( "--pccoef-in",     f_pccoef_in,     USE="input PC coefficient file name")

  CALL getoption( "--coef-out",       f_coef_out,       USE="output coefficient file name")
  CALL getoption( "--scaer-out",      f_scaer_out,      USE="output aerosol coefficient file name")
  CALL getoption( "--sccld-out",      f_sccld_out,      USE="output cloud coefficient file name")
!   CALL getoption( "--mfasis_aer-out", f_mfasis_aer_out, USE="output MFASIS aerosol file name")
  CALL getoption( "--mfasis_cld-out", f_mfasis_cld_out, USE="output MFASIS cloud file name")
  CALL getoption( "--pccoef-out",     f_pccoef_out,     USE="output PC coefficient file name")

  CALL checkoptions()

  force_double = .NOT. force_single

  IF (TRIM(coef_format_out) .EQ. "HDF5") THEN
    IF (f_coef_out       .EQ. "" .AND. f_coef_in       .NE. "") f_coef_out = TRIM(f_coef_in)//'.H5'
    IF (f_scaer_out      .EQ. "" .AND. f_scaer_in      .NE. "") f_scaer_out = TRIM(f_scaer_in)//'.H5'
    IF (f_sccld_out      .EQ. "" .AND. f_sccld_in      .NE. "") f_sccld_out = TRIM(f_sccld_in)//'.H5'
!     IF (f_mfasis_aer_out .EQ. "" .AND. f_mfasis_aer_in .NE. "") f_mfasis_aer_out = TRIM(f_mfasis_aer_in)//'.H5'
    IF (f_mfasis_cld_out .EQ. "" .AND. f_mfasis_cld_in .NE. "") f_mfasis_cld_out = TRIM(f_mfasis_cld_in)//'.H5'
    IF (f_pccoef_out     .EQ. "" .AND. f_pccoef_in     .NE. "") f_pccoef_out = TRIM(f_pccoef_in)//'.H5'
    IF (AllInOne) THEN
      IF (f_scaer_in      .NE. "") f_scaer_out = f_coef_out
      IF (f_sccld_in      .NE. "") f_sccld_out = f_coef_out
!       IF (f_mfasis_aer_in .NE. "") f_mfasis_aer_out = f_coef_out
      IF (f_mfasis_cld_in .NE. "") f_mfasis_cld_out = f_coef_out
      IF (f_pccoef_in     .NE. "") f_pccoef_out = f_coef_out
    EndIf
  ELSE
    IF (f_coef_out       .EQ. "" .AND. f_coef_in       .NE. "") f_coef_out = TRIM(f_coef_in)//'.bin'
    IF (f_scaer_out      .EQ. "" .AND. f_scaer_in      .NE. "") f_scaer_out = TRIM(f_scaer_in)//'.bin'
    IF (f_sccld_out      .EQ. "" .AND. f_sccld_in      .NE. "") f_sccld_out = TRIM(f_sccld_in)//'.bin'
!     IF (f_mfasis_aer_out .EQ. "" .AND. f_mfasis_aer_in .NE. "") f_mfasis_aer_out = TRIM(f_mfasis_aer_in)//'.bin'
    IF (f_mfasis_cld_out .EQ. "" .AND. f_mfasis_cld_in .NE. "") f_mfasis_cld_out = TRIM(f_mfasis_cld_in)//'.bin'
    IF (f_pccoef_out     .EQ. "" .AND. f_pccoef_in     .NE. "") f_pccoef_out = TRIM(f_pccoef_in)//'.bin'
  ENDIF

  opts%rt_ir%addaerosl    = f_scaer_out .NE. ""
  opts%rt_ir%addclouds    = f_sccld_out .NE. ""
  opts%rt_ir%pc%addpc     = f_pccoef_in .NE. ""
  IF (f_mfasis_cld_out .NE. "") opts%rt_ir%vis_scatt_model = vis_scatt_mfasis
!   IF (f_mfasis_cld_out .NE. "" .OR. f_mfasis_aer_out .NE. "") opts%rt_ir%vis_scatt_model = vis_scatt_mfasis
  opts%rt_ir%pc%addradrec = opts%rt_ir%pc%addpc .AND. ASSOCIATED(channels_rec)

#ifdef _RTTOV_HDF
  HDF5REALS64 = .NOT. HDF5REALS32   ! 64 bits (JPRB) is the default
  CALL open_hdf(HDF5REALS64, err)
  THROW(err.NE.0)
#endif

!
! Do not use channels_rec argument. Reducing the rtcoef file only,
! not the pcscore coef file, allows to process any reconstructed radiances.
!
  IF (ASSOCIATED(channels)) THEN
    CALL rttov_read_coefs(err, coefs, opts,                       &
                          channels        = channels,             &
                          form_coef       = TRIM(coef_format_in), &
                          form_scaer      = TRIM(coef_format_in), &
                          form_sccld      = TRIM(coef_format_in), &
!                           form_mfasis_aer = TRIM(coef_format_in), &
                          form_mfasis_cld = TRIM(coef_format_in), &
                          form_pccoef     = TRIM(coef_format_in), &
                          file_coef       = f_coef_in,            &
                          file_scaer      = f_scaer_in,           &
                          file_sccld      = f_sccld_in,           &
!                           file_mfasis_aer = f_mfasis_aer_in,      &
                          file_mfasis_cld = f_mfasis_cld_in,      &
                          file_pccoef     = f_pccoef_in)
    THROW(err.NE.0)
  ELSE
    CALL rttov_read_coefs(err, coefs, opts,                       &
                          form_coef       = TRIM(coef_format_in), &
                          form_scaer      = TRIM(coef_format_in), &
                          form_sccld      = TRIM(coef_format_in), &
!                           form_mfasis_aer = TRIM(coef_format_in), &
                          form_mfasis_cld = TRIM(coef_format_in), &
                          form_pccoef     = TRIM(coef_format_in), &
                          file_coef       = f_coef_in,            &
                          file_scaer      = f_scaer_in,           &
                          file_sccld      = f_sccld_in,           &
!                           file_mfasis_aer = f_mfasis_aer_in,      &
                          file_mfasis_cld = f_mfasis_cld_in,      &
                          file_pccoef     = f_pccoef_in)
    THROW(err.NE.0)
  ENDIF

  CALL rttov_write_coefs(err, coefs, opts,                        &
                         form_coef       = TRIM(coef_format_out), &
                         form_scaer      = TRIM(coef_format_out), &
                         form_sccld      = TRIM(coef_format_out), &
!                          form_mfasis_aer = TRIM(coef_format_out), &
                         form_mfasis_cld = TRIM(coef_format_out), &
                         form_pccoef     = TRIM(coef_format_out), &
                         file_coef       = f_coef_out,            &
                         file_scaer      = f_scaer_out,           &
                         file_sccld      = f_sccld_out,           &
!                          file_mfasis_aer = f_mfasis_aer_out,      &
                         file_mfasis_cld = f_mfasis_cld_out,      &
                         file_pccoef     = f_pccoef_out,          &
                         compress        = compress,              &
                         force_double    = force_double)
  THROW(err.NE.0)

  IF (ASSOCIATED(channels))     DEALLOCATE(channels)
  IF (ASSOCIATED(channels_rec)) DEALLOCATE(channels_rec)

  CALL rttov_dealloc_coefs(err, coefs)
  THROWM(err.NE.0,"deallocation of coefs ")

#ifdef _RTTOV_HDF
  CALL close_hdf(err)
  THROW(err.NE.0)
#endif

  PCATCH
END PROGRAM rttov_conv_coef
