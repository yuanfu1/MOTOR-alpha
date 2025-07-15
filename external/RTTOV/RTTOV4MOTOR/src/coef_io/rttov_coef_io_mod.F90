! Description:
!> @file
!!   Useful subroutines for coefficient file I/O
!
!> @brief
!!   Useful subroutines for coefficient file I/O
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
!    Copyright 2015, EUMETSAT, All Rights Reserved.
!
MODULE rttov_coef_io_mod

  USE parkind1, ONLY : jpim, jplm

#include "throw.h"

  IMPLICIT NONE

#include "rttov_errorreport.interface"

  PRIVATE
  PUBLIC :: getlun, closelun

CONTAINS

!> Open an RTTOV coefficient file (rtcoef, scaer/cldcoef, pccoef) for I/O. If a
!! logical unit lun is passed in this is used; otherwise if a filename f is
!! passed in, this file is opened; otherwise the coefficient filename is
!! created from the instrument triplet and opened. The instrument triplet IDs
!! are defined in src/main/rttov_const.F90. Only one lun, f or instrument
!! needs to be supplied. If the file format is not specified in form the code
!! determines the format automatically for file input. The form argument is
!! mandatory for output files. The path argument can be used with the instrument
!! argument to specify the directory containing the coefficient file (if path is
!! omitted the file - or a symbolic link to it - must be in the current
!! directory). Formatted/unformatted files are opened, but HDF5 files are not.
!! The code which determines a free logical unit is not thread-safe.
!! @param[out]     err             status on exit
!! @param[out]     lun1            logical unit file was opened with
!! @param[out]     f1              filename of the opened file
!! @param[out]     form1           format (unformatted/formatted/hdf5) of the opened file
!! @param[in]      loutput         FALSE => open for input, TRUE => open for output
!! @param[in]      filetype        should be one of 'rtcoef', 'scaercoef', 'sccldcoef', 'pccoef'
!! @param[in]      lun             if supplied this is copied to lun1 (only use with formatted/unformatted files), optional
!! @param[in]      f               filename (including path) of file to open, optional
!! @param[in]      form            format (unformatted/formatted/hdf5) of file to open, optional
!! @param[in]      instrument      (platform_id, sat_id, instrument_id) for coef file to open, optional
!! @param[in]      path            directory containing coefficient file, for use with instrument argument, optional
!! @param[in]      suffix          optional suffix to append to filename, for use with instrument argument, optional
SUBROUTINE getlun(err, lun1, f1, form1, loutput, filetype, lun, f, form, instrument, path, suffix)

  USE rttov_const, ONLY : rttov_magic_string

  INTEGER(KIND=jpim), INTENT(OUT)          :: err
  INTEGER(KIND=jpim), INTENT(OUT)          :: lun1
  CHARACTER(LEN=*),   INTENT(OUT)          :: f1
  CHARACTER(LEN=*),   INTENT(OUT)          :: form1
  LOGICAL(KIND=jplm), INTENT(IN)           :: loutput
  CHARACTER(LEN=*),   INTENT(IN)           :: filetype
  INTEGER(KIND=jpim), INTENT(IN), OPTIONAL :: lun
  CHARACTER(LEN=*),   INTENT(IN), OPTIONAL :: f
  CHARACTER(LEN=*),   INTENT(IN), OPTIONAL :: form
  INTEGER(KIND=jpim), INTENT(IN), OPTIONAL :: instrument(3)
  CHARACTER(LEN=*),   INTENT(IN), OPTIONAL :: path
  CHARACTER(LEN=*),   INTENT(IN), OPTIONAL :: suffix

#include "rttov_cmpuc.interface"
#include "rttov_coeffname.interface"

  CHARACTER(LEN=8)   :: hdf5_format_signature
  ! DAR - removed kind=jplm from line below to comply with 2003 std which states
  !       that a logical must be of default type for an intrinsic (paraphrasing)
  LOGICAL            :: lopen, lexist         ! Use native KIND
  CHARACTER(LEN=4)   :: int4ch
  CHARACTER(LEN=8)   :: int8ch
  INTEGER            :: recl                  ! Use native KIND
  CHARACTER(LEN=16)  :: bin_check_string
  CHARACTER(LEN=8)   :: hdf5_check_signature
  CHARACTER(LEN=300) :: filestem

  TRY

  hdf5_format_signature = CHAR(137)//'HDF'//CHAR(13)//CHAR(10)//CHAR(26)//CHAR(10)

  form1 = ''
  IF (PRESENT(form)) form1 = form

  ! If lun was passed as argument, then take it and return
  IF (PRESENT(lun)) THEN
    lun1 = lun

    ! If form was passed as argument then take it, otherwise make a guess
    IF (PRESENT(form)) THEN
      IF (form /= '') RETURN ! form1 = form (see above)
    ENDIF
    INQUIRE(lun1, form = form1)
    RETURN
  ENDIF

  ! Try to guess the filename
  IF (PRESENT(f)) THEN

    ! If f was passed as argument, then take it
    f1 = f
    IF (loutput .AND. (form1 == '')) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0,'Format argument missing for output file')
    ENDIF

  ELSE IF (PRESENT(instrument)) THEN
    ! Guess the filename from the instrument triplet

    ! Get the filename from instrument triplet (excluding file extension)
    CALL rttov_coeffname(err, instrument, filetype, filestem)
    THROW(err.NE.0)
    IF (PRESENT(path)) filestem = TRIM(path)//'/'//TRIM(filestem)
    IF (PRESENT(suffix)) filestem = TRIM(filestem)//TRIM(suffix)

    IF (loutput .OR. form1 /= '') THEN

      ! Formatted output is the default
      IF (rttov_cmpuc('unformatted', form1)) THEN
        f1 = TRIM(filestem)//'.bin'
      ELSE IF (rttov_cmpuc('hdf5', form1)) THEN
        f1 = TRIM(filestem)//'.H5'
      ELSE
        form1 = 'formatted'
        f1 = TRIM(filestem)//'.dat'
      ENDIF

    ELSE

      ! Check for file existence using common file extensions
      ! File format is determined below (this way the code is more robust
      ! against users accidentally getting file extensions wrong)
      form1 = ''
      f1 = TRIM(filestem)//'.bin'
      INQUIRE(file = f1, exist = lexist)

      IF (.NOT. lexist) THEN
        f1 = TRIM(filestem)//'.h5'
        INQUIRE(file = f1, exist = lexist)
      ENDIF

      IF (.NOT. lexist) THEN
        f1 = TRIM(filestem)//'.H5'
        INQUIRE(file = f1, exist = lexist)
      ENDIF

      IF (.NOT. lexist) THEN
        f1 = TRIM(filestem)//'.dat'
        INQUIRE(file = f1, exist = lexist)
      ENDIF

      IF (.NOT. lexist) THEN
        err = errorstatus_fatal
        THROWM(err.NE.0,'Cannot find any coefficient file named: '//TRIM(filestem)//'.dat/.bin/.h5/.H5')
      ENDIF

    ENDIF
  ELSE
    err = errorstatus_fatal
    THROWM(err.NE.0,'Opening the coefficient file requires the filename or the instrument ID')
  ENDIF

  ! Look for a free logical unit (not thread-safe)
  DO lun1 = 9, 99
    INQUIRE(lun1, opened = lopen)
    IF (.NOT. lopen) EXIT
  ENDDO
  IF (lopen) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0,'Cannot find a free lun')
  ENDIF

  ! HDF5 files are not opened by this routine: just return file name,
  !   open will be done by load/save routines from rttov_hdf modules
  IF (loutput) THEN

    IF ( .NOT. rttov_cmpuc('hdf5', form1)) THEN
      OPEN(lun1, file = f1, form = form1, status = 'replace', action = 'write', iostat = err)
      THROWM(err.NE.0,'Cannot open '//TRIM(f1))
    ENDIF

  ELSE

    ! This prevents the test suite hanging if a scattering coefficient file is not specified
    ! (actually compiler-dependent: some compilers return exist=True for directories...)
    INQUIRE(file = f1, exist = lexist)
    IF (.NOT. lexist) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0,'File of filetype '//TRIM(filetype)//' does not exist: '//TRIM(f1))
    ENDIF

    IF (form1 == '') THEN
      ! Guess the input file format with the embedded magic string

      ! Do not just open as unformatted and try reading the magic string because ifort takes the
      ! first 4/8 bytes as the record length and may allocate a huge amount of memory for the read
      ! if the file is actually formatted.
      form1 = 'formatted'

      INQUIRE(iolength = recl) int4ch, bin_check_string
      OPEN(lun1, file = f1, form = 'unformatted', access = 'direct',recl = recl, &
                  status = 'old', action = 'read', iostat = err)
      THROWM(err.NE.0,'Cannot open '//TRIM(f1))
      READ(lun1, rec = 1, iostat = err) int4ch, bin_check_string
      CLOSE(lun1, iostat = err)
      THROW(err.NE.0)
      IF (bin_check_string == rttov_magic_string) form1 = 'unformatted'

      INQUIRE(iolength = recl) int8ch, bin_check_string
      OPEN(lun1, file = f1, form = 'unformatted', access = 'direct',recl = recl, &
                  status = 'old', action = 'read', iostat = err)
      THROWM(err.NE.0,'Cannot open '//TRIM(f1))
      READ(lun1, rec = 1, iostat = err) int8ch, bin_check_string
      CLOSE(lun1, iostat = err)
      THROW(err.NE.0)
      IF (bin_check_string == rttov_magic_string) form1 = 'unformatted'

      INQUIRE(iolength = recl) hdf5_check_signature
      OPEN(lun1, file = f1, form = 'unformatted', access = 'direct',recl = recl, &
                  status = 'old', action = 'read', iostat = err)
      THROWM(err.NE.0,'Cannot open '//TRIM(f1))
      READ(lun1, rec = 1, iostat = err) hdf5_check_signature
      CLOSE(lun1, iostat = err)
      THROW(err.NE.0)
      IF (hdf5_check_signature == hdf5_format_signature) form1 = 'hdf5'

    ENDIF

    IF (.NOT. rttov_cmpuc('hdf5', form1)) THEN
      OPEN(lun1, file = f1, form = form1, status = 'old', action = 'read', iostat = err)
      THROWM(err.NE.0,'Cannot open '//TRIM(f1))
    ENDIF

  ENDIF

  CATCH
END SUBROUTINE getlun

!> Close logical unit lun1. It is not closed if lun is not present: this is
!! useful to avoid closing a logical unit that was opened by a user.
!! param[out]   err       status on exit
!! param[in]    lun1      logical unit to close
!! param[in]    lun       lun1 is not closed if this argument is present, optional
SUBROUTINE closelun(err, lun1, lun)
  INTEGER(KIND=jpim), INTENT(OUT)          :: err
  INTEGER(KIND=jpim), INTENT(IN)           :: lun1
  INTEGER(KIND=jpim), INTENT(IN), OPTIONAL :: lun

  TRY
  IF (.NOT. PRESENT(lun)) THEN
    CLOSE(lun1, iostat = err)
    THROW(err.NE.0)
  ENDIF
  CATCH
END SUBROUTINE closelun

END MODULE
