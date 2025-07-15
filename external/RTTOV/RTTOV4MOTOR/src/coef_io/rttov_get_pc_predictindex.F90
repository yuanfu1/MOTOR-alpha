! Description:
!> @file
!!   Obtain the channel indices of the PC-RTTOV predictor channels
!!   for a given PC coefficient file and the specified PC regression
!!   set.
!
!> @brief
!!   Obtain the channel indices of the PC-RTTOV predictor channels.
!!
!! @details
!!   PC-RTTOV uses standard RTTOV radiances for a subset of instrument
!!   channels to calculate the PC scores for the whole spectrum. The
!!   chanprof array must be populated with the indices of these specific
!!   predictor channels. This subroutine can be used to obtain the
!!   channel indices from a PC coefficient file given the selected
!!   regression set and band specified in opts%rt_ir%pc%ipcreg and
!!   opts%rt_ir%pc%ipcbnd respectively.
!!
!!   The predictindex argument should be passed in as an unallocated
!!   pointer. If the subroutine runs successfully then on exit the array
!!   will have been allocated and will contain the channel indices.
!!
!!   The pccoef file can be read either by specifying filenames explicitly
!!   (the recommended method), by opening the file(s) and passing the
!!   relevant logical unit number(s), or by specifying the instrument
!!   ID triplet so that RTTOV can construct the filenames (see
!!   rttov_const.F90 for the platform and instrument IDs).
!!
!!   RTTOV can normally detect the format of coefficient files so
!!   the file format argument is not usually required.
!!
!!
!! @param[out]    err             status on exit
!! @param[in]     opts            options to configure the simulations
!! @param[out]    predictindex    unallocated pointer, allocated on exit containing the predictor channel indices
!! @param[in]     form_pccoef     format of PC-RTTOV pccoef file, optional
!! @param[in]     file_pccoef     file name of PC-RTTOV pccoef file, optional
!! @param[in]     file_id_pccoef  file name of PC-RTTOV pccoef file, optional
!! @param[in]     instrument      (platform,satellite,instrument) ID triplet, optional
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
SUBROUTINE rttov_get_pc_predictindex( &
              err,           &
              opts,          &
              predictindex,  &
              form_pccoef,   &
              file_pccoef,   &
              file_id_pccoef,&
              instrument)
!INTF_OFF
#include "throw.h"
!INTF_ON
  USE rttov_types, ONLY : rttov_options
  USE parkind1, ONLY : jpim
!INTF_OFF
  USE rttov_const, ONLY :    &
         rttov_magic_string, &
         rttov_magic_number, &
         lensection,         &
         errorstatus_fatal  ! PGI v2013 complains if this is not explicitly USEd
  USE parkind1, ONLY : jprb, jplm
  USE rttov_coef_io_mod, ONLY : getlun, closelun
#ifdef _RTTOV_HDF
  USE rttov_hdf_mod, ONLY : open_hdf, close_hdf
#endif
!INTF_ON
  IMPLICIT NONE
  INTEGER(KIND=jpim),  INTENT(OUT)          :: err
  TYPE(rttov_options), INTENT(IN)           :: opts
  INTEGER(KIND=jpim),  POINTER              :: predictindex(:)
  CHARACTER(LEN=*),    INTENT(IN), OPTIONAL :: form_pccoef
  CHARACTER(LEN=*),    INTENT(IN), OPTIONAL :: file_pccoef
  INTEGER(KIND=jpim),  INTENT(IN), OPTIONAL :: file_id_pccoef
  INTEGER(KIND=jpim),  INTENT(IN), OPTIONAL :: instrument(3)
!INTF_END

#include "rttov_errorreport.interface"
#include "rttov_skipcommentline.interface"
#include "rttov_cmpuc.interface"
#include "rttov_findnextsection.interface"

#ifdef _RTTOV_HDF
#include "rttov_hdf_load.interface"
#endif

  INTEGER(KIND=jpim) :: file_id_pccoef1
  CHARACTER(LEN=32)  :: form_pccoef1
  CHARACTER(LEN=256) :: file_pccoef1

  INTEGER(KIND=jpim) :: io_status
  INTEGER(KIND=jpim) :: i, m, n

  CHARACTER(LEN=16)  :: bin_check_string
  REAL(KIND=jprb)    :: bin_check_number
  REAL(KIND=jprb)    :: bin_check_value

  INTEGER(KIND=jpim) :: fmv_pc_comp_pc, fmv_pc_cld
  INTEGER(KIND=jpim) :: fmv_pc_bands, fmv_pc_msets
  INTEGER(KIND=jpim) :: fmv_pc_mnum, fmv_pc_mchn
  INTEGER(KIND=jpim) :: fmv_pc_sets, fmv_pc_npred

  CHARACTER(LEN=64)  :: dname
  CHARACTER(LEN=lensection) :: section
!- End of header --------------------------------------------------------

  TRY

  IF (opts%rt_ir%addclouds) THEN
    CALL getlun(err, file_id_pccoef1, file_pccoef1, form_pccoef1, .FALSE._jplm, "pccoefcld",&
        file_id_pccoef, file_pccoef, form_pccoef, instrument )
  ELSE
    CALL getlun(err, file_id_pccoef1, file_pccoef1, form_pccoef1, .FALSE._jplm, "pccoef",&
        file_id_pccoef, file_pccoef, form_pccoef, instrument )
  ENDIF
  THROW(err .NE. 0)

! Read the file
  IF (rttov_cmpuc(form_pccoef1, "formatted")) THEN
    readascii : DO
      CALL rttov_findnextsection(file_id_pccoef1, io_status, section)
      IF (io_status < 0) EXIT ! end-of-file

      SELECT CASE (TRIM(section))
      CASE ('PRINCOMP_PREDICTORS')
        CALL rttov_skipcommentline(file_id_pccoef1, err)
        THROWM(err .NE. 0, 'io status while reading section '//section)

        READ (file_id_pccoef1, *, iostat=err)fmv_pc_comp_pc
        THROWM(err .NE. 0, 'io status while reading section '//section)

        READ (file_id_pccoef1, *, iostat=err)fmv_pc_cld
        THROWM(err .NE. 0, 'io status while reading section '//section)

        READ (file_id_pccoef1, *, iostat=err)fmv_pc_bands
        THROWM(err .NE. 0, 'io status while reading section '//section)

        READ (file_id_pccoef1, *, iostat=err)fmv_pc_msets
        THROWM(err .NE. 0, 'io status while reading section '//section)

        READ (file_id_pccoef1, *, iostat=err)fmv_pc_mnum
        THROWM(err .NE. 0, 'io status while reading section '//section)

        READ (file_id_pccoef1, *, iostat=err)fmv_pc_mchn
        THROWM(err .NE. 0, 'io status while reading section '//section)

        ! Loop on bands
        DO m = 1, fmv_pc_bands
          READ (file_id_pccoef1, *, iostat=err)fmv_pc_sets
          THROWM(err .NE. 0, 'io status while reading section '//section)

          ! Loop on predictor sets
          DO n = 1,  fmv_pc_sets
            READ (file_id_pccoef1, *, iostat=err)fmv_pc_npred
            THROWM(err .NE. 0, 'io status while reading section '//section)

            CALL rttov_skipcommentline(file_id_pccoef1, err)
            THROWM(err .NE. 0, 'io status while reading section '//section)

            ALLOCATE (predictindex(fmv_pc_npred), STAT = err)
            THROWM(err .NE. 0, "allocation of predictindex array")

            READ (file_id_pccoef1, *, iostat=err)(predictindex(i), i = 1, fmv_pc_npred)
            THROWM(err .NE. 0, "io status while reading predictindex array")

            IF (opts%rt_ir%pc%ipcreg == n .AND. opts%rt_ir%pc%ipcbnd == m) EXIT readascii

            DEALLOCATE (predictindex, STAT = err)
            THROW(err .NE. 0)
          ENDDO
        ENDDO
        EXIT readascii

      CASE ('END')
        EXIT readascii
      CASE DEFAULT
        CYCLE readascii
      END SELECT

    ENDDO readascii

  ELSE IF (rttov_cmpuc(form_pccoef1,"unformatted")) THEN
    READ (file_id_pccoef1, iostat=err)bin_check_string, bin_check_number
    THROW(err .NE. 0)

    ! Verification of header string
    IF (bin_check_string /= rttov_magic_string) err = errorstatus_fatal
    THROWM(err .NE. 0, 'Wrong header string in file')

    ! Verification of single/double precision using a 5 digit number
    ! with exponent 12, which is always OK for single precision
    bin_check_value = 1._jprb - ABS(bin_check_number - rttov_magic_number)
    IF (bin_check_value > 1.01_jprb .OR. bin_check_value < 0.99_jprb) err = errorstatus_fatal
    THROWM(err .NE. 0, 'File created with a different real precision (R4<->R8)')

    READ (file_id_pccoef1, iostat=err)fmv_pc_comp_pc
    THROWM(err .NE. 0, "reading fmv_pc_comp_pc")

    READ (file_id_pccoef1, iostat=err)fmv_pc_cld
    THROWM(err .NE. 0, "reading fmv_pc_cld")

    READ (file_id_pccoef1, iostat=err)fmv_pc_bands
    THROWM(err .NE. 0, "reading fmv_pc_bands")

    READ (file_id_pccoef1, iostat=err)fmv_pc_msets
    THROWM(err .NE. 0, "reading fmv_pc_msets")

    READ (file_id_pccoef1, iostat=err)fmv_pc_mnum
    THROWM(err .NE. 0, "reading fmv_pc_mnum")

    READ (file_id_pccoef1, iostat=err)fmv_pc_mchn
    THROWM(err .NE. 0, "reading fmv_pc_mchn")

    ! Loop on bands
    readbin: DO m = 1, fmv_pc_bands
      READ (file_id_pccoef1, iostat=err)fmv_pc_sets
      THROWM(err .NE. 0, "reading fmv_pc_sets")

      ! Loop on predictor sets
      DO n = 1, fmv_pc_sets
        READ (file_id_pccoef1, iostat=err)fmv_pc_npred
        THROWM(err .NE. 0, "reading fmv_pc_npred")

        ALLOCATE (predictindex(fmv_pc_npred), STAT = err)
        THROWM(err .NE. 0, "allocation of predictindex array")

        READ (file_id_pccoef1, iostat=err)(predictindex(i), i = 1, fmv_pc_npred)
        THROWM(err .NE. 0, "reading predictindex")

        IF (opts%rt_ir%pc%ipcreg == n .AND. opts%rt_ir%pc%ipcbnd == m) EXIT readbin

        DEALLOCATE (predictindex, STAT = err)
        THROW(err .NE. 0)
      ENDDO
    ENDDO readbin

  ELSE IF (rttov_cmpuc(form_pccoef1,"hdf5")) THEN
#ifndef _RTTOV_HDF
    err = errorstatus_fatal
    THROWM(err .NE. 0, "This program is not compiled with HDF5 capability; use RTTOV_HDF=1 with Makefile.PL")
#else
    CALL open_hdf(.FALSE., err) ! 32 bits - doesn't matter here
    THROW(err .NE. 0)
    WRITE(dname,"('/PC/PCREG/',i2.2,'/',i2.2)") opts%rt_ir%pc%ipcbnd, opts%rt_ir%pc%ipcreg
    CALL rttov_hdf_load(err, file_pccoef1, TRIM(dname), PI1 = predictindex, SNAME = 'PREDICTINDEX')
    THROW(err .NE. 0)
    CALL close_hdf( err )
    THROW(err .NE. 0)
#endif
  ELSE
    err = errorstatus_fatal
    THROWM(err .NE. 0, "Unknown format "//TRIM(form_pccoef1))
  ENDIF

  CALL closelun(err, file_id_pccoef1, file_id_pccoef)
  THROW(err .NE. 0)

  IF (.NOT. ASSOCIATED(predictindex)) err = errorstatus_fatal
  THROWM(err .NE. 0, "PC regression set not found")

  CATCH
END SUBROUTINE rttov_get_pc_predictindex
