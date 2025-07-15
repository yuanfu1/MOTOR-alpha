! Description:
!> @file
!!   Read a binary MFASIS LUT file, optionally extracting a subset of channels.
!
!> @brief
!!   Read a binary MFASIS LUT file, optionally extracting a subset of channels.
!!
!! @details
!!   The file unit must be open when this subroutine is called.
!!
!!   Note that after reading a subset of channels RTTOV will identify them by
!!   indexes 1...SIZE(channels), not by the original channel numbers.
!!
!! @param[out]    err             status on exit
!! @param[in]     coef            RTTOV optical depth coefficient structure
!! @param[in,out] coef_mfasis     MFASIS cloud or aerosol coefficient structure
!! @param[in]     file_id         logical unit for input MFASIS file
!! @param[in]     channels        list of channels to read, optional
!! @param[in]     file_type       if present a check is done to ensure MFASIS file is of
!!                                  the specified type (1=>cloud, 2=>aerosol), optional
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
SUBROUTINE rttov_read_binary_mfasis_file( &
              err,           &
              coef,          &
              coef_mfasis,   &
              file_id,       &
              channels,      &
              file_type)
!INTF_OFF
#include "throw.h"
!INTF_ON
  USE rttov_types, ONLY : rttov_coef, rttov_coef_mfasis
  USE parkind1, ONLY : jpim
!INTF_OFF
  USE rttov_const, ONLY : rttov_magic_string, rttov_magic_number, &
                          mfasislut_version_compatible_min,       &
                          mfasislut_version_compatible_max,       &
                          mfasis_maxzenangle
  USE parkind1, ONLY : jprb, jplm
!INTF_ON
  IMPLICIT NONE

  INTEGER(KIND=jpim),      INTENT(OUT)          :: err
  TYPE(rttov_coef),        INTENT(IN)           :: coef
  TYPE(rttov_coef_mfasis), INTENT(INOUT)        :: coef_mfasis
  INTEGER(KIND=jpim),      INTENT(IN)           :: file_id
  INTEGER(KIND=jpim),      INTENT(IN), OPTIONAL :: channels(:)
  INTEGER(KIND=jpim),      INTENT(IN), OPTIONAL :: file_type
!INTF_END
#include "rttov_errorreport.interface"
#include "rttov_nullify_coef_mfasis.interface"
#include "rttov_channel_extract_sublist.interface"

  INTEGER(KIND=jpim) :: nchannels_mfasis
  LOGICAL(KIND=jplm) :: all_channels
  INTEGER(KIND=jpim) :: i, j, lut_size, nluts

  INTEGER(KIND=jpim), ALLOCATABLE :: channels_rtcoef(:)       ! List of rtcoef channels
  INTEGER(KIND=jpim), ALLOCATABLE :: channels_mfasis(:)       ! List of MFASIS channels
  INTEGER(KIND=jpim), ALLOCATABLE :: mfasis_ext_index(:)      ! Indexes of LUTs to extract
  REAL(KIND=jprb),    ALLOCATABLE :: temp(:)

  CHARACTER(LEN=16) :: bin_check_string
  REAL(KIND=jprb)   :: bin_check_number
  REAL(KIND=jprb)   :: bin_check_value
!- End of header --------------------------------------------------------
  TRY

  all_channels = .NOT. PRESENT(channels)

  CALL rttov_nullify_coef_mfasis(coef_mfasis)

  READ (file_id, iostat=err) bin_check_string, bin_check_number
  THROWM(err.NE.0,'io status while reading header')

  ! Verification of header string
  IF (bin_check_string /= rttov_magic_string) err = errorstatus_fatal
  THROWM(err.NE.0,'Wrong header string in file')

  ! Verification of single/double precision using a 5 digit number
  ! with exponent 12, which is always Ok for single precision
  bin_check_value = 1._jprb - ABS(bin_check_number - rttov_magic_number)
  IF (bin_check_value > 1.01_jprb .OR. bin_check_value < 0.99_jprb) err = errorstatus_fatal
  THROWM(err.NE.0,'File created with a different real precision (R4<->R8)')


  READ (file_id, iostat=err)coef_mfasis%file_type,      &
                            coef_mfasis%version,        &
                            nchannels_mfasis,           &
                            coef_mfasis%nchannels_coef, &
                            coef_mfasis%nparticles,     &
                            coef_mfasis%ndims
  THROWM(err.NE.0, 'reading MFASIS file dimensions')

  IF (PRESENT(file_type)) THEN
    IF (coef_mfasis%file_type /= file_type) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0, 'Error: MFASIS file type mismatch (cloud/aerosol)')
    ENDIF
  ENDIF

  IF (coef_mfasis%version < mfasislut_version_compatible_min .OR. &
      coef_mfasis%version > mfasislut_version_compatible_max) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, "Version of MFASIS LUT file is incompatible with RTTOV library")
  ENDIF

  IF (coef%fmv_ori_nchn /= coef_mfasis%nchannels_coef) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, "Incompatible channels between rtcoef and MFASIS files")
  ENDIF

  IF (.NOT. all_channels) THEN
    ALLOCATE(channels_rtcoef(SIZE(channels)))
    channels_rtcoef = channels
  ELSE
    ALLOCATE(channels_rtcoef(coef%fmv_chn))
    channels_rtcoef = (/(i, i = 1, coef%fmv_chn)/)
  ENDIF

  ! Take care of the user list of channels
  ! coef_mfasis%nchannels_coef is the total number of channels that the user requests

  IF (.NOT. all_channels) coef_mfasis%nchannels_coef = SIZE(channels)

  ALLOCATE(channels_mfasis(nchannels_mfasis))
  READ (file_id, iostat=err)channels_mfasis(:)
  THROWM(err.NE.0, 'reading MFASIS channel list')

  ! Determine MFASIS channels to be extracted
  ALLOCATE(mfasis_ext_index(nchannels_mfasis))
  CALL rttov_channel_extract_sublist( &
          err,                           &
          channels_mfasis,               &
          channels_rtcoef,               &
          coef_mfasis%nchannels,         &
          coef_mfasis%channel_list,      &
          coef_mfasis%channel_lut_index, &
          mfasis_ext_index)
  THROW(err.NE.0)

  IF (coef_mfasis%nchannels == 0) THEN
    DEALLOCATE(channels_rtcoef, channels_mfasis, mfasis_ext_index)
    err = errorstatus_fatal
    THROWM(err.NE.0, 'No MFASIS channels to extract')
  ENDIF

  ! See rttov_read_ascii_mfasis_file.F90 for a description of what the various arrays contain.

  IF (coef_mfasis%file_type == 2) THEN
    coef_mfasis%clw_scheme = 0
    coef_mfasis%ice_scheme = 0

    ALLOCATE(coef_mfasis%aer_types(coef_mfasis%nparticles), STAT = err)
    THROWM(err.NE.0, "allocation of coef_mfasis%aer_types array")

    READ (file_id, iostat=err)coef_mfasis%aer_types(:)
    THROWM(err.NE.0, 'while reading MFASIS aer_types')
  ELSE
    NULLIFY(coef_mfasis%aer_types)

    READ (file_id, iostat=err)coef_mfasis%clw_scheme, coef_mfasis%ice_scheme
    THROWM(err.NE.0, 'while reading MFASIS clw_scheme and ice_scheme')
  ENDIF

  IF (coef_mfasis%version > 0) THEN
    READ (file_id, iostat=err)coef_mfasis%maxzenangle
    THROWM(err.NE.0, 'while reading MFASIS maxzenangle')
  ELSE ! coef_mfasis%version = 0
    coef_mfasis%maxzenangle = mfasis_maxzenangle
  ENDIF

  READ (file_id, iostat=err)coef_mfasis%readme_lut

  ALLOCATE(coef_mfasis%lut_axes(coef_mfasis%ndims), STAT = err)
  THROWM(err.NE.0, "allocation of coef_mfasis%lut_axes array")

  DO j = 1, coef_mfasis%ndims
    READ (file_id, iostat=err)coef_mfasis%lut_axes(j)%name
    THROWM(err.NE.0, 'while reading LUT axis name')

    READ (file_id,  iostat=err)coef_mfasis%lut_axes(j)%dim_type
    THROWM(err.NE.0, 'while reading LUT axis dim_type')

    READ (file_id, iostat=err)coef_mfasis%lut_axes(j)%nvalues
    THROWM(err.NE.0, 'while reading LUT axis nvalues')

    ALLOCATE(coef_mfasis%lut_axes(j)%values(coef_mfasis%lut_axes(j)%nvalues), STAT = err)
    THROWM(err.NE.0, "allocation of coef_mfasis%lut_axes%values array")

    READ (file_id, iostat=err)coef_mfasis%lut_axes(j)%values(:)
    THROWM(err.NE.0, 'while reading LUT axis values')
  ENDDO


  ALLOCATE(coef_mfasis%lut(coef_mfasis%nchannels), STAT = err)
  THROWM(err.NE.0, "allocation of coef_mfasis%lut array")

  lut_size = PRODUCT(coef_mfasis%lut_axes(:)%nvalues)

  i = 1
  DO j = 1, nchannels_mfasis
    IF (j == mfasis_ext_index(i)) THEN
      READ (file_id, iostat=err)coef_mfasis%lut(i)%nluts
      THROWM(err.NE.0, 'while reading LUT data')

      ALLOCATE(coef_mfasis%lut(i)%qint(2,coef_mfasis%lut(i)%nluts), STAT = err)
      THROWM(err.NE.0, "allocation of coef_mfasis%lut%qint array")

      READ (file_id, iostat=err)coef_mfasis%lut(i)%qint(:,:)
      THROWM(err.NE.0, 'while reading LUT qint values')

      ALLOCATE(coef_mfasis%lut(i)%data(lut_size,coef_mfasis%lut(i)%nluts), STAT = err)
      THROWM(err.NE.0, "allocation of coef_mfasis%lut%data array")

      READ (file_id, iostat=err)coef_mfasis%lut(i)%data(:,:)
      THROWM(err.NE.0, 'while reading LUT data')
      i = i + 1
    ELSE
      READ (file_id, iostat=err)nluts
      THROWM(err.NE.0, 'while reading LUT data')

      ALLOCATE(temp(lut_size * nluts))
      READ (file_id, iostat=err)temp(1:nluts)  ! Read qint array
      THROWM(err.NE.0, 'while reading LUT data')

      READ (file_id, iostat=err)temp(:)        ! Read data array
      THROWM(err.NE.0, 'while reading LUT data')
      DEALLOCATE(temp)
    ENDIF
  ENDDO

  IF (ALLOCATED(mfasis_ext_index)) DEALLOCATE(mfasis_ext_index)
  IF (ALLOCATED(channels_mfasis))  DEALLOCATE(channels_mfasis)

  CATCH
END SUBROUTINE rttov_read_binary_mfasis_file
