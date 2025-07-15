! Description:
!> @file
!!   Write a binary MFASIS LUT file.
!
!> @brief
!!   Write a binary MFASIS LUT file.
!!
!! @details
!!   The file unit must be open when this subroutine is called.
!!
!! @param[out]    err             status on exit
!! @param[in,out] coef_mfasis     MFASIS cloud or aerosol coefficient structure
!! @param[in]     file_id         logical unit for input MFASIS file
!! @param[in]     verbose         flag to switch verbose output on/off (default TRUE), optional
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
SUBROUTINE rttov_write_binary_mfasis_file( &
              err,           &
              coef_mfasis,   &
              file_id,       &
              verbose)
!INTF_OFF
#include "throw.h"
!INTF_ON
  USE rttov_types, ONLY : rttov_coef_mfasis
  USE parkind1, ONLY : jpim, jplm
!INTF_OFF
  USE rttov_const, ONLY : rttov_magic_string, rttov_magic_number
!INTF_ON
  IMPLICIT NONE

  INTEGER(KIND=jpim),      INTENT(OUT)          :: err
  TYPE(rttov_coef_mfasis), INTENT(IN)           :: coef_mfasis
  INTEGER(KIND=jpim),      INTENT(IN)           :: file_id
  LOGICAL(KIND=jplm),      INTENT(IN), OPTIONAL :: verbose
!INTF_END
#include "rttov_errorreport.interface"

  INTEGER(KIND=jpim) :: j
  LOGICAL(KIND=jplm) :: lverbose
  CHARACTER(LEN=80)  :: errMessage
!- End of header --------------------------------------------------------
  TRY
  IF (PRESENT(verbose)) THEN
    lverbose = verbose
  ELSE
    lverbose = .TRUE._jplm
  ENDIF

  IF (lverbose) THEN
    WRITE (errMessage, '( "write coefficient to file_id ", i2, " in binary format")') file_id
    INFO(errMessage)
  ENDIF

  WRITE (file_id, iostat=err) rttov_magic_string, rttov_magic_number
  THROW(err.NE.0)

  WRITE (file_id, iostat=err)coef_mfasis%file_type,      &
                             coef_mfasis%version,        &
                             coef_mfasis%nchannels,      &
                             coef_mfasis%nchannels_coef, &
                             coef_mfasis%nparticles,     &
                             coef_mfasis%ndims
  THROW(err.NE.0)

  WRITE (file_id, iostat=err)coef_mfasis%channel_list(:)
  THROW(err.NE.0)

  IF (coef_mfasis%file_type == 2) THEN
    WRITE (file_id, iostat=err)coef_mfasis%aer_types(:)
    THROW(err.NE.0)
  ELSE
    WRITE (file_id, iostat=err)coef_mfasis%clw_scheme, coef_mfasis%ice_scheme
    THROW(err.NE.0)
  ENDIF

  IF (coef_mfasis%version > 0) THEN
    WRITE(file_id, iostat=err) coef_mfasis%maxzenangle
    THROW(err.NE.0)
  ENDIF

  WRITE (file_id, iostat=err)coef_mfasis%readme_lut

  DO j = 1, coef_mfasis%ndims
    WRITE (file_id, iostat=err)coef_mfasis%lut_axes(j)%name
    THROW(err.NE.0)

    WRITE (file_id, iostat=err)coef_mfasis%lut_axes(j)%dim_type
    THROW(err.NE.0)

    WRITE (file_id, iostat=err)coef_mfasis%lut_axes(j)%nvalues
    THROW(err.NE.0)

    WRITE (file_id, iostat=err)coef_mfasis%lut_axes(j)%values(:)
    THROW(err.NE.0)
  ENDDO

  DO j = 1, coef_mfasis%nchannels
    WRITE (file_id, iostat=err)coef_mfasis%lut(j)%nluts
    THROW(err.NE.0)

    WRITE (file_id, iostat=err)coef_mfasis%lut(j)%qint(:,:)
    THROW(err.NE.0)

    WRITE (file_id, iostat=err)coef_mfasis%lut(j)%data(:,:)
    THROW(err.NE.0)
  ENDDO

  CATCH
END SUBROUTINE rttov_write_binary_mfasis_file
