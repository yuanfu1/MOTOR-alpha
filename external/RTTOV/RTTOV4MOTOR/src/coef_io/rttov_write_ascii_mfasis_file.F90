! Description:
!> @file
!!   Write an ASCII MFASIS LUT file.
!
!> @brief
!!   Write an ASCII MFASIS LUT file.
!!
!! @details
!!   The file unit must be open when this subroutine is called.
!!
!! @param[out]    err             status on exit
!! @param[in]     coef            RTTOV optical depth coefficient structure
!! @param[in,out] coef_mfasis     MFASIS cloud or aerosol coefficient structure
!! @param[in]     file_id         logical unit for output MFASIS file
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
SUBROUTINE rttov_write_ascii_mfasis_file( &
              err,           &
              coef,          &
              coef_mfasis,   &
              file_id,       &
              verbose)
!INTF_OFF
#include "throw.h"
!INTF_ON
  USE rttov_types, ONLY : rttov_coef, rttov_coef_mfasis
  USE parkind1, ONLY : jpim, jplm
!INTF_ON
  IMPLICIT NONE

  INTEGER(KIND=jpim),      INTENT(OUT)          :: err
  TYPE(rttov_coef),        INTENT(IN)           :: coef
  TYPE(rttov_coef_mfasis), INTENT(IN)           :: coef_mfasis
  INTEGER(KIND=jpim),      INTENT(IN)           :: file_id
  LOGICAL(KIND=jplm),      INTENT(IN), OPTIONAL :: verbose
!INTF_END
#include "rttov_errorreport.interface"

  INTEGER(KIND=jpim) :: j
  LOGICAL(KIND=jplm) :: lverbose
  CHARACTER(LEN=32)  :: section
  CHARACTER(LEN=80)  :: errMessage
  CHARACTER(LEN=*), PARAMETER :: routinename = 'rttov_write_ascii_mfasis_file'
!- End of header --------------------------------------------------------
  TRY
  IF (PRESENT(verbose)) THEN
    lverbose = verbose
  ELSE
    lverbose = .TRUE._jplm
  END IF

  IF (lverbose) THEN
    WRITE (errMessage, '( "write coefficient to file_id ", i2, " in ASCII format")') file_id
    INFO(errMessage)
  END IF

  WRITE (file_id, '(a)', iostat=err) ' ! RTTOV coefficient file '//TRIM(coef%id_common_name)
  THROW(err.NE.0)

  WRITE (file_id, '(a)', iostat=err) ' ! Automatic creation by subroutine '//routinename
  THROW(err.NE.0)

  WRITE (file_id, '(a)', iostat=err) ' ! ------------------------------------------------------'
  THROW(err.NE.0)


  section = 'MFASIS_GENERAL'
  WRITE (file_id, '(a)', iostat=err) ' ! ------------------------------------------------------'
  THROW(err.NE.0)

  WRITE (file_id, '(a)', iostat=err) TRIM(section)
  THROW(err.NE.0)

  WRITE (file_id, '(a)', iostat=err) ' !'
  THROW(err.NE.0)

  WRITE(file_id,'(1x,i8,t20,a)', iostat=err) &
    coef_mfasis%file_type, '! MFASIS file type: 1 => clouds, 2 => aerosols'
  THROW(err.NE.0)

  WRITE(file_id,'(1x,i8,t20,a)', iostat=err) &
    coef_mfasis%version, '! MFASIS file version'
  THROW(err.NE.0)

  WRITE(file_id,'(1x,i8,t20,a)', iostat=err) &
    coef_mfasis%nchannels_coef, '! Number of channels in associated rtcoef file'
  THROW(err.NE.0)

  WRITE(file_id,'(1x,i8,t20,a)', iostat=err) &
    coef_mfasis%nchannels, '! Number of channels supported by MFASIS'
  THROW(err.NE.0)

  WRITE(file_id,'(1x,a)', iostat=err) '! Channels for which MFASIS LUTs are stored'
  THROW(err.NE.0)
  WRITE(file_id,'(10i6)', iostat=err) coef_mfasis%channel_list(:)
  THROW(err.NE.0)

  WRITE(file_id,'(1x,i8,t20,a)', iostat=err) &
    coef_mfasis%nparticles, '! Number of particle types'
  THROW(err.NE.0)

  IF (coef_mfasis%file_type == 2) THEN
    WRITE(file_id,'(1x,a)', iostat=err) '! Aerosol particle types'
    THROW(err.NE.0)
    WRITE(file_id,'(10i6)', iostat=err) coef_mfasis%aer_types(:)
    THROW(err.NE.0)
  ELSE
    WRITE(file_id,'(1x,a)', iostat=err) '! Cloud liquid and ice water schemes used for training'
    THROW(err.NE.0)
    WRITE(file_id,'(10i6)', iostat=err) coef_mfasis%clw_scheme, coef_mfasis%ice_scheme
    THROW(err.NE.0)
  ENDIF

  IF (coef_mfasis%version > 0) THEN
    WRITE(file_id,'(1x,a)', iostat=err) '! Maximum zenith angle used for training'
    THROW(err.NE.0)
    WRITE(file_id,'(f10.2)', iostat=err) coef_mfasis%maxzenangle
    THROW(err.NE.0)
  ENDIF

  IF (coef_mfasis%readme_lut(1) .NE. 'xxxx') THEN
    section = 'README_MFASIS'
    WRITE (file_id, '(a)', iostat=err) ' ! ------------------------------------------------------'
    THROW(err.NE.0)

    WRITE (file_id, '(a)', iostat=err) TRIM(section)
    THROW(err.NE.0)

    WRITE (file_id, '(a)', iostat=err) ' !'
    THROW(err.NE.0)

    DO j = 1, SIZE(coef_mfasis%readme_lut)
      IF (coef_mfasis%readme_lut(j) .EQ. 'xxxx') EXIT
      WRITE (file_id, '(a)', iostat=err)TRIM(coef_mfasis%readme_lut(j))
      THROW(err.NE.0)
    ENDDO
  ENDIF


  section = 'MFASIS_DIMS'
  WRITE (file_id, '(a)', iostat=err) ' ! ------------------------------------------------------'
  THROW(err.NE.0)

  WRITE (file_id, '(a)', iostat=err) TRIM(section)
  THROW(err.NE.0)

  WRITE (file_id, '(a)', iostat=err) ' !'
  THROW(err.NE.0)

  WRITE(file_id,'(1x,i8,t20,a)', iostat=err) &
    coef_mfasis%ndims, '! Number of dimensions in each LUT'
  THROW(err.NE.0)

  WRITE(file_id,'(1x,a)', iostat=err) '! Type of each dimension'
  THROW(err.NE.0)
  WRITE(file_id,'(1x,a)', iostat=err) '! 0 : theta+ Fourier index (k)'
  THROW(err.NE.0)
  WRITE(file_id,'(1x,a)', iostat=err) '! 1 : theta- Fourier index (l)'
  THROW(err.NE.0)
  WRITE(file_id,'(1x,a)', iostat=err) '! 2 : Albedo (no interp., use Jonkheid 2012)'
  THROW(err.NE.0)
  WRITE(file_id,'(1x,a)', iostat=err) '! 3 : tau-like (interpolate in logarithm)'
  THROW(err.NE.0)
  WRITE(file_id,'(1x,a)', iostat=err) '! 4 : Reff-like (interpolate linearly)'
  THROW(err.NE.0)
  WRITE(file_id,'(1x,a)', iostat=err) '! 5 : scattering angle (constr. lin. interp.)'
  THROW(err.NE.0)

  WRITE(file_id,'(1x,a)', iostat=err) '! For each dimension: name, dim_type, number of values, value list'
  THROW(err.NE.0)
  DO j = 1, coef_mfasis%ndims
    WRITE(file_id,'(a)', iostat=err) TRIM(coef_mfasis%lut_axes(j)%name)
    THROW(err.NE.0)

    WRITE(file_id,'(1x,i8)', iostat=err) coef_mfasis%lut_axes(j)%dim_type
    THROW(err.NE.0)

    WRITE(file_id,'(1x,i8)', iostat=err)coef_mfasis%lut_axes(j)%nvalues
    THROW(err.NE.0)

    WRITE(file_id,'(5e16.8)', iostat=err) coef_mfasis%lut_axes(j)%values(:)
    THROW(err.NE.0)
  ENDDO


  section = 'MFASIS_LUTS'
  WRITE (file_id, '(a)', iostat=err) ' ! ------------------------------------------------------'
  THROW(err.NE.0)

  WRITE (file_id, '(a)', iostat=err) TRIM(section)
  THROW(err.NE.0)

  WRITE (file_id, '(a)', iostat=err) ' !'
  THROW(err.NE.0)

  DO j = 1, coef_mfasis%nchannels
    WRITE (file_id, '(a,i6)', iostat=err) '! nluts, q values, LUT(s) for channel ', coef_mfasis%channel_list(j)

    WRITE(file_id,'(1x,i6)', iostat=err) coef_mfasis%lut(j)%nluts
    THROW(err.NE.0)

    WRITE(file_id,'(10e16.8)', iostat=err) coef_mfasis%lut(j)%qint(:,:)
    THROW(err.NE.0)

    WRITE(file_id,'(10e16.8)', iostat=err) coef_mfasis%lut(j)%data(:,:)
    THROW(err.NE.0)
  ENDDO

  IF (lverbose) INFO("end of write coefficient")

  CATCH
END SUBROUTINE rttov_write_ascii_mfasis_file
