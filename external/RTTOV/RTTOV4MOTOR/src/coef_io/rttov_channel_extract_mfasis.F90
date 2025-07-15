! Description:
!> @file
!!   Extract data for given channel list from an MFASIS LUT structure.
!
!> @brief
!!   Extract data for given channel list from an MFASIS LUT structure.
!!
!! @details
!!   This is used by HDF5 I/O code to read in a subset of channels from a
!!   coefficient file. The first coef argument contains the coefficients
!!   from the file. The second argument is an uninitialised structure
!!   which contains the extracted coefficients on exit.
!!
!!   NB The channel list MUST include at least one channel supported by
!!   the input MFASIS structure.
!!
!! @param[out]     err            status on exit
!! @param[in]      coef_mfasis1   input coefficients read from file
!! @param[in,out]  coef_mfasis2   output coefficients, uninitialised on entry
!! @param[in]      channels       list of channels to extract
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
!    Copyright 2017, EUMETSAT, All Rights Reserved.
!
SUBROUTINE rttov_channel_extract_mfasis(err, coef_mfasis1, coef_mfasis2, channels)
!INTF_OFF
#include "throw.h"
!INTF_ON
  USE parkind1, ONLY : jpim
  USE rttov_types, ONLY : rttov_coef_mfasis
!INTF_ON
  IMPLICIT NONE

  INTEGER(jpim),           INTENT(OUT)   :: err
  TYPE(rttov_coef_mfasis), INTENT(IN)    :: coef_mfasis1
  TYPE(rttov_coef_mfasis), INTENT(INOUT) :: coef_mfasis2
  INTEGER(jpim),           INTENT(IN)    :: channels(:)
!INTF_END

#include "rttov_errorreport.interface"
#include "rttov_channel_extract_sublist.interface"

  INTEGER(jpim)              :: j, lut_size
  INTEGER(jpim), ALLOCATABLE :: channels_mfasis(:)
  INTEGER(jpim), ALLOCATABLE :: mfasis_ext_index(:)
! ----------------------------------------------------------------------------

  TRY

  ! Work out the new channel list

  coef_mfasis2%nchannels_coef = SIZE(channels)

  ALLOCATE(channels_mfasis(coef_mfasis1%nchannels))

  ! Determine the solar channel numbers
  channels_mfasis(:) = coef_mfasis1%channel_list(:)

  ! Determine solar channels/phase functions to be extracted
  ALLOCATE(mfasis_ext_index(coef_mfasis1%nchannels))
  CALL rttov_channel_extract_sublist( &
        err,                            &
        channels_mfasis,                &
        channels,                       &
        coef_mfasis2%nchannels,         &
        coef_mfasis2%channel_list,      &
        coef_mfasis2%channel_lut_index, &
        mfasis_ext_index)
  THROW(err.NE.0)

  IF (coef_mfasis2%nchannels == 0) THEN
    DEALLOCATE(channels_mfasis, mfasis_ext_index)
    err = errorstatus_fatal
    THROWM(err.NE.0, 'No MFASIS channels to extract')
  ENDIF

  ! See rttov_read_ascii_mfasis_file.F90 for a description of what the various arrays contain.

  lut_size = PRODUCT(coef_mfasis1%lut_axes(:)%nvalues)

  ALLOCATE(coef_mfasis2%lut(coef_mfasis2%nchannels), stat=err)
  THROWM(err.NE.0, 'allocation of lut')
  DO j = 1, coef_mfasis2%nchannels
    coef_mfasis2%lut(j)%nluts = coef_mfasis1%lut(mfasis_ext_index(j))%nluts
    ALLOCATE(coef_mfasis2%lut(j)%qint(2,coef_mfasis2%lut(j)%nluts), stat=err)
    THROWM(err.NE.0, 'allocation of lut%qint')
    coef_mfasis2%lut(j)%qint = coef_mfasis1%lut(mfasis_ext_index(j))%qint
    ALLOCATE(coef_mfasis2%lut(j)%data(lut_size,coef_mfasis2%lut(j)%nluts), stat=err)
    THROWM(err.NE.0, 'allocation of lut%data')
    coef_mfasis2%lut(j)%data = coef_mfasis1%lut(mfasis_ext_index(j))%data
  ENDDO

  IF (ALLOCATED(channels_mfasis))  DEALLOCATE(channels_mfasis)
  IF (ALLOCATED(mfasis_ext_index)) DEALLOCATE(mfasis_ext_index)


  ! Everything else is the same...

  coef_mfasis2%file_type    = coef_mfasis1%file_type
  coef_mfasis2%version      = coef_mfasis1%version
  coef_mfasis2%readme_lut   = coef_mfasis1%readme_lut

  coef_mfasis2%clw_scheme   = coef_mfasis1%clw_scheme
  coef_mfasis2%ice_scheme   = coef_mfasis1%ice_scheme

  coef_mfasis2%maxzenangle  = coef_mfasis1%maxzenangle

  coef_mfasis2%nparticles   = coef_mfasis1%nparticles

  IF (coef_mfasis1%file_type == 2) THEN
    ALLOCATE(coef_mfasis2%aer_types(coef_mfasis2%nparticles), stat=err)
    THROWM(err.NE.0, 'allocation of aer_types')
    coef_mfasis2%aer_types  = coef_mfasis1%aer_types
  ELSE
    NULLIFY(coef_mfasis2%aer_types)
  ENDIF

  coef_mfasis2%ndims        = coef_mfasis1%ndims

  ALLOCATE(coef_mfasis2%lut_axes(coef_mfasis2%ndims), stat=err)
  THROWM(err.NE.0, 'allocation of lut_axes')
  DO j = 1, coef_mfasis1%ndims
    coef_mfasis2%lut_axes(j)%name     = coef_mfasis1%lut_axes(j)%name
    coef_mfasis2%lut_axes(j)%dim_type = coef_mfasis1%lut_axes(j)%dim_type
    coef_mfasis2%lut_axes(j)%nvalues  = coef_mfasis1%lut_axes(j)%nvalues
    ALLOCATE(coef_mfasis2%lut_axes(j)%values(coef_mfasis2%lut_axes(j)%nvalues), stat=err)
    THROWM(err.NE.0, 'allocation of lut_axes%values')
    coef_mfasis2%lut_axes(j)%values   = coef_mfasis1%lut_axes(j)%values
  ENDDO

  CATCH
END SUBROUTINE rttov_channel_extract_mfasis
