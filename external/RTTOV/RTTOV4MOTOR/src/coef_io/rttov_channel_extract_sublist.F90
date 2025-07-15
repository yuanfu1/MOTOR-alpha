! Description
!> @file
!!   Given a sub-list of a larger channel list and a list of channels to
!!   extract from the larger list, returns the new sub-list containing the
!!   intersection of the original sub-list and the extracted channel list.
!
!> @brief
!!   Given a sub-list of a larger channel list and a list of channels to
!!   extract from the larger list, returns the new sub-list containing the
!!   intersection of the original sub-list and the extracted channel list.
!!
!! @details
!!   The channels_newsub(:) and index_newsub_ext(:) arguments should be
!!   unallocated pointers. They are allocated on exit of the subroutine unless
!!   nchannels_newsub is zero (i.e. there are no channels to extract from the
!!   sub-list).
!!
!!   The index_oldsub_newsub(:) array should have SIZE(channels_oldsub)
!!   elements. On exit the first nchannels_newsub elements contain the indices
!!   into the old sub-list of the new sub-list.
!!
!! @param[out]    err                   status on exit
!! @param[in]     channels_oldsub       channels in original sub-list
!! @param[in]     channels_ext          list of channels to extract
!! @param[in,out] nchannels_newsub      number of channels in new sub-list
!! @param         channels_newsub       channels in new sub-list (indices into the extracted channel list)
!! @param         index_newsub_ext      indices in new sub-list for each extracted channel
!! @param[in,out] index_oldsub_newsub   indices in original sub-list for each channel in new sub-list (i.e. indices to extract)
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
SUBROUTINE rttov_channel_extract_sublist( &
              err,                &
              channels_oldsub,    &
              channels_ext,       &
              nchannels_newsub,   &
              channels_newsub,    &
              index_newsub_ext,   &
              index_oldsub_newsub)
#include "throw.h"

  USE parkind1, ONLY : jpim

  IMPLICIT NONE

  INTEGER(jpim),  INTENT(OUT)   :: err
  INTEGER(jpim),  INTENT(IN)    :: channels_oldsub(:)
  INTEGER(jpim),  INTENT(IN)    :: channels_ext(:)
  INTEGER(jpim),  INTENT(INOUT) :: nchannels_newsub
  INTEGER(jpim),  POINTER       :: channels_newsub(:)
  INTEGER(jpim),  POINTER       :: index_newsub_ext(:)
  INTEGER(jpim),  INTENT(INOUT) :: index_oldsub_newsub(:)
!INTF_END

#include "rttov_errorreport.interface"

  INTEGER(jpim) :: i, j, nchannels_oldsub, nchannels_ext
  INTEGER(jpim) :: channels_newsub_tmp(SIZE(channels_ext))
  INTEGER(jpim) :: index_oldsub_newsub_tmp(SIZE(channels_ext))

  TRY

  nchannels_oldsub = SIZE(channels_oldsub)
  nchannels_ext = SIZE(channels_ext)

  ! For each extracted channel in channels_ext, see if it is in the old sub-list (channels_oldsub):
  !   if so:
  !   - increment the count of channels to extract from the old sub-list
  !   - record the channel number (this is the index of the channel in the extracted channel list)
  !   - record the index into the old sub-list

  nchannels_newsub = 0
  DO i = 1, nchannels_ext
    DO j = 1, nchannels_oldsub
      IF (channels_ext(i) == channels_oldsub(j)) THEN
        nchannels_newsub = nchannels_newsub + 1         ! Count of channels in new sub-list
        channels_newsub_tmp(nchannels_newsub) = i       ! New channel sub-list: indices are into extracted channel list
        index_oldsub_newsub_tmp(nchannels_newsub) = j   ! Indices in old sub-list of new sub-list (i.e. indices to extract)
        EXIT
      ENDIF
    ENDDO
  ENDDO

  index_oldsub_newsub = 0
  IF (nchannels_newsub > 0) THEN
    ! New channel sub-list for extracted channels
    ALLOCATE(channels_newsub(nchannels_newsub), STAT=err)
    THROWM(err.NE.0,'Error allocating channels_newsub')
    channels_newsub(:) = channels_newsub_tmp(1:nchannels_newsub)

    ! Indices into new sub-list for each extracted channel
    ALLOCATE(index_newsub_ext(nchannels_ext), STAT=err)
    THROWM(err.NE.0,'Error allocating index_newsub_ext')
    index_newsub_ext = 0
    index_newsub_ext(channels_newsub(1:nchannels_newsub)) = &
          (/ (i, i = 1, nchannels_newsub) /)

    ! Indices into old sub-list for each channel in new sub-list
    index_oldsub_newsub(1:nchannels_newsub) = index_oldsub_newsub_tmp(1:nchannels_newsub)
  ENDIF

  CATCH
END SUBROUTINE
