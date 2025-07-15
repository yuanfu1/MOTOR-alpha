! Description:
!> @file
!!   Deallocate an MFASIS LUT structure.
!!   This should usually be called via rttov_dealloc_coefs.
!
!> @brief
!!   Deallocate an MFASIS LUT structure.
!!   This should usually be called via rttov_dealloc_coefs.
!!
!! @param[out]     err           status on exit
!! @param[in,out]  coef_mfasis   the MFASIS LUT structure to deallocate
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
SUBROUTINE rttov_dealloc_coef_mfasis(err, coef_mfasis)
!INTF_OFF
#include "throw.h"
!INTF_ON
  USE rttov_types, ONLY : rttov_coef_mfasis
  USE parkind1, ONLY : jpim
  IMPLICIT NONE
  INTEGER(KIND=jpim),      INTENT(OUT)   :: err
  TYPE(rttov_coef_mfasis), INTENT(INOUT) :: coef_mfasis
!INTF_END
#include "rttov_errorreport.interface"
#include "rttov_nullify_coef_mfasis.interface"
  INTEGER(KIND=jpim) :: i
!- End of header --------------------------------------------------------
  TRY

  IF (ASSOCIATED(coef_mfasis%aer_types)) &
      DEALLOCATE(coef_mfasis%aer_types, STAT = err)
  THROW(err.NE.0)

  IF (ASSOCIATED(coef_mfasis%channel_list)) &
      DEALLOCATE(coef_mfasis%channel_list, STAT = err)
  THROW(err.NE.0)

  IF (ASSOCIATED(coef_mfasis%channel_lut_index)) &
      DEALLOCATE(coef_mfasis%channel_lut_index, STAT = err)
  THROW(err.NE.0)

  IF (ASSOCIATED(coef_mfasis%lut_axes)) THEN
    DO i = 1, SIZE(coef_mfasis%lut_axes)
      IF (ASSOCIATED(coef_mfasis%lut_axes(i)%values)) &
          DEALLOCATE(coef_mfasis%lut_axes(i)%values, STAT = err)
      THROW(err.NE.0)
    ENDDO
    DEALLOCATE(coef_mfasis%lut_axes, STAT = err)
    THROW(err.NE.0)
  ENDIF

  IF (ASSOCIATED(coef_mfasis%lut)) THEN
    DO i = 1, SIZE(coef_mfasis%lut)
      IF (ASSOCIATED(coef_mfasis%lut(i)%qint)) &
          DEALLOCATE(coef_mfasis%lut(i)%qint, STAT = err)
      THROW(err.NE.0)
      IF (ASSOCIATED(coef_mfasis%lut(i)%data)) &
          DEALLOCATE(coef_mfasis%lut(i)%data, STAT = err)
      THROW(err.NE.0)
    ENDDO
    DEALLOCATE(coef_mfasis%lut, STAT = err)
    THROW(err.NE.0)
  ENDIF

  CALL rttov_nullify_coef_mfasis(coef_mfasis)

  CATCH
END SUBROUTINE rttov_dealloc_coef_mfasis
