! Description:
!> @file
!!   Calculate a weighting function given an optical depth profile
!
!> @brief
!!   Calculate a weighting function given an optical depth profile
!!
!! @details
!!   Weighting function is calculated as
!!      d(transmittance)/d(-log(p))
!!
!!   The optical depths could be taken directly from
!!      traj\%opdp_path\%atm_level(:,ichan) for example.
!!
!!   Note that the RTTOV convention is to store optical depths as
!!   negative values and this subroutine follows this convention
!!   so the transmittance is calculated as EXP(opdep) (rather than
!!   EXP(-opdep)).
!!
!! @param[out]  err           status on exit
!! @param[in]   p             pressure levels
!! @param[in]   opdep         optical depths from each level to space
!! @param[out]  weighting_fn  calculated weighting function
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
SUBROUTINE rttov_calc_weighting_fn(err, p, opdep, weighting_fn)

  USE parkind1, ONLY : jprb, jpim

!INTF_OFF
#include "throw.h"
!INTF_ON

  IMPLICIT NONE

  INTEGER(KIND=jpim), INTENT(OUT) :: err
  REAL(KIND=jprb),    INTENT(IN)  :: p(:)
  REAL(KIND=jprb),    INTENT(IN)  :: opdep(SIZE(p))
  REAL(KIND=jprb),    INTENT(OUT) :: weighting_fn(SIZE(p)-1)
!INTF_END

#include "rttov_errorreport.interface"

  INTEGER(KIND=jpim) :: i
  INTEGER(KIND=jpim) :: nlevels

  TRY
  nlevels = SIZE(p)

  DO i = 2, nlevels

    IF (p(i-1) >= p(i)) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0, "Error in pressure levels")
    ELSE
      weighting_fn(i-1) = (EXP(opdep(i)) - EXP(opdep(i-1))) / (LOG(p(i-1)) - LOG(p(i)))
    ENDIF
  ENDDO

  CATCH
END SUBROUTINE rttov_calc_weighting_fn
