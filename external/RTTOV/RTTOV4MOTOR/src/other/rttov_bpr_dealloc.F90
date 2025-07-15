! Description:
!> @file
!!   Deallocate look-up tables for bpr calculations.
!
!> @brief
!!   Deallocate look-up tables for bpr (back-scattering
!!   parameter) calculations.
!!
!! @details
!!   This should be called once at the end of processing
!!   when deallocating other memory.
!!
!! @param[out]  err       status on exit
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
SUBROUTINE RTTOV_BPR_DEALLOC( ERR )

!INTF_OFF
USE RTTOV_BPR_MOD
#include "throw.h"
!INTF_ON
USE PARKIND1,    ONLY : JPIM

IMPLICIT NONE
    INTEGER(KIND=JPIM), INTENT(OUT) :: ERR
!INTF_END

#include "rttov_errorreport.interface"

TRY

    IF (PHASE_INIT) THEN
      IF (ASSOCIATED(ARX)) DEALLOCATE(ARX, STAT = ERR)
      THROWM(err.NE.0,"Deallocation of ARX failed")
      IF (ASSOCIATED(XARR0)) DEALLOCATE(XARR0, STAT = ERR)
      THROWM(err.NE.0,"Deallocation of XARR0 failed")
      IF (ASSOCIATED(CXARR0)) DEALLOCATE(CXARR0, STAT = ERR)
      THROWM(err.NE.0,"Deallocation of CXARR0 failed")
      IF (ASSOCIATED(MUX)) DEALLOCATE(MUX, STAT = ERR)
      THROWM(err.NE.0,"Deallocation of MUX failed")

      PHASE_INIT = .FALSE.
    ENDIF

CATCH
END SUBROUTINE RTTOV_BPR_DEALLOC
