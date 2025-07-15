! Description:
!> @file
!!   Initialise auxiliary phase function arrays.
!
!> @brief
!!   Initialise auxiliary phase function arrays.
!!
!! @details
!!   This subroutine is called internally by RTTOV to pre-calculate
!!   some data to allow fast interpolation of phase functions.
!!
!! @param[out]    err            status on exit
!! @param[in]     phangle        phase function angle grid array
!! @param[in,out] phfn_int       phase function interpolation data structure
!! @param[in]     asw            1_jpim => allocate; 0_jpim => deallocate
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
SUBROUTINE rttov_alloc_phfn_int(err, phangle, phfn_int, asw)

!INTF_OFF
#include "throw.h"
!INTF_ON
  USE rttov_types, ONLY : rttov_phasefn_int
  USE parkind1, ONLY : jpim, jprb
!INTF_OFF
  USE rttov_const, ONLY : deg2rad
!INTF_ON
  IMPLICIT NONE
  INTEGER(jpim),           INTENT(OUT)   :: err
  REAL(jprb),              INTENT(IN)    :: phangle(:)
  TYPE(rttov_phasefn_int), INTENT(INOUT) :: phfn_int
  INTEGER(jpim),           INTENT(IN)    :: asw
!INTF_END
  INTEGER(jpim) :: nphangle, i, k, icount
  REAL(jprb)    :: minphadiff

#include "rttov_errorreport.interface"

  TRY

  IF (asw == 1_jpim) THEN
    NULLIFY(phfn_int%cosphangle, phfn_int%iphangle)
    nphangle = SIZE(phangle)
    IF (nphangle > 0) THEN
      ALLOCATE(phfn_int%cosphangle(nphangle), STAT=err)
      THROWM(err.NE.0, "allocation of cosphangle")

      phfn_int%cosphangle(:) = COS(phangle(:) * deg2rad)
      minphadiff = 0.5_jprb * MINVAL(phangle(2:nphangle) - phangle(1:nphangle-1))
      phfn_int%zminphadiff = 1._jprb / minphadiff

      icount = phangle(nphangle) * phfn_int%zminphadiff
      ALLOCATE (phfn_int%iphangle(icount), STAT=err)
      THROWM(err.NE.0, "allocation of iphangle")

      k = 1
      DO i = 1, icount
        IF (phangle(k) >= (i + 1) * minphadiff) THEN
          phfn_int%iphangle(i) = k
        ELSE
          k = k + 1
          phfn_int%iphangle(i) = k
        ENDIF
      ENDDO
    ENDIF
  ENDIF

  IF (asw == 0_jpim) THEN
    IF (ASSOCIATED(phfn_int%cosphangle)) DEALLOCATE(phfn_int%cosphangle, STAT=err)
    THROWM(err.NE.0, "deallocation of cosphangle")
    IF (ASSOCIATED(phfn_int%iphangle))   DEALLOCATE(phfn_int%iphangle, STAT=err)
    THROWM(err.NE.0, "deallocation of iphangle")
    NULLIFY(phfn_int%cosphangle, phfn_int%iphangle)
  ENDIF

  CATCH

END SUBROUTINE rttov_alloc_phfn_int
