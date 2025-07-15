! Description:
!> @file
!!   Allocate/deallocate internal mfasis_refl structure.
!
!> @brief
!!   Allocate/deallocate internal mfasis_refl structure.
!!
!! @details
!!   This structure holds quantities computed by direct and needed for TL/AD/K
!!
!!   Some arrays are only required by the direct model instance
!!   of this structure: these are allocated if the optional "direct"
!!   argument is TRUE. If omitted, it is assumed to be FALSE.
!!
!! @param[out]    err                    status on exit
!! @param[in,out] mfasis_refl            structure to (de)allocate
!! @param[in]     nchanprof              total number of channels being simulated (SIZE(chanprof))
!! @param[in]     asw                    1_jpim => allocate; 0_jpim => deallocate
!! @param[in]     ncolumns               number of cloud columns as calculated by rttov_cloud_overlap
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
!    Copyright 2018, EUMETSAT, All Rights Reserved.
!
SUBROUTINE rttov_alloc_mfasis_refl( &
              err,            &
              mfasis_refl,    &
              nchanprof,      &
              asw,            &
              ncolumns)
!INTF_OFF
#include "throw.h"
!INTF_ON
  USE rttov_types, ONLY : rttov_mfasis_refl
  USE parkind1, ONLY : jpim
!INTF_OFF
!INTF_ON
  IMPLICIT NONE

  INTEGER(KIND=jpim),                INTENT(OUT)          :: err
  TYPE(rttov_mfasis_refl),           POINTER              :: mfasis_refl(:,:)
  INTEGER(KIND=jpim),                INTENT(IN)           :: nchanprof
  INTEGER(KIND=jpim),                INTENT(IN)           :: asw
  INTEGER(KIND=jpim),                INTENT(IN)           :: ncolumns
!INTF_END
#include "rttov_errorreport.interface"

  INTEGER(KIND=jpim) :: i, col
!- End of header --------------------------------------------------------
  TRY

  IF (asw == 1_jpim) THEN
    ALLOCATE(mfasis_refl(0:ncolumns,nchanprof), STAT=err)
    THROWM(err.NE.0,"allocation of mfasis_refl")
    DO i = 1, nchanprof
      DO col = 0, ncolumns
        NULLIFY(mfasis_refl(col,i)%refl_lin_coef)
      ENDDO
    ENDDO
  ENDIF

  IF (asw == 0_jpim) THEN
    DO i = 1, nchanprof
      DO col = 0, ncolumns
        IF (ASSOCIATED(mfasis_refl(col,i)%refl_lin_coef)) THEN
          DEALLOCATE(mfasis_refl(col,i)%refl_lin_coef, STAT = err)
          THROWM(err.NE.0,"deallocation of mfasis_refl%refl_lin_coef")
        ENDIF
        NULLIFY(mfasis_refl(col,i)%refl_lin_coef)
      ENDDO
    ENDDO

    DEALLOCATE(mfasis_refl, STAT=err)
    THROWM(err.NE.0,"deallocation of mfasis_refl")
  ENDIF

  CATCH


END SUBROUTINE rttov_alloc_mfasis_refl
