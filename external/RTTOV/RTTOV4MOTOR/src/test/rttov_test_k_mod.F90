! Description:
!> @file
!!   Helper subroutines for K_TL and K_BF calculation
!
!> @brief
!!   Helper subroutines for K_TL and K_BF calculation
!!
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
MODULE rttov_test_k_mod

#include "throw.h"

  USE rttov_chain, ONLY : chain, pchain, chain_rttov_profile, chain_rttov_profile_cloud, &
                          free_chain, advance_chain, get_pointer_chain, &
                          zero_chain
  USE parkind1, ONLY : jpim, jprb, jplm
  USE rttov_types, ONLY : rttov_profile, rttov_profile_cloud

  IMPLICIT NONE

  INTEGER(jpim), PARAMETER :: radar_k_lev = 40

  PRIVATE
  PUBLIC :: make_chain_profile, make_chain_cld_profile, &
            free_chain_profile, assign_chain_profile, radar_k_lev

#include "rttov_errorreport.interface"

CONTAINS

  SUBROUTINE make_chain_profile(err, chain_profiles, c, t, profiles, lzero)
!
    INTEGER(KIND=jpim),  INTENT(OUT) :: err
    TYPE(chain),         POINTER     :: chain_profiles(:)
    TYPE(pchain),        POINTER     :: c(:)
    CHARACTER(LEN=*),    INTENT(IN)  :: t
    TYPE(rttov_profile), INTENT(IN), TARGET :: profiles(:)
    LOGICAL(KIND=jplm),  INTENT(IN)  :: lzero
!
    INTEGER(KIND=jpim) :: iprof, nprof
!
    TRY
    nprof = SIZE(profiles)

    ALLOCATE (chain_profiles(nprof), c(nprof), STAT = err)
    THROWM(err.NE.0,"Allocation of chain_profiles failed")

    DO iprof = 1, nprof

      CALL chain_rttov_profile( &
              err,                    &
              chain_profiles(iprof),  &
              TRIM(t),                &
              a0 = profiles(iprof))
      THROW(err.NE.0)

      c(iprof)%p => chain_profiles(iprof)

      IF (lzero) CALL zero_chain(c(iprof)%p)

    ENDDO
    CATCH
  END SUBROUTINE

  SUBROUTINE make_chain_cld_profile(err, chain_profiles, c, t, profiles, lzero)
!
    INTEGER(KIND=jpim),        INTENT(OUT) :: err
    TYPE(chain),               POINTER     :: chain_profiles(:)
    TYPE(pchain),              POINTER     :: c(:)
    CHARACTER(LEN=*),          INTENT(IN)  :: t
    TYPE(rttov_profile_cloud), INTENT(IN), TARGET :: profiles(:)
    LOGICAL(KIND=jplm),        INTENT(IN)  :: lzero
!
    INTEGER(KIND=jpim) :: iprof, nprof
!
    TRY
    nprof = SIZE(profiles)

    ALLOCATE (chain_profiles(nprof), c(nprof), STAT = err)
    THROWM(err.NE.0,"Allocation of chain_profiles failed")

    DO iprof = 1, nprof

      CALL chain_rttov_profile_cloud( &
              err,                    &
              chain_profiles(iprof),  &
              TRIM(t),                &
              a0 = profiles(iprof))
      THROW(err.NE.0)

      c(iprof)%p => chain_profiles(iprof)

      IF (lzero) CALL zero_chain(c(iprof)%p)

    ENDDO
    CATCH
  END SUBROUTINE

  SUBROUTINE free_chain_profile(err, chain_profiles, c)
!
    INTEGER(KIND=jpim), INTENT(OUT) :: err
    TYPE(chain),        POINTER     :: chain_profiles(:)
    TYPE(pchain),       POINTER     :: c(:)
!
    INTEGER(KIND=jpim) :: iprof, nprof
!
    TRY
    IF (ASSOCIATED(chain_profiles)) THEN
      nprof = SIZE(chain_profiles)
      DO iprof = 1, nprof
        CALL free_chain(chain_profiles(iprof))
      ENDDO
      DEALLOCATE (chain_profiles, c, STAT = err)
      THROWM(err.NE.0,"Deallocation of chain_profiles failed")
    ENDIF
    CATCH
  END SUBROUTINE


  SUBROUTINE assign_chain_profile(c, w)
!
    TYPE(pchain), INTENT(INOUT) :: c(:)
    REAL(KIND=jprb), INTENT(IN) :: w(:)
!
    INTEGER(KIND=jpim) :: nprof, iprof
    REAL(KIND=jprb), POINTER :: x
!
    nprof = SIZE(w)
    DO iprof = 1, nprof
      CALL get_pointer_chain(c(iprof)%p, x)
      x = w(iprof)
      CALL advance_chain(c(iprof)%p)
    ENDDO

  END SUBROUTINE

END MODULE
