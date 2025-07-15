! Description:
!> @file
!!   Explicit interfaces for all LAPACK subroutines called by RTTOV.
!
!> @brief
!!   Explicit interfaces for all LAPACK subroutines called by RTTOV.
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
!    Copyright 2016, EUMETSAT, All Rights Reserved.
!
MODULE rttov_lapack_mod

  IMPLICIT NONE

  INTERFACE
    SUBROUTINE DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )

      IMPLICIT NONE
      INTEGER          :: INFO, LDA, LDB, N, NRHS
      INTEGER          :: IPIV(*)
      DOUBLE PRECISION :: A(LDA,*), B(LDB,*)

    END SUBROUTINE DGESV
  END INTERFACE

  INTERFACE
    SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )

      IMPLICIT NONE
      INTEGER          :: INFO, LDA, M, N
      INTEGER          :: IPIV(*)
      DOUBLE PRECISION :: A(LDA,*)

    END SUBROUTINE DGETRF
  END INTERFACE

  INTERFACE
    SUBROUTINE DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )

      IMPLICIT NONE
      CHARACTER        :: TRANS
      INTEGER          :: INFO, LDA, LDB, N, NRHS
      INTEGER          :: IPIV(*)
      DOUBLE PRECISION :: A(LDA,*), B(LDB,*)

    END SUBROUTINE DGETRS
  END INTERFACE

  INTERFACE
    SUBROUTINE DGBTRF( M, N, KL, KU, AB, LDAB, IPIV, INFO )

      IMPLICIT NONE
      INTEGER          :: INFO, KL, KU, LDAB, M, N
      INTEGER          :: IPIV(*)
      DOUBLE PRECISION :: AB(LDAB,*)

    END SUBROUTINE DGBTRF
  END INTERFACE

  INTERFACE
    SUBROUTINE DGBTRS( TRANS, N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )

      IMPLICIT NONE
      CHARACTER        :: TRANS
      INTEGER          :: INFO, KL, KU, LDAB, LDB, N, NRHS
      INTEGER          :: IPIV(*)
      DOUBLE PRECISION :: AB(LDAB,*), B(LDB,*)

    END SUBROUTINE DGBTRS
  END INTERFACE

  INTERFACE
    SUBROUTINE DGELS(TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO)

      IMPLICIT NONE
      CHARACTER        :: TRANS
      INTEGER          :: INFO, LDA, LDB, LWORK, M, N, NRHS
      DOUBLE PRECISION :: A(LDA,*), B(LDB,*), WORK(*)
    END SUBROUTINE
  END INTERFACE

END MODULE rttov_lapack_mod
