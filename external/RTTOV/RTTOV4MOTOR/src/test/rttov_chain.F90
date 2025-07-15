! Description:
!> @file
!!   Subroutines for creating and managing linked lists of data in a profile
!!   structure.
!
!> @brief
!!   Subroutines for creating and managing linked lists of data in a profile
!!   structure.
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
MODULE RTTOV_CHAIN

USE PARKIND1, ONLY: JPIM, JPRB, JPLM
USE RTTOV_TYPES, ONLY: RTTOV_PROFILE, RTTOV_S2M, RTTOV_SKIN, RTTOV_PROFILE_CLOUD
USE RTTOV_CONST, ONLY : ERRORSTATUS_FATAL

IMPLICIT NONE

PRIVATE
PUBLIC :: CHAIN, PCHAIN, SIZE_CHAIN, CHAIN_RTTOV_PROFILE, CHAIN_RTTOV_PROFILE_CLOUD, &
          FREE_CHAIN, ADVANCE_CHAIN, GET_POINTER_CHAIN, &
          ZERO_CHAIN, PRINT_CHAIN, PRINT_ARRAY

TYPE CHAIN
  TYPE(CHAIN),  POINTER :: NEXT => NULL()
  REAL(KIND=JPRB), POINTER :: R_8_0 => NULL()
  REAL(KIND=JPRB), POINTER :: R_8_1(:) => NULL()
  REAL(KIND=JPRB), POINTER :: R_8_2(:,:) => NULL()
  REAL(KIND=JPRB), POINTER :: R_8_3(:,:,:) => NULL()
  REAL(KIND=JPRB), POINTER :: R_8_4(:,:,:,:) => NULL()
  REAL(KIND=JPRB), POINTER :: R_8_5(:,:,:,:,:) => NULL()
  REAL(KIND=JPRB), POINTER :: R_8_6(:,:,:,:,:,:) => NULL()
  REAL(KIND=JPRB), POINTER :: R_8_7(:,:,:,:,:,:,:) => NULL()
  INTEGER(KIND=JPIM) :: II(7) = (/ 1, 1, 1, 1, 1, 1, 1 /)
  CHARACTER(LEN=128) :: NAME  = ""
END TYPE

TYPE PCHAIN
  TYPE(CHAIN), POINTER :: P
END TYPE

CONTAINS

SUBROUTINE CLEAR_CHAIN( C )
TYPE(CHAIN), INTENT(INOUT) :: C
TYPE(CHAIN) :: C1

 C = C1

END SUBROUTINE

SUBROUTINE ADVANCE_II( ADVANCED, II, UBOUNDA )
LOGICAL(KIND=JPLM) :: ADVANCED
INTEGER(KIND=JPIM), DIMENSION(7) :: II, UBOUNDA
INTEGER(KIND=JPIM) :: I, J

ADVANCED = .TRUE.

DO I = 7, 1, -1
  IF( II(I) .LT. UBOUNDA(I) ) THEN
    II(I) = II(I) + 1
    DO J = I+1, 7
      II(J) = 1
    ENDDO
    RETURN
  ENDIF
ENDDO

ADVANCED = .FALSE.

END SUBROUTINE

SUBROUTINE SET_VALUE_CHAIN( C, X )
TYPE(CHAIN),  INTENT(INOUT) :: C
REAL(KIND=JPRB), INTENT(IN) :: X

IF( ASSOCIATED( C%R_8_0 ) ) C%R_8_0                                                          = X
IF( ASSOCIATED( C%R_8_1 ) ) C%R_8_1(C%II(1))                                                 = X
IF( ASSOCIATED( C%R_8_2 ) ) C%R_8_2(C%II(1),C%II(2))                                         = X
IF( ASSOCIATED( C%R_8_3 ) ) C%R_8_3(C%II(1),C%II(2),C%II(3))                                 = X
IF( ASSOCIATED( C%R_8_4 ) ) C%R_8_4(C%II(1),C%II(2),C%II(3),C%II(4))                         = X
IF( ASSOCIATED( C%R_8_5 ) ) C%R_8_5(C%II(1),C%II(2),C%II(3),C%II(4),C%II(5))                 = X
IF( ASSOCIATED( C%R_8_6 ) ) C%R_8_6(C%II(1),C%II(2),C%II(3),C%II(4),C%II(5),C%II(6))         = X
IF( ASSOCIATED( C%R_8_7 ) ) C%R_8_7(C%II(1),C%II(2),C%II(3),C%II(4),C%II(5),C%II(6),C%II(7)) = X

END SUBROUTINE

SUBROUTINE GET_VALUE_CHAIN( C, X )
TYPE(CHAIN),  INTENT(IN)  :: C
REAL(KIND=JPRB), INTENT(OUT) :: X

IF( ASSOCIATED( C%R_8_0 ) ) X = C%R_8_0                                                     
IF( ASSOCIATED( C%R_8_1 ) ) X = C%R_8_1(C%II(1))                                                
IF( ASSOCIATED( C%R_8_2 ) ) X = C%R_8_2(C%II(1),C%II(2))                                        
IF( ASSOCIATED( C%R_8_3 ) ) X = C%R_8_3(C%II(1),C%II(2),C%II(3))                                
IF( ASSOCIATED( C%R_8_4 ) ) X = C%R_8_4(C%II(1),C%II(2),C%II(3),C%II(4))                        
IF( ASSOCIATED( C%R_8_5 ) ) X = C%R_8_5(C%II(1),C%II(2),C%II(3),C%II(4),C%II(5))                
IF( ASSOCIATED( C%R_8_6 ) ) X = C%R_8_6(C%II(1),C%II(2),C%II(3),C%II(4),C%II(5),C%II(6))        
IF( ASSOCIATED( C%R_8_7 ) ) X = C%R_8_7(C%II(1),C%II(2),C%II(3),C%II(4),C%II(5),C%II(6),C%II(7))

END SUBROUTINE

SUBROUTINE GET_POINTER_CHAIN( C, X )
TYPE(CHAIN), POINTER :: C
REAL(KIND=JPRB), POINTER :: X


X => NULL()

IF( .NOT. ASSOCIATED( C ) ) RETURN

IF( ASSOCIATED( C%R_8_0 ) ) X => C%R_8_0                                                     
IF( ASSOCIATED( C%R_8_1 ) ) X => C%R_8_1(C%II(1))                                                 
IF( ASSOCIATED( C%R_8_2 ) ) X => C%R_8_2(C%II(1),C%II(2))                                         
IF( ASSOCIATED( C%R_8_3 ) ) X => C%R_8_3(C%II(1),C%II(2),C%II(3))                                 
IF( ASSOCIATED( C%R_8_4 ) ) X => C%R_8_4(C%II(1),C%II(2),C%II(3),C%II(4))                         
IF( ASSOCIATED( C%R_8_5 ) ) X => C%R_8_5(C%II(1),C%II(2),C%II(3),C%II(4),C%II(5))                 
IF( ASSOCIATED( C%R_8_6 ) ) X => C%R_8_6(C%II(1),C%II(2),C%II(3),C%II(4),C%II(5),C%II(6))         
IF( ASSOCIATED( C%R_8_7 ) ) X => C%R_8_7(C%II(1),C%II(2),C%II(3),C%II(4),C%II(5),C%II(6),C%II(7)) 

END SUBROUTINE

SUBROUTINE ADVANCE_CHAIN( C )
TYPE(CHAIN), POINTER :: C
INTEGER(KIND=JPIM) :: UBOUNDA(7)
LOGICAL(KIND=JPLM) :: ADV

UBOUNDA = 1

IF( ASSOCIATED( C%R_8_1 ) ) THEN
  UBOUNDA(1:1) = UBOUND( C%R_8_1 )
END IF
IF( ASSOCIATED( C%R_8_2 ) ) THEN
  UBOUNDA(1:2) = UBOUND( C%R_8_2 )
END IF
IF( ASSOCIATED( C%R_8_3 ) ) THEN
  UBOUNDA(1:3) = UBOUND( C%R_8_3 )
END IF
IF( ASSOCIATED( C%R_8_4 ) ) THEN
  UBOUNDA(1:4) = UBOUND( C%R_8_4 )
END IF
IF( ASSOCIATED( C%R_8_5 ) ) THEN
  UBOUNDA(1:5) = UBOUND( C%R_8_5 )
END IF
IF( ASSOCIATED( C%R_8_6 ) ) THEN
  UBOUNDA(1:6) = UBOUND( C%R_8_6 )
END IF
IF( ASSOCIATED( C%R_8_7 ) ) THEN
  UBOUNDA(1:7) = UBOUND( C%R_8_7 )
END IF

CALL ADVANCE_II( ADV, C%II, UBOUNDA )

IF( .NOT. ADV ) THEN
  DO WHILE( .TRUE. )
    C => C%NEXT
    IF( .NOT. ASSOCIATED( C ) ) EXIT
    IF( ASSOCIATED( C%R_8_0 ) ) EXIT
    IF( ASSOCIATED( C%R_8_1 ) ) EXIT
    IF( ASSOCIATED( C%R_8_2 ) ) EXIT
    IF( ASSOCIATED( C%R_8_3 ) ) EXIT
    IF( ASSOCIATED( C%R_8_4 ) ) EXIT
    IF( ASSOCIATED( C%R_8_5 ) ) EXIT
    IF( ASSOCIATED( C%R_8_6 ) ) EXIT
    IF( ASSOCIATED( C%R_8_7 ) ) EXIT
  ENDDO
ENDIF

END SUBROUTINE

SUBROUTINE FREE_CHAIN( C )
TYPE(CHAIN), INTENT(INOUT) :: C
TYPE(CHAIN), POINTER :: C1, C2

C2 => C%NEXT

DO WHILE( ASSOCIATED( C2 ) )

  C1 => C2%NEXT
  DEALLOCATE( C2 )
  C2 => C1

END DO

CALL CLEAR_CHAIN( C )

END SUBROUTINE

SUBROUTINE ZERO_CHAIN( C ) 
TYPE(CHAIN), TARGET, INTENT(IN) :: C
TYPE(CHAIN), POINTER :: C1

C1 => C

DO WHILE( ASSOCIATED( C1 ) )

  IF( ASSOCIATED( C1%R_8_0 ) ) C1%R_8_0 = 0.
  IF( ASSOCIATED( C1%R_8_1 ) ) C1%R_8_1 = 0.
  IF( ASSOCIATED( C1%R_8_2 ) ) C1%R_8_2 = 0.
  IF( ASSOCIATED( C1%R_8_3 ) ) C1%R_8_3 = 0.
  IF( ASSOCIATED( C1%R_8_4 ) ) C1%R_8_4 = 0.
  IF( ASSOCIATED( C1%R_8_5 ) ) C1%R_8_5 = 0.
  IF( ASSOCIATED( C1%R_8_6 ) ) C1%R_8_6 = 0.
  IF( ASSOCIATED( C1%R_8_7 ) ) C1%R_8_7 = 0.
  
  C1 => C1%NEXT

ENDDO

END SUBROUTINE

SUBROUTINE SIZE_CHAIN( N, C ) 
INTEGER(KIND=JPIM), INTENT(OUT) :: N
TYPE(CHAIN), TARGET, INTENT(IN) :: C
TYPE(CHAIN), POINTER :: C1

C1 => C
N = 0

DO WHILE( ASSOCIATED( C1 ) )

  IF( ASSOCIATED( C1%R_8_0 ) ) N = N + 1
  IF( ASSOCIATED( C1%R_8_1 ) ) N = N + SIZE( C1%R_8_1 )
  IF( ASSOCIATED( C1%R_8_2 ) ) N = N + SIZE( C1%R_8_2 )
  IF( ASSOCIATED( C1%R_8_3 ) ) N = N + SIZE( C1%R_8_3 )
  IF( ASSOCIATED( C1%R_8_4 ) ) N = N + SIZE( C1%R_8_4 )
  IF( ASSOCIATED( C1%R_8_5 ) ) N = N + SIZE( C1%R_8_5 )
  IF( ASSOCIATED( C1%R_8_6 ) ) N = N + SIZE( C1%R_8_6 )
  IF( ASSOCIATED( C1%R_8_7 ) ) N = N + SIZE( C1%R_8_7 )
  
  C1 => C1%NEXT

ENDDO

END SUBROUTINE

RECURSIVE SUBROUTINE PRINT_ARRAY( LUN, FORMAT, A0, A1, A2, A3, A4, A5, A6, A7, AI1 )
INTEGER(KIND=JPIM), INTENT(IN)           :: LUN
CHARACTER(LEN=*),   INTENT(IN)           :: FORMAT
REAL(KIND=JPRB),    INTENT(IN), OPTIONAL :: A0 
REAL(KIND=JPRB),    INTENT(IN), OPTIONAL :: A1(:) 
REAL(KIND=JPRB),    INTENT(IN), OPTIONAL :: A2(:,:)
REAL(KIND=JPRB),    INTENT(IN), OPTIONAL :: A3(:,:,:)
REAL(KIND=JPRB),    INTENT(IN), OPTIONAL :: A4(:,:,:,:)
REAL(KIND=JPRB),    INTENT(IN), OPTIONAL :: A5(:,:,:,:,:)
REAL(KIND=JPRB),    INTENT(IN), OPTIONAL :: A6(:,:,:,:,:,:) 
REAL(KIND=JPRB),    INTENT(IN), OPTIONAL :: A7(:,:,:,:,:,:,:)
INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: AI1(:)
INTEGER(KIND=JPIM) :: I, J, L1, U1
INTEGER(KIND=JPIM), PARAMETER :: N = 10

IF( PRESENT( A0 ) ) THEN
  WRITE( LUN, FMT = FORMAT ) A0
ELSE IF( PRESENT( A1 ) ) THEN
  L1 = LBOUND( A1, 1 )
  U1 = UBOUND( A1, 1 )
  DO I = L1, U1, N
    J = I + N - 1
    IF( J .GT. U1 ) J = U1
    WRITE( LUN, FMT = FORMAT ) A1(I:J)
  ENDDO
ELSE IF( PRESENT( A2 ) ) THEN
  DO I = 1, SIZE( A2, 1 )
    CALL PRINT_ARRAY( LUN, FORMAT, A1 = A2(I,:) )
    WRITE( LUN, * )
  ENDDO
ELSE IF( PRESENT( A3 ) ) THEN
  DO I = 1, SIZE( A3, 1 )
    CALL PRINT_ARRAY( LUN, FORMAT, A2 = A3(I,:,:) )
    WRITE( LUN, * )
  ENDDO
ELSE IF( PRESENT( A4 ) ) THEN
  DO I = 1, SIZE( A4, 1 )
    CALL PRINT_ARRAY( LUN, FORMAT, A3 = A4(I,:,:,:) )
    WRITE( LUN, * )
  ENDDO
ELSE IF( PRESENT( A5 ) ) THEN
  DO I = 1, SIZE( A5, 1 )
    CALL PRINT_ARRAY( LUN, FORMAT, A4 = A5(I,:,:,:,:) )
    WRITE( LUN, * )
  ENDDO
ELSE IF( PRESENT( A6 ) ) THEN
  DO I = 1, SIZE( A6, 1 )
    CALL PRINT_ARRAY( LUN, FORMAT, A5 = A6(I,:,:,:,:,:) )
    WRITE( LUN, * )
  ENDDO
ELSE IF( PRESENT( A7 ) ) THEN
  DO I = 1, SIZE( A7, 1 )
    CALL PRINT_ARRAY( LUN, FORMAT, A6 = A7(I,:,:,:,:,:,:) )
    WRITE( LUN, * )
  ENDDO
ELSE IF( PRESENT( AI1 ) ) THEN
  L1 = LBOUND( AI1, 1 )
  U1 = UBOUND( AI1, 1 )
  DO I = L1, U1, N
    J = I + N - 1
    IF( J .GT. U1 ) J = U1
    WRITE( LUN, FMT = FORMAT ) AI1(I:J)
  ENDDO
ENDIF

END SUBROUTINE


SUBROUTINE PRINT_CHAIN( LUN, C, FORMAT )
INTEGER(KIND=JPIM), INTENT(IN) :: LUN
TYPE(CHAIN),   TARGET,   INTENT(IN) :: C
TYPE(CHAIN),   POINTER              :: C1
CHARACTER(*),  OPTIONAL, INTENT(IN) :: FORMAT
INTEGER(KIND=JPIM) :: LEN_MAX 
CHARACTER(LEN=32) :: FORMAT_N
CHARACTER(LEN=32) :: FORMAT1

IF( PRESENT( FORMAT ) ) THEN
  FORMAT1 = '(10000(' // TRIM( FORMAT ) // '))'
ELSE
  FORMAT1 = '(10000(G16.7E3))'
ENDIF

LEN_MAX = 0
C1 => C
DO WHILE( ASSOCIATED( C1 ) )
  IF( LEN_MAX .LT. LEN( TRIM( C1%NAME ) ) ) THEN
    LEN_MAX = LEN( TRIM( C1%NAME ) )
  ENDIF
  C1 => C1%NEXT
ENDDO


WRITE( FORMAT_N, '("(A",I5.5,","" = ("")")' ) LEN_MAX

C1 => C
DO WHILE( ASSOCIATED( C1 ) )
  WRITE( LUN, FORMAT_N ) C1%NAME
  IF( ASSOCIATED( C1%R_8_0 ) ) CALL PRINT_ARRAY( LUN, FORMAT1, A0 = C1%R_8_0 )
  IF( ASSOCIATED( C1%R_8_1 ) ) CALL PRINT_ARRAY( LUN, FORMAT1, A1 = C1%R_8_1 )
  IF( ASSOCIATED( C1%R_8_2 ) ) CALL PRINT_ARRAY( LUN, FORMAT1, A2 = C1%R_8_2 )
  IF( ASSOCIATED( C1%R_8_3 ) ) CALL PRINT_ARRAY( LUN, FORMAT1, A3 = C1%R_8_3 )
  IF( ASSOCIATED( C1%R_8_4 ) ) CALL PRINT_ARRAY( LUN, FORMAT1, A4 = C1%R_8_4 )
  IF( ASSOCIATED( C1%R_8_5 ) ) CALL PRINT_ARRAY( LUN, FORMAT1, A5 = C1%R_8_5 )
  IF( ASSOCIATED( C1%R_8_6 ) ) CALL PRINT_ARRAY( LUN, FORMAT1, A6 = C1%R_8_6 )
  IF( ASSOCIATED( C1%R_8_7 ) ) CALL PRINT_ARRAY( LUN, FORMAT1, A7 = C1%R_8_7 )
  WRITE( LUN, * ) ')'
  C1 => C1%NEXT
ENDDO


END SUBROUTINE


SUBROUTINE chain_rttov_profile( err, chn, name, a0, a1, a2, a3, a4, a5, a6, a7, p0, p1, p2, p3, p4, p5, p6, p7 )
  INTEGER(KIND=JPIM),                    INTENT(OUT)   :: err
  TYPE(chain),                           INTENT(INOUT) :: chn
  CHARACTER(LEN=*),                      INTENT(IN)    :: name
  TYPE(rttov_profile), OPTIONAL, TARGET, INTENT(IN)    :: a0
  TYPE(rttov_profile), OPTIONAL, TARGET, INTENT(IN)    :: a1(:)
  TYPE(rttov_profile), OPTIONAL, TARGET, INTENT(IN)    :: a2(:,:)
  TYPE(rttov_profile), OPTIONAL, TARGET, INTENT(IN)    :: a3(:,:,:)
  TYPE(rttov_profile), OPTIONAL, TARGET, INTENT(IN)    :: a4(:,:,:,:)
  TYPE(rttov_profile), OPTIONAL, TARGET, INTENT(IN)    :: a5(:,:,:,:,:)
  TYPE(rttov_profile), OPTIONAL, TARGET, INTENT(IN)    :: a6(:,:,:,:,:,:)
  TYPE(rttov_profile), OPTIONAL, TARGET, INTENT(IN)    :: a7(:,:,:,:,:,:,:)
  TYPE(rttov_profile), OPTIONAL, POINTER               :: p0
  TYPE(rttov_profile), OPTIONAL, POINTER               :: p1(:)
  TYPE(rttov_profile), OPTIONAL, POINTER               :: p2(:,:)
  TYPE(rttov_profile), OPTIONAL, POINTER               :: p3(:,:,:)
  TYPE(rttov_profile), OPTIONAL, POINTER               :: p4(:,:,:,:)
  TYPE(rttov_profile), OPTIONAL, POINTER               :: p5(:,:,:,:,:)
  TYPE(rttov_profile), OPTIONAL, POINTER               :: p6(:,:,:,:,:,:)
  TYPE(rttov_profile), OPTIONAL, POINTER               :: p7(:,:,:,:,:,:,:)
  INTEGER(KIND=JPIM)                                   :: nargs
  INTEGER(KIND=JPIM)                                   :: lbounda(7)
  INTEGER(KIND=JPIM)                                   :: ubounda(7)
  INTEGER(KIND=JPIM)                                   :: i1
  INTEGER(KIND=JPIM)                                   :: i2
  INTEGER(KIND=JPIM)                                   :: i3
  INTEGER(KIND=JPIM)                                   :: i4
  INTEGER(KIND=JPIM)                                   :: i5
  INTEGER(KIND=JPIM)                                   :: i6
  INTEGER(KIND=JPIM)                                   :: i7
  TYPE(rttov_profile),           POINTER               :: p
  CHARACTER(LEN=128)                                   :: str
  err   = 0
  nargs = 0
  IF( present( a0 ) ) nargs = nargs + 1
  IF( present( a1 ) ) nargs = nargs + 1
  IF( present( a2 ) ) nargs = nargs + 1
  IF( present( a3 ) ) nargs = nargs + 1
  IF( present( a4 ) ) nargs = nargs + 1
  IF( present( a5 ) ) nargs = nargs + 1
  IF( present( a6 ) ) nargs = nargs + 1
  IF( present( a7 ) ) nargs = nargs + 1
  IF( present( p0 ) ) nargs = nargs + 1
  IF( present( p1 ) ) nargs = nargs + 1
  IF( present( p2 ) ) nargs = nargs + 1
  IF( present( p3 ) ) nargs = nargs + 1
  IF( present( p4 ) ) nargs = nargs + 1
  IF( present( p5 ) ) nargs = nargs + 1
  IF( present( p6 ) ) nargs = nargs + 1
  IF( present( p7 ) ) nargs = nargs + 1
  IF( nargs .NE. 1 ) GOTO 999
  lbounda(:) = 0
  ubounda(:) = 0
  IF( present( a1 ) ) THEN
    lbounda(1:1) = lbound( a1 )
    ubounda(1:1) = ubound( a1 )
  ENDIF
  IF( present( a2 ) ) THEN
    lbounda(1:2) = lbound( a2 )
    ubounda(1:2) = ubound( a2 )
  ENDIF
  IF( present( a3 ) ) THEN
    lbounda(1:3) = lbound( a3 )
    ubounda(1:3) = ubound( a3 )
  ENDIF
  IF( present( a4 ) ) THEN
    lbounda(1:4) = lbound( a4 )
    ubounda(1:4) = ubound( a4 )
  ENDIF
  IF( present( a5 ) ) THEN
    lbounda(1:5) = lbound( a5 )
    ubounda(1:5) = ubound( a5 )
  ENDIF
  IF( present( a6 ) ) THEN
    lbounda(1:6) = lbound( a6 )
    ubounda(1:6) = ubound( a6 )
  ENDIF
  IF( present( a7 ) ) THEN
    lbounda(1:7) = lbound( a7 )
    ubounda(1:7) = ubound( a7 )
  ENDIF
  IF( present( p1 ) ) THEN
    IF( associated( p1 ) ) THEN
      lbounda(1:1) = lbound( p1 )
      ubounda(1:1) = ubound( p1 )
    ENDIF
  ENDIF
  IF( present( p2 ) ) THEN
    IF( associated( p2 ) ) THEN
      lbounda(1:2) = lbound( p2 )
      ubounda(1:2) = ubound( p2 )
    ENDIF
  ENDIF
  IF( present( p3 ) ) THEN
    IF( associated( p3 ) ) THEN
      lbounda(1:3) = lbound( p3 )
      ubounda(1:3) = ubound( p3 )
    ENDIF
  ENDIF
  IF( present( p4 ) ) THEN
    IF( associated( p4 ) ) THEN
      lbounda(1:4) = lbound( p4 )
      ubounda(1:4) = ubound( p4 )
    ENDIF
  ENDIF
  IF( present( p5 ) ) THEN
    IF( associated( p5 ) ) THEN
      lbounda(1:5) = lbound( p5 )
      ubounda(1:5) = ubound( p5 )
    ENDIF
  ENDIF
  IF( present( p6 ) ) THEN
    IF( associated( p6 ) ) THEN
      lbounda(1:6) = lbound( p6 )
      ubounda(1:6) = ubound( p6 )
    ENDIF
  ENDIF
  IF( present( p7 ) ) THEN
    IF( associated( p7 ) ) THEN
      lbounda(1:7) = lbound( p7 )
      ubounda(1:7) = ubound( p7 )
    ENDIF
  ENDIF
  IF( present( a0 ) ) p => a0
  IF( present( p0 ) ) p => p0
  IF( present( a0 ) .OR. present( p0 ) ) THEN
    WRITE( UNIT = str, FMT = "(A)" ) trim( name )
  ENDIF
  DO i1 = ubounda(1), lbounda(1), -1
    IF( present( a1 ) ) p => a1(i1)
    IF( present( p1 ) ) p => p1(i1)
    IF( present( a1 ) .OR. present( p1 ) ) THEN
      WRITE( UNIT = str, FMT = "(A,'(',I4,')')" ) trim( name ), i1
    ENDIF
    DO i2 = ubounda(2), lbounda(2), -1
      IF( present( a2 ) ) p => a2(i1,i2)
      IF( present( p2 ) ) p => p2(i1,i2)
      IF( present( a2 ) .OR. present( p2 ) ) THEN
        WRITE( UNIT = str, FMT = "(A,'(',I4,',',I4,')')" ) trim( name ), i1, i2
      ENDIF
      DO i3 = ubounda(3), lbounda(3), -1
        IF( present( a3 ) ) p => a3(i1,i2,i3)
        IF( present( p3 ) ) p => p3(i1,i2,i3)
        IF( present( a3 ) .OR. present( p3 ) ) THEN
          WRITE( UNIT = str, FMT = "(A,'(',I4,',',I4,',',I4,')')" ) trim( name ), i1, i2, i3
        ENDIF
        DO i4 = ubounda(4), lbounda(4), -1
          IF( present( a4 ) ) p => a4(i1,i2,i3,i4)
          IF( present( p4 ) ) p => p4(i1,i2,i3,i4)
          IF( present( a4 ) .OR. present( p4 ) ) THEN
            WRITE( UNIT = str, FMT = "(A,'(',I4,',',I4,',',I4,',',I4,')')" ) trim( name ), i1, i2, i3, i4
          ENDIF
          DO i5 = ubounda(5), lbounda(5), -1
            IF( present( a5 ) ) p => a5(i1,i2,i3,i4,i5)
            IF( present( p5 ) ) p => p5(i1,i2,i3,i4,i5)
            IF( present( a5 ) .OR. present( p5 ) ) THEN
              WRITE( UNIT = str, FMT = "(A,'(',I4,',',I4,',',I4,',',I4,',',I4,')')" ) trim( name ), i1, i2, i3, i4, i5
            ENDIF
            DO i6 = ubounda(6), lbounda(6), -1
              IF( present( a6 ) ) p => a6(i1,i2,i3,i4,i5,i6)
              IF( present( p6 ) ) p => p6(i1,i2,i3,i4,i5,i6)
              IF( present( a6 ) .OR. present( p6 ) ) THEN
                WRITE( UNIT = str, FMT = "(A,'(',I4,',',I4,',',I4,',',I4,',',I4,',',I4,')')" ) trim( name ), i1, i2, i3, i4, i5, i6
              ENDIF
              DO i7 = ubounda(7), lbounda(7), -1
                IF( present( a7 ) ) p => a7(i1,i2,i3,i4,i5,i6,i7)
                IF( present( p7 ) ) p => p7(i1,i2,i3,i4,i5,i6,i7)
                IF( present( a7 ) .OR. present( p7 ) ) THEN
                  WRITE( UNIT = str, FMT = "(A,'(',I4,',',I4,',',I4,',',I4,',',I4,',',I4,',',I4,')')" ) trim( name ),& 
                   & i1, i2, i3, i4, i5, i6, i7
                ENDIF
                CALL chain_real_JPRB( err, chn, trim( str ) // "%CFRACTION", a0 = p%cfraction )
                IF( err .NE. 0 ) RETURN
                CALL chain_real_JPRB( err, chn, trim( str ) // "%CTP", a0 = p%ctp )
                IF( err .NE. 0 ) RETURN
                CALL chain_rttov_s2m( err, chn, trim( str ) // "%S2M", a0 = p%s2m )
                IF( err .NE. 0 ) RETURN
                CALL chain_rttov_skin( err, chn, trim( str ) // "%SKIN", a0 = p%skin )
                IF( err .NE. 0 ) RETURN
                CALL chain_real_JPRB( err, chn, trim( str ) // "%CFRAC", p1 = p%cfrac )
                IF( err .NE. 0 ) RETURN
                CALL chain_real_JPRB( err, chn, trim( str ) // "%CLOUD", p2 = p%cloud )
                IF( err .NE. 0 ) RETURN
                CALL chain_real_JPRB( err, chn, trim( str ) // "%CLWDE", p1 = p%clwde )
                IF( err .NE. 0 ) RETURN
                CALL chain_real_JPRB( err, chn, trim( str ) // "%ICEDE", p1 = p%icede )
                IF( err .NE. 0 ) RETURN
                CALL chain_real_JPRB( err, chn, trim( str ) // "%AEROSOLS", p2 = p%aerosols )
                IF( err .NE. 0 ) RETURN
                CALL chain_real_JPRB( err, chn, trim( str ) // "%CLW", p1 = p%clw )
                IF( err .NE. 0 ) RETURN
                CALL chain_real_JPRB( err, chn, trim( str ) // "%CH4", p1 = p%ch4 )
                IF( err .NE. 0 ) RETURN
                CALL chain_real_JPRB( err, chn, trim( str ) // "%SO2", p1 = p%so2 )
                IF( err .NE. 0 ) RETURN
                CALL chain_real_JPRB( err, chn, trim( str ) // "%CO", p1 = p%co )
                IF( err .NE. 0 ) RETURN
                CALL chain_real_JPRB( err, chn, trim( str ) // "%N2O", p1 = p%n2o )
                IF( err .NE. 0 ) RETURN
                CALL chain_real_JPRB( err, chn, trim( str ) // "%CO2", p1 = p%co2 )
                IF( err .NE. 0 ) RETURN
                CALL chain_real_JPRB( err, chn, trim( str ) // "%O3", p1 = p%o3 )
                IF( err .NE. 0 ) RETURN
                CALL chain_real_JPRB( err, chn, trim( str ) // "%Q", p1 = p%q )
                IF( err .NE. 0 ) RETURN
                CALL chain_real_JPRB( err, chn, trim( str ) // "%T", p1 = p%t )
                IF( err .NE. 0 ) RETURN
                CALL chain_real_JPRB( err, chn, trim( str ) // "%P", p1 = p%p )
                IF( err .NE. 0 ) RETURN
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  RETURN
  999 CONTINUE
  err = errorstatus_fatal
END SUBROUTINE chain_rttov_profile
SUBROUTINE chain_real_JPRB( err, chn, name, a0, a1, a2, a3, a4, a5, a6, a7, p0, p1, p2, p3, p4, p5, p6, p7 )
  INTEGER(KIND=JPIM),                 INTENT(OUT)   :: err
  TYPE(chain),                        INTENT(INOUT) :: chn
  CHARACTER(LEN=*),                   INTENT(IN)    :: name
  REAL(KIND=JPRB), OPTIONAL, TARGET,  INTENT(IN)    :: a0
  REAL(KIND=JPRB), OPTIONAL, TARGET,  INTENT(IN)    :: a1(:)
  REAL(KIND=JPRB), OPTIONAL, TARGET,  INTENT(IN)    :: a2(:,:)
  REAL(KIND=JPRB), OPTIONAL, TARGET,  INTENT(IN)    :: a3(:,:,:)
  REAL(KIND=JPRB), OPTIONAL, TARGET,  INTENT(IN)    :: a4(:,:,:,:)
  REAL(KIND=JPRB), OPTIONAL, TARGET,  INTENT(IN)    :: a5(:,:,:,:,:)
  REAL(KIND=JPRB), OPTIONAL, TARGET,  INTENT(IN)    :: a6(:,:,:,:,:,:)
  REAL(KIND=JPRB), OPTIONAL, TARGET,  INTENT(IN)    :: a7(:,:,:,:,:,:,:)
  REAL(KIND=JPRB), OPTIONAL, POINTER                :: p0
  REAL(KIND=JPRB), OPTIONAL, POINTER                :: p1(:)
  REAL(KIND=JPRB), OPTIONAL, POINTER                :: p2(:,:)
  REAL(KIND=JPRB), OPTIONAL, POINTER                :: p3(:,:,:)
  REAL(KIND=JPRB), OPTIONAL, POINTER                :: p4(:,:,:,:)
  REAL(KIND=JPRB), OPTIONAL, POINTER                :: p5(:,:,:,:,:)
  REAL(KIND=JPRB), OPTIONAL, POINTER                :: p6(:,:,:,:,:,:)
  REAL(KIND=JPRB), OPTIONAL, POINTER                :: p7(:,:,:,:,:,:,:)
  INTEGER(KIND=JPIM)                                :: nargs
  TYPE(chain),               POINTER                :: next
  err   = 0
  nargs = 0
  IF( present( a0 ) ) nargs = nargs + 1
  IF( present( a1 ) ) nargs = nargs + 1
  IF( present( a2 ) ) nargs = nargs + 1
  IF( present( a3 ) ) nargs = nargs + 1
  IF( present( a4 ) ) nargs = nargs + 1
  IF( present( a5 ) ) nargs = nargs + 1
  IF( present( a6 ) ) nargs = nargs + 1
  IF( present( a7 ) ) nargs = nargs + 1
  IF( present( p0 ) ) nargs = nargs + 1
  IF( present( p1 ) ) nargs = nargs + 1
  IF( present( p2 ) ) nargs = nargs + 1
  IF( present( p3 ) ) nargs = nargs + 1
  IF( present( p4 ) ) nargs = nargs + 1
  IF( present( p5 ) ) nargs = nargs + 1
  IF( present( p6 ) ) nargs = nargs + 1
  IF( present( p7 ) ) nargs = nargs + 1
  IF( nargs .NE. 1 ) GOTO 999
  IF( present( p0 ) ) THEN
    IF( .NOT. associated( p0 ) ) RETURN
  ENDIF
  IF( present( p1 ) ) THEN
    IF( .NOT. associated( p1 ) ) RETURN
  ENDIF
  IF( present( p2 ) ) THEN
    IF( .NOT. associated( p2 ) ) RETURN
  ENDIF
  IF( present( p3 ) ) THEN
    IF( .NOT. associated( p3 ) ) RETURN
  ENDIF
  IF( present( p4 ) ) THEN
    IF( .NOT. associated( p4 ) ) RETURN
  ENDIF
  IF( present( p5 ) ) THEN
    IF( .NOT. associated( p5 ) ) RETURN
  ENDIF
  IF( present( p6 ) ) THEN
    IF( .NOT. associated( p6 ) ) RETURN
  ENDIF
  IF( present( p7 ) ) THEN
    IF( .NOT. associated( p7 ) ) RETURN
  ENDIF
  ALLOCATE( next )
  next = chn
  CALL clear_chain( chn )
  chn%next => next
  chn%name = name
  IF( present( a0 ) ) chn%r_8_0 => a0
  IF( present( a1 ) ) chn%r_8_1 => a1
  IF( present( a2 ) ) chn%r_8_2 => a2
  IF( present( a3 ) ) chn%r_8_3 => a3
  IF( present( a4 ) ) chn%r_8_4 => a4
  IF( present( a5 ) ) chn%r_8_5 => a5
  IF( present( a6 ) ) chn%r_8_6 => a6
  IF( present( a7 ) ) chn%r_8_7 => a7
  IF( present( p0 ) ) chn%r_8_0 => p0
  IF( present( p1 ) ) chn%r_8_1 => p1
  IF( present( p2 ) ) chn%r_8_2 => p2
  IF( present( p3 ) ) chn%r_8_3 => p3
  IF( present( p4 ) ) chn%r_8_4 => p4
  IF( present( p5 ) ) chn%r_8_5 => p5
  IF( present( p6 ) ) chn%r_8_6 => p6
  IF( present( p7 ) ) chn%r_8_7 => p7
  RETURN
  999 CONTINUE
  err = errorstatus_fatal
END SUBROUTINE chain_real_JPRB
SUBROUTINE chain_rttov_skin( err, chn, name, a0, a1, a2, a3, a4, a5, a6, a7, p0, p1, p2, p3, p4, p5, p6, p7 )
  INTEGER(KIND=JPIM),                  INTENT(OUT)   :: err
  TYPE(chain),                         INTENT(INOUT) :: chn
  CHARACTER(LEN=*),                    INTENT(IN)    :: name
  TYPE(rttov_skin), OPTIONAL, TARGET,  INTENT(IN)    :: a0
  TYPE(rttov_skin), OPTIONAL, TARGET,  INTENT(IN)    :: a1(:)
  TYPE(rttov_skin), OPTIONAL, TARGET,  INTENT(IN)    :: a2(:,:)
  TYPE(rttov_skin), OPTIONAL, TARGET,  INTENT(IN)    :: a3(:,:,:)
  TYPE(rttov_skin), OPTIONAL, TARGET,  INTENT(IN)    :: a4(:,:,:,:)
  TYPE(rttov_skin), OPTIONAL, TARGET,  INTENT(IN)    :: a5(:,:,:,:,:)
  TYPE(rttov_skin), OPTIONAL, TARGET,  INTENT(IN)    :: a6(:,:,:,:,:,:)
  TYPE(rttov_skin), OPTIONAL, TARGET,  INTENT(IN)    :: a7(:,:,:,:,:,:,:)
  TYPE(rttov_skin), OPTIONAL, POINTER                :: p0
  TYPE(rttov_skin), OPTIONAL, POINTER                :: p1(:)
  TYPE(rttov_skin), OPTIONAL, POINTER                :: p2(:,:)
  TYPE(rttov_skin), OPTIONAL, POINTER                :: p3(:,:,:)
  TYPE(rttov_skin), OPTIONAL, POINTER                :: p4(:,:,:,:)
  TYPE(rttov_skin), OPTIONAL, POINTER                :: p5(:,:,:,:,:)
  TYPE(rttov_skin), OPTIONAL, POINTER                :: p6(:,:,:,:,:,:)
  TYPE(rttov_skin), OPTIONAL, POINTER                :: p7(:,:,:,:,:,:,:)
  INTEGER(KIND=JPIM)                                 :: nargs
  INTEGER(KIND=JPIM)                                 :: lbounda(7)
  INTEGER(KIND=JPIM)                                 :: ubounda(7)
  INTEGER(KIND=JPIM)                                 :: i1
  INTEGER(KIND=JPIM)                                 :: i2
  INTEGER(KIND=JPIM)                                 :: i3
  INTEGER(KIND=JPIM)                                 :: i4
  INTEGER(KIND=JPIM)                                 :: i5
  INTEGER(KIND=JPIM)                                 :: i6
  INTEGER(KIND=JPIM)                                 :: i7
  TYPE(rttov_skin),           POINTER                :: p
  CHARACTER(LEN=128)                                 :: str
  err   = 0
  nargs = 0
  IF( present( a0 ) ) nargs = nargs + 1
  IF( present( a1 ) ) nargs = nargs + 1
  IF( present( a2 ) ) nargs = nargs + 1
  IF( present( a3 ) ) nargs = nargs + 1
  IF( present( a4 ) ) nargs = nargs + 1
  IF( present( a5 ) ) nargs = nargs + 1
  IF( present( a6 ) ) nargs = nargs + 1
  IF( present( a7 ) ) nargs = nargs + 1
  IF( present( p0 ) ) nargs = nargs + 1
  IF( present( p1 ) ) nargs = nargs + 1
  IF( present( p2 ) ) nargs = nargs + 1
  IF( present( p3 ) ) nargs = nargs + 1
  IF( present( p4 ) ) nargs = nargs + 1
  IF( present( p5 ) ) nargs = nargs + 1
  IF( present( p6 ) ) nargs = nargs + 1
  IF( present( p7 ) ) nargs = nargs + 1
  IF( nargs .NE. 1 ) GOTO 999
  lbounda(:) = 0
  ubounda(:) = 0
  IF( present( a1 ) ) THEN
    lbounda(1:1) = lbound( a1 )
    ubounda(1:1) = ubound( a1 )
  ENDIF
  IF( present( a2 ) ) THEN
    lbounda(1:2) = lbound( a2 )
    ubounda(1:2) = ubound( a2 )
  ENDIF
  IF( present( a3 ) ) THEN
    lbounda(1:3) = lbound( a3 )
    ubounda(1:3) = ubound( a3 )
  ENDIF
  IF( present( a4 ) ) THEN
    lbounda(1:4) = lbound( a4 )
    ubounda(1:4) = ubound( a4 )
  ENDIF
  IF( present( a5 ) ) THEN
    lbounda(1:5) = lbound( a5 )
    ubounda(1:5) = ubound( a5 )
  ENDIF
  IF( present( a6 ) ) THEN
    lbounda(1:6) = lbound( a6 )
    ubounda(1:6) = ubound( a6 )
  ENDIF
  IF( present( a7 ) ) THEN
    lbounda(1:7) = lbound( a7 )
    ubounda(1:7) = ubound( a7 )
  ENDIF
  IF( present( p1 ) ) THEN
    IF( associated( p1 ) ) THEN
      lbounda(1:1) = lbound( p1 )
      ubounda(1:1) = ubound( p1 )
    ENDIF
  ENDIF
  IF( present( p2 ) ) THEN
    IF( associated( p2 ) ) THEN
      lbounda(1:2) = lbound( p2 )
      ubounda(1:2) = ubound( p2 )
    ENDIF
  ENDIF
  IF( present( p3 ) ) THEN
    IF( associated( p3 ) ) THEN
      lbounda(1:3) = lbound( p3 )
      ubounda(1:3) = ubound( p3 )
    ENDIF
  ENDIF
  IF( present( p4 ) ) THEN
    IF( associated( p4 ) ) THEN
      lbounda(1:4) = lbound( p4 )
      ubounda(1:4) = ubound( p4 )
    ENDIF
  ENDIF
  IF( present( p5 ) ) THEN
    IF( associated( p5 ) ) THEN
      lbounda(1:5) = lbound( p5 )
      ubounda(1:5) = ubound( p5 )
    ENDIF
  ENDIF
  IF( present( p6 ) ) THEN
    IF( associated( p6 ) ) THEN
      lbounda(1:6) = lbound( p6 )
      ubounda(1:6) = ubound( p6 )
    ENDIF
  ENDIF
  IF( present( p7 ) ) THEN
    IF( associated( p7 ) ) THEN
      lbounda(1:7) = lbound( p7 )
      ubounda(1:7) = ubound( p7 )
    ENDIF
  ENDIF
  IF( present( a0 ) ) p => a0
  IF( present( p0 ) ) p => p0
  IF( present( a0 ) .OR. present( p0 ) ) THEN
    WRITE( UNIT = str, FMT = "(A)" ) trim( name )
  ENDIF
  DO i1 = ubounda(1), lbounda(1), -1
    IF( present( a1 ) ) p => a1(i1)
    IF( present( p1 ) ) p => p1(i1)
    IF( present( a1 ) .OR. present( p1 ) ) THEN
      WRITE( UNIT = str, FMT = "(A,'(',I4,')')" ) trim( name ), i1
    ENDIF
    DO i2 = ubounda(2), lbounda(2), -1
      IF( present( a2 ) ) p => a2(i1,i2)
      IF( present( p2 ) ) p => p2(i1,i2)
      IF( present( a2 ) .OR. present( p2 ) ) THEN
        WRITE( UNIT = str, FMT = "(A,'(',I4,',',I4,')')" ) trim( name ), i1, i2
      ENDIF
      DO i3 = ubounda(3), lbounda(3), -1
        IF( present( a3 ) ) p => a3(i1,i2,i3)
        IF( present( p3 ) ) p => p3(i1,i2,i3)
        IF( present( a3 ) .OR. present( p3 ) ) THEN
          WRITE( UNIT = str, FMT = "(A,'(',I4,',',I4,',',I4,')')" ) trim( name ), i1, i2, i3
        ENDIF
        DO i4 = ubounda(4), lbounda(4), -1
          IF( present( a4 ) ) p => a4(i1,i2,i3,i4)
          IF( present( p4 ) ) p => p4(i1,i2,i3,i4)
          IF( present( a4 ) .OR. present( p4 ) ) THEN
            WRITE( UNIT = str, FMT = "(A,'(',I4,',',I4,',',I4,',',I4,')')" ) trim( name ), i1, i2, i3, i4
          ENDIF
          DO i5 = ubounda(5), lbounda(5), -1
            IF( present( a5 ) ) p => a5(i1,i2,i3,i4,i5)
            IF( present( p5 ) ) p => p5(i1,i2,i3,i4,i5)
            IF( present( a5 ) .OR. present( p5 ) ) THEN
              WRITE( UNIT = str, FMT = "(A,'(',I4,',',I4,',',I4,',',I4,',',I4,')')" ) trim( name ), i1, i2, i3, i4, i5
            ENDIF
            DO i6 = ubounda(6), lbounda(6), -1
              IF( present( a6 ) ) p => a6(i1,i2,i3,i4,i5,i6)
              IF( present( p6 ) ) p => p6(i1,i2,i3,i4,i5,i6)
              IF( present( a6 ) .OR. present( p6 ) ) THEN
                WRITE( UNIT = str, FMT = "(A,'(',I4,',',I4,',',I4,',',I4,',',I4,',',I4,')')" ) trim( name ), i1, i2, i3, i4, i5, i6
              ENDIF
              DO i7 = ubounda(7), lbounda(7), -1
                IF( present( a7 ) ) p => a7(i1,i2,i3,i4,i5,i6,i7)
                IF( present( p7 ) ) p => p7(i1,i2,i3,i4,i5,i6,i7)
                IF( present( a7 ) .OR. present( p7 ) ) THEN
                  WRITE( UNIT = str, FMT = "(A,'(',I4,',',I4,',',I4,',',I4,',',I4,',',I4,',',I4,')')" ) trim( name ),& 
                   & i1, i2, i3, i4, i5, i6, i7
                ENDIF
                CALL chain_real_JPRB( err, chn, trim( str ) // "%FASTEM", a1 = p%fastem )
                IF( err .NE. 0 ) RETURN
                CALL chain_real_JPRB( err, chn, trim( str ) // "%FOAM_FRACTION", a0 = p%foam_fraction )
                IF( err .NE. 0 ) RETURN
                CALL chain_real_JPRB( err, chn, trim( str ) // "%SALINITY", a0 = p%salinity )
                IF( err .NE. 0 ) RETURN
                CALL chain_real_JPRB( err, chn, trim( str ) // "%T", a0 = p%t )
                IF( err .NE. 0 ) RETURN
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  RETURN
  999 CONTINUE
  err = errorstatus_fatal
END SUBROUTINE chain_rttov_skin
SUBROUTINE chain_rttov_s2m( err, chn, name, a0, a1, a2, a3, a4, a5, a6, a7, p0, p1, p2, p3, p4, p5, p6, p7 )
  INTEGER(KIND=JPIM),                INTENT(OUT)   :: err
  TYPE(chain),                       INTENT(INOUT) :: chn
  CHARACTER(LEN=*),                  INTENT(IN)    :: name
  TYPE(rttov_s2m), OPTIONAL, TARGET,  INTENT(IN)   :: a0
  TYPE(rttov_s2m), OPTIONAL, TARGET,  INTENT(IN)   :: a1(:)
  TYPE(rttov_s2m), OPTIONAL, TARGET,  INTENT(IN)   :: a2(:,:)
  TYPE(rttov_s2m), OPTIONAL, TARGET,  INTENT(IN)   :: a3(:,:,:)
  TYPE(rttov_s2m), OPTIONAL, TARGET,  INTENT(IN)   :: a4(:,:,:,:)
  TYPE(rttov_s2m), OPTIONAL, TARGET,  INTENT(IN)   :: a5(:,:,:,:,:)
  TYPE(rttov_s2m), OPTIONAL, TARGET,  INTENT(IN)   :: a6(:,:,:,:,:,:)
  TYPE(rttov_s2m), OPTIONAL, TARGET,  INTENT(IN)   :: a7(:,:,:,:,:,:,:)
  TYPE(rttov_s2m), OPTIONAL, POINTER               :: p0
  TYPE(rttov_s2m), OPTIONAL, POINTER               :: p1(:)
  TYPE(rttov_s2m), OPTIONAL, POINTER               :: p2(:,:)
  TYPE(rttov_s2m), OPTIONAL, POINTER               :: p3(:,:,:)
  TYPE(rttov_s2m), OPTIONAL, POINTER               :: p4(:,:,:,:)
  TYPE(rttov_s2m), OPTIONAL, POINTER               :: p5(:,:,:,:,:)
  TYPE(rttov_s2m), OPTIONAL, POINTER               :: p6(:,:,:,:,:,:)
  TYPE(rttov_s2m), OPTIONAL, POINTER               :: p7(:,:,:,:,:,:,:)
  INTEGER(KIND=JPIM)                               :: nargs
  INTEGER(KIND=JPIM)                               :: lbounda(7)
  INTEGER(KIND=JPIM)                               :: ubounda(7)
  INTEGER(KIND=JPIM)                               :: i1
  INTEGER(KIND=JPIM)                               :: i2
  INTEGER(KIND=JPIM)                               :: i3
  INTEGER(KIND=JPIM)                               :: i4
  INTEGER(KIND=JPIM)                               :: i5
  INTEGER(KIND=JPIM)                               :: i6
  INTEGER(KIND=JPIM)                               :: i7
  TYPE(rttov_s2m),           POINTER               :: p
  CHARACTER(LEN=128)                               :: str
  err   = 0
  nargs = 0
  IF( present( a0 ) ) nargs = nargs + 1
  IF( present( a1 ) ) nargs = nargs + 1
  IF( present( a2 ) ) nargs = nargs + 1
  IF( present( a3 ) ) nargs = nargs + 1
  IF( present( a4 ) ) nargs = nargs + 1
  IF( present( a5 ) ) nargs = nargs + 1
  IF( present( a6 ) ) nargs = nargs + 1
  IF( present( a7 ) ) nargs = nargs + 1
  IF( present( p0 ) ) nargs = nargs + 1
  IF( present( p1 ) ) nargs = nargs + 1
  IF( present( p2 ) ) nargs = nargs + 1
  IF( present( p3 ) ) nargs = nargs + 1
  IF( present( p4 ) ) nargs = nargs + 1
  IF( present( p5 ) ) nargs = nargs + 1
  IF( present( p6 ) ) nargs = nargs + 1
  IF( present( p7 ) ) nargs = nargs + 1
  IF( nargs .NE. 1 ) GOTO 999
  lbounda(:) = 0
  ubounda(:) = 0
  IF( present( a1 ) ) THEN
    lbounda(1:1) = lbound( a1 )
    ubounda(1:1) = ubound( a1 )
  ENDIF
  IF( present( a2 ) ) THEN
    lbounda(1:2) = lbound( a2 )
    ubounda(1:2) = ubound( a2 )
  ENDIF
  IF( present( a3 ) ) THEN
    lbounda(1:3) = lbound( a3 )
    ubounda(1:3) = ubound( a3 )
  ENDIF
  IF( present( a4 ) ) THEN
    lbounda(1:4) = lbound( a4 )
    ubounda(1:4) = ubound( a4 )
  ENDIF
  IF( present( a5 ) ) THEN
    lbounda(1:5) = lbound( a5 )
    ubounda(1:5) = ubound( a5 )
  ENDIF
  IF( present( a6 ) ) THEN
    lbounda(1:6) = lbound( a6 )
    ubounda(1:6) = ubound( a6 )
  ENDIF
  IF( present( a7 ) ) THEN
    lbounda(1:7) = lbound( a7 )
    ubounda(1:7) = ubound( a7 )
  ENDIF
  IF( present( p1 ) ) THEN
    IF( associated( p1 ) ) THEN
      lbounda(1:1) = lbound( p1 )
      ubounda(1:1) = ubound( p1 )
    ENDIF
  ENDIF
  IF( present( p2 ) ) THEN
    IF( associated( p2 ) ) THEN
      lbounda(1:2) = lbound( p2 )
      ubounda(1:2) = ubound( p2 )
    ENDIF
  ENDIF
  IF( present( p3 ) ) THEN
    IF( associated( p3 ) ) THEN
      lbounda(1:3) = lbound( p3 )
      ubounda(1:3) = ubound( p3 )
    ENDIF
  ENDIF
  IF( present( p4 ) ) THEN
    IF( associated( p4 ) ) THEN
      lbounda(1:4) = lbound( p4 )
      ubounda(1:4) = ubound( p4 )
    ENDIF
  ENDIF
  IF( present( p5 ) ) THEN
    IF( associated( p5 ) ) THEN
      lbounda(1:5) = lbound( p5 )
      ubounda(1:5) = ubound( p5 )
    ENDIF
  ENDIF
  IF( present( p6 ) ) THEN
    IF( associated( p6 ) ) THEN
      lbounda(1:6) = lbound( p6 )
      ubounda(1:6) = ubound( p6 )
    ENDIF
  ENDIF
  IF( present( p7 ) ) THEN
    IF( associated( p7 ) ) THEN
      lbounda(1:7) = lbound( p7 )
      ubounda(1:7) = ubound( p7 )
    ENDIF
  ENDIF
  IF( present( a0 ) ) p => a0
  IF( present( p0 ) ) p => p0
  IF( present( a0 ) .OR. present( p0 ) ) THEN
    WRITE( UNIT = str, FMT = "(A)" ) trim( name )
  ENDIF
  DO i1 = ubounda(1), lbounda(1), -1
    IF( present( a1 ) ) p => a1(i1)
    IF( present( p1 ) ) p => p1(i1)
    IF( present( a1 ) .OR. present( p1 ) ) THEN
      WRITE( UNIT = str, FMT = "(A,'(',I4,')')" ) trim( name ), i1
    ENDIF
    DO i2 = ubounda(2), lbounda(2), -1
      IF( present( a2 ) ) p => a2(i1,i2)
      IF( present( p2 ) ) p => p2(i1,i2)
      IF( present( a2 ) .OR. present( p2 ) ) THEN
        WRITE( UNIT = str, FMT = "(A,'(',I4,',',I4,')')" ) trim( name ), i1, i2
      ENDIF
      DO i3 = ubounda(3), lbounda(3), -1
        IF( present( a3 ) ) p => a3(i1,i2,i3)
        IF( present( p3 ) ) p => p3(i1,i2,i3)
        IF( present( a3 ) .OR. present( p3 ) ) THEN
          WRITE( UNIT = str, FMT = "(A,'(',I4,',',I4,',',I4,')')" ) trim( name ), i1, i2, i3
        ENDIF
        DO i4 = ubounda(4), lbounda(4), -1
          IF( present( a4 ) ) p => a4(i1,i2,i3,i4)
          IF( present( p4 ) ) p => p4(i1,i2,i3,i4)
          IF( present( a4 ) .OR. present( p4 ) ) THEN
            WRITE( UNIT = str, FMT = "(A,'(',I4,',',I4,',',I4,',',I4,')')" ) trim( name ), i1, i2, i3, i4
          ENDIF
          DO i5 = ubounda(5), lbounda(5), -1
            IF( present( a5 ) ) p => a5(i1,i2,i3,i4,i5)
            IF( present( p5 ) ) p => p5(i1,i2,i3,i4,i5)
            IF( present( a5 ) .OR. present( p5 ) ) THEN
              WRITE( UNIT = str, FMT = "(A,'(',I4,',',I4,',',I4,',',I4,',',I4,')')" ) trim( name ), i1, i2, i3, i4, i5
            ENDIF
            DO i6 = ubounda(6), lbounda(6), -1
              IF( present( a6 ) ) p => a6(i1,i2,i3,i4,i5,i6)
              IF( present( p6 ) ) p => p6(i1,i2,i3,i4,i5,i6)
              IF( present( a6 ) .OR. present( p6 ) ) THEN
                WRITE( UNIT = str, FMT = "(A,'(',I4,',',I4,',',I4,',',I4,',',I4,',',I4,')')" ) trim( name ), i1, i2, i3, i4, i5, i6
              ENDIF
              DO i7 = ubounda(7), lbounda(7), -1
                IF( present( a7 ) ) p => a7(i1,i2,i3,i4,i5,i6,i7)
                IF( present( p7 ) ) p => p7(i1,i2,i3,i4,i5,i6,i7)
                IF( present( a7 ) .OR. present( p7 ) ) THEN
                  WRITE( UNIT = str, FMT = "(A,'(',I4,',',I4,',',I4,',',I4,',',I4,',',I4,',',I4,')')" ) trim( name ),&
                   & i1, i2, i3, i4, i5, i6, i7
                ENDIF
                CALL chain_real_JPRB( err, chn, trim( str ) // "%WFETC", a0 = p%wfetc )
                IF( err .NE. 0 ) RETURN
                CALL chain_real_JPRB( err, chn, trim( str ) // "%V", a0 = p%v )
                IF( err .NE. 0 ) RETURN
                CALL chain_real_JPRB( err, chn, trim( str ) // "%U", a0 = p%u )
                IF( err .NE. 0 ) RETURN
                CALL chain_real_JPRB( err, chn, trim( str ) // "%P", a0 = p%p )
                IF( err .NE. 0 ) RETURN
                CALL chain_real_JPRB( err, chn, trim( str ) // "%O", a0 = p%o )
                IF( err .NE. 0 ) RETURN
                CALL chain_real_JPRB( err, chn, trim( str ) // "%Q", a0 = p%q )
                IF( err .NE. 0 ) RETURN
                CALL chain_real_JPRB( err, chn, trim( str ) // "%T", a0 = p%t )
                IF( err .NE. 0 ) RETURN
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  RETURN
  999 CONTINUE
  err = errorstatus_fatal
END SUBROUTINE chain_rttov_s2m
SUBROUTINE chain_rttov_profile_cloud( err, chn, name, a0, a1, a2, a3, a4, a5, a6, a7, p0, p1, p2, p3, p4, p5, p6, p7 )
  INTEGER(KIND=JPIM),                          INTENT(OUT)   :: err
  TYPE(chain),                                 INTENT(INOUT) :: chn
  CHARACTER(LEN=*),                            INTENT(IN)    :: name
  TYPE(rttov_profile_cloud), OPTIONAL, TARGET, INTENT(IN)    :: a0
  TYPE(rttov_profile_cloud), OPTIONAL, TARGET, INTENT(IN)    :: a1(:)
  TYPE(rttov_profile_cloud), OPTIONAL, TARGET, INTENT(IN)    :: a2(:,:)
  TYPE(rttov_profile_cloud), OPTIONAL, TARGET, INTENT(IN)    :: a3(:,:,:)
  TYPE(rttov_profile_cloud), OPTIONAL, TARGET, INTENT(IN)    :: a4(:,:,:,:)
  TYPE(rttov_profile_cloud), OPTIONAL, TARGET, INTENT(IN)    :: a5(:,:,:,:,:)
  TYPE(rttov_profile_cloud), OPTIONAL, TARGET, INTENT(IN)    :: a6(:,:,:,:,:,:)
  TYPE(rttov_profile_cloud), OPTIONAL, TARGET, INTENT(IN)    :: a7(:,:,:,:,:,:,:)
  TYPE(rttov_profile_cloud), OPTIONAL, POINTER               :: p0
  TYPE(rttov_profile_cloud), OPTIONAL, POINTER               :: p1(:)
  TYPE(rttov_profile_cloud), OPTIONAL, POINTER               :: p2(:,:)
  TYPE(rttov_profile_cloud), OPTIONAL, POINTER               :: p3(:,:,:)
  TYPE(rttov_profile_cloud), OPTIONAL, POINTER               :: p4(:,:,:,:)
  TYPE(rttov_profile_cloud), OPTIONAL, POINTER               :: p5(:,:,:,:,:)
  TYPE(rttov_profile_cloud), OPTIONAL, POINTER               :: p6(:,:,:,:,:,:)
  TYPE(rttov_profile_cloud), OPTIONAL, POINTER               :: p7(:,:,:,:,:,:,:)
  INTEGER(KIND=JPIM)                                         :: nargs
  INTEGER(KIND=JPIM)                                         :: lbounda(7)
  INTEGER(KIND=JPIM)                                         :: ubounda(7)
  INTEGER(KIND=JPIM)                                         :: i1
  INTEGER(KIND=JPIM)                                         :: i2
  INTEGER(KIND=JPIM)                                         :: i3
  INTEGER(KIND=JPIM)                                         :: i4
  INTEGER(KIND=JPIM)                                         :: i5
  INTEGER(KIND=JPIM)                                         :: i6
  INTEGER(KIND=JPIM)                                         :: i7
  TYPE(rttov_profile_cloud),           POINTER               :: p
  CHARACTER(LEN=128)                                         :: str
  err   = 0
  nargs = 0
  IF( present( a0 ) ) nargs = nargs + 1
  IF( present( a1 ) ) nargs = nargs + 1
  IF( present( a2 ) ) nargs = nargs + 1
  IF( present( a3 ) ) nargs = nargs + 1
  IF( present( a4 ) ) nargs = nargs + 1
  IF( present( a5 ) ) nargs = nargs + 1
  IF( present( a6 ) ) nargs = nargs + 1
  IF( present( a7 ) ) nargs = nargs + 1
  IF( present( p0 ) ) nargs = nargs + 1
  IF( present( p1 ) ) nargs = nargs + 1
  IF( present( p2 ) ) nargs = nargs + 1
  IF( present( p3 ) ) nargs = nargs + 1
  IF( present( p4 ) ) nargs = nargs + 1
  IF( present( p5 ) ) nargs = nargs + 1
  IF( present( p6 ) ) nargs = nargs + 1
  IF( present( p7 ) ) nargs = nargs + 1
  IF( nargs .NE. 1 ) GOTO 999
  lbounda(:) = 0
  ubounda(:) = 0
  IF( present( a1 ) ) THEN
    lbounda(1:1) = lbound( a1 )
    ubounda(1:1) = ubound( a1 )
  ENDIF
  IF( present( a2 ) ) THEN
    lbounda(1:2) = lbound( a2 )
    ubounda(1:2) = ubound( a2 )
  ENDIF
  IF( present( a3 ) ) THEN
    lbounda(1:3) = lbound( a3 )
    ubounda(1:3) = ubound( a3 )
  ENDIF
  IF( present( a4 ) ) THEN
    lbounda(1:4) = lbound( a4 )
    ubounda(1:4) = ubound( a4 )
  ENDIF
  IF( present( a5 ) ) THEN
    lbounda(1:5) = lbound( a5 )
    ubounda(1:5) = ubound( a5 )
  ENDIF
  IF( present( a6 ) ) THEN
    lbounda(1:6) = lbound( a6 )
    ubounda(1:6) = ubound( a6 )
  ENDIF
  IF( present( a7 ) ) THEN
    lbounda(1:7) = lbound( a7 )
    ubounda(1:7) = ubound( a7 )
  ENDIF
  IF( present( p1 ) ) THEN
    IF( associated( p1 ) ) THEN
      lbounda(1:1) = lbound( p1 )
      ubounda(1:1) = ubound( p1 )
    ENDIF
  ENDIF
  IF( present( p2 ) ) THEN
    IF( associated( p2 ) ) THEN
      lbounda(1:2) = lbound( p2 )
      ubounda(1:2) = ubound( p2 )
    ENDIF
  ENDIF
  IF( present( p3 ) ) THEN
    IF( associated( p3 ) ) THEN
      lbounda(1:3) = lbound( p3 )
      ubounda(1:3) = ubound( p3 )
    ENDIF
  ENDIF
  IF( present( p4 ) ) THEN
    IF( associated( p4 ) ) THEN
      lbounda(1:4) = lbound( p4 )
      ubounda(1:4) = ubound( p4 )
    ENDIF
  ENDIF
  IF( present( p5 ) ) THEN
    IF( associated( p5 ) ) THEN
      lbounda(1:5) = lbound( p5 )
      ubounda(1:5) = ubound( p5 )
    ENDIF
  ENDIF
  IF( present( p6 ) ) THEN
    IF( associated( p6 ) ) THEN
      lbounda(1:6) = lbound( p6 )
      ubounda(1:6) = ubound( p6 )
    ENDIF
  ENDIF
  IF( present( p7 ) ) THEN
    IF( associated( p7 ) ) THEN
      lbounda(1:7) = lbound( p7 )
      ubounda(1:7) = ubound( p7 )
    ENDIF
  ENDIF
  IF( present( a0 ) ) p => a0
  IF( present( p0 ) ) p => p0
  IF( present( a0 ) .OR. present( p0 ) ) THEN
    WRITE( UNIT = str, FMT = "(A)" ) trim( name )
  ENDIF
  DO i1 = ubounda(1), lbounda(1), -1
    IF( present( a1 ) ) p => a1(i1)
    IF( present( p1 ) ) p => p1(i1)
    IF( present( a1 ) .OR. present( p1 ) ) THEN
      WRITE( UNIT = str, FMT = "(A,'(',I4,')')" ) trim( name ), i1
    ENDIF
    DO i2 = ubounda(2), lbounda(2), -1
      IF( present( a2 ) ) p => a2(i1,i2)
      IF( present( p2 ) ) p => p2(i1,i2)
      IF( present( a2 ) .OR. present( p2 ) ) THEN
        WRITE( UNIT = str, FMT = "(A,'(',I4,',',I4,')')" ) trim( name ), i1, i2
      ENDIF
      DO i3 = ubounda(3), lbounda(3), -1
        IF( present( a3 ) ) p => a3(i1,i2,i3)
        IF( present( p3 ) ) p => p3(i1,i2,i3)
        IF( present( a3 ) .OR. present( p3 ) ) THEN
          WRITE( UNIT = str, FMT = "(A,'(',I4,',',I4,',',I4,')')" ) trim( name ), i1, i2, i3
        ENDIF
        DO i4 = ubounda(4), lbounda(4), -1
          IF( present( a4 ) ) p => a4(i1,i2,i3,i4)
          IF( present( p4 ) ) p => p4(i1,i2,i3,i4)
          IF( present( a4 ) .OR. present( p4 ) ) THEN
            WRITE( UNIT = str, FMT = "(A,'(',I4,',',I4,',',I4,',',I4,')')" ) trim( name ), i1, i2, i3, i4
          ENDIF
          DO i5 = ubounda(5), lbounda(5), -1
            IF( present( a5 ) ) p => a5(i1,i2,i3,i4,i5)
            IF( present( p5 ) ) p => p5(i1,i2,i3,i4,i5)
            IF( present( a5 ) .OR. present( p5 ) ) THEN
              WRITE( UNIT = str, FMT = "(A,'(',I4,',',I4,',',I4,',',I4,',',I4,')')" ) trim( name ), i1, i2, i3, i4, i5
            ENDIF
            DO i6 = ubounda(6), lbounda(6), -1
              IF( present( a6 ) ) p => a6(i1,i2,i3,i4,i5,i6)
              IF( present( p6 ) ) p => p6(i1,i2,i3,i4,i5,i6)
              IF( present( a6 ) .OR. present( p6 ) ) THEN
                WRITE( UNIT = str, FMT = "(A,'(',I4,',',I4,',',I4,',',I4,',',I4,',',I4,')')" ) trim( name ), i1, i2, i3, i4, i5, i6
              ENDIF
              DO i7 = ubounda(7), lbounda(7), -1
                IF( present( a7 ) ) p => a7(i1,i2,i3,i4,i5,i6,i7)
                IF( present( p7 ) ) p => p7(i1,i2,i3,i4,i5,i6,i7)
                IF( present( a7 ) .OR. present( p7 ) ) THEN
                  WRITE( UNIT = str, FMT = "(A,'(',I4,',',I4,',',I4,',',I4,',',I4,',',I4,',',I4,')')" ) trim( name ),&
                   & i1, i2, i3, i4, i5, i6, i7
                ENDIF
                CALL chain_real_JPRB( err, chn, trim( str ) // "%PH", p1 = p%ph )
                IF( err .NE. 0 ) RETURN
                CALL chain_real_JPRB( err, chn, trim( str ) // "%HYDRO", p2 = p%hydro )
                IF( err .NE. 0 ) RETURN
                CALL chain_real_JPRB( err, chn, trim( str ) // "%HYDRO_FRAC", p2 = p%hydro_frac )
                IF( err .NE. 0 ) RETURN
                CALL chain_real_JPRB( err, chn, trim( str ) // "%CFRAC", a0 = p%cfrac )
                IF( err .NE. 0 ) RETURN
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  RETURN
  999 CONTINUE
  err = errorstatus_fatal
END SUBROUTINE chain_rttov_profile_cloud

END MODULE
