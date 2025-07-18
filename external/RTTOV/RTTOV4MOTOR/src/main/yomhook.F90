MODULE YOMHOOK

USE PARKIND1  ,ONLY : JPIM, JPRB, JPRM, JPRD

IMPLICIT NONE

! Used by "hook" function
! LHOOK = true implies "hook" function will be called
! Altough initialized to TRUE it will be reset by first call to 
! DR_HOOK unless we really want to use the hook function

SAVE
PUBLIC

LOGICAL :: LHOOK=.TRUE.

!#include "dr_hook_util.h"
!#include "dr_hook_util_multi.h"

INTERFACE DR_HOOK  
! We want compile time mapping of DR_HOOK-arguments and not
! to test OPTIONAL-arguments with PRESENT()-function, since
! it costs more.
! However, this "unrolling" approach cannot be streched much more
! than this without making number of member-functions too large
! (i.e. all the possible permutations of these "optional" args;
!  arguments that are not present in the DR_HOOK_DEFAULT -version)

MODULE PROCEDURE &
  DR_HOOK_DEFAULT4, &
  DR_HOOK_DEFAULT8, &
  DR_HOOK_FILE, &
  DR_HOOK_SIZE, &
  DR_HOOK_FILE_SIZE, &
  DR_HOOK_MULTI_DEFAULT, &
  DR_HOOK_MULTI_FILE, &
  DR_HOOK_MULTI_SIZE, &
  DR_HOOK_MULTI_FILE_SIZE
END INTERFACE

CONTAINS 

SUBROUTINE DR_HOOK_DEFAULT4(CDNAME,KSWITCH,PKEY)
CHARACTER(LEN=*), INTENT(IN) :: CDNAME
INTEGER(KIND=JPIM),        INTENT(IN) :: KSWITCH
REAL(KIND=JPRM),        INTENT(INOUT) :: PKEY
REAL(KIND=JPRB) :: ZKEY
ZKEY = TRANSFER(PKEY,ZKEY)
CALL DR_HOOK_UTIL(LHOOK,CDNAME,KSWITCH,ZKEY,'',0_JPIM)
PKEY = TRANSFER(ZKEY,PKEY)
END SUBROUTINE DR_HOOK_DEFAULT4

SUBROUTINE DR_HOOK_DEFAULT8(CDNAME,KSWITCH,PKEY)
CHARACTER(LEN=*), INTENT(IN) :: CDNAME
INTEGER(KIND=JPIM),        INTENT(IN) :: KSWITCH
REAL(KIND=JPRD),        INTENT(INOUT) :: PKEY
REAL(KIND=JPRB) :: ZKEY
ZKEY = TRANSFER(PKEY,ZKEY)
CALL DR_HOOK_UTIL(LHOOK,CDNAME,KSWITCH,ZKEY,'',0_JPIM)
PKEY = TRANSFER(ZKEY,PKEY)
END SUBROUTINE DR_HOOK_DEFAULT8



SUBROUTINE DR_HOOK_MULTI_DEFAULT(CDNAME,KSWITCH,PKEY)
CHARACTER(LEN=*), INTENT(IN) :: CDNAME
INTEGER(KIND=JPIM),        INTENT(IN) :: KSWITCH
REAL(KIND=JPRB),        INTENT(INOUT) :: PKEY(:)
CALL DR_HOOK_UTIL_MULTI(LHOOK,CDNAME,KSWITCH,PKEY,INT(SIZE(PKEY)),'',0_JPIM)
END SUBROUTINE DR_HOOK_MULTI_DEFAULT



SUBROUTINE DR_HOOK_FILE(CDNAME,KSWITCH,PKEY,CDFILE)
CHARACTER(LEN=*), INTENT(IN) :: CDNAME,CDFILE
INTEGER(KIND=JPIM),        INTENT(IN) :: KSWITCH
REAL(KIND=JPRB),        INTENT(INOUT) :: PKEY
CALL DR_HOOK_UTIL(LHOOK,CDNAME,KSWITCH,PKEY,CDFILE,0_JPIM)
END SUBROUTINE DR_HOOK_FILE

SUBROUTINE DR_HOOK_MULTI_FILE(CDNAME,KSWITCH,PKEY,CDFILE)
CHARACTER(LEN=*), INTENT(IN) :: CDNAME,CDFILE
INTEGER(KIND=JPIM),        INTENT(IN) :: KSWITCH
REAL(KIND=JPRB),        INTENT(INOUT) :: PKEY(:)
CALL DR_HOOK_UTIL_MULTI(LHOOK,CDNAME,KSWITCH,PKEY,INT(SIZE(PKEY)),CDFILE,0_JPIM)
END SUBROUTINE DR_HOOK_MULTI_FILE



SUBROUTINE DR_HOOK_SIZE(CDNAME,KSWITCH,PKEY,KSIZEINFO)
CHARACTER(LEN=*), INTENT(IN) :: CDNAME
INTEGER(KIND=JPIM),        INTENT(IN) :: KSWITCH,KSIZEINFO
REAL(KIND=JPRB),        INTENT(INOUT) :: PKEY
CALL DR_HOOK_UTIL(LHOOK,CDNAME,KSWITCH,PKEY,'',KSIZEINFO)
END SUBROUTINE DR_HOOK_SIZE

SUBROUTINE DR_HOOK_MULTI_SIZE(CDNAME,KSWITCH,PKEY,KSIZEINFO)
CHARACTER(LEN=*), INTENT(IN) :: CDNAME
INTEGER(KIND=JPIM),        INTENT(IN) :: KSWITCH,KSIZEINFO
REAL(KIND=JPRB),        INTENT(INOUT) :: PKEY(:)
CALL DR_HOOK_UTIL_MULTI(LHOOK,CDNAME,KSWITCH,PKEY,INT(SIZE(PKEY)),'',KSIZEINFO)
END SUBROUTINE DR_HOOK_MULTI_SIZE



SUBROUTINE DR_HOOK_FILE_SIZE(CDNAME,KSWITCH,PKEY,CDFILE,KSIZEINFO)
CHARACTER(LEN=*), INTENT(IN) :: CDNAME,CDFILE
INTEGER(KIND=JPIM),        INTENT(IN) :: KSWITCH,KSIZEINFO
REAL(KIND=JPRB),        INTENT(INOUT) :: PKEY
CALL DR_HOOK_UTIL(LHOOK,CDNAME,KSWITCH,PKEY,CDFILE,KSIZEINFO)
END SUBROUTINE DR_HOOK_FILE_SIZE

SUBROUTINE DR_HOOK_MULTI_FILE_SIZE(CDNAME,KSWITCH,PKEY,CDFILE,KSIZEINFO)
CHARACTER(LEN=*), INTENT(IN) :: CDNAME,CDFILE
INTEGER(KIND=JPIM),        INTENT(IN) :: KSWITCH,KSIZEINFO
REAL(KIND=JPRB),        INTENT(INOUT) :: PKEY(:)
CALL DR_HOOK_UTIL_MULTI(LHOOK,CDNAME,KSWITCH,PKEY,INT(SIZE(PKEY)),CDFILE,KSIZEINFO)
END SUBROUTINE DR_HOOK_MULTI_FILE_SIZE

SUBROUTINE DR_HOOK_UTIL(LHOOK,CDNAME,KCASE,PKEY,CDFILENAME,KSIZEINFO)
USE PARKIND1  ,ONLY : JPIM     ,JPRB
IMPLICIT NONE
LOGICAL, INTENT(IN) :: LHOOK
CHARACTER(LEN=*),INTENT(IN) :: CDNAME,CDFILENAME
INTEGER(KIND=JPIM),INTENT(IN) :: KCASE,KSIZEINFO
REAL(KIND=JPRB),INTENT(INOUT) :: PKEY
PKEY=0._jprb
END SUBROUTINE DR_HOOK_UTIL

SUBROUTINE DR_HOOK_UTIL_MULTI(LHOOK,CDNAME,KCASE,PKEY,KPKEY,CDFILENAME,KSIZEINFO)
USE PARKIND1  ,ONLY : JPIM     ,JPRB
IMPLICIT NONE
LOGICAL, INTENT(IN) :: LHOOK
CHARACTER(LEN=*),INTENT(IN) :: CDNAME,CDFILENAME
INTEGER(KIND=JPIM),INTENT(IN) :: KPKEY, KCASE,KSIZEINFO
REAL(KIND=JPRB),INTENT(INOUT) :: PKEY(KPKEY)
PKEY=0._jprb
END SUBROUTINE DR_HOOK_UTIL_MULTI

END MODULE YOMHOOK
