!----------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Centers
!                     for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Jiongming Pang
! VERSION           : V 0.0
! HISTORY           :
!   Created by Jiongming Pang (pang.j.m@hotmail.com), 2022-03-31, @GBA-MWF, Shenzhen
!----------------------------------------------------------------------------------------

MODULE Widgets_m
  USE kinds_m
  IMPLICIT NONE

CONTAINS

  FUNCTION getStrIndex(varStrs, usrStr)
    CHARACTER(LEN=*), INTENT(IN) :: varStrs(:)
    CHARACTER(LEN=*), INTENT(IN) :: usrStr
    INTEGER(i_kind) :: i
    INTEGER(i_kind) :: getStrIndex

    getStrIndex = -1
    DO i = 1, SIZE(varStrs)
      IF (TRIM(usrStr) == TRIM(varStrs(i))) THEN
        getStrIndex = i
        EXIT
      END IF
    END DO

  END FUNCTION getStrIndex

  FUNCTION getStrIndexs(varStrs, usrStrs)
    CHARACTER(LEN=*), INTENT(IN) :: varStrs(:)
    CHARACTER(LEN=*), INTENT(IN) :: usrStrs(:)
    INTEGER(i_kind) :: i, j
    INTEGER(i_kind), ALLOCATABLE :: getStrIndexs(:)

    ALLOCATE (getStrIndexs(SIZE(usrStrs)))
    getStrIndexs = -1
    DO j = 1, SIZE(usrStrs)
      DO i = 1, SIZE(varStrs)
        IF (TRIM(usrStrs(j)) == TRIM(varStrs(i))) THEN
          getStrIndexs(j) = i
          CYCLE
        END IF
      END DO
    END DO

  END FUNCTION getStrIndexs

END MODULE Widgets_m
