!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA.IO.GrapesIO
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
! AUTOHR(S)         : Sanshan Tu
! VERSION           : V 0.0
! HISTORY           :
!   Created by Sanshan Tu (tss71618@163.com), 2021/10/29, @SZSC, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module contains the data type for model_states
!! method.
!! @copyright (C) 2020 GBA-MWF, All rights reserved.
! @note
! @warning
! @attention

MODULE EndianConverter_m
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !       SYNOPSIS: Converts a 32 bit, 4 byte, REAL or Integer from
  !       big
  !                 Endian to little Endian, or conversely from little
  !                 Endian
  !                 to big Endian.
  !                  --------    --------    --------    --------
  !                 |    D   |  |    C   |  |    B   |  |    A   |  4
  !                 Bytes
  !                  --------    --------    --------    --------
  !                                                             |
  !                                                              -> 1
  !                                                              bit
  !                                       ||
  !                                     MVBITS
  !                                       ||
  !                                       \/
  !
  !                  --------    --------    --------    --------
  !                 |    A   |  |    B   |  |    C   |  |    D   |  4
  !                 Bytes
  !                  --------    --------    --------    --------
  !                         |           |           |           |
  !                         24          16          8           0   <-
  !                         bit
  !                                                                 position
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTERFACE EndianConvert
    MODULE PROCEDURE EndianConvert_Integer4
    MODULE PROCEDURE EndianConvert_Integer4_Array1D
    MODULE PROCEDURE EndianConvert_Integer4_Array2D
    MODULE PROCEDURE EndianConvert_Real4
    MODULE PROCEDURE EndianConvert_Real4_Array1D
    MODULE PROCEDURE EndianConvert_Real4_Array2D
    MODULE PROCEDURE EndianConvert_Real4_Array3D
  END INTERFACE

CONTAINS

  FUNCTION Native_Integer4(IntegerIn) RESULT(IntegerOut)
    IMPLICIT NONE
    INTEGER(4), INTENT(IN)                              :: IntegerIn
    INTEGER(4)                                          :: IntegerOut
    CALL MVBITS(IntegerIn, 24, 8, IntegerOut, 0)
    CALL MVBITS(IntegerIn, 16, 8, IntegerOut, 8)
    CALL MVBITS(IntegerIn, 8, 8, IntegerOut, 16)
    CALL MVBITS(IntegerIn, 0, 8, IntegerOut, 24)

  END FUNCTION

  FUNCTION Native_Real4(realIn) RESULT(realOut)
    IMPLICIT NONE
    REAL, INTENT(IN)                              :: realIn
    REAL                                          :: realOut
    INTEGER(4)                                       :: i_element
    INTEGER(4)                                       :: i_element_br

    i_element = TRANSFER(realIn, 0)

    i_element_br = Native_Integer4(i_element)

    realOut = TRANSFER(i_element_br, 0.0)
  END FUNCTION

  SUBROUTINE EndianConvert_Integer4(IntegerVal)
    IMPLICIT NONE
    INTEGER(4), INTENT(INOUT)                 :: IntegerVal

    IntegerVal = Native_Integer4(IntegerVal)

  END SUBROUTINE

  SUBROUTINE EndianConvert_Integer4_Array1D(IntegerVals)
    IMPLICIT NONE
    INTEGER(4), INTENT(INOUT)                 :: IntegerVals(:)
    INTEGER(4) ::i
    DO i = LBOUND(IntegerVals, 1), UBOUND(IntegerVals, 1)
      CALL EndianConvert_Integer4(IntegerVals(i))
    END DO
  END SUBROUTINE

  SUBROUTINE EndianConvert_Integer4_Array2D(IntegerVals)
    IMPLICIT NONE
    INTEGER(4), INTENT(INOUT)                 :: IntegerVals(:, :)
    INTEGER(4) ::i
    DO i = LBOUND(IntegerVals, 2), UBOUND(IntegerVals, 2)
      CALL EndianConvert_Integer4_Array1D(IntegerVals(:, i))
    END DO
  END SUBROUTINE

  SUBROUTINE EndianConvert_Real4(RealVal)
    IMPLICIT NONE
    REAL(4), INTENT(INOUT)                 :: RealVal

    RealVal = Native_Real4(RealVal)

  END SUBROUTINE

  SUBROUTINE EndianConvert_Real4_Array1D(RealVals)
    IMPLICIT NONE
    REAL(4), INTENT(INOUT)                 :: RealVals(:)
    INTEGER(4) ::i
    DO i = LBOUND(RealVals, 1), UBOUND(RealVals, 1)
      CALL EndianConvert_Real4(RealVals(i))
    END DO
  END SUBROUTINE

  SUBROUTINE EndianConvert_Real4_Array2D(RealVals)
    IMPLICIT NONE
    REAL(4), INTENT(INOUT)                 :: RealVals(:, :)
    INTEGER(4) ::i
    DO i = LBOUND(RealVals, 2), UBOUND(RealVals, 2)
      CALL EndianConvert_Real4_Array1D(RealVals(:, i))
    END DO
  END SUBROUTINE

  SUBROUTINE EndianConvert_Real4_Array3D(RealVals)
    IMPLICIT NONE
    REAL(4), INTENT(INOUT)                 :: RealVals(:, :, :)
    INTEGER(4) ::i
    DO i = LBOUND(RealVals, 3), UBOUND(RealVals, 3)
      CALL EndianConvert_Real4_Array2D(RealVals(:, :, i))
    END DO
  END SUBROUTINE

END MODULE EndianConverter_m

