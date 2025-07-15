!!--------------------------------------------------------------------------------------------------
! PROJECT           : GRAPES IO
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Sanshan Tu
! VERSION           : Beta 0.0
! HISTORY           :
!   Created by Sanshan Tu (tss71618@163.com), 2021/3/12, @SZSC, Shenzhen
! Modified by Zilong Qin (zilong.qin@gmail.com), 2021/3/18, @GBA-MWF, Shenzhen
! Modified by Zilong Qin (zilong.qin@gmail.com), 2021/3/29, @GBA-MWF, Shenzhen
    !!--------------------------------------------------------------------------------------------------

    !!===================================================================
!> @brief
    !! # GRAPES IO Module
    !!
    !!  *This module defines data structures for GRAPES input namelist*
    !! @author Sanshan Tu
    !! @copyright (C) 2020 GBA-MWF, All rights reserved.
    !!===================================================================
MODULE NMLRead_m

  INTERFACE namelist_read
    MODULE PROCEDURE nml_read_int4
    MODULE PROCEDURE nml_read_int8
    MODULE PROCEDURE nml_read_real4
    MODULE PROCEDURE nml_read_real8
    MODULE PROCEDURE nml_read_bool
    MODULE PROCEDURE nml_read_string
    MODULE PROCEDURE nml_read_string_array
  END INTERFACE

  ! interface nml_get
  !   module procedure nml_get_int4
  !   module procedure nml_get_string
  ! end interface nml_get

CONTAINS

  SUBROUTINE nml_read_int4(nmlst_fn, varName, var)
    IMPLICIT NONE
    INTEGER*4              :: var
    CHARACTER(*)           :: nmlst_fn
    CHARACTER(*)           :: varName
    CHARACTER*256          :: rdVarStr

    INCLUDE 'namelist_read.inc'
    ! The following print is useless as it checks the variable value to be read: Yuanfu Xie turns it off for now
    ! PRINT *, 'nml_read_int4 - var: ', var
    IF (hasPara) READ (rdVarStr, *) var

  END SUBROUTINE nml_read_int4

  SUBROUTINE nml_read_int8(nmlst_fn, varName, var)
    IMPLICIT NONE
    INTEGER*8              :: var
    CHARACTER(*)           :: nmlst_fn
    CHARACTER(*)           :: varName
    CHARACTER*256          :: rdVarStr

    INCLUDE 'namelist_read.inc'
    READ (rdVarStr, *) var
  END SUBROUTINE nml_read_int8

  SUBROUTINE nml_read_real4(nmlst_fn, varName, var)
    IMPLICIT NONE
    REAL*4                 :: var
    CHARACTER(*)           :: nmlst_fn
    CHARACTER(*)           :: varName
    CHARACTER*256          :: rdVarStr

    INCLUDE 'namelist_read.inc'
    IF (hasPara) READ (rdVarStr, *) var
  END SUBROUTINE nml_read_real4

  SUBROUTINE nml_read_real8(nmlst_fn, varName, var)
    IMPLICIT NONE
    REAL*8                 :: var
    CHARACTER(*)           :: nmlst_fn
    CHARACTER(*)           :: varName
    CHARACTER*256          :: rdVarStr

    INCLUDE 'namelist_read.inc'
    IF (hasPara) READ (rdVarStr, *) var
  END SUBROUTINE nml_read_real8

  SUBROUTINE nml_read_bool(nmlst_fn, varName, var)
    IMPLICIT NONE
    LOGICAL              :: var
    CHARACTER(*)           :: nmlst_fn
    CHARACTER(*)           :: varName
    CHARACTER*256          :: rdVarStr

    INCLUDE 'namelist_read.inc'
    IF (hasPara) READ (rdVarStr, *) var
  END SUBROUTINE nml_read_bool

  SUBROUTINE nml_read_string(nmlst_fn, varName, var)
    IMPLICIT NONE
    CHARACTER(*) :: var
    CHARACTER(*)           :: nmlst_fn
    CHARACTER(*)           :: varName
    CHARACTER*1024          :: rdVarStr

    INCLUDE 'namelist_read.inc'
    IF (hasPara) READ (rdVarStr, *) var
  END SUBROUTINE nml_read_string

  FUNCTION nml_get_string(nmlst_fn, varName) RESULT(var)
    IMPLICIT NONE
    CHARACTER(LEN=1024) :: var
    CHARACTER(*)           :: nmlst_fn
    CHARACTER(*)           :: varName

    CALL nml_read_string(nmlst_fn, varName, var)
  END FUNCTION nml_get_string

  FUNCTION nml_get_int4(nmlst_fn, varName) RESULT(var)
    IMPLICIT NONE
    INTEGER*4              :: var
    CHARACTER(*)           :: nmlst_fn
    CHARACTER(*)           :: varName

    CALL nml_read_int4(nmlst_fn, varName, var)
  END FUNCTION nml_get_int4

  SUBROUTINE nml_read_string_array(nmlst_fn, varName, var)
    IMPLICIT NONE
    CHARACTER(*), ALLOCATABLE           :: var(:)
    CHARACTER(*)           :: nmlst_fn
    CHARACTER(*)           :: varName
    CHARACTER*1024          :: rdVarStr

    INCLUDE 'namelist_read.inc'

    IF (hasPara) THEN
      ! Read (rdVarStr, *) var
      BLOCK
        CHARACTER*1024          :: sgVarStr
        INTEGER :: edPos, varIndex, i, numVars
        edPos = 1
        varIndex = 1
        numVars = 1
        IF (ALLOCATED(var)) DEALLOCATE (var)

        DO i = 1, LEN(TRIM(rdVarStr))
          IF (rdVarStr(i:i) .EQ. ',') numVars = numVars + 1
        END DO
        ALLOCATE (var(numVars))
        DO
          edPos = INDEX(rdVarStr, ",")
          IF (edPos .EQ. 0) EXIT
          sgVarStr = TRIM(rdVarStr(1:edPos - 1))
          READ (sgVarStr, *) var(varIndex)
          varIndex = varIndex + 1
          rdVarStr = rdVarStr(edPos + 1:LEN(TRIM(rdVarStr)))
        END DO

        READ (rdVarStr, *) var(varIndex)
      END BLOCK
    END IF
  END SUBROUTINE nml_read_string_array

END MODULE NMLRead_m
