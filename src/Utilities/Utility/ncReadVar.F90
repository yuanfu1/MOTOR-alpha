!----------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Centers
!                     for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Jiongming Pang
! VERSION           : V 0.0
! HISTORY           :
!   Created by Jiongming Pang (pang.j.m@hotmail.com), 2021-9-14, @GBA-MWF, Shenzhen
!----------------------------------------------------------------------------------------

MODULE ncReadVar_m
  USE kinds_m
  USE netcdf
  IMPLICIT NONE

  TYPE ncDimInfo_t
    CHARACTER(LEN=1024) :: inFileName
    CHARACTER(LEN=20), ALLOCATABLE :: file_dimName(:)
    INTEGER(i_kind), ALLOCATABLE :: file_dimLen(:)

  CONTAINS
    PROCEDURE, PUBLIC  :: ncFileDimsLen => getFileDimsLen_sub
  END TYPE ncDimInfo_t

  TYPE ncVarInfo_t
    CHARACTER(LEN=1024) :: inFileName
    CHARACTER(LEN=10) :: varName
    INTEGER(i_kind) :: num_dims = 0
    INTEGER(i_kind) :: var_types = 0
    INTEGER(i_kind), DIMENSION(:), ALLOCATABLE :: var_dimsLen

  CONTAINS
    ! PROCEDURE, PUBLIC  :: ncInInitialize => initialize_sub
    PROCEDURE, PUBLIC  :: ncVarDimsLen => getVarDimsLen_sub
  END TYPE ncVarInfo_t

  !! FOR REAL/INTEGER DATA TYPE
  TYPE, EXTENDS(ncVarInfo_t) :: ncData1D_t
    REAL(r_kind), DIMENSION(:), ALLOCATABLE :: ncVar

  CONTAINS
    PROCEDURE, PUBLIC  :: ncGetData1d => getdata1d_sub
  END TYPE ncData1D_t

  TYPE, EXTENDS(ncVarInfo_t) :: ncData2D_t
    REAL(r_kind), DIMENSION(:, :), ALLOCATABLE :: ncVar

  CONTAINS
    PROCEDURE, PUBLIC  :: ncGetData2d => getdata2d_sub
  END TYPE ncData2D_t

  TYPE, EXTENDS(ncVarInfo_t) :: ncData3D_t
    REAL(r_kind), DIMENSION(:, :, :), ALLOCATABLE :: ncVar

  CONTAINS
    PROCEDURE, PUBLIC  :: ncGetData3d => getdata3d_sub
  END TYPE ncData3D_t

  TYPE, EXTENDS(ncVarInfo_t) :: ncData4D_t
    REAL(r_kind), DIMENSION(:, :, :, :), ALLOCATABLE :: ncVar

  CONTAINS
    PROCEDURE, PUBLIC  :: ncGetData4d => getdata4d_sub
  END TYPE ncData4D_t

  TYPE, EXTENDS(ncVarInfo_t) :: ncData5D_t
    REAL(r_kind), DIMENSION(:, :, :, :, :), ALLOCATABLE :: ncVar

  CONTAINS
    PROCEDURE, PUBLIC  :: ncGetData5d => getdata5d_sub
  END TYPE ncData5D_t

  ! !! FOR CHARACTER DATA TYPE
  TYPE, EXTENDS(ncVarInfo_t) :: ncStr_t
    CHARACTER(len=20), DIMENSION(:), ALLOCATABLE :: ncVar

  CONTAINS
    PROCEDURE, PUBLIC  :: ncGetStr => getstr_sub
  END TYPE ncStr_t

  ! PRIVATE :: initialize_sub
  PRIVATE :: getVarDimsLen_sub, getFileDimsLen_sub
  PRIVATE :: getdata1d_sub, getdata2d_sub, getdata3d_sub
  PRIVATE :: getdata4d_sub, getdata5d_sub
  PRIVATE :: getstr_sub

CONTAINS

  SUBROUTINE getFileDimsLen_sub(this, inFileName, dimName)
    IMPLICIT NONE
    CLASS(ncDimInfo_t) :: this
    CHARACTER(LEN=*), INTENT(IN) :: inFileName
    CHARACTER(LEN=20), INTENT(IN) :: dimName(:)
    INTEGER(i_kind) :: ncFileID
    INTEGER(i_kind) :: ncDimID
    INTEGER(i_kind) :: num_dims
    ! INTEGER(i_kind), DIMENSION(:), ALLOCATABLE :: dimsLen
    INTEGER(i_kind) :: i

    num_dims = UBOUND(dimName, 1)
    ALLOCATE (this%file_dimName(num_dims), this%file_dimLen(num_dims))
    this%inFileName = TRIM(inFileName)
    this%file_dimName = dimName

    CALL check(nf90_open(TRIM(this%inFileName), nf90_nowrite, ncFileID))
    DO i = 1, num_dims
      CALL check(nf90_inq_dimid(ncFileID, TRIM(dimName(i)), ncDimID))
      CALL check(nf90_inquire_dimension(ncFileID, ncDimID, len=this%file_dimLen(i)))
    END DO
    CALL check(nf90_close(ncFileID))
  END SUBROUTINE getFileDimsLen_sub

  SUBROUTINE getVarDimsLen_sub(this, inFileName, varName)
    IMPLICIT NONE
    CLASS(ncVarInfo_t) :: this
    CHARACTER(LEN=*), INTENT(IN) :: inFileName
    CHARACTER(LEN=*), INTENT(IN) :: varName
    INTEGER(i_kind) :: ncFileID
    INTEGER(i_kind) :: ncVarID
    INTEGER(i_kind) :: num_dims
    INTEGER(i_kind), DIMENSION(:), ALLOCATABLE :: var_dimsID
    INTEGER(i_kind), DIMENSION(:), ALLOCATABLE :: var_dimsLen
    INTEGER(i_kind) :: i

    this%inFileName = TRIM(inFileName)
    this%varName = TRIM(varName)

    CALL check(nf90_open(TRIM(this%inFileName), nf90_nowrite, ncFileID))
    CALL check(nf90_inq_varid(ncFileID, TRIM(this%varName), ncVarID))
    CALL check(nf90_inquire_variable(ncFileID, ncVarID, xtype=this%var_types, ndims=this%num_dims))

    ALLOCATE (var_dimsID(this%num_dims), this%var_dimsLen(this%num_dims))
    CALL check(nf90_inquire_variable(ncFileID, ncVarID, dimids=var_dimsID))
    DO i = 1, this%num_dims
      CALL check(nf90_inquire_dimension(ncFileID, var_dimsID(i), len=this%var_dimsLen(i)))
    END DO
    CALL check(nf90_close(ncFileID))
  END SUBROUTINE getVarDimsLen_sub

  SUBROUTINE getdata1d_sub(this, inFileName, varName)
    IMPLICIT NONE
    CLASS(ncData1D_t) :: this
    INTEGER(i_kind) :: ncFileID
    INTEGER(i_kind) :: ncVarID
    CHARACTER(len=*), INTENT(IN) :: inFileName
    CHARACTER(len=*), INTENT(IN) :: varName

    CALL this%ncVarDimsLen(TRIM(inFileName), TRIM(varName))

    IF (this%num_dims /= 1) THEN
      STOP "Wrong dimensions of input data: ONE-D of data"
    END IF

    ALLOCATE (this%ncVar(this%var_dimsLen(1)))
    CALL check(nf90_open(TRIM(this%inFileName), nf90_nowrite, ncFileID))
    CALL check(nf90_inq_varid(ncFileID, TRIM(this%varName), ncVarID))

    IF (this%var_types == NF90_FLOAT) THEN
      BLOCK
        REAL(r_single), DIMENSION(:), ALLOCATABLE :: tmp
        ALLOCATE (tmp(this%var_dimsLen(1)))
        CALL check(nf90_get_var(ncFileID, ncVarID, tmp))
        this%ncVar = tmp
        DEALLOCATE (tmp)
      END BLOCK
    ELSE IF (this%var_types == NF90_DOUBLE) THEN
      BLOCK
        REAL(r_double), DIMENSION(:), ALLOCATABLE :: tmp
        ALLOCATE (tmp(this%var_dimsLen(1)))
        CALL check(nf90_get_var(ncFileID, ncVarID, tmp))
        this%ncVar = tmp
        DEALLOCATE (tmp)
      END BLOCK
    ELSE IF (this%var_types == NF90_SHORT) THEN
      BLOCK
        INTEGER(i_short), DIMENSION(:), ALLOCATABLE :: tmp
        ALLOCATE (tmp(this%var_dimsLen(1)))
        CALL check(nf90_get_var(ncFileID, ncVarID, tmp))
        this%ncVar = tmp
        DEALLOCATE (tmp)
      END BLOCK
    ELSE IF (this%var_types == NF90_INT) THEN
      BLOCK
        INTEGER(i_long), DIMENSION(:), ALLOCATABLE :: tmp
        ALLOCATE (tmp(this%var_dimsLen(1)))
        CALL check(nf90_get_var(ncFileID, ncVarID, tmp))
        this%ncVar = tmp
        DEALLOCATE (tmp)
      END BLOCK
    ELSE
      STOP "Wrong data type of input data: 1d-(REAL/INTEGER)"
    END IF
    CALL check(nf90_close(ncFileID))
  END SUBROUTINE getdata1d_sub

  SUBROUTINE getdata2d_sub(this, inFileName, varName)
    IMPLICIT NONE
    CLASS(ncData2D_t) :: this
    INTEGER(i_kind) :: ncFileID
    INTEGER(i_kind) :: ncVarID
    CHARACTER(len=*), INTENT(IN) :: inFileName
    CHARACTER(len=*), INTENT(IN) :: varName

    CALL this%ncVarDimsLen(TRIM(inFileName), TRIM(varName))

    IF (this%num_dims /= 2) THEN
      STOP "Wrong dimensions of input data: TWO-D of data"
    END IF

    ALLOCATE (this%ncVar(this%var_dimsLen(1), this%var_dimsLen(2)))
    CALL check(nf90_open(TRIM(this%inFileName), nf90_nowrite, ncFileID))
    CALL check(nf90_inq_varid(ncFileID, TRIM(this%varName), ncVarID))

    IF (this%var_types == NF90_FLOAT) THEN
      BLOCK
        REAL(r_single), DIMENSION(:, :), ALLOCATABLE :: tmp
        ALLOCATE (tmp(this%var_dimsLen(1), this%var_dimsLen(2)))
        CALL check(nf90_get_var(ncFileID, ncVarID, tmp))
        this%ncVar = tmp
        DEALLOCATE (tmp)
      END BLOCK
    ELSE IF (this%var_types == NF90_DOUBLE) THEN
      BLOCK
        REAL(r_double), DIMENSION(:, :), ALLOCATABLE :: tmp
        ALLOCATE (tmp(this%var_dimsLen(1), this%var_dimsLen(2)))
        CALL check(nf90_get_var(ncFileID, ncVarID, tmp))
        this%ncVar = tmp
        DEALLOCATE (tmp)
      END BLOCK
    ELSE IF (this%var_types == NF90_SHORT) THEN
      BLOCK
        INTEGER(i_short), DIMENSION(:, :), ALLOCATABLE :: tmp
        ALLOCATE (tmp(this%var_dimsLen(1), this%var_dimsLen(2)))
        CALL check(nf90_get_var(ncFileID, ncVarID, tmp))
        this%ncVar = tmp
        DEALLOCATE (tmp)
      END BLOCK
    ELSE IF (this%var_types == NF90_INT) THEN
      BLOCK
        INTEGER(i_long), DIMENSION(:, :), ALLOCATABLE :: tmp
        ALLOCATE (tmp(this%var_dimsLen(1), this%var_dimsLen(2)))
        CALL check(nf90_get_var(ncFileID, ncVarID, tmp))
        this%ncVar = tmp
        DEALLOCATE (tmp)
      END BLOCK
    ELSE
      STOP "Wrong data type of input data: 2d-(REAL/INTEGER)"
    END IF
    CALL check(nf90_close(ncFileID))
  END SUBROUTINE getdata2d_sub

  SUBROUTINE getdata3d_sub(this, inFileName, varName)
    IMPLICIT NONE
    CLASS(ncData3D_t) :: this
    INTEGER(i_kind) :: ncFileID
    INTEGER(i_kind) :: ncVarID
    CHARACTER(len=*), INTENT(IN) :: inFileName
    CHARACTER(len=*), INTENT(IN) :: varName

    CALL this%ncVarDimsLen(TRIM(inFileName), TRIM(varName))

    IF (this%num_dims /= 3) THEN
      STOP "Wrong dimensions of input data: THREE-D of data"
    END IF

    ALLOCATE (this%ncVar(this%var_dimsLen(1), this%var_dimsLen(2), this%var_dimsLen(3)))
    CALL check(nf90_open(TRIM(this%inFileName), nf90_nowrite, ncFileID))
    CALL check(nf90_inq_varid(ncFileID, TRIM(this%varName), ncVarID))

    IF (this%var_types == NF90_FLOAT) THEN
      BLOCK
        REAL(r_single), DIMENSION(:, :, :), ALLOCATABLE :: tmp
        ALLOCATE (tmp(this%var_dimsLen(1), this%var_dimsLen(2), this%var_dimsLen(3)))
        CALL check(nf90_get_var(ncFileID, ncVarID, tmp))
        this%ncVar = tmp
        DEALLOCATE (tmp)
      END BLOCK
    ELSE IF (this%var_types == NF90_DOUBLE) THEN
      BLOCK
        REAL(r_double), DIMENSION(:, :, :), ALLOCATABLE :: tmp
        ALLOCATE (tmp(this%var_dimsLen(1), this%var_dimsLen(2), this%var_dimsLen(3)))
        CALL check(nf90_get_var(ncFileID, ncVarID, tmp))
        this%ncVar = tmp
        DEALLOCATE (tmp)
      END BLOCK
    ELSE IF (this%var_types == NF90_SHORT) THEN
      BLOCK
        INTEGER(i_short), DIMENSION(:, :, :), ALLOCATABLE :: tmp
        ALLOCATE (tmp(this%var_dimsLen(1), this%var_dimsLen(2), this%var_dimsLen(3)))
        CALL check(nf90_get_var(ncFileID, ncVarID, tmp))
        this%ncVar = tmp
        DEALLOCATE (tmp)
      END BLOCK
    ELSE IF (this%var_types == NF90_INT) THEN
      BLOCK
        INTEGER(i_long), DIMENSION(:, :, :), ALLOCATABLE :: tmp
        ALLOCATE (tmp(this%var_dimsLen(1), this%var_dimsLen(2), this%var_dimsLen(3)))
        CALL check(nf90_get_var(ncFileID, ncVarID, tmp))
        this%ncVar = tmp
        DEALLOCATE (tmp)
      END BLOCK
    ELSE
      STOP "Wrong data type of input data: 3d-(REAL/INTEGER)"
    END IF
    CALL check(nf90_close(ncFileID))
  END SUBROUTINE getdata3d_sub

  SUBROUTINE getdata4d_sub(this, inFileName, varName)
    IMPLICIT NONE
    CLASS(ncData4D_t) :: this
    INTEGER(i_kind) :: ncFileID
    INTEGER(i_kind) :: ncVarID
    CHARACTER(len=*), INTENT(IN) :: inFileName
    CHARACTER(len=*), INTENT(IN) :: varName

    CALL this%ncVarDimsLen(TRIM(inFileName), TRIM(varName))

    IF (this%num_dims /= 4) THEN
      STOP "Wrong dimensions of input data: FOUR-D of data"
    END IF

    ALLOCATE (this%ncVar(this%var_dimsLen(1), this%var_dimsLen(2), this%var_dimsLen(3), this%var_dimsLen(4)))
    CALL check(nf90_open(TRIM(this%inFileName), nf90_nowrite, ncFileID))
    CALL check(nf90_inq_varid(ncFileID, TRIM(this%varName), ncVarID))

    IF (this%var_types == NF90_FLOAT) THEN
      BLOCK
        REAL(r_single), DIMENSION(:, :, :, :), ALLOCATABLE :: tmp
        ALLOCATE (tmp(this%var_dimsLen(1), this%var_dimsLen(2), this%var_dimsLen(3), this%var_dimsLen(4)))
        CALL check(nf90_get_var(ncFileID, ncVarID, tmp))
        this%ncVar = tmp
        DEALLOCATE (tmp)
      END BLOCK
    ELSE IF (this%var_types == NF90_DOUBLE) THEN
      BLOCK
        REAL(r_double), DIMENSION(:, :, :, :), ALLOCATABLE :: tmp
        ALLOCATE (tmp(this%var_dimsLen(1), this%var_dimsLen(2), this%var_dimsLen(3), this%var_dimsLen(4)))
        CALL check(nf90_get_var(ncFileID, ncVarID, tmp))
        this%ncVar = tmp
        DEALLOCATE (tmp)
      END BLOCK
    ELSE IF (this%var_types == NF90_SHORT) THEN
      BLOCK
        INTEGER(i_short), DIMENSION(:, :, :, :), ALLOCATABLE :: tmp
        ALLOCATE (tmp(this%var_dimsLen(1), this%var_dimsLen(2), this%var_dimsLen(3), this%var_dimsLen(4)))
        CALL check(nf90_get_var(ncFileID, ncVarID, tmp))
        this%ncVar = tmp
        DEALLOCATE (tmp)
      END BLOCK
    ELSE IF (this%var_types == NF90_INT) THEN
      BLOCK
        INTEGER(i_long), DIMENSION(:, :, :, :), ALLOCATABLE :: tmp
        ALLOCATE (tmp(this%var_dimsLen(1), this%var_dimsLen(2), this%var_dimsLen(3), this%var_dimsLen(4)))
        CALL check(nf90_get_var(ncFileID, ncVarID, tmp))
        this%ncVar = tmp
        DEALLOCATE (tmp)
      END BLOCK
    ELSE
      STOP "Wrong data type of input data: 4d-(REAL/INTEGER)"
    END IF
    CALL check(nf90_close(ncFileID))
  END SUBROUTINE getdata4d_sub

  SUBROUTINE getdata5d_sub(this, inFileName, varName)
    IMPLICIT NONE
    CLASS(ncData5D_t) :: this
    INTEGER(i_kind) :: ncFileID
    INTEGER(i_kind) :: ncVarID
    CHARACTER(len=*), INTENT(IN) :: inFileName
    CHARACTER(len=*), INTENT(IN) :: varName

    CALL this%ncVarDimsLen(TRIM(inFileName), TRIM(varName))

    IF (this%num_dims /= 5) THEN
      STOP "Wrong dimensions of input data: FIVE-D of data"
    END IF

    ALLOCATE (this%ncVar(this%var_dimsLen(1), this%var_dimsLen(2), this%var_dimsLen(3), this%var_dimsLen(4), this%var_dimsLen(5)))
    CALL check(nf90_open(TRIM(this%inFileName), nf90_nowrite, ncFileID))
    CALL check(nf90_inq_varid(ncFileID, TRIM(this%varName), ncVarID))

    IF (this%var_types == NF90_FLOAT) THEN
      BLOCK
        REAL(r_single), DIMENSION(:, :, :, :, :), ALLOCATABLE :: tmp
        ALLOCATE (tmp(this%var_dimsLen(1), this%var_dimsLen(2), this%var_dimsLen(3), this%var_dimsLen(4), this%var_dimsLen(5)))
        CALL check(nf90_get_var(ncFileID, ncVarID, tmp))
        this%ncVar = tmp
        DEALLOCATE (tmp)
      END BLOCK
    ELSE IF (this%var_types == NF90_DOUBLE) THEN
      BLOCK
        REAL(r_double), DIMENSION(:, :, :, :, :), ALLOCATABLE :: tmp
        ALLOCATE (tmp(this%var_dimsLen(1), this%var_dimsLen(2), this%var_dimsLen(3), this%var_dimsLen(4), this%var_dimsLen(5)))
        CALL check(nf90_get_var(ncFileID, ncVarID, tmp))
        this%ncVar = tmp
        DEALLOCATE (tmp)
      END BLOCK
    ELSE IF (this%var_types == NF90_SHORT) THEN
      BLOCK
        INTEGER(i_short), DIMENSION(:, :, :, :, :), ALLOCATABLE :: tmp
        ALLOCATE (tmp(this%var_dimsLen(1), this%var_dimsLen(2), this%var_dimsLen(3), this%var_dimsLen(4), this%var_dimsLen(5)))
        CALL check(nf90_get_var(ncFileID, ncVarID, tmp))
        this%ncVar = tmp
        DEALLOCATE (tmp)
      END BLOCK
    ELSE IF (this%var_types == NF90_INT) THEN
      BLOCK
        INTEGER(i_long), DIMENSION(:, :, :, :, :), ALLOCATABLE :: tmp
        ALLOCATE (tmp(this%var_dimsLen(1), this%var_dimsLen(2), this%var_dimsLen(3), this%var_dimsLen(4), this%var_dimsLen(5)))
        CALL check(nf90_get_var(ncFileID, ncVarID, tmp))
        this%ncVar = tmp
        DEALLOCATE (tmp)
      END BLOCK
    ELSE
      STOP "Wrong data type of input data: 5d-(REAL/INTEGER)"
    END IF
    CALL check(nf90_close(ncFileID))
  END SUBROUTINE getdata5d_sub

  !! FOR CHARACTER DATA TYPE
  SUBROUTINE getstr_sub(this, inFileName, varName)
    IMPLICIT NONE
    CLASS(ncStr_t) :: this
    INTEGER(i_kind) :: ncFileID
    INTEGER(i_kind) :: ncVarID
    CHARACTER(len=*), INTENT(IN) :: inFileName
    CHARACTER(len=*), INTENT(IN) :: varName
    INTEGER(i_kind) :: istart(2), iend(2)
    INTEGER(i_kind) :: i, istate

    CALL this%ncVarDimsLen(TRIM(inFileName), TRIM(varName))

    IF (this%num_dims /= 2) THEN
      STOP "Missmatch dimensions of the character variable in input data"
    END IF

    ALLOCATE (this%ncVar(this%var_dimsLen(2)))
    CALL check(nf90_open(TRIM(this%inFileName), nf90_nowrite, ncFileID))
    CALL check(nf90_inq_varid(ncFileID, TRIM(this%varName), ncVarID))

    IF (this%var_types == NF90_CHAR) THEN
      ! nchars=this%var_dimsLen(1)
      DO i = 1, this%var_dimsLen(2)
        istart(1) = 1
        iend(1) = this%var_dimsLen(1)
        istart(2) = i
        iend(2) = 1
        ! istate = nf_get_vara_text(ncFileID, ncVarID, istart, iend, this%ncVar(i))
        CALL check(nf90_get_var(ncFileID, ncVarID, this%ncVar(i), start=istart, count=iend))
      END DO
      ! CALL check( nf90_get_var(ncFileID, ncVarID, this%ncVar) )
    ELSE
      STOP "Wrong data type of input data: CHARACTER"
    END IF
    CALL check(nf90_close(ncFileID))
  END SUBROUTINE getstr_sub

  SUBROUTINE check(status)
    INTEGER, INTENT(IN) :: status

    IF (status /= nf90_noerr) THEN
      PRINT *, TRIM(nf90_strerror(status))
      STOP "Stopped"
    END IF
  END SUBROUTINE check

END MODULE ncReadVar_m
