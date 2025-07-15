!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2021/1/26, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!
MODULE mpddBase_m
  USE kinds_m, ONLY: i_kind, r_kind, r_Double, i_llong
  INCLUDE "mpif.h"

  TYPE mpddBase_t
    INTEGER(i_kind) :: nProc = 0                  !< Number of total processes running in the group.
    INTEGER(i_kind) :: rank = -1                  !< Rank number of current process.
    INTEGER(i_kind) :: group = 0                  !< Group number of current process
    INTEGER(i_kind) :: comm = MPI_COMM_NULL              !< MPI_COMM_WORLD
    INTEGER(i_kind) :: rankBase = -1                     !< Rank of base proc
    INTEGER(i_kind) :: myrank = -1                       !< myrank = rank + 1
    INTEGER(i_kind) :: ierr = 0                          !< Error hanlder of MPI
    INTEGER(i_kind) :: STATUS(MPI_STATUS_SIZE)           !< Status info of MPI
  CONTAINS
    PROCEDURE, PUBLIC, PASS :: barrier               !< MPI barrier in global
    PROCEDURE, PUBLIC, PASS :: isActiveProc
    PROCEDURE, PUBLIC, PASS :: isBaseProc

    PROCEDURE, PUBLIC, PASS :: AllReduceSumReal
    PROCEDURE, PUBLIC, PASS :: AllReduceSumInt

    PROCEDURE, PUBLIC, PASS :: AllReduceMaxReal
    PROCEDURE, PUBLIC, PASS :: AllReduceMinReal
    PROCEDURE, PUBLIC, PASS :: AllRedMaxReal
    PROCEDURE, PUBLIC, PASS :: AllRedMinReal
    PROCEDURE, PUBLIC, PASS :: AllRedLOR

    PROCEDURE, PUBLIC, PASS :: bCastInt
    PROCEDURE, PUBLIC, PASS :: bCastReal
    PROCEDURE, PUBLIC, PASS :: bCastInt1D
    PROCEDURE, PUBLIC, PASS :: bCastString
    PROCEDURE, PUBLIC, PASS :: bCastInt81D
    PROCEDURE, PUBLIC, PASS :: bCastReal1D
    PROCEDURE, PUBLIC, PASS :: bCastReal2D

    GENERIC, PUBLIC :: bCast &                        !< Boardcast the varibales in single int, single real and 1D int
      => bCastInt, bCastReal, bCastInt1D, bCastReal1D, bCastReal2D

  END TYPE mpddBase_t

  INTERFACE mpddBase_t
    PROCEDURE :: constructor
  END INTERFACE mpddBase_t

CONTAINS

  FUNCTION constructor() RESULT(this)
    TYPE(mpddBase_t) :: this
  END FUNCTION constructor

  SUBROUTINE barrier(this)
    IMPLICIT NONE
    CLASS(mpddBase_t) :: this
    CALL mpi_barrier(this%comm, this%ierr)
  END SUBROUTINE barrier

  LOGICAL FUNCTION isBaseProc(this)
    CLASS(mpddBase_t) :: this
    IF (this%rank .EQ. this%rankBase) THEN
      isBaseProc = .TRUE.
    ELSE
      isBaseProc = .FALSE.
    END IF
  END FUNCTION isBaseProc

  LOGICAL FUNCTION isActiveProc(this)
    CLASS(mpddBase_t) :: this
    IF (this%comm .NE. MPI_COMM_NULL) THEN
      isActiveProc = .TRUE.
    ELSE
      isActiveProc = .FALSE.
    END IF
  END FUNCTION isActiveProc

  SUBROUTINE AllReduceSumReal(this, send_data, recv_data)
    CLASS(mpddBase_t) :: this
    REAL(r_kind), INTENT(IN) :: send_data
    REAL(r_kind), INTENT(INOUT) :: recv_data

    CALL MPI_ALLREDUCE(send_data, recv_data, 1, MPI_DOUBLE_PRECISION, MPI_SUM, this%comm, this%ierr)
  END SUBROUTINE AllReduceSumReal

  SUBROUTINE AllReduceSumInt(this, send_data, recv_data)
    CLASS(mpddBase_t) :: this
    INTEGER(i_kind), INTENT(IN) :: send_data
    INTEGER(i_kind), INTENT(INOUT) :: recv_data

    CALL MPI_ALLREDUCE(send_data, recv_data, 1, MPI_INTEGER, MPI_SUM, this%comm, this%ierr)
  END SUBROUTINE AllReduceSumInt

  SUBROUTINE AllReduceMaxReal(this, send_data, recv_data)
    CLASS(mpddBase_t) :: this
    REAL(r_kind), INTENT(IN) :: send_data
    REAL(r_kind), INTENT(INOUT) :: recv_data

    CALL MPI_ALLREDUCE(send_data, recv_data, 1, MPI_DOUBLE_PRECISION, MPI_MAX, this%comm, this%ierr)
  END SUBROUTINE

  FUNCTION AllRedMaxReal(this, send_data) RESULT(recv_data)
    CLASS(mpddBase_t) :: this
    REAL(r_kind), INTENT(IN) :: send_data
    REAL(r_kind) :: recv_data

    CALL MPI_ALLREDUCE(send_data, recv_data, 1, MPI_DOUBLE_PRECISION, MPI_MAX, this%comm, this%ierr)
  END FUNCTION

  FUNCTION AllRedMinReal(this, send_data) RESULT(recv_data)
    CLASS(mpddBase_t) :: this
    REAL(r_kind), INTENT(IN) :: send_data
    REAL(r_kind) :: recv_data

    CALL MPI_ALLREDUCE(send_data, recv_data, 1, MPI_DOUBLE_PRECISION, MPI_MIN, this%comm, this%ierr)
  END FUNCTION

  SUBROUTINE AllReduceMinReal(this, send_data, recv_data)
    CLASS(mpddBase_t) :: this
    REAL(r_kind), INTENT(IN) :: send_data
    REAL(r_kind), INTENT(INOUT) :: recv_data

    CALL MPI_ALLREDUCE(send_data, recv_data, 1, MPI_DOUBLE_PRECISION, MPI_MIN, this%comm, this%ierr)
  END SUBROUTINE

  SUBROUTINE AllRedLOR(this, send_data, recv_data)
    CLASS(mpddBase_t) :: this
    LOGICAL, INTENT(IN) :: send_data
    LOGICAL, INTENT(INOUT) :: recv_data

    CALL MPI_ALLREDUCE(send_data, recv_data, 1, MPI_LOGICAL, MPI_LOR, this%comm, this%ierr)
  END SUBROUTINE

  SUBROUTINE bCastInt(this, addr)
    CLASS(mpddBase_t) :: this
    INTEGER(i_kind), INTENT(INOUT) :: addr

    CALL mpi_bcast(addr, 1, MPI_INTEGER, this%rankBase, this%comm, this%ierr)
  END SUBROUTINE bCastInt

  SUBROUTINE bCastInt1D(this, addr)
    CLASS(mpddBase_t) :: this
    INTEGER(i_kind), INTENT(INOUT) :: addr(:)
    INTEGER(i_kind) :: sizeAddr

    sizeAddr = SIZE(addr)
    CALL mpi_bcast(addr, sizeAddr, MPI_INTEGER, this%rankBase, this%comm, this%ierr)
  END SUBROUTINE bCastInt1D

  SUBROUTINE bCastInt81D(this, addr)
    CLASS(mpddBase_t) :: this
    INTEGER*8, INTENT(INOUT) :: addr(:)
    INTEGER(i_kind) :: sizeAddr

    sizeAddr = SIZE(addr) * sizeof(INTEGER * 8)
    CALL mpi_bcast(addr, sizeAddr, MPI_CHAR, this%rankBase, this%comm, this%ierr)
  END SUBROUTINE bCastInt81D

  SUBROUTINE bCastString(this, addr)
    CLASS(mpddBase_t) :: this
    CHARACTER(len=*), INTENT(INOUT) :: addr
    INTEGER(i_kind) :: sizeAddr

    sizeAddr = LEN(addr)
    CALL mpi_bcast(addr, sizeAddr, MPI_CHAR, this%rankBase, this%comm, this%ierr)
  END SUBROUTINE bCastString

  SUBROUTINE bCastReal(this, addr)
    CLASS(mpddBase_t) :: this
    REAL(r_kind), INTENT(INOUT) :: addr

    CALL mpi_bcast(addr, 1, MPI_DOUBLE_PRECISION, this%rankBase, this%comm, this%ierr)
  END SUBROUTINE bCastReal

  SUBROUTINE bCastReal1D(this, addr)
    CLASS(mpddBase_t) :: this
    REAL(r_kind), INTENT(INOUT) :: addr(:)

    CALL mpi_bcast(addr, SIZE(addr), MPI_DOUBLE_PRECISION, this%rankBase, this%comm, this%ierr)
  END SUBROUTINE bCastReal1D

  SUBROUTINE bCastReal2D(this, addr)
    CLASS(mpddBase_t) :: this
    REAL(r_kind), INTENT(INOUT) :: addr(:, :)

    CALL mpi_bcast(addr, SIZE(addr), MPI_DOUBLE_PRECISION, this%rankBase, this%comm, this%ierr)
  END SUBROUTINE bCastReal2D

END MODULE mpddBase_m
