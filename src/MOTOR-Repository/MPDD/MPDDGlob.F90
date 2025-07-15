!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA.mpdd.mpddInfo_sg
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2020/11/25, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module contains the data type for model_states
!! method.
!! @author Zilong Qin
!! @copyright (C) 2020 GBA-MWF, All rights reserved.
! @note
! @warning
! @attention
MODULE mpddGlob_m
  USE kinds_m, ONLY: i_kind, r_kind, r_Double, i_llong
  USE mpddBase_m, ONLY: mpddBase_t

!  USE mpi_f08
  INCLUDE "mpif.h"

!> @brief
!!
  TYPE, EXTENDS(MPDDBase_t) :: mpddGlob_t
    ! INTEGER(i_kind) :: nProc = 0                  !< Number of total processes running.
    ! INTEGER(i_kind) :: rank = -1                  !< Rank number of current process.
    ! INTEGER(i_kind) :: group = 0                  !< World group number of current process
    ! INTEGER(i_kind) :: comm = -9999               !< MPI_COMM_WORLD
    ! INTEGER(i_kind) :: rankBase = -1                     !< Rank of base proc
    ! INTEGER(i_kind) :: myrank = -1                       !< myrank = rank + 1
    ! INTEGER(i_kind) :: ierr = 0                          !< Error hanlder of MPI
    ! INTEGER(i_kind) :: STATUS(MPI_STATUS_SIZE)           !< Status info of MPI

  CONTAINS
    ! FINAL :: destructor               !< Deconstructor
    PROCEDURE, PUBLIC, PASS :: isNProcMatch          !< Reserved
    PROCEDURE, PUBLIC, PASS :: finalize
    PROCEDURE, PUBLIC :: initialize
  END TYPE mpddGlob_t

CONTAINS

  SUBROUTINE isNProcMatch(this, nProc)
    CLASS(mpddGlob_t) :: this
    INTEGER(i_kind), INTENT(IN) :: nProc

    IF (this%nProc .NE. nProc) THEN
      ERROR STOP 'Running proc number does not match with configuration file.'
    END IF
  END SUBROUTINE isNProcMatch

!> @brief
!! Initialize of mpdd
! @see
! @note
! @warning
! @attention
  SUBROUTINE initialize(this, rankBase)
    CLASS(mpddGlob_t) :: this

    INTEGER(i_kind), INTENT(IN), OPTIONAL :: rankBase              !< Whether prolongate the cells on the boundary

    IF (PRESENT(rankBase)) THEN
      this%rankBase = rankBase
    ELSE
      this%rankBase = 0
    END IF

    CALL mpi_init(this%ierr)
    CALL mpi_comm_rank(MPI_COMM_WORLD, this%rank, this%ierr)
    CALL mpi_comm_size(MPI_COMM_WORLD, this%nProc, this%ierr)
    CALL MPI_Comm_group(MPI_COMM_WORLD, this%group, this%ierr)

    this%comm = MPI_COMM_WORLD
    this%myrank = this%rank + 1
  END SUBROUTINE

  ! SUBROUTINE destructor(this)
  !   TYPE(mpddGlob_t) :: this
  !   ! IF (this%comm .ne. MPI_COMM_NULL) THEN
  !   !   PRINT *, "Proc", this%myrank, "is on finalize barrier."
  !   !   CALL this%barrier
  !   !   CALL mpi_finalize(this%ierr)
  !   ! END IF
  ! END SUBROUTINE destructor

  SUBROUTINE finalize(this)
    CLASS(mpddGlob_t) :: this
    IF (this%comm .NE. MPI_COMM_NULL) THEN
      PRINT *, "Proc", this%myrank, "is on finalize barrier."
      CALL this%barrier
      CALL mpi_finalize(this%ierr)
    END IF
  END SUBROUTINE finalize
END MODULE mpddGlob_m
