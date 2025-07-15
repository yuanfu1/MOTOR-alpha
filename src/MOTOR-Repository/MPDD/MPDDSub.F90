!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA.mpdd.mpddSub
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2021/1/26, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!!
!! method.
!! @copyright (C) 2020 GBA-MWF, All rights reserved.
! @note
! @warning
! @attention
MODULE mpddSub_m
  USE kinds_m, ONLY: i_kind, r_kind, r_Double, i_llong
  USE mpddBase_m, ONLY: mpddBase_t
!  USE mpi_f08
  INCLUDE "mpif.h"

  TYPE, EXTENDS(MPDDBase_t) :: mpddSub_t
    ! INTEGER(i_kind) :: group = 0
    ! INTEGER(i_kind) :: ierr = 0
    ! INTEGER(i_kind) :: comm = 0
    ! INTEGER(i_kind) :: rankBase = -1                  !< Rank of base proc
    ! INTEGER(i_kind) :: STATUS(MPI_STATUS_SIZE)

    ! INTEGER(i_kind) :: nProc = 0
    ! INTEGER(i_kind) :: rank = -1
    ! INTEGER(i_kind) :: myrank = -1
  CONTAINS
    ! PROCEDURE, PUBLIC, PASS :: isActiveProc
    ! PROCEDURE, PUBLIC, PASS :: isBaseProc
    PROCEDURE, PUBLIC :: initialize
    FINAL :: destructor               !< Deconstructor
  END TYPE mpddSub_t
  ! INTERFACE mpddSub_t
  !   PROCEDURE :: constructor
  ! END interface mpddSub_t

CONTAINS

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    TYPE(mpddSub_t), INTENT(INOUT) :: this
    IF (this%isBaseProc()) PRINT *, 'mpddSub_t destructor is invoked!'
    IF (this%isActiveProc()) THEN
      IF (this%comm .NE. MPI_COMM_NULL) THEN
        CALL MPI_Group_free(this%group, this%ierr)
        CALL MPI_Comm_free(this%comm, this%ierr)
      END IF
    END IF
  END SUBROUTINE destructor

  SUBROUTINE initialize(this, nProc, group, rank, rankBase)
    CLASS(mpddSub_t) :: this
    INTEGER(i_kind), INTENT(IN) :: nProc, rank
    INTEGER(i_kind), INTENT(IN) :: group
    INTEGER(i_kind), INTENT(IN), OPTIONAL :: rankBase              !< Whether prolongate the cells on the boundary
    INTEGER(i_kind), ALLOCATABLE :: ranksSub(:)
    INTEGER(i_kind) :: i

    IF (PRESENT(rankBase)) THEN
      this%rankBase = rankBase
    ELSE
      this%rankBase = 0
    END IF

    this%nProc = nProc

    ALLOCATE (ranksSub(this%nProc))

    DO i = 1, this%nProc
      ranksSub(i) = i - 1
    END DO

    IF (rank >= nProc) this%nProc = 0

    ! Construct a group containing all of the prime ranks in global_group
    CALL MPI_Group_incl(group, this%nProc, ranksSub, this%group, this%ierr)

    ! Create a new communicator based on the group
    CALL MPI_Comm_create(MPI_COMM_WORLD, this%group, this%comm, this%ierr)

    IF (this%isActiveProc()) CALL MPI_Comm_rank(this%comm, this%rank, this%ierr)
    IF (this%isActiveProc()) CALL MPI_Comm_size(this%comm, this%nProc, this%ierr)
    IF (this%isActiveProc()) this%myrank = this%rank + 1

    DEALLOCATE (ranksSub)
  END SUBROUTINE

END MODULE mpddSub_m
