!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA.MPObs
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2021/1/26, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!!
!! @copyright (C) 2020 GBA-MWF, All rights reserved.
MODULE MPObs_m
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE mpddSub_m, ONLY: mpddSub_t

  TYPE, EXTENDS(mpddSub_t) :: MPObs_t
    TYPE(SingleGrid_t), POINTER :: sg

  CONTAINS
    FINAL :: destructor
    PROCEDURE, PUBLIC :: aggr_to_base_proc
    PROCEDURE, PUBLIC :: initializeMPObs
  END TYPE MPObs_t

CONTAINS

  SUBROUTINE initializeMPObs(this, sg)
    IMPLICIT NONE
    CLASS(MPObs_t) :: this
    TYPE(SingleGrid_t), TARGET, INTENT(IN) :: sg

    this%sg => sg

    CALL this%initialize(sg%mpddInfo_sg%nProc, sg%mpddInfo_sg%group, &
                         sg%mpddInfo_sg%rank, sg%mpddInfo_sg%rankBase)

    ! this%nProc = sg%mpddInfo_sg%nProc                  !< Number of total processes running in the group.
    ! this%rank = sg%mpddInfo_sg%rank                    !< Rank number of current process.
    ! this%group = sg%mpddInfo_sg%group                  !< Group number of current process
    ! this%comm = sg%mpddInfo_sg%comm                    !< MPI_COMM_WORLD
    ! this%rankBase = sg%mpddInfo_sg%rankBase            !< Rank of base proc
    ! this%myrank = sg%mpddInfo_sg%myrank                !< myrank = rank + 1  END FUNCTION cons_replication
  END SUBROUTINE

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(MPObs_t), INTENT(INOUT) :: this
  END SUBROUTINE destructor

  SUBROUTINE aggr_to_base_proc(this, a)
    CLASS(MPObs_t) :: this
    TYPE(*) :: a(:)

  END SUBROUTINE aggr_to_base_proc

END MODULE MPObs_m
