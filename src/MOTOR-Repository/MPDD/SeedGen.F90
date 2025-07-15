!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2021/5/7, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module is design for id for each object.
MODULE SeedGen_m
  INTEGER(kind=8) :: seed = 1000000
  INCLUDE "mpif.h"

CONTAINS

  INTEGER(kind=8) FUNCTION GenSeed()
    INTEGER(kind=8) :: ierr = 0                          !< Error hanlder of MPI
    seed = seed + 1
    GenSeed = seed
    CALL mpi_bcast(seed, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  END FUNCTION

END MODULE SeedGen_m
