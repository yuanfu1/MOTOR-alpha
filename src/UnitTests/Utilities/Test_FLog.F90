!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA.flog
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Chris MacMackin <cmacmackin@gmail.com> (https://github.com/cmacmackin/flogging)
! VERSION           : V 0.1
! HISTORY           :
!   Modified by Zilong Qin (zilong.qin@gmail.com), 2021/3/23, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
! @note
! @warning
! @attention
PROGRAM Test_FLog
  USE FLog_m, ONLY: logger, Enable_Log

  INCLUDE "mpif.h"
  INTEGER :: ierr                          !< Error hanlder of MPI
  INTEGER :: rank
  CHARACTER(LEN=1024) :: configFile, logFile

  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", configFile)
  configFile = TRIM(configFile)//"/UnitTest/Test_FLog.yaml"

  CALL mpi_init(ierr)

  ! Initialise the logger prior to use
  CALL Enable_Log(configFile, .TRUE.)

  ! Write some debugging information
  CALL logger%debug('logger_example', 'Starting program logger_example')

  ! Perform some calculation
  ! ...
  CALL logger%info('logger_example', 'Found result of calculation')
  ! Perform another calculation

  ! ...
  ! Oh no, an error has occurred
  CALL logger%error('logger_example', 'Calculation failed due to error')

  CALL logger%warning('logger_example', 'Ending program logger_example')

  CALL mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)

  IF (rank .EQ. 0) PRINT *, 'Test passed!'
  CALL logger%destroy

  CALL mpi_finalize(ierr)

END PROGRAM Test_FLog
