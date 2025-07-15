! Created by Yali Wu, 2022/07/25
! Usage:
! Step 01: prepare the following two files by yourself:
! ###
! ### Prepare 3 data files + namelist.input before your run:
! ### the easiest way to do this is by setting both  and  as /public/home/simi
! ### /namelist.input (provide grid info)
! ### /grapesinput  (from SI)
! ### /grapesbdy (from SI)
! ### /grapesinput-MOTORDA
! ###
! ### The updated bdy file will be generatead as /grapesbdy_update
! ###
!
! Step 02: Submit the following script for conversion:
!
! #!/bin/bash
! #SBATCH -J DA
! #SBATCH -N 1
! #SBATCH --ntasks-per-node 1
! #SBATCH -p MOTOR
! #SBATCH -o log
! #SBATCH -e log

! export INPUT_DIR=/public/home/simi
! export OUTPUT_DIR=/public/home/simi
! /public/software/openmpi-4.0.2/bin/mpirun -n 1 -mca btl_tcp_if_include ib0 ./Test_updatebdy.exe

PROGRAM Test_updatebdy
  USE GrapesIO_m, ONLY: GrapesIO_t
  USE module_configure, ONLY: grid_config_rec_type
  USE module_variables
  USE module_utility
  USE kinds_m
  USE NMLRead_m
  USE YAMLRead_m
  USE UpdateBdy
  IMPLICIT NONE

  CHARACTER(LEN=1024) :: bdyPathName, anaPathName, nmlPathName

  ! Get the configFile
!   CALL GET_ENVIRONMENT_VARIABLE("INPUT_BDY", inputDir)
!   CALL GET_ENVIRONMENT_VARIABLE("OUTPUT_BDY", output_Dir)
!   CALL GET_ENVIRONMENT_VARIABLE("NEW_INPUT", DA_input_name)
  ! OutputFileName: grapesinput-MOTORDA

  CALL getarg(1, bdyPathName)
  CALL getarg(2, anaPathName)
  CALL getarg(3, nmlPathName)

  PRINT *, 'bdyPathName: ', TRIM(bdyPathName)
  PRINT *, 'anaPathName: ', TRIM(anaPathName)
  PRINT *, 'nmlPathName: ', TRIM(nmlPathName)

  IF (TRIM(bdyPathName) .EQ. '' .OR. TRIM(anaPathName) .EQ. '' &
      .OR. TRIM(nmlPathName) .EQ. '') THEN
    PRINT *, 'Wrong input arguements, test failed!'
    STOP
  END IF

!   IF (TRIM(configFile) .EQ. '') THEN
!     PRINT *, TRIM(configFile)

!     CALL getarg(0, configFile)
!     PRINT *, TRIM(configFile)
!     WRITE (*, 1)
! 1   FORMAT('Usage of this driver: mpirun -n <n> Debug/App_CMA_GD.exe configFile', /, &
!            ' Check your configure file and rerun!')
!     configFile = './App_3DVarVerification.yaml'
!   ELSE
!     WRITE (*, 2) TRIM(configFile)
! 2   FORMAT('ConfigFile is: ', A)
!   END IF

  CALL update_bdy(bdyPathName, anaPathName, nmlPathName)

  PRINT *, "Test passed!"

END PROGRAM Test_updatebdy
