!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA.IO.GrapesIO
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
! AUTOHR(S)         : Sanshan Tu, Zilong Qin
! VERSION           : V 0.0
! HISTORY           :
!   Created by Sanshan Tu (tss71618@163.com), 2020/12/31, @SZSC, Shenzhen
!   Modified by Sanshan Tu (tss71618@163.com), 2021/11/01, @SZSC, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!! GrapesIO module contains the data type for model_states
!! method.
!! @copyright (C) 2020 GBA-MWF, All rights reserved.
! @note
! @warning
! @attention
program GRAPES_INPUT_TEST
  USE ncWriteGrid_m
  USE GrapesIO_m, ONLY: GrapesIO_t
  USE kinds_m, ONLY: i_kind, r_kind, r_single
  ! USE Filter_m, ONLY: smoothField, guidedfilter

  IMPLICIT NONE
  REAL(r_kind)      :: bgn_time, end_time, value_tmp
  REAL(r_kind), DIMENSION(:, :, :), ALLOCATABLE :: value_sin, value_u, value_cmp
  INTEGER           :: ierr, i, j, k
  type(GrapesIO_t)  :: GrapesIO
  INTEGER*8         :: offset_1d_i, offset_1d_j, offset_1d_k
  CHARACTER(LEN=1024) :: ncOutDir, nlFileName, giFileName, qcqrFileName

  CALL GET_ENVIRONMENT_VARIABLE("INPUT_DIR", nlFileName)
  CALL GET_ENVIRONMENT_VARIABLE("OUTPUT_DIR", ncOutDir)
  CALL GET_ENVIRONMENT_VARIABLE("INPUT_DIR", giFileName)

  ! yaml_get_var(configFile, 'IO', 'input_dir', nlFileName)
  ! yaml_get_var(configFile, 'IO', 'output_dir', ncOutDir)
  ! yaml_get_var(configFile, 'IO', 'input_dir', giFileName)

  giFileName = "/Users/yaliwu/Desktop/2022052700case/model-cma-gd/grapesinput2022052700"
  qcqrFileName = "/Users/yaliwu/Desktop/2022052700case/model-cma-gd/qcqr2022052700.dat"
  nlFileName = "/Users/yaliwu/Desktop/MOTOR/MOTOR/input/namelist.input"

  ! nlFileName = ''
  ! giFileName = ''
  ! ncOutDir = ''

  ! CALL getarg(1, giFileName)
  ! CALL getarg(2, nlFileName)
  ! CALL getarg(3, ncOutDir)

  PRINT *, 'giFileName: ', TRIM(giFileName)
  PRINT *, 'nlFileName: ', TRIM(nlFileName)
  PRINT *, 'ncOutDir: ', TRIM(ncOutDir)

  IF (TRIM(giFileName) .EQ. '' .OR. TRIM(nlFileName) .EQ. '' &
      .OR. TRIM(ncOutDir) .EQ. '') THEN
    PRINT *, 'Wrong input arguements, test failed!'
    STOP
  END IF

  print *, nlFileName

  CALL CPU_TIME(bgn_time)
  GrapesIO = GrapesIO_t(nlFileName, giFileName, qcqrFileName=qcqrFileName, isTest=.true.)
  CALL CPU_TIME(end_time)

  PRINT *, 'SHAPE', SHAPE(GrapesIO%hgrid%u)

  BLOCK
    REAL(r_single), ALLOCATABLE :: lat1D(:), lon1D(:)             !< Degree

    ALLOCATE (lat1D(GrapesIO%hgrid%jds:GrapesIO%hgrid%jde), lon1D(GrapesIO%hgrid%ids:GrapesIO%hgrid%ide))
    FORALL (i=GrapesIO%hgrid%ids:GrapesIO%hgrid%ide) lon1D(i) = &
      GrapesIO%hgrid%config%xs_we + (i - 1)*GrapesIO%hgrid%config%xd
    FORALL (j=GrapesIO%hgrid%jds:GrapesIO%hgrid%jde) lat1D(j) = &
      GrapesIO%hgrid%config%ys_sn + (j - 1)*GrapesIO%hgrid%config%yd

    PRINT *, 'Going here.'

    GrapesIO%zSigma_s(0) = GrapesIO%zSigma_s(1) - GrapesIO%zSigma_s(2)
    CALL writegrid3D_s(TRIM(ncOutDir)//"/grapesinput_qi.nc", GrapesIO%hgrid%ide - GrapesIO%hgrid%ids + 1, &
                       GrapesIO%hgrid%jde - GrapesIO%hgrid%jds + 1, GrapesIO%hgrid%kde - GrapesIO%hgrid%kds + 1, &
                       lon1D, lat1d, GrapesIO%zSigma_s, "qi", GrapesIO%hgrid%qi)
    PRINT *, 'qi is successfully written to the output file'

    DEALLOCATE (lat1D, lon1D)
  END BLOCK

end program GRAPES_INPUT_TEST
