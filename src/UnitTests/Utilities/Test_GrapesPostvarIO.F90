!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA.IO.GrapesPostvarIO.test
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
! AUTOHR(S)         : Jiongming Pang
! VERSION           : V 0.0
! HISTORY           :
!   Created by Jiongming Pang (pang.j.m@hotmail.com), 2023/08/09, @SIMI, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!! GrapesPostvarIO module contains the data type for model_states
!! method.
!! @copyright (C) 2020 GBA-MWF, All rights reserved.
! @note
! @warning
! @attention
PROGRAM MODELVAR_INPUT_TEST
  USE GrapesPostvarIO_m, ONLY: GrapesPostvarIO_t
  USE kinds_m, ONLY: i_kind, r_kind, r_single
  USE parameters_m, ONLY: degree2radian
  USE ncWriteGrid_m
  USE WriteVar_m

  IMPLICIT NONE
  TYPE(GrapesPostvarIO_t) :: GrapesPostvarIO
  REAL(r_kind)            :: beg_time, end_time, value_tmp
  INTEGER(i_kind)         :: ids, ide, jds, jde, kds, kde, idn, jdn, kdn, i
  INTEGER(i_kind)         :: mdate(6, 1)
  CHARACTER(LEN=1024)     :: postvarFn(1), postvarCtlFn, OutDir

  REAL(r_kind), ALLOCATABLE, DIMENSION(:, :) :: lon2D, lat2D
  REAL(r_kind), ALLOCATABLE :: vardata(:, :, :)

  mdate(1, 1) = 2023
  mdate(2, 1) = 8
  mdate(3, 1) = 8
  mdate(4, 1) = 6
  mdate(5, 1) = 0
  mdate(6, 1) = 0

  postvarFn(1) = "/Users/pang/Documents/SIMI/Data/MOTOR/DATA_GRAPES/CMA-GD_postvar/2023080800/postvar006"
  postvarCtlFn = "/Users/pang/Documents/SIMI/Data/MOTOR/DATA_GRAPES/CMA-GD_postvar/2023080800/post.ctl"
  OutDir = "/Users/pang/Documents/SIMI/Data/MOTOR/MOTOR/output"

  PRINT *, 'mdate: ', mdate
  PRINT *, 'postvarFn: ', postvarFn
  PRINT *, 'postvarCtlFn: ', TRIM(postvarCtlFn)
  PRINT *, 'OutDir: ', TRIM(OutDir)

  CALL CPU_TIME(beg_time)
  GrapesPostvarIO = GrapesPostvarIO_t(postvarFn, postvarCtlFn, mdate)
  CALL CPU_TIME(end_time)
  PRINT *, 'cost time GrapesPostvarIO_t: ', end_time - beg_time

  idn = GrapesPostvarIO%idn
  jdn = GrapesPostvarIO%jdn
  kdn = GrapesPostvarIO%kdn

  PRINT *, 'idn, jdn, kdn:', idn, jdn, kdn

  ALLOCATE (lon2D(idn, jdn), lat2D(idn, jdn))
  lon2D = GrapesPostvarIO%lon2D / degree2radian
  lat2D = GrapesPostvarIO%lat2D / degree2radian
  PRINT *, 'min and max of lat:', MINVAL(lat2D), MAXVAL(lat2D)
  PRINT *, 'min and max of lon:', MINVAL(lon2D), MAXVAL(lon2D)

  CALL writeVar2Dnc(TRIM(OutDir)//"/postvar_t2m.nc", idn, jdn, lon2D, lat2D, "t2m", GrapesPostvarIO%t2m * 1.0D0)
  CALL writeVar2Dnc(TRIM(OutDir)//"/postvar_ts.nc", idn, jdn, lon2D, lat2D, "ts", GrapesPostvarIO%ts * 1.0D0)
  CALL writeVar2Dnc(TRIM(OutDir)//"/postvar_q2m.nc", idn, jdn, lon2D, lat2D, "q2m", GrapesPostvarIO%q2m * 1.0D0)
  CALL writeVar2Dnc(TRIM(OutDir)//"/postvar_psl.nc", idn, jdn, lon2D, lat2D, "psl", GrapesPostvarIO%psl * 1.0D0)
  CALL writeVar2Dnc(TRIM(OutDir)//"/postvar_u10m.nc", idn, jdn, lon2D, lat2D, "u10m", GrapesPostvarIO%u10m * 1.0D0)
  CALL writeVar2Dnc(TRIM(OutDir)//"/postvar_v10m.nc", idn, jdn, lon2D, lat2D, "v10m", GrapesPostvarIO%v10m * 1.0D0)
  CALL writeVar2Dnc(TRIM(OutDir)//"/postvar_psfc.nc", idn, jdn, lon2D, lat2D, "ps", GrapesPostvarIO%psfc * 1.0D0)

  ALLOCATE (vardata(idn, jdn, 1))
  vardata(:, :, 1) = GrapesPostvarIO%zs * 1.0D0
  CALL writeVar2Dnc(TRIM(OutDir)//"/postvar_zs.nc", idn, jdn, lon2D, lat2D, "zs", vardata)

  DEALLOCATE (vardata)

END PROGRAM MODELVAR_INPUT_TEST
