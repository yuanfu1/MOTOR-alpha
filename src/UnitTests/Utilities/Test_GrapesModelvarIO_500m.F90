!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA.IO.GrapesModelvarIO.test
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
! AUTOHR(S)         : Jiongming Pang
! VERSION           : V 0.0
! HISTORY           :
!   Created by Jiongming Pang (pang.j.m@hotmail.com), 2023/06/28, @SIMI, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!! GrapesModelvarIO module contains the data type for model_states to creating BEC
!! method.
!! @copyright (C) 2020 GBA-MWF, All rights reserved.
! @note
! @warning
! @attention
PROGRAM MODELVAR_INPUT_TEST
  USE GrapesModelvarIO_500m_m, ONLY: GrapesModelvarIO_500m_t
  USE kinds_m, ONLY: i_kind, r_kind, r_single
  USE ncWriteGrid_m
  USE WriteVar_m

  IMPLICIT NONE
  TYPE(GrapesModelvarIO_500m_t)  :: GrapesModelvarIO_500m
  REAL(r_kind)              :: beg_time, end_time, value_tmp
  INTEGER(i_kind)           :: ids, ide, jds, jde, kds, kde, idn, jdn, kdn, i
  INTEGER(i_kind)           :: mdate(6, 1)
  CHARACTER(LEN=1024)       :: modvarFn(1), modvarCtlFn, hfFn, OutDir

  ! REAL(r_kind), ALLOCATABLE, DIMENSION(:) :: ztmp, lon1D, lat1D
  REAL(r_kind), ALLOCATABLE, DIMENSION(:, :) :: lon2D, lat2D
  REAL(r_kind), ALLOCATABLE :: vardata(:, :, :, :)

  mdate(1, 1) = 2022
  mdate(2, 1) = 5
  mdate(3, 1) = 27
  mdate(4, 1) = 0
  mdate(5, 1) = 0
  mdate(6, 1) = 0

  ! modvarFn(1) = "/Users/pang/Documents/SIMI/Data/MOTOR/DATA_GRAPES/CMA-GD_modelvar_v02/2022052500/modelvar042"
  modvarFn(1) = "/home/mgq/proj/MOTOR/input/ana_500m/input/model/modelvar001"
  modvarCtlFn = "/home/mgq/proj/MOTOR/input/ana_500m/input/model/post.ctl"
  hfFn = "/home/mgq/proj/MOTOR/input/ana_500m/input/model/zHight.nc"
  OutDir = "/home/mgq/proj/MOTOR/input/ana_500m/output"

  PRINT *, 'mdate: ', mdate
  PRINT *, 'modvarFn: ', modvarFn
  PRINT *, 'modvarCtlFn: ', TRIM(modvarCtlFn)
  PRINT *, 'hfFn: ', TRIM(hfFn)
  PRINT *, 'OutDir: ', TRIM(OutDir)

  CALL CPU_TIME(beg_time)
  ! GrapesModelvarIO = GrapesModelvarIO_t(TRIM(modvarFn), TRIM(modvarCtlFn), TRIM(hfFn), mdate)
  GrapesModelvarIO_500m = GrapesModelvarIO_500m_t(modvarFn, modvarCtlFn, mdate, 1)
  CALL CPU_TIME(end_time)
  PRINT *, 'cost time GrapesModelvarIO_500m_t: ', end_time - beg_time

  idn = GrapesModelvarIO_500m%idn
  jdn = GrapesModelvarIO_500m%jdn
  kdn = GrapesModelvarIO_500m%kdn - 1

  PRINT *, 'idn, jdn, kdn:', idn, jdn, kdn

  lon2D = GrapesModelvarIO_500m%lon2D
  lat2D = GrapesModelvarIO_500m%lat2D

  CALL writeVar3Dnc(TRIM(OutDir)//"/modelvar_temp.nc", idn, jdn, kdn, & !kdn-1
                    GrapesModelvarIO_500m%lon2D, GrapesModelvarIO_500m%lat2D, "temp", GrapesModelvarIO_500m%temp) !(1:kdn-1)
  CALL writeVar3Dnc(TRIM(OutDir)//"/modelvar_pres.nc", idn, jdn, kdn, &
                    GrapesModelvarIO_500m%lon2D, GrapesModelvarIO_500m%lat2D, "pres", GrapesModelvarIO_500m%pres)
  CALL writeVar3Dnc(TRIM(OutDir)//"/modelvar_uwnd.nc", idn, jdn, kdn, &
                    GrapesModelvarIO_500m%lon2D, GrapesModelvarIO_500m%lat2D, "uwnd", GrapesModelvarIO_500m%uwnd)
  CALL writeVar3Dnc(TRIM(OutDir)//"/modelvar_vwnd.nc", idn, jdn, kdn, &
                    GrapesModelvarIO_500m%lon2D, GrapesModelvarIO_500m%lat2D, "vwnd", GrapesModelvarIO_500m%vwnd)
  CALL writeVar3Dnc(TRIM(OutDir)//"/modelvar_qvapor.nc", idn, jdn, kdn, &
                    GrapesModelvarIO_500m%lon2D, GrapesModelvarIO_500m%lat2D, "qvapor", GrapesModelvarIO_500m%qvapor)

  ALLOCATE (vardata(idn, jdn, kdn + 1, 1))
  vardata(:, :, :, 1) = GrapesModelvarIO_500m%zRHght_u * 1.0D0
  CALL writeVar3Dnc(TRIM(OutDir)//"/modelvar_zRHght_u.nc", idn, jdn, kdn + 1, &
                    GrapesModelvarIO_500m%lon2D, GrapesModelvarIO_500m%lat2D, "zRHght_u", vardata)

  DEALLOCATE (vardata)
  ALLOCATE (vardata(idn, jdn, kdn, 1))
  vardata(:, :, :, 1) = GrapesModelvarIO_500m%zRHght_s * 1.0D0
  CALL writeVar3Dnc(TRIM(OutDir)//"/modelvar_zRHght_s.nc", idn, jdn, kdn, &
                    GrapesModelvarIO_500m%lon2D, GrapesModelvarIO_500m%lat2D, "zRHght_s", vardata)

END PROGRAM MODELVAR_INPUT_TEST
