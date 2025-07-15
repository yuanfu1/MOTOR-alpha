!!--------------------------------------------------------------------------------------------------
! PROJECT           : Grid generation
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : Beta 0.0
! HISTORY           :
!   Created by Yuanfu Xie (yuanfu_xie@yahoo.com), 2020/06/11, @GBA-MWF, Shenzhen
!   Modified by ..... (???@gmail.com), YYYY/MM/DD, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!!===================================================================
!> @brief
!! # Icosahedral triangle grid generation routine collection
!!
!!  *This file defines a set of routines implementing icosahedral
!!   triangle grid*
!!
!! ## Routine set:
!!
!! ### getDomain_icos
!!  *The routine reads domain parameters*
!!
!! @author Yuanfu Xie
!! @copyright (C) 2020 GBA-MWF, All rights reserved.
!!===================================================================
SUBROUTINE getDomain_icosTr(configFile, gl, ls, le)
  USE kinds_m, ONLY: i_kind, r_kind
  USE YAMLRead_m
  IMPLICIT NONE

  CHARACTER(LEN=1024), INTENT(IN) :: configFile
  INTEGER(i_kind)                 :: rc, istatus

  ! Numbers of vertices and cells for all levels:
  ! G lvl, lvl start, lvl end, order of accuracy
  INTEGER(i_kind), INTENT(OUT) :: gl, ls, le

  ! Local variables:
  CHARACTER(LEN=256) :: path2static
  ! glvl: icos G level; slvl: start level;
  ! elvl: end level; ordr: order of accuracy
  INTEGER(i_kind) :: i, glvl, slvl
  !NAMELIST /icosTr_grid/glvl,slvl

  ! Check whether file exists.
  INQUIRE (file=TRIM(configFile), iostat=rc)

  ! Read the namelist file
  IF (rc .EQ. 0) THEN
    !OPEN (10, FILE=TRIM(configFile))
    !READ (10, NML=icosTr_grid)
    !CLOSE (10)
    istatus = yaml_get_var(TRIM(configFile), 'icosTr_grid', 'glvl', glvl)
    istatus = yaml_get_var(TRIM(configFile), 'icosTr_grid', 'slvl', slvl)
  ELSE
    ERROR STOP "Namelist file of *getDomain_icosTr* does not exist!"
  END IF

  gl = glvl
  ls = slvl
  le = glvl

END SUBROUTINE getDomain_icosTr
