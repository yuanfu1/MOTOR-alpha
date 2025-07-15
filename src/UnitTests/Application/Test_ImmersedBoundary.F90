!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR.MOTOR-Repository.Geometry
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by YUANFU XIE  (yuanfu_xie@yahoo.com), 2024/08, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This is a test program for the immersed boundary module:
PROGRAM test_immersedBoundary
  USE geometry_m
  USE singleGrid_m, ONLY: singleGrid_t
  USE immersedBoundary_m, ONLY: immersedBoundary_t
  USE mpddGlob_m
  USE bkgMG_m
  !USE linkedList_m, ONLY : linkedList_t,insert,to1DArray

  IMPLICIT NONE

  TYPE(immersedBoundary_t) :: ibm
  TYPE(singleGrid_t) :: sg
  TYPE(mpddGlob_t) :: mpddGlob
  TYPE(geometry_t) :: geometry
  TYPE(bkgMG_t) :: bkg
  !TYPE(linkedList_t), pointer :: head => null()

  CHARACTER(LEN=1024) :: yamlFile
  INTEGER(i_kind) :: numSurfPts, i

  ! Initialize the multiprocessors:
  CALL mpddGlob%initialize()

  ! Get the configuration file
  CALL getarg(1, yamlFile)
  ! Initialize geometry
  CALL geometry%initialize(yamlFile, mpddGlob)

  CALL bkg%initialize_s(geometry, yamlFile)
  CALL bkg%getBackgrd_s(yamlFile)

  PRINT *, 'geometry%mg%mg_finest: ', geometry%mg%mg_finest, geometry%mg%sg(geometry%mg%mg_finest)%num_cell
  CALL ibm%initialization(geometry%mg%sg(geometry%mg%mg_finest))

END PROGRAM test_immersedBoundary
