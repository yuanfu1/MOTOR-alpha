!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA.TestMixedSolver.F90
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring 
!                     Warning and Forecasting (GBA-MWF) Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Yuanfu Xie on 2024-12-09
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This test program tests a mixed solver of Psi and Chi with U and V boundary conditions.
!! @copyright (C) 2020 GBA-MWF, All rights reserved.
! @note
! @warning
! @attention
PROGRAM testPsiChi2UV
  USE kinds_m, ONLY: i_kind, r_kind
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE geometry_m, ONLY: geometry_t

  IMPLICIT NONE

  CHARACTER(LEN=1024) :: yamlFile
  INTEGER(i_kind) :: ig
  LOGICAL :: pass
  REAL(r_kind), ALLOCATABLE :: errMax(:)
  REAL(r_kind), ALLOCATABLE :: psi(:,:),chi(:,:),uan(:,:),van(:,:),vor(:,:),div(:,:)
  TYPE(mpddGlob_t) :: mpddGlob
  TYPE(geometry_t) :: geo

  ! Fix the yaml file for this test:
  ! yamlFile = '../static/unitests/psiChi2UV.yaml'
  CALL getarg(1, yamlFile)
  IF (TRIM(yamlFile) .EQ. '') THEN
    PRINT*,'Usage: <this executable> <yaml file>, please and find and supply your yamlfile and rerun'
    STOP
  ELSE
    PRINT*,'yaml: ',TRIM(yamlFile)
  END IF

  ! Set up MPDD and Geometry:
  CALL mpddGlob%initialize()
  CALL geo%initialize(yamlFile,mpddGlob)

  pass = .TRUE.
  ALLOCATE(errMax(geo%mg%mg_coarsest:geo%mg%mg_finest))

  DO ig=geo%mg%mg_coarsest,geo%mg%mg_finest
    ASSOCIATE(sg => geo%mg%sg(ig), vlvl => geo%mg%sg(ig)%vLevel, num_cell => geo%mg%sg(ig)%num_cell)
      ALLOCATE(psi(vlvl,num_cell),chi(vlvl,num_cell),uan(vlvl,num_cell),van(vlvl,num_cell), &
               vor(vlvl,num_cell),div(vlvl,num_cell))
      ! Get the analytic functions:
      CALL LinearLatCase(sg%num_cell,sg%vLevel,psi,chi,vor,div,uan,van,sg%cell_cntr)

      ! Calculate the UV based on the analytic psi and chi:
      CALL checkPsiChi2UV_s(psi,chi,uan,van,sg,errMax(ig))

      DEALLOCATE(psi,chi,uan,van,vor,div)
    END ASSOCIATE
  END DO

  DO ig=geo%mg%mg_coarsest+1,geo%mg%mg_finest
    IF (errMax(ig-1)/errMax(ig) .LT. 3.7) pass = .FALSE.
    PRINT*,'Error reduction in PsiChi to UV: ',errMax(ig-1)/errMax(ig)
  END DO

  IF (pass) THEN
    PRINT*,"Test passed!"
  ELSE
    PRINT*,"Test failed!"
  END IF

  ! Deallocate memory:
  DEALLOCATE(errMax)
END PROGRAM testPsiChi2UV