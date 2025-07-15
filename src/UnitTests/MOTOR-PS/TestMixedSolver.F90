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
PROGRAM testMixedSolver
  USE kinds_m, ONLY: i_kind, r_kind
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE geometry_m, ONLY: geometry_t
  USE mixedPsiChiSolver_m, ONLY: mixedPsiChiSolver_t
  USE State_m, ONLY: State_t
  USE State2NC_m
  USE YAMLRead_m

  ! Plot:
  USE gzm_m, ONLY: gzm_t

  IMPLICIT NONE

  TYPE(mpddGlob_t) :: mpddGlob
  TYPE(geometry_t) :: geo
  TYPE(mixedPsiChiSolver_t) :: mix

  ! For plot:
  TYPE(State_t) :: X
  TYPE(gzm_t) :: gzm

  CHARACTER(LEN=1024) :: yamlFile, ncOutputFile
  CHARACTER(LEN=20) :: task
  INTEGER(i_kind) :: mgStart, mgEnd   ! mgEnd: Glevel solution is solved
  INTEGER(i_kind) :: istatus,vlevel,ngauss,nedges,ncells
  REAL(r_kind), ALLOCATABLE :: vor(:,:),div(:,:),u(:,:),v(:,:), &
                               psi(:,:),chi(:,:),dpsi(:,:),dchi(:,:)

  ! Get the configuration file
  CALL getarg(1, yamlFile)
  

  ! Set up MPDD and Geometry:
  CALL mpddGlob%initialize()

  print*,'TestMixed yaml 1: ',TRIM(yamlFile)
  CALL geo%initialize(yamlFile,mpddGlob)
  print*,'TestMixed yaml 2: ',TRIM(yamlFile),' pc: ',geo%mg%sg(2)%mpddInfo_sg%myrank
  !STOP

  ! Initialize the mix solver:
  mgStart = 2; mgEnd = 2  ! mgEnd is the level where the solution is calculated not neccesary the finest mg level
  CALL mix%initialize_s(mgStart,mgEnd,geo)

  ! Plotting calculated vorticity:
  gzm = gzm_t(geo%mg%sg(mgEnd))
  
  ! Vorticity, divergence, u and v:
  vlevel = geo%mg%sg(mgEnd)%vLevel
  ngauss = geo%mg%sg(mgEnd)%numQuadPerEdge
  nedges = geo%mg%sg(mgEnd)%num_edge
  ncells = geo%mg%sg(mgEnd)%num_cell
  ALLOCATE(vor(vlevel,ncells),div(vlevel,ncells))
  ALLOCATE(psi(vlevel,ncells),chi(vlevel,ncells),dpsi(vlevel,ncells),dchi(vlevel,ncells))
  ALLOCATE(u(vlevel,ncells),v(vlevel,ncells))

  ! Initial values of the required fields:
  vor = 0.0D0; div = 0.0D0; u = 0.0D0; v = 0.0D0

  ! Plot the solutions:
  istatus = yaml_get_var(TRIM(yamlFile), 'IO', 'output_dir', ncOutputFile)
  istatus = yaml_get_var(TRIM(yamlFile), 'RunMode', 'Task', task)
  CALL X%initialize(yamlFile,geo%mg%sg(mgEnd))

  ! Get analytic functions:
  CALL analytics(psi,chi,vor,div,u,v,geo%mg%sg(mgEnd))
  X%fields(1)%DATA(:,:,2) = dpsi
  X%fields(2)%DATA(:,:,2) = dchi
  X%fields(3)%DATA(:,:,2) = u
  X%fields(4)%DATA(:,:,2) = v
  X%fields(5)%DATA(:,:,2) = vor
  X%fields(6)%DATA(:,:,2) = div

  ! Test an ideal initial guess: the actual initial guess will be determined in the mixed solver
  dpsi = psi
  dchi = chi

  ! Solve:
  CALL mix%mixedSolver_s(geo,vor,div,u,v,dpsi,dchi)

  ! Plot calculcated vorticity:
  CALL gzm%Laplacia(dpsi,vor)
  CALL gzm%Laplacia(dchi,div)

  CALL uvVelocityOnInterior(dpsi, dchi, u, v, geo%mg%sg(mgEnd))

  X%fields(1)%DATA(:,:,1) = dpsi
  X%fields(2)%DATA(:,:,1) = dchi
  X%fields(3)%DATA(:,:,1) = u
  X%fields(4)%DATA(:,:,1) = v
  X%fields(5)%DATA(:,:,1) = vor
  X%fields(6)%DATA(:,:,1) = div
  CALL Output_NC_State_AV(X, TRIM(ncOutputFile), TRIM(task)//"_mixSoln", .TRUE., .TRUE.)

  ! Deallocate memory:
  DEALLOCATE(vor,div,u,v,psi,chi,dpsi,dchi)

  PRINT*,'TestMixedSolver: Finished!'

END PROGRAM testMixedSolver
