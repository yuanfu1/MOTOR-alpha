!!--------------------------------------------------------------------------------------------------
! PROJECT           : Unit test of a shallow water equation
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : V 1.0
! HISTORY           : 2024-01-17, created by Yuanfu Xie.
!
!   Created by Yuanfu Xie (yuanfu_xie@yahoo.com), 2024/01/17, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------
!!
!> @brief
!! This is a MOTOR-DA 4DVAR with analysis
PROGRAM testShallowWater
  USE kinds_m, ONLY: i_kind, r_kind
  USE geometry_m, ONLY: geometry_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE rhsBase_m, ONLY: rhsBase_t
  USE rhsShallowWater_m, ONLY: rhsShallowWater_t
  USE rhsZGrid_SW_m, ONLY: rhsZGrid_SW_t
  USE State_m, ONLY: State_t
  USE C2MBase_m, ONLY: C2MBase_t
  USE cDefault_m, ONLY: cDefault_t
  USE obsMG_m, ONLY: obsMG_t
  USE C2O_m, ONLY: C2O_t
  USE BMatrix_m, ONLY: BMatrix_t
  USE TimeIntegrationBase_m, ONLY: TimeIntegrationBase_t
  USE TimeIntegrationRK4_m, ONLY: TimeIntegrationRK4_t
  USE TimeIntegrationAB3_m, ONLY: TimeIntegrationAB3_t
  USE dyCoresBase_m, ONLY: dyCoresBase_t
  USE dyCoreGZM_m, ONLY: dyCoreGZM_t
  USE C2MBase_m, ONLY: C2MBase_t
  USE c4DVar_m, ONLY: c4DVar_t
  USE poissonSolver_m, ONLY: poissonSolver_t
  USE bkgMG_m, ONLY: bkgMG_t
  USE GenContainers_m, ONLY: GenContainers_t
  USE State2NC_m
  USE YAMLRead_m

  ! For unit test using Rossby-Haurwitz function:
  USE RossbyHaurwitzSphere_m, ONLY: RossbyHaurwitzSphere_t ! rossbyhaurwitz function from utility

  IMPLICIT NONE

  CHARACTER(LEN=1024) :: yamlFile
  INTEGER(i_kind) :: nic
  REAL(r_kind) :: dt
  ! Multi-processors handler:
  TYPE(mpddGlob_t) :: mpddGlob
  ! Geometry:
  TYPE(geometry_t) :: geometry
  TYPE(GenContainers_t) :: GenContainers

  ! TYPE(rhsShallowWater_t) :: rhs
  TYPE(rhsZGrid_SW_t), ALLOCATABLE :: rhs(:)
  CLASS(TimeIntegrationBase_t), POINTER :: p
  ! TYPE(TimeIntegrationRK4_t) :: tim
  TYPE(TimeIntegrationAB3_t), TARGET :: tim
  TYPE(dyCoreGZM_t) :: gzm
  TYPE(State_t) :: X
  TYPE(State_t), ALLOCATABLE :: ic(:), y(:)   !< ic saves model initial condition, RK4 has 1 time frame; AB3 has 3;
  TYPE(cDefault_t) :: C
  TYPE(obsMG_t) :: O
  TYPE(C2O_t) :: H
  TYPE(BMatrix_t) :: L  ! Temporarily hold for Laplace covariance
  TYPE(BMatrix_t) :: B

  TYPE(c4DVar_t) :: fdvar
  TYPE(bkgMG_t) :: bkg

  TYPE(poissonSolver_t) :: poissonSolver

  ! Unit test of Rossby-Haurwitz function:
  TYPE(RossbyHaurwitzSphere_t), ALLOCATABLE :: rossby(:)

  CHARACTER(LEN=1024) :: ncOutputFile
  INTEGER(i_kind) :: mgStart, mgEnd, ig, i, istatus
  REAL(r_kind), ALLOCATABLE :: anaRHS(:,:), errmx(:,:), valmx(:,:)

  ! CALL mpddGlob%initialize()

  ! Get the configuration file
  CALL getarg(1, yamlFile)

  ! Initialize geometry
  GenContainers = GenContainers_t(TRIM(yamlFile))
  PRINT*,'GenContainers...'
  CALL GenContainers%GenGeometry(geometry)
  PRINT*,'Hello...'
  ! CALL geometry%initialize(yamlFile, mpddGlob)
  mgStart = geometry%mg%mg_coarsest
  mgEnd = geometry%mg%mg_finest
  PRINT*,'Geometry sg bdy_type: ' !,geometry%mg%sg(mgEnd)%bdy_type(100)

  CALL bkg%initialize_s(geometry, yamlFile)
  ! CALL bkg%getBackgrd_s(yamlFile) ! Use of the operation background fields

  CALL X%initialize(yamlFile, geometry%mg%sg(geometry%mg%mg_finest))

  ! Initialize the Poisson solver
  poissonSolver = poissonSolver_t(yamlFile, geometry)

  p => tim
  SELECT TYPE (p)
  TYPE IS (TimeIntegrationAB3_t)
    nic = 3
  TYPE IS (TimeIntegrationRK4_t)
    nic = 1
  CLASS DEFAULT
    WRITE (*, 1)
1   FORMAT('testShallowWater - unsupported time integration scheme!')
    STOP
  END SELECT

  dt = 900.0D0

  ! Unit test of RHS:
  ALLOCATE(rhs(mgStart:mgEnd),ic(mgStart:mgEnd),y(mgStart:mgEnd), &
    rossby(mgStart:mgEnd), errmx(3,mgStart:mgEnd), valmx(3,mgStart:mgEnd))
  istatus = yaml_get_var(yamlFile, 'IO', 'output_dir', ncOutputFile)
  print*,'NC outputfile: ',TRIM(ncOutputFile)
  DO ig=mgStart,mgEnd
    CALL ic(ig)%initialize(yamlFile, geometry%mg%sg(ig), nic)
    CALL y(ig)%initialize(yamlFile, geometry%mg%sg(ig), nic)
    ic(ig) = ic(ig)%zeroCopy()
    y(ig) = y(ig)%zeroCopy()

    ALLOCATE(anaRHS(ic(ig)%sg%num_icell,3))
    CALL rossby(ig)%initializ(ic(ig)%sg%num_icell)
    CALL rossby(ig)%UpdateVar(ic(ig)%sg%num_icell, ic(ig)%sg%cell_cntr, 0.0D0, &
      ic(ig)%fields(4)%DATA(1, :, 3), ic(ig)%fields(1)%DATA(1, :, 3), &
      ic(ig)%fields(3)%DATA(1, :, 3), anaRHS)

    CALL rhs(ig)%initialize_s(geometry%mg%sg(ig), poissonSolver)

    PRINT *, 'Model integration...',nic, ic(ig)%sg%num_icell, UBOUND(ic(ig)%sg%cell_cntr,2)

    CALL rhs(ig)%rightHandSide(3,ic(ig),y(ig))

    ! Plot the fields:
    ic(ig)%fields(1)%DATA(1,:,1) = anaRHS(:,1)
    ic(ig)%fields(2)%DATA(1,:,1) = anaRHS(:,2)
    ic(ig)%fields(3)%DATA(1,:,1) = anaRHS(:,3)
    print*,'MinMax ic 1: ',MINVAL(ic(ig)%fields(1)%DATA(1,:,3)), &
      MAXVAL(ic(ig)%fields(1)%DATA(1,:,3)),ig,ic(ig)%sg%gLevel
    print*,'MinMax ic 3: ',MINVAL(ic(ig)%fields(3)%DATA(1,:,3)), &
      MAXVAL(ic(ig)%fields(3)%DATA(1,:,3))
    print*,'MinMax ic 4: ',MINVAL(ic(ig)%fields(4)%DATA(1,:,3)), &
      MAXVAL(ic(ig)%fields(4)%DATA(1,:,3))
    CALL Output_NC_State_AV(ic(ig), ncOutputFile, &
      "4DVAR"//"_IC", .TRUE., .TRUE.)
    ! Plot the numerical RHS:
    CALL Output_NC_State_AV(y(ig), ncOutputFile, &
      "4DVAR"//"_RHS", .TRUE., .TRUE.)

    ! Maximum errors:
    errmx(1,ig) = MAXVAL(ABS(y(ig)%fields(1)%DATA(1,:,3)-anaRHS(:,1)))
    errmx(2,ig) = MAXVAL(ABS(y(ig)%fields(2)%DATA(1,:,3)-anaRHS(:,2)))
    errmx(3,ig) = MAXVAL(ABS(y(ig)%fields(3)%DATA(1,:,3)-anaRHS(:,3)))
    valmx(1,ig) = MAXVAL(ABS(anaRHS(:,1)))
    valmx(2,ig) = MAXVAL(ABS(anaRHS(:,2)))
    valmx(3,ig) = MAXVAL(ABS(anaRHS(:,3)))

    DEALLOCATE(anaRHS)
  END DO

  ! Output max errs:
  PRINT*,'Output RHS errors: ',errmx(1:3,5)
  PRINT*,'Analytic RHS max: ',valmx(1:3,5)
  DO ig=mgStart+1,mgEnd
    WRITE(*,2) (errmx(i,ig-1)/errmx(i,ig), i=1,3)
2   FORMAT('Errors in RHS reduction ratio: ',3E10.2)
  END DO

  CALL tim%initialization_s(rhs(mgEnd), X, yamlFile)
  CALL gzm%initialize_s(yamlFile, tim, rhs(mgEnd), geometry%mg%sg(mgEnd), dt)

  ! CALL gzm%dyCore_s(dt, ic, X)
  ! Yuanfu Xie added 1 to the second argument for testing only, 2024-10-10
  ! CALL tim%TimeIntegrationBase_fwd(6.0D1, 1, ic, X, y)
  ! CALL tim%TimeIntegrationBase_adj(1.0D-4, X, O, L, B)

  PRINT *, 'Model integration is completed'

  ! ! Test the adjoint of gzm:
  BLOCK

    ! CALL fdvar%initialize_s(mgStart, mgEnd, bkg%XbMG, gzm)
    ! y = fdvar%adjointOpr_s(bkg%XbMG(mgEnd), dt)
    ! CALL gzm%dyAdjoint_s(dt, 3, ic, bkg%XbMG(mgEnd))

  END BLOCK

  ! Deallocate memory:
  DEALLOCATE(rhs,ic,y,rossby)

END PROGRAM testShallowWater
