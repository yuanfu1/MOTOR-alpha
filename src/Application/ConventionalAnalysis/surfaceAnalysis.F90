!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           : 2021-10-27, created by Yuanfu Xie
!                     2022-01-20, unified the modules and types of OBS by Jiongming Pang
!
!   Created by Yuanfu Xie (yuanfu_xie@yahoo.com.com), 2021/10/27, @GBA-MWF, Shenzhen
!   Modified by Jiongming Pang (pang.j.m@hotmail.com), 2022/1/20, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This is a unit test of MOTOR-DA surface analysis
!!   Uses a module of a multiscale surface analytic function to test the multiscale capabity of MOTOR-DA
!
PROGRAM uniTest_sfc
  USE SolverLBFGS_m, ONLY: SolverLBFGS_t
  USE SolverFRCG_m, ONLY: SolverFRCG_t
  USE MiniSolver_m, ONLY: MiniSolver_t
  USE geometry_m, ONLY: geometry_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE MPObs_m, ONLY: MPObs_t
  USE State_m, ONLY: State_t
  USE ObsSet_m, ONLY: ObsSet_t
  USE C2O_m, ONLY: C2O_t
  USE BMatrix_m, ONLY: BMatrix_t
  USE JFunc_m, ONLY: JFunc_t
  USE State2NC_m
  USE Mock_m

  ! Yuanfu Xie adds a real data test:
  ! Unified the modules and types of OBS. -Jiongming Pang 2022-01-20
  USE ObsBase_m, ONLY: ObsBase_t
  USE ObsSurface_m, ONLY: ObsSurface_t

  USE unitTest_sfc_m, ONLY: unitTest_sfc_t

  ! Define types
  TYPE(mpddGlob_t), TARGET :: mpddGlob
  TYPE(geometry_t), TARGET :: geometry
  TYPE(MPObs_t), TARGET :: mpObs
  TYPE(State_t) :: X
  TYPE(ObsSet_t) :: Y
  TYPE(C2O_t) :: H
  TYPE(BMatrix_t) :: B
  TYPE(JFunc_t) :: JFunc
  REAL(r_kind) :: t1, t2

  CLASS(MiniSolver_t), ALLOCATABLE :: miniSolver

  CHARACTER(LEN=1024) :: configFile, ncOutputFile !< Config and output file name.

  ! Yuanfu Xie test:
  CHARACTER(LEN=1024) :: filename
  CHARACTER(LEN=20), ALLOCATABLE :: varname(:)
  INTEGER(i_kind) :: numVars, i
  TYPE(ObsSurface_t) :: obs_in
  TYPE(unitTest_sfc_t) :: unitTest

  REAL(r_kind), ALLOCATABLE :: xyt(:, :)

  ! Get the configFile
  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", configFile)
  configFile = TRIM(configFile)//"/Application/App_3DVAR_SimpleCase.yaml"

  IF (yaml_get_var(configFile, 'IO', 'output_dir', outputDir) /= 0) STOP
  ncOutputFile = TRIM(ncOutputFile)//"/test.nc"

  ! Initializer
  ! Auxtypes
  CALL mpddGlob%initialize()                                   ! Initialize the mpdd

  CALL mpddGlob%barrier
  CALL CPU_TIME(t1)

  ! Initialize geometry
  CALL geometry%initialize(configFile, mpddGlob)               ! Initialize the geometry

  ! Initialize solver
  ALLOCATE (miniSolver)
!  CALL miniSolver%initialize(configFile)
  CALL miniSolver%initialize(configFile)

  ASSOCIATE (sg => geometry%mg%sg(8))
    CALL mpObs%initializeMPObs(sg)    ! Initialize the observation parallel processing proc

    ! Initialize X, Y, B, H
    CALL X%initialize(configFile, sg)       ! Initialize the state
    CALL B%initialize(configFile, sg)     ! Initialize the B matrix
    ! PRINT *, X%fields(1)%data
    X%fields(1)%DATA = 0.0D0
    ! CALL Output_NC_State_SV(X, ncOutputFile, "ua")

    ! Mock observation data:
    IF (.FALSE.) THEN ! Fake data test from Zilong Qin
      CALL Set_Mock_Single_Pt(configFile, Y, mpObs, 1.0D0, sg)     ! Initialize Y
    ELSE              ! Real data test from Yuanfu Xie
      filename = '../data/20210906_1210'
      numVars = 1
      ALLOCATE (varname(numVars))
      varname(1) = 'temperature'

      ! Initialize obs:
      CALL obs_in%initial(numVars, configFile, mpObs)

      ! Ingest obs:
      CALL obs_in%obsIngest(filename, numVars, varname)

      ! Use of unit test analytic observation data:
      ALLOCATE (xyt(3, obs_in%numObs))
      xyt(1, :) = (obs_in%olatlon(1, 1:obs_in%numObs) - MINVAL(sg%cell_cntr(1, :))) / &
                  (MAXVAL(sg%cell_cntr(1, :)) - MINVAL(sg%cell_cntr(1, :)))
      xyt(2, :) = (obs_in%olatlon(2, 1:obs_in%numObs) - MINVAL(sg%cell_cntr(2, :))) / &
                  (MAXVAL(sg%cell_cntr(2, :)) - MINVAL(sg%cell_cntr(2, :)))
      xyt(3, :) = obs_in%obsTime(1:obs_in%numObs) - sg%tt(1)
      CALL unitTest%analytic(obs_in%numObs, xyt, obs_in%obsData)

      ! Thinning observations:
      CALL obs_in%superObs(X, Y)
    END IF
    CALL H%initialize(configFile, X, Y)
    CALL Output_NC_State_SV(X, './test.nc', "temp")

    PRINT *, 'temp: ', Y%ObsFields(1)%values
    CALL Y%ObsFields(1)%Set_Name('temp')
    ! Initialize J Function
    CALL JFunc%initialize(configFile, X, Y, H, B, sg)                                ! Initialize H
    ! goto 1

    ! Run minimization
    CALL miniSolver%run(X, JFunc, sg)

  END ASSOCIATE
1 CONTINUE
  CALL mpddGlob%barrier

  ! Output state file to NC for view.
  CALL Output_NC_State_SV(X, ncOutputFile, "temp")
  CALL CPU_TIME(t2)
  IF (mpddGlob%isBaseProc()) PRINT *, 'Time cost:', t2 - t1, 's'

  DEALLOCATE (xyt)
  DEALLOCATE (miniSolver)

  ! Finalize
  CALL mpddGlob%finalize

END PROGRAM uniTest_sfc
