!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2021/1/26, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!!
PROGRAM Test_SolverFRCG
  USE SolverFRCG_m, ONLY: SolverFRCG_t
  USE geometry_m, ONLY: geometry_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE MPObs_m, ONLY: MPObs_t
  USE State_m, ONLY: State_t
  USE ObsSet_m, ONLY: ObsSet_t
  USE C2O_m, ONLY: C2O_t
  USE BMatrix_m, ONLY: BMatrix_t
  USE RMatrix_m, ONLY: RMatrix_t
  USE JFunc_m, ONLY: JFunc_t
  USE State2NC_m
  USE Mock_m
  USE YAMLRead_m
  TYPE(mpddGlob_t), TARGET :: mpddGlob

  ! Define types
  BLOCK

    TYPE(geometry_t), TARGET :: geometry
    TYPE(MPObs_t), TARGET :: mpObs
    TYPE(State_t) :: X
    TYPE(ObsSet_t) :: Y
    TYPE(C2O_t) :: H
    TYPE(BMatrix_t) :: B
    TYPE(RMatrix_t) :: R
    TYPE(JFunc_t) :: JFunc
    TYPE(SolverFRCG_t) :: miniSolver
    INTEGER(i_kind) :: fileStat

    CHARACTER(LEN=1024) :: configFile, outputDir !< Config and output file name.

    ! Get the configFile
    CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", configFile)
    configFile = TRIM(configFile)//"/UnitTest/testMiniSolver.yaml"

    fileStat = yaml_get_var(configFile, 'IO', 'output_dir', outputDir)
    ! Initializer
    ! Auxtypes
    CALL mpddGlob%initialize()
    CALL mpddGlob%barrier
    CALL CPU_TIME(t1)                               ! Initialize the mpdd

    CALL geometry%initialize(configFile, mpddGlob)               ! Initialize the geometry

    ASSOCIATE (sg => geometry%mg%sg(8))
      CALL mpObs%initializeMPObs(sg)    ! Initialize the observation parallel processing proc

      ! Initialize X, Y, B, H
      CALL X%initialize(configFile, sg)       ! Initialize the state
      CALL B%initialize(configFile, sg, 'Laplace')     ! Initialize the B matrix
      ! PRINT *, X%fields(1)%data
      X%fields(1)%DATA = 0.0D0

      ! Mock observation data:
      CALL Set_Mock_Single_Pt(configFile, Y, mpObs, 1.0D0, sg, 'temp')     ! Initialize Y
      CALL H%initialize(configFile, X, Y)                                  ! Initialize H

      CALL R%initialize(configFile, Y, X%sg) ! Initialize R

      ! Initialize J Function
      CALL JFunc%initialize(configFile, X, Y, H, B, R, sg, X)

      ! Initialize solver
      CALL miniSolver%initialize(configFile)

      ! Run minimization
      CALL miniSolver%run(X, JFunc, sg)

    END ASSOCIATE
    CALL mpddGlob%barrier

    ! Output state file to NC for view.
    CALL Output_NC_State_SV(X, outputDir, 'test', "temp")
    CALL CPU_TIME(t2)
    IF (mpddGlob%isBaseProc()) PRINT *, 'Time cost:', t2 - t1, 's'
    IF (mpddGlob%isBaseProc()) THEN
      PRINT *, 'J: ', miniSolver%J
      IF (J < 527.0D0) THEN
        PRINT *, 'Test passed!'
      ELSE
        PRINT *, 'Test failed!'
      END IF
    END IF
    ! Finalize
  END BLOCK
  CALL mpddGlob%finalize

END PROGRAM Test_SolverFRCG
