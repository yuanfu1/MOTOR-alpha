!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather
! Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Jilong CHEN, Zilong Qin
! VERSION           : V 0.0
! HISTORY           :
!   Created by Jilong CHEN (jchen@link.cuhk.edu.hk), 2021/1/26, @GBA-MWF,
!   Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!!
PROGRAM Cum_Jfun_test
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
  USE RMatrix_m, ONLY: RMatrix_t
  USE JFunc_m, ONLY: JFunc_t
  USE MGOpts_m
  USE State2NC_m
  USE Mock_m
  USE ObsField_m, ONLY: ObsField_t
  USE kinds_m, ONLY: i_kind, r_kind
  !USE NMLRead_m, ONLY: namelist_read
  USE YAMLRead_m
  USE parameters_m
  USE Obs2State_m
  USE UnitTestData_m

  IMPLICIT NONE

  CHARACTER(LEN=1024) :: configFile
  TYPE(mpddGlob_t), TARGET :: mpddGlob
  TYPE(geometry_t), TARGET :: geometry
  TYPE(State_t) :: X, delX, Xtmp, Xb, grad_J

  TYPE(BMatrix_t) :: B
  TYPE(MPObs_t), TARGET :: mpObs

  TYPE(ObsSet_t) :: Y
  TYPE(C2O_t) :: H
  TYPE(RMatrix_t) :: R
  TYPE(JFunc_t) :: JFunc
  ! INTEGER(i_kind) :: iobs, it, ivlevel
  CHARACTER(len=2) :: tmp
  REAL(r_kind) :: dtemp, dqv, J0, J1, J2, delta_J_temp, delta_J_qv, gradJ_temp, gradJ_qv

  REAL(r_kind) :: t1, t2

  PRINT *, 'TEST START'
  CALL CPU_TIME(t1)

  ! Get the configFile

  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", configFile)
  configFile = TRIM(configFile)//"/UnitTest/testCum.yaml"

  ! Initialize
  CALL mpddGlob%initialize()
  CALL geometry%initialize(configFile, mpddGlob)          ! Initialize the geometry

  ASSOCIATE (sg => geometry%mg%sg(1))
    CALL X%initialize(configFile, sg)
    delCALL X%initialize(configFile, sg)
    CALL XB%initialize(configFile, sg)

    CALL mpObs%initializeMPObs(geometry%mg%sg(geometry%mg%mg_finest)) ! Initialize the observation parallel processing proc

    ! CALL Set_MP_ObsFields(configFile, Y, mpObs, 0.0D0, sg, 2)

    CALL UnitTestData(X)

    X%Fields(1)%DATA(:, :, 2) = 1.01D0 * X%Fields(1)%DATA(:, :, 1)
    X%Fields(2)%DATA(:, :, 2) = 1.01D0 * X%Fields(2)%DATA(:, :, 1)
    X%Fields(3)%DATA(:, :, 2) = 1.0D0 * X%Fields(3)%DATA(:, :, 1)
    X%Fields(4)%DATA(:, :, 2) = 1.0D0 * X%Fields(4)%DATA(:, :, 1)
    X%Fields(5)%DATA(:, :, 2) = 1.0D0 * X%Fields(5)%DATA(:, :, 1)
    X%Fields(6)%DATA(:, :, 2) = 1.0D0 * X%Fields(6)%DATA(:, :, 1)

    CALL B%initialize(configFile, sg)     ! Initialize the B matrix

    ! BLOCK
    Y = ObsSet_t(configFile, mpObs)

    IF (ALLOCATED(Y%ObsFields)) STOP

    ALLOCATE (Y%ObsFields(2))

    Y%ObsFields(1) = ObsField_t(configFile, mpObs)
    Y%ObsFields(2) = ObsField_t(configFile, mpObs)

    CALL Y%ObsFields(1)%Set_Name("theta")
    CALL Y%ObsFields(2)%Set_Name("qvapor")
    ! CALL Y%ObsFields(6)%Set_Name("pres")

    ALLOCATE (Y%ObsFields(1)%idx(2))
    ALLOCATE (Y%ObsFields(1)%values(2))

    ALLOCATE (Y%ObsFields(2)%idx(2))
    ALLOCATE (Y%ObsFields(2)%values(2))

    Y%ObsFields(1)%idx(1)%tIdx = 1
    Y%ObsFields(1)%idx(1)%hIdx = 5
    Y%ObsFields(1)%idx(1)%vIdx = 20
    Y%ObsFields(1)%values(1) = 312.0D0

    Y%ObsFields(1)%idx(2)%tIdx = 2
    Y%ObsFields(1)%idx(2)%hIdx = 5
    Y%ObsFields(1)%idx(2)%vIdx = 20
    Y%ObsFields(1)%values(2) = 312.0D0

    Y%ObsFields(2)%idx(1)%tIdx = 1
    Y%ObsFields(2)%idx(1)%hIdx = 5
    Y%ObsFields(2)%idx(1)%vIdx = 20
    Y%ObsFields(2)%values(1) = 0.001D0

    Y%ObsFields(2)%idx(2)%tIdx = 2
    Y%ObsFields(2)%idx(2)%hIdx = 5
    Y%ObsFields(2)%idx(2)%vIdx = 20
    Y%ObsFields(2)%values(2) = 0.001D0

    Xb = X
    delX = X
    grad_J = X
    CALL delX%setAllFieldData(0.0D0)
    CALL grad_J%setAllFieldData(0.0D0)

    CALL R%initialize(configFile, Y, X%sg) ! Initialize R

    H = C2O_t(configFile, X, Y)

    BLOCK
      TYPE(ObsSet_t) :: D, DD
      DD = H%fwdNL_opr(X)
      D = H%fwdNL_opr(X) - Y

    END BLOCK

    CALL JFunc%initialize(configFile, X, Y, H, B, R, sg)

    CALL JFunc%get_J_value(X, J0)
    CALL JFunc%get_grad_J_vec(X, grad_J)
    PRINT *, 'J0 = ', J0
    ! ivlevel = level_pert
    gradJ_temp = grad_J%fields(4)%DATA(20, 5, 1)
    gradJ_qv = grad_J%fields(5)%DATA(20, 5, 1)

    PRINT *, 'grad_J = ', gradJ_temp, gradJ_qv

    dtemp = 5.0D-2
    dqv = 1.0D-3
    delX%fields(4)%DATA(20, 5, 1) = dtemp

    Xtmp = Xb + delX

    CALL JFunc%get_J_value(Xtmp, J1)
    PRINT *, 'J1 = ', J1

    Xtmp = Xb - delX

    CALL JFunc%get_J_value(Xtmp, J2)
    PRINT *, 'J2 = ', J2

    delta_J_temp = (J1 - J2) / (2.0 * dtemp)
    PRINT *, 'delta_J_temp = ', delta_J_temp

    CALL delX%setAllFieldData(0.0D0)
    delX%fields(5)%DATA(20, 5, 1) = dqv
    Xtmp = Xb + delX

    CALL JFunc%get_J_value(Xtmp, J1)
    PRINT *, 'J1 = ', J1

    Xtmp = Xb - delX

    CALL JFunc%get_J_value(Xtmp, J2)
    PRINT *, 'J2 = ', J2

    delta_J_qv = (J1 - J2) / (2.0 * dqv)
    PRINT *, 'delta_J_qv = ', delta_J_qv

    CALL mpddGlob%barrier

    PRINT *, 'T check: ', delta_J_temp - gradJ_temp, delta_J_temp / gradJ_temp
    PRINT *, 'QV check: ', delta_J_qv - gradJ_qv, delta_J_qv / gradJ_qv

    IF (delta_J_temp / gradJ_temp - 1.0D0 < 1E-3 .AND. delta_J_qv / gradJ_qv - 1.0D0 < 1E-3) THEN
      PRINT *, "tested passed!"
    ELSE
      PRINT *, "tested failed!"
    END IF

    ! END BLOCK

  END ASSOCIATE

  ! END DO

  CALL CPU_TIME(t2)
  PRINT *, 'Time cost:', t2 - t1, 's'

  PRINT *, 'DEALLOCATE OVER'

  ! CALL mpddGlob%finalize

  PRINT *, 'finalize OVER'

  CALL mpddGlob%finalize

END PROGRAM Cum_Jfun_test

