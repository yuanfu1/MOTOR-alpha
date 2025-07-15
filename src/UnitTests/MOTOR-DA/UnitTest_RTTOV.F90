!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-QC.FY4_inst_App
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Yali Wu
! VERSION           : V 0.1
! HISTORY           : 2022-01-30, created by Yali Wu, for FY4 inst DA with WRF.
!
!   Created by Yali Wu (wuyali@gbamwf.com), 2022/01/30, @GBA-MWF, Shenzhen
!
!!--------------------------------------------------------------------------------------------------

!> @brief
!!
PROGRAM UnitTest_RTTOV
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
  USE YAMLRead_m
  USE parameters_m
  USE Obs2State_m
  USE RTTOV_utils_m, ONLY: RTTOV_utils_t
  USE Ctl2State_m, ONLY: Ctl2State_t

  IMPLICIT NONE
  INTEGER(i_kind), PARAMETER :: tlevel01 = 51, tlevel02 = 6
  TYPE(mpddGlob_t), TARGET :: mpddGlob
  TYPE(geometry_t), TARGET :: geometry
  REAL(r_kind) :: t1, t2
  INTEGER(i_kind) :: curr_sg, iv
  LOGICAL :: IncFlag = .TRUE.
  TYPE(MPObs_t), TARGET :: mpObs

  TYPE(State_t) :: X
  TYPE(RTTOV_utils_t) :: RTTOV_utils
  TYPE(Ctl2State_t) :: Ctl2State
  CHARACTER(LEN=1024) :: configFile, ObsFile, ProfFile, StaticDir
  CHARACTER(LEN=10), ALLOCATABLE, DIMENSION(:) :: varList
  INTEGER(i_kind)     :: vLevel, numChans, numVars

  ! Local: for tests only. Delete later.
  INTEGER(i_kind) :: i, j, k, ivar, numObs, numObsTotal
  INTEGER(i_kind) :: n_insts, i_inst, i_chan, ilevel, nlevs
  INTEGER(i_kind) :: nchans, ifile
  CHARACTER(LEN=1024), ALLOCATABLE :: inst_name(:)

  ! Local variables
  REAL(r_kind) :: zenangle = 0.0, azangle = 0.0 , sunzenangle = 0.78539, sunazangle = 0.52359
  REAL(r_kind) :: temp_pert = 0.1D0, qv_pert = 0.001D0
  REAL(4) :: temp_scaling, qv_scaling
  INTEGER(i_kind) :: level_pert
  LOGICAL :: Thin_w_bkg = .TRUE., Use_real_obs = .TRUE.

  CALL CPU_TIME(t1)

  ! Get the configFile
  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", StaticDir)
  configFile = TRIM(StaticDir)//"/UnitTest/UnitTest_RTTOV.yaml"
  ProfFile = TRIM(StaticDir)//"/Satellite/prof.dat"
  PRINT *, 'Test_Sat_fy4_inst Config: ', TRIM(configFile)
  ! ============ Get var info ============
  !CALL namelist_read(TRIM(configFile), 'vLevel', vLevel)
  ifile = yaml_get_var(TRIM(configFile), 'modelState', 'vLevel', vLevel)
  PRINT *, 'vLevel is: ', vLevel
  !CALL namelist_read(TRIM(configFile), "varList", varList)
  ifile = yaml_get_var(TRIM(configFile), 'modelState', 'varList', varList)
  numVars = SIZE(varList, 1)
  PRINT *, 'numVars = ', numVars

  ifile = yaml_get_var(TRIM(configFile), 'RTTOV', 'level_pert', level_pert)
  ifile = yaml_get_var(TRIM(configFile), 'RTTOV', 'temp_scaling', temp_scaling)
  ifile = yaml_get_var(TRIM(configFile), 'RTTOV', "qv_scaling", qv_scaling)
  ifile = yaml_get_var(TRIM(configFile), 'RTTOV', 'inst_name', inst_name) 

  ! Initializer
  CALL mpddGlob%initialize()

  CALL geometry%initialize(configFile, mpddGlob)

  ifile = yaml_get_var(TRIM(configFile), 'FY4-AGRI', 'curr_sg', curr_sg)
  ifile = yaml_get_var(TRIM(configFile), 'FY4-AGRI', 'Thin_w_bkg', Thin_w_bkg)
  ifile = yaml_get_var(TRIM(configFile), 'FY4-AGRI', 'Use_real_obs', Use_real_obs)

  ! Initialize the Scaling implementation
  CALL Ctl2State%initialize(configFile)

  ASSOCIATE (sg => geometry%mg%sg(curr_sg))
    CALL X%initialize(configFile, sg)

    PRINT *, 'Read RTTOV profiles'
    ProfFile = TRIM(StaticDir)//"/Satellite/prof.dat"
    CALL RTTOV_utils%initialize(configFile, X)
    CALL RTTOV_utils%rttovProfile2Bkgd(TRIM(ProfFile), X)
    ! CALL Ctl2State%transBackward(X)

    BLOCK
      TYPE(State_t) :: DD, dX, Xtmp, Xb, grad_J
      TYPE(MPObs_t), TARGET :: mpObs
      TYPE(BMatrix_t) :: B
      TYPE(ObsSet_t) :: Y1, YAngles
      TYPE(ObsSet_t) :: Y2
      TYPE(ObsSet_t) :: Y
      TYPE(C2O_t) :: H
      TYPE(RMatrix_t) :: R
      TYPE(JFunc_t) :: JFunc
      INTEGER(i_kind) :: iobs, it, ivlevel, i_inst
      CHARACTER(len=2) :: tmp
      REAL(r_kind) :: dtemp, dqv, J0, J1, J2, delta_J_temp, delta_J_qv, gradJ_temp, gradJ_qv

      CALL mpObs%initializeMPObs(sg)    ! Initialize the observation parallel processing proc

      IF (.NOT. Use_real_obs) THEN
        CALL Set_MP_ObsFields(configFile, Y1, mpObs, 0.0D0, zenangle, azangle, sunzenangle, sunazangle, sg, 'tbb', 7)
        DO i_inst = 1, SIZE(inst_name)
          SELECT CASE (inst_name(i_inst))
          CASE('agri')
            CALL Y1%ObsFields(1)%Set_ObsType('fy4_1_agri')
            CALL Y1%ObsFields(2)%Set_ObsType('fy4_2_agri')
            CALL Y1%ObsFields(3)%Set_ObsType('fy4_3_agri')
            CALL Y1%ObsFields(4)%Set_ObsType('fy4_4_agri')
            CALL Y1%ObsFields(5)%Set_ObsType('fy4_5_agri')
            CALL Y1%ObsFields(6)%Set_ObsType('fy4_6_agri')
            CALL Y1%ObsFields(7)%Set_ObsType('fy4_7_agri')
          CASE('giirs')
            CALL Y1%ObsFields(1)%Set_ObsType('fy4_1_giirs')
            CALL Y1%ObsFields(2)%Set_ObsType('fy4_2_giirs')
            CALL Y1%ObsFields(3)%Set_ObsType('fy4_3_giirs')
            CALL Y1%ObsFields(4)%Set_ObsType('fy4_4_giirs')
            CALL Y1%ObsFields(5)%Set_ObsType('fy4_5_giirs')
            CALL Y1%ObsFields(6)%Set_ObsType('fy4_6_giirs')
            CALL Y1%ObsFields(7)%Set_ObsType('fy4_7_giirs')
          END SELECT
        END DO
      END IF
      PRINT *, 'FINISH ObsSuper_satellite'

      Y = Y1
      Xb = X
      dX = X
      grad_J = X
      PRINT *, 'Check pres profile: ', X%sg%FGPres(:, 1, 1)
      CALL dX%setAllFieldData(0.0D0)
      CALL grad_J%setAllFieldData(0.0D0)

      CALL R%initialize(configFile, Y, X%sg) ! Initialize R
      CALL B%initialize(configFile, sg, 'Laplace')     ! Initialize the B matrix

      CALL H%initialize(configFile, X, Y)

      ! ! Initialize J Function
      ! PRINT *, Y%ObsFields(1)%values
      IF (sg%isActiveProc()) THEN
        ! Y%ObsFields(1)%values = 0.0
        ! Y%ObsFields(2)%values = 0.0

        CALL JFunc%initialize(configFile, X, Y, H, B, R, sg)

        CALL JFunc%get_J_and_grad_J_vec(X, J0, grad_J)
        PRINT *, 'J0 = ', J0
        ivlevel = level_pert
        ! nlevs = SIZE(X%fields(X%getVarIdx("temp"))%data,1)
        nlevs = SIZE(X%fields(X%getVarIdx("temp"))%DATA, 1)

        gradJ_temp = SUM(grad_J%fields(grad_J%getVarIdx("temp"))%data(ivlevel,  :, :))
        gradJ_qv = SUM(grad_J%fields(grad_J%getVarIdx("qvapor"))%data(ivlevel, :, :))
        ! gradJ_temp = SUM(grad_J%fields(grad_J%getVarIdx("temp"))%data(ivlevel,  :, :))
        ! gradJ_qv = SUM(grad_J%fields(grad_J%getVarIdx("qvapor_ctl"))%data(ivlevel, :, :))
        PRINT *, 'grad_J = ', gradJ_temp, gradJ_qv

        dtemp = temp_pert * temp_scaling
        ! dtemp = dtemp / Xb%sg%STemp(ivlevel,1) ! for temp_ctl
        dqv = qv_pert*qv_scaling
        ! dqv = dqv / Xb%sg%SVapor(ivlevel,1)  ! for qvapor_ctl
        ! dX%fields(dX%getVarIdx("temp"))%data(ivlevel, :, :) = dtemp
        dX%fields(dX%getVarIdx("temp"))%DATA(ivlevel, :, :) = dtemp
        Xtmp = Xb + dX

        CALL JFunc%get_J_value(Xtmp, J1)
        PRINT *, 'J1 = ', J1

        Xtmp = Xb - dX
        CALL JFunc%get_J_value(Xtmp, J2)
        PRINT *, 'J2 = ', J2

        delta_J_temp = (J1 - J2) / (2.0D0 * dtemp)
        PRINT *, 'J2-J1=, delta_J_temp = ', J2 - J1, delta_J_temp

        CALL dX%setAllFieldData(0.0D0)
        dX%fields(dX%getVarIdx("qvapor"))%data(ivlevel, :, :) = dqv
        ! dX%fields(dX%getVarIdx("qvapor_ctl"))%data(ivlevel, :, :) = dqv
        Xtmp = Xb + dX
        CALL JFunc%get_J_value(Xtmp, J1)
        PRINT *, 'J1 = ', J1

        Xtmp = Xb - dX

        CALL JFunc%get_J_value(Xtmp, J2)
        PRINT *, 'J2 = ', J2

        delta_J_qv = (J1 - J2) / (2.0D0 * dqv)
        PRINT *, 'J2-J1=, delta_J_qv = ', J2 - J1, delta_J_qv

        PRINT *, 'T check: ', delta_J_temp - gradJ_temp, delta_J_temp / gradJ_temp
        PRINT *, 'QV check: ', delta_J_qv - gradJ_qv, delta_J_qv / gradJ_qv

        IF (mpddGlob%isBaseProc()) THEN
          ! IF (ABS(delta_J_temp/gradJ_temp - 1.0D0) < 1E-3 .AND. ABS(delta_J_qv/gradJ_qv - 1.0D0) < 1E-3) THEN
          IF (ABS(delta_J_qv / gradJ_qv - 1.0D0) < 1E-7) THEN
            PRINT *, mpddGlob%myrank, "Test passed!"
          ELSE
            PRINT *, "Test failed! result is ", ABS(delta_J_temp/gradJ_temp - 1.0D0), ' AND ', ABS(delta_J_qv/gradJ_qv - 1.0D0)
            PRINT *, 'Temp field is: ', dtemp, X%fields(X%getVarIdx("temp"))%data(ivlevel, 1, 1)
            PRINT *, 'qvapor field is: ', dqv, X%fields(X%getVarIdx("qvapor"))%data(ivlevel, 1, 1)
          END IF
        END IF

        CALL sg%mpddInfo_sg%barrier

      END IF

    END BLOCK

  END ASSOCIATE

  DEALLOCATE (varList)
  CALL mpddGlob%finalize

END PROGRAM UnitTest_RTTOV
