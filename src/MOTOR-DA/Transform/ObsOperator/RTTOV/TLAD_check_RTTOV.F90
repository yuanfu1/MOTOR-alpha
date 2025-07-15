!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-QC.FY4_AGRI_App
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Yali Wu
! VERSION           : V 0.1
! HISTORY           : 2022-01-30, created by Yali Wu, for FY4 AGRI DA with WRF.
!
!   Created by Yali Wu (wuyali@gbamwf.com), 2022/01/30, @GBA-MWF, Shenzhen
!
!!--------------------------------------------------------------------------------------------------

!> @brief
!!
PROGRAM TLAD_check_RTTOV
  USE geometry_m, ONLY: geometry_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE MPObs_m, ONLY: MPObs_t
  USE State_m, ONLY: State_t
  USE ObsSet_m, ONLY: ObsSet_t
  USE C2O_m, ONLY: C2O_t
  USE MGOpts_m
  USE State2NC_m
  USE Mock_m
  USE ObsField_m, ONLY: ObsField_t
  USE kinds_m, ONLY: i_kind, r_kind
  !USE NMLRead_m, ONLY: namelist_read
  USE YAMLRead_m
  USE parameters_m
  USE Obs2State_m
  USE RTTOV_utils_m, ONLY: RTTOV_utils_t
  USE Ctl2State_m, ONLY: Ctl2State_t
  ! USE ObsUtilities_m

  IMPLICIT NONE
  INTEGER(i_kind), PARAMETER :: tlevel01 = 51, tlevel02 = 6
  TYPE(mpddGlob_t), TARGET :: mpddGlob
  TYPE(geometry_t), TARGET :: geometry
  TYPE(Ctl2State_t) :: Ctl2State
  TYPE(State_t) :: XRes
  REAL(r_kind) :: t1, t2
  INTEGER(i_kind) :: curr_sg, iv
  LOGICAL :: IncFlag = .TRUE.
  TYPE(MPObs_t), TARGET :: mpObs

  TYPE(State_t) :: X
  TYPE(RTTOV_utils_t) :: RTTOV_utils
  CHARACTER(LEN=1024) :: configFile, ObsFile, ProfFile, StaticDir
  CHARACTER(LEN=10), ALLOCATABLE, DIMENSION(:) :: varList
  INTEGER(i_kind)     :: vLevel, numChans, numVars

  ! Local: for tests only. Delete later.
  INTEGER(i_kind) :: i, j, k, ivar, numObs, numObsTotal
  INTEGER(i_kind) :: n_insts, i_inst, i_chan, ilevel, ifile, level_pert
  INTEGER(i_kind) :: nchans

  ! Local variables
  REAL(r_kind) :: zenangle = 0.0, azangle = 0.0 , sunzenangle = 0.78539, sunazangle = 0.52359
  REAL(r_kind) :: temp_pert = 0.01D0, qv_pert = 0.0001D0, dqice
  INTEGER(i_kind) :: qice_pert = 33
  LOGICAL :: Thin_w_bkg = .True., Use_real_obs = .True.
  REAL(4) :: temp_scaling, qv_scaling

  REAL(r_kind) :: rhs, lhs, ad_check, tl_check_min, tl_check_max, tl_check_avg

  ! Get the configFile
  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", StaticDir)
  configFile = TRIM(StaticDir)//"/UnitTest/UnitTest_RTTOV.yaml"
  ProfFile = TRIM(StaticDir)//"/Satellite/prof.dat"
  PRINT *, 'Test_Sat_fy4_agri Config: ', TRIM(configFile)
  ! ============ Get var info ============
  !CALL namelist_read(TRIM(configFile), 'vLevel', vLevel)
  ifile = yaml_get_var(TRIM(configFile), 'modelState', 'vLevel', vLevel)
  PRINT *, 'vLevel is: ', vLevel
  !CALL namelist_read(TRIM(configFile), "varList", varList)
  ifile = yaml_get_var(TRIM(configFile), 'modelState', 'varList', varList)
  numVars = SIZE(varList, 1)
  PRINT *, 'numVars = ', numVars

  ! Initializer
  CALL mpddGlob%initialize()

  CALL geometry%initialize(configFile, mpddGlob)

  ifile = yaml_get_var(TRIM(configFile), 'FY4-AGRI', 'curr_sg', curr_sg)
  ifile = yaml_get_var(TRIM(configFile), 'FY4-AGRI', 'Thin_w_bkg', Thin_w_bkg)
  ifile = yaml_get_var(TRIM(configFile), 'FY4-AGRI', 'Use_real_obs', Use_real_obs)
  !ifile = yaml_get_var(TRIM(configFile), 'FY4-GIIRS', 'curr_sg', curr_sg)
  !ifile = yaml_get_var(TRIM(configFile), 'FY4-GIIRS', 'Thin_w_bkg', Thin_w_bkg)
  !ifile = yaml_get_var(TRIM(configFile), 'FY4-GIIRS', 'Use_real_obs', Use_real_obs)
  ifile = yaml_get_var(TRIM(configFile), 'RTTOV', 'temp_scaling', temp_scaling)
  ifile = yaml_get_var(TRIM(configFile), 'RTTOV', "qv_scaling", qv_scaling)
  ifile = yaml_get_var(TRIM(configFile), 'RTTOV', 'level_pert', level_pert)

  ! Initialize the Scaling implementation
  CALL Ctl2State%initialize(configFile)

  ASSOCIATE (sg => geometry%mg%sg(curr_sg))
    CALL X%initialize(configFile, sg)

    PRINT *, 'Read RTTOV profiles'
    ProfFile = TRIM(StaticDir)//"/Satellite/prof.dat"
    CALL RTTOV_utils%initialize(configFile, X)
    CALL RTTOV_utils%rttovProfile2Bkgd(TRIM(ProfFile), X)
    ! CALL Ctl2State%transBackward(X)

  END ASSOCIATE

  ASSOCIATE (sg => geometry%mg%sg(curr_sg))
  ! Run 3DVAR in each single grid.
  BLOCK
    TYPE(State_t) :: DD
    TYPE(MPObs_t), TARGET :: mpObs
    TYPE(ObsSet_t) :: Y1, YAngles
    TYPE(ObsSet_t) :: Y2
    TYPE(ObsSet_t) :: Y
    TYPE(C2O_t) :: H
    INTEGER(i_kind) :: iobs, it, ivlevel
    
    CALL mpObs%initializeMPObs(sg)    ! Initialize the observation parallel processing proc
    
    IF (.NOT. Use_real_obs) THEN
      ! CALL Set_MP_ObsFields(configFile, Y1, mpObs, 253.0D0, zenangle, azangle, sunzenangle, sunazangle, sg, 'tbb', 7)
      CALL Set_MP_ObsFields(configFile, Y1, mpObs, 240.0D0, zenangle, azangle, sunzenangle, sunazangle, sg, 'tbb', 7)
      CALL Y1%ObsFields(1)%Set_ObsType('fy4_1_agri')
      CALL Y1%ObsFields(2)%Set_ObsType('fy4_2_agri')
      CALL Y1%ObsFields(3)%Set_ObsType('fy4_3_agri')
      CALL Y1%ObsFields(4)%Set_ObsType('fy4_4_agri')
      CALL Y1%ObsFields(5)%Set_ObsType('fy4_5_agri')
      CALL Y1%ObsFields(6)%Set_ObsType('fy4_6_agri')
      CALL Y1%ObsFields(7)%Set_ObsType('fy4_7_agri')
    END IF
    PRINT *, 'FINISH ObsSuper_satellite'
    
    Y = Y1

    ! X%Fields(numvars)%data = ZERO
    CALL H%initialize(configFile, X, Y)

    ! For TL check
    ! delta_Hx/dHx ~1
    BLOCK
      TYPE(State_t) :: dX
      TYPE(ObsSet_t) :: Y_tl_check, tl_00, tl_01, tl_02, tl_03

      CALL dX%initialize(configFile, sg)
      CALL dX%setAllFieldData(0.0D0)
      dqice = qv_scaling * qv_pert * 0.000001
      dX%fields(dX%getVarIdx("qice"))%data(qice_pert, :, :) = dqice
      ! dX%fields(dX%getVarIdx("qcloud"))%data(level_pert, :, :) = qv_scaling*qv_pert  ! near 400hPa
      ! dX%fields(dX%getVarIdx("qvapor"))%data(level_pert, :, :) = qv_scaling*qv_pert  ! near 400hPa
      ! dX%fields(dX%getVarIdx("temp"))%data(level_pert, :, :) = temp_scaling*temp_pert  ! near 400hPa
      ! CALL Ctl2State%transBackward(dX)

      ! --- Attention: do NOT multiply dx. H'(x)*delta_x = H'(delta_X)!
      print *, 'here: ', dX%fields(dX%getVarIdx("qice"))%data(qice_pert, 1, 1), X%fields(X%getVarIdx("qice"))%data(qice_pert, 1, 1)

      tl_00 = H%fwdNL_opr(X + dX)
      tl_01 = H%fwdNL_opr(X - dX)
      tl_02 = H%fwdNL_opr(X + dX) - H%fwdNL_opr(X - dX)
      tl_03 = (H%fwdTL_opr(dX, X)) * 2.0D0

      print *, 'here 0: ', MAXVAL(tl_00%ObsFields(2)%values), MINVAL(tl_00%ObsFields(2)%values)
      print *, 'here 1: ', MAXVAL(tl_01%ObsFields(2)%values), MINVAL(tl_01%ObsFields(2)%values)
      print *, 'here 2: ', MAXVAL(tl_02%ObsFields(2)%values), MINVAL(tl_02%ObsFields(2)%values)
      print *, 'here 3: ', MAXVAL(tl_03%ObsFields(2)%values), MINVAL(tl_03%ObsFields(2)%values)
      Y_tl_check = ((H%fwdNL_opr(X + dX)) - (H%fwdNL_opr(X - dX)))/( (H%fwdTL_opr(dX, X)) * 2.0D0 )

      tl_check_max = MAXVAL(Y_tl_check%ObsFields(2)%values)
      tl_check_min = MINVAL(Y_tl_check%ObsFields(2)%values)
      tl_check_avg = SUM(Y_tl_check%ObsFields(2)%values)/REAL(SIZE(Y_tl_check%ObsFields(2)%values))
      
      IF (mpddGlob%isBaseProc()) PRINT *, 'tl_check = : ', tl_check_max, tl_check_min, tl_check_avg

    END BLOCK

    BLOCK
    ! For AD check
    ! <y, Lx> = <L*y, x>
    ! y = L(x), L - linear model
    ! L* -adjoint model
      TYPE(State_t) :: dX, dX_ad
      CALL dX%initialize(configFile, sg)
      CALL dX_ad%initialize(configFile, sg)
      CALL dX%setAllFieldData(0.0D0)
      CALL dX_ad%setAllFieldData(0.0D0)
      dqice = qv_scaling * qv_pert * 0.0001
      dX%fields(dX%getVarIdx("qice"))%data(qice_pert, :, :) = dqice
      dX%fields(dX%getVarIdx("qcloud"))%data(level_pert, :, :) = qv_scaling*qv_pert
      ! dX%fields(dX%getVarIdx("qvapor"))%data(level_pert, :, :) = qv_scaling*qv_pert  ! near 400hPa
      dX%fields(dX%getVarIdx("temp"))%data(level_pert, :, :) = temp_scaling*temp_pert  ! near 400hPa
      ! CALL Ctl2State%transBackward(dX)

      lhs = Y .DOT. (H%fwdTL_opr(dX, X))
      rhs = (H%adjMul_opr(Y, X)) .DOT. dX
      ad_check = lhs - rhs

      IF (mpddGlob%isBaseProc()) THEN
        PRINT *, "ad_check: ", ad_check
      END IF
      PRINT *, 'FINISH ObsSuper_satellite'

      Y = Y1

      ! X%Fields(numvars)%data = ZERO
      CALL H%initialize(configFile, X, Y)

      ! For TL check
      ! delta_Hx/dHx ~1
      BLOCK
        TYPE(State_t) :: dX
        TYPE(ObsSet_t) :: Y_tl_check, tl_01, tl_02, tl_03

        CALL dX%initialize(configFile, sg)
        CALL dX%setAllFieldData(0.0D0)
        dX%fields(dX%getVarIdx("qvapor"))%DATA(level_pert, :, :) = qv_scaling * qv_pert  ! near 400hPa
        ! dX%fields(dX%getVarIdx("temp"))%data = temp_scaling*temp_pert  ! near 400hPa
        CALL Ctl2State%transBackward(dX)

        ! --- Attention: do NOT multiply dx. H'(x)*delta_x = H'(delta_X)!
        Y_tl_check = ((H%fwdNL_opr(X + dX)) - (H%fwdNL_opr(X - dX))) / ((H%fwdTL_opr(dX, X)) * 2.0D0)

        tl_check_max = MAXVAL(Y_tl_check%ObsFields(1)%values)
        tl_check_min = MINVAL(Y_tl_check%ObsFields(1)%values)
        tl_check_avg = SUM(Y_tl_check%ObsFields(1)%values) / REAL(SIZE(Y_tl_check%ObsFields(1)%values))

        IF (mpddGlob%isBaseProc()) PRINT *, 'tl_check = : ', tl_check_max, tl_check_min, tl_check_avg

      END BLOCK

      BLOCK
        ! For AD check
        ! <y, Lx> = <L*y, x>
        ! y = L(x), L - linear model
        ! L* -adjoint model
        TYPE(State_t) :: dX, dX_ad
        CALL dX%initialize(configFile, sg)
        CALL dX_ad%initialize(configFile, sg)
        CALL dX%setAllFieldData(0.0D0)
        CALL dX_ad%setAllFieldData(0.0D0)
        dX%fields(dX%getVarIdx("qvapor"))%DATA(level_pert, :, :) = qv_scaling * qv_pert  ! near 400hPa
        ! dX%fields(dX%getVarIdx("temp"))%data = temp_scaling*temp_pert  ! near 400hPa
        CALL Ctl2State%transBackward(dX)

        lhs = Y.DOT. (H%fwdTL_opr(dX, X))
        rhs = (H%adjMul_opr(Y, X)) .DOT.dX
        ad_check = lhs - rhs

        IF (mpddGlob%isBaseProc()) THEN
          PRINT *, "ad_check: ", ad_check
        END IF

      END BLOCK

    END BLOCK
    END BLOCK
  END ASSOCIATE

  IF (mpddGlob%isBaseProc()) THEN
    IF ( ABS(tl_check_avg - 1.0D0) .LE. 1e-3 .AND. ABS(ad_check) .LE. 1e-7) THEN
      PRINT *, 'Test passed!'
      PRINT *, 'TL check accuracy: ', ABS(tl_check_avg - 1.0D0)
      PRINT *, 'AD check accuracy: ', ABS(ad_check)
    ELSE
      PRINT *, 'Test failed!'
    END IF
  END IF

! Finalize
  DEALLOCATE (varList)
  CALL mpddGlob%finalize

END PROGRAM TLAD_check_RTTOV
