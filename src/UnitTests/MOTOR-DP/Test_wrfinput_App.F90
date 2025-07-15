!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-QC.FY4_AGRI_App
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Yali Wu
! VERSION           : V 0.1
! HISTORY           : 2022-01-30, created by Yali Wu, for FY4 AGRI DA with WRF.
!
!   Created by Yali Wu (wuyali@gbamwf.com), 2022/01/30, @GBA-MWF, Shenzhen
!   Modified by Yuanfu Xie, 2022/05/16 replacing RTTOV default background using
!   RTTOV_utils.F90 procedure, rttovProfile2Bkgd
!
!!--------------------------------------------------------------------------------------------------

!> @brief
!!
PROGRAM Test_wrfinput_App
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
  USE YAMLRead_m
  USE parameters_m
  USE Obs2State_m
  USE ObsUtilities_m
  USE ObsAttr_m, ONLY: ObsAttrSAT_t
  USE RTTOV_diag_out_m
  USE RTTOV_utils_m, ONLY: RTTOV_utils_t
  USE IOGrapes_m, ONLY: IOGrapes_t
  USE Ctl2State_m, ONLY: Ctl2State_t

  IMPLICIT NONE
  INTEGER(i_kind), PARAMETER :: tlevel01 = 51, tlevel02 = 6
  TYPE(mpddGlob_t), TARGET :: mpddGlob
  TYPE(geometry_t), TARGET :: geometry
  TYPE(State_t) :: xANA
  TYPE(State_t), ALLOCATABLE  :: XbMG(:)
  TYPE(ObsSet_t), ALLOCATABLE :: HxMG(:)
  REAL(r_kind) :: t1, t2
  INTEGER(i_kind) :: mgStart, mgEnd, iv
  LOGICAL :: IncFlag = .TRUE.
  CLASS(MiniSolver_t), ALLOCATABLE :: miniSolver
  TYPE(MPObs_t), TARGET :: mpObs

  TYPE(State_t) :: X
  TYPE(RTTOV_utils_t) :: RTTOV_utils
  TYPE(Ctl2State_t) :: Ctl2State
  CHARACTER(LEN=1024) :: configFile, ObsFile, ProfFile, StaticDir
  CHARACTER(LEN=10), ALLOCATABLE, DIMENSION(:) :: varList
  INTEGER(i_kind)     :: vLevel, numChans, numVars
  CHARACTER(LEN=70) :: outputDir, varName, tmpstring
  REAL(r_kind), ALLOCATABLE :: theta_wrf(:), qvapor_wrf(:), pres_wrf(:), temp_wrf(:), rho_wrf(:)
  CHARACTER(LEN=20) :: rhov_scale_scheme = 'DomainAvg'
  CHARACTER(LEN=20) :: temp_scale_scheme = 'DomainAvg'
  CHARACTER(LEN=20) :: qvapor_scale_scheme = 'DomainAvg'

  INTEGER(i_kind) :: i, j, k, ivar, numObs, numObsTotal
  INTEGER(i_kind) :: n_insts, i_inst, i_chan, ilevel
  INTEGER(i_kind) :: nchans, istatus
  TYPE(ObsAttrSAT_t) :: ObsAttrSAT_AGRI
  CHARACTER(len=10) :: minimizer

  ! Local variables
  REAL(r_kind) :: zenangle = ZERO, azangle = ZERO, sunzenangle = 0.78539, sunazangle = 0.52359
  LOGICAL :: Thin_w_bkg = .TRUE., Use_real_obs = .TRUE., Use_real_bkg = .TRUE., write_diag = .FALSE.
  CHARACTER(LEN=128) :: tempchar, test_name

  ! Initializer
  CALL mpddGlob%initialize()

  PRINT *, 'TEST START'
  CALL CPU_TIME(t1)

  ALLOCATE (theta_wrf(57))
  ALLOCATE (qvapor_wrf(57))
  ALLOCATE (pres_wrf(57))
  ALLOCATE (temp_wrf(57))
  ALLOCATE (rho_wrf(57))
  ! Read WRF profiles
  OPEN (unit=11, file="/MOTOR/MOTOR-TEST-DATA/20210504/CumUniTest/theta.txt", status="old")
  OPEN (unit=12, file="/MOTOR/MOTOR-TEST-DATA/20210504/CumUniTest/qvapor.txt", status="old")
  OPEN (unit=13, file="/MOTOR/MOTOR-TEST-DATA/20210504/CumUniTest/pres.txt", status="old")

  DO i = 1, 57
    READ (11, *) theta_wrf(i)
    READ (12, *) qvapor_wrf(i)
    READ (13, *) pres_wrf(i)
  END DO
  CLOSE (11)
  CLOSE (12)
  CLOSE (13)
  temp_wrf = theta_wrf / (1.0D05 / pres_wrf)**k_d
  rho_wrf = pres_wrf / (r_d * temp_wrf * (1.0D0 + qvapor_wrf))
  PRINT *, 'theta_wrf = ', theta_wrf
  PRINT *, ''
  PRINT *, 'qvapor_wrf = ', qvapor_wrf
  PRINT *, ''
  PRINT *, 'pres_wrf = ', pres_wrf
  PRINT *, ''
  PRINT *, 'temp_wrf = ', temp_wrf
  PRINT *, ''
  PRINT *, 'rho_wrf = ', rho_wrf
  PRINT *, ''
  CALL mpddGlob%barrier
  ! End of WRF profiles

  ! Get the configFile
  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", StaticDir)

  configFile = TRIM(StaticDir)//"/UnitTest/Test_wrfinput_App.yaml"
  PRINT *, 'Test_Sat_fy4_agri Config: ', TRIM(configFile)
  ! ============ Get var info ============
  !CALL namelist_read(TRIM(configFile), 'vLevel', vLevel)
  istatus = yaml_get_var(TRIM(configFile), 'modelState', 'vLevel', vLevel)
  PRINT *, 'vLevel is: ', vLevel
  !CALL namelist_read(TRIM(configFile), "varList", varList)
  istatus = yaml_get_var(TRIM(configFile), 'modelState', 'varList', varList)
  numVars = SIZE(varList, 1)
  PRINT *, 'numVars = ', numVars

  CALL geometry%initialize(configFile, mpddGlob)
  ALLOCATE (XbMG(geometry%mg%mg_coarsest:geometry%mg%mg_finest))
  ALLOCATE (HxMG(geometry%mg%mg_coarsest:geometry%mg%mg_finest))

  IF (yaml_get_var(configFile, 'IO', 'output_dir', outputDir) /= 0) STOP
  istatus = yaml_get_var(TRIM(configFile), 'geometry', 'mgStart', mgStart)
  istatus = yaml_get_var(TRIM(configFile), 'geometry', 'mgEnd', mgEnd)
  istatus = yaml_get_var(TRIM(configFile), 'FY4-AGRI', 'Use_real_obs', Use_real_obs)
  istatus = yaml_get_var(TRIM(configFile), 'FY4-AGRI', 'Use_real_bkg', Use_real_bkg)
  istatus = yaml_get_var(TRIM(configFile), 'FY4-AGRI', 'test_name', test_name)
  istatus = yaml_get_var(TRIM(configFile), 'FY4-AGRI', 'write_diag', write_diag)
  istatus = yaml_get_var(TRIM(configFile), 'FY4-AGRI', 'minimizer', minimizer)
  istatus = yaml_get_var(TRIM(configFile), 'CV_Transform', 'rhov_scale_scheme', rhov_scale_scheme)
  istatus = yaml_get_var(TRIM(configFile), 'CV_Transform', 'temp_scale_scheme', temp_scale_scheme)
  istatus = yaml_get_var(TRIM(configFile), 'CV_Transform', 'temp_scale_scheme', qvapor_scale_scheme)

  ! Initialize the Scaling implementation
  CALL Ctl2State%initialize(configFile)

  ALLOCATE (SolverFRCG_t::miniSolver)
  SELECT TYPE (miniSolver)
  TYPE IS (SolverFRCG_t)
    CALL miniSolver%initialize(configFile)
  TYPE IS (SolverLBFGS_t)
    CALL miniSolver%initialize(configFile)
  END SELECT

  DO i = mgEnd, mgStart, -1
    PRINT *, 'IO at sg ', i
    IF (i .EQ. mgEnd) THEN
      ! Initialize a zeros background field at finest grid.
      ASSOCIATE (sg => geometry%mg%sg(mgEnd))

        CALL XbMG(mgEnd)%initialize(configFile, sg)

        IF (use_real_bkg) THEN
          CALL ioGrapes%initialize(configFile, geometry)
          CALL ioGrapes%m_read_bcg_into_Xm(XbMG(mgEnd), sg)
        ELSE
          BLOCK
            PRINT *, 'Read RTTOV profiles'
            ProfFile = TRIM(StaticDir)//"/Satellite/prof.dat"
            CALL RTTOV_utils%initialize(configFile, X)
            CALL RTTOV_utils%rttovProfile2Bkgd(TRIM(ProfFile), XbMG(mgEnd))
          END BLOCK
        END IF

        ! ! Read observations at the finest grid
        CALL OBS_AGRI%ObsInitial(configFile)
        CALL OBS_AGRI%ObsIngest(XbMG(mgEnd))

        numObs_AGRI = OBS_AGRI%numObs
        PRINT *, 'numObs_AGRI = ', numObs_AGRI
        CALL mpddGlob%barrier

      END ASSOCIATE

    ELSE
      ! Restrict to each coarser grid
      ASSOCIATE (sgFiner => geometry%mg%sg(i + 1), sgCoarser => geometry%mg%sg(i))
        CALL XbMG(i)%initialize(configFile, sgCoarser)
        CALL restrictionMG(XbMG(i), XbMG(i + 1), geometry%mg)
        CALL geometry%mg%restrictionOfStatics(sgFiner, sgCoarser)
        CALL update_restrictionOfStatics(configFile, XbMG(i), sgCoarser)
      END ASSOCIATE
    END IF

  END DO

  DO i = mgEnd, mgStart, -1
    ! Add the Variable replacement
    CALL Ctl2State%transBackward(XbMG(i))
  END DO

  ! ! Give the initial value at the coarsest grid.
  xANA = XbMG(mgStart)
  IF (xANA%sg%isActiveProc()) THEN
    PRINT *, 'check outside pres: mgStart ', MAXVAL(XbMG(mgStart)%sg%FGPres), MINVAL(XbMG(mgStart)%sg%FGPres)
    PRINT *, 'check outside pres: mgEnd ', MAXVAL(XbMG(mgEnd)%sg%FGPres), MINVAL(XbMG(mgEnd)%sg%FGPres)
  END IF

  ! MultiGrid
  DO i = mgStart, mgEnd
    PRINT *, 'Iteration at sg ', i
    ! Read observations at the finest grid

    ASSOCIATE (sg => geometry%mg%sg(i))
      ! Run 3DVAR in each single grid.
      BLOCK
        TYPE(State_t) :: X
        TYPE(BMatrix_t) :: B
        TYPE(MPObs_t), TARGET :: mpObs
        TYPE(ObsSet_t) :: Y1, YAngles
        TYPE(ObsSet_t) :: Y, Y_bkg
        TYPE(C2O_t) :: H
        TYPE(RMatrix_t) :: R
        TYPE(JFunc_t) :: JFunc
        INTEGER(i_kind) :: iobs, it, ivlevel, nvars
        CHARACTER(len=2) :: sgstring

        CALL mpObs%initializeMPObs(sg)    ! Initialize the observation parallel processing proc

        ! Initialize X, Y, B, H
        X = xANA
        XbRef = Ctl2State%fwdNL_opr(XbMG(i))

        CALL B%initialize(configFile, sg)     ! Initialize the B matrix

        ! Need to loop for instruments when assimilating multiple satellites together

        IF (.NOT. Use_real_obs) THEN
          CALL Set_MP_ObsFields(configFile, Y1, mpObs, 2.0D0, zenangle, azangle, sunzenangle, sunazangle, sg, 'tbb', 7)
          CALL Y1%ObsFields(1)%Set_ObsType('fy4_1_agri')
          CALL Y1%ObsFields(2)%Set_ObsType('fy4_1_agri')
          CALL Y1%ObsFields(3)%Set_ObsType('fy4_1_agri')
          CALL Y1%ObsFields(4)%Set_ObsType('fy4_1_agri')
          CALL Y1%ObsFields(5)%Set_ObsType('fy4_1_agri')
          CALL Y1%ObsFields(6)%Set_ObsType('fy4_1_agri')
          CALL Y1%ObsFields(7)%Set_ObsType('fy4_1_agri')
        ELSE

          IF (sg%isActiveProc()) THEN

            CALL OBS_AGRI%ObsPrepareForSg(XbRef)

            PRINT *, 'START ObsSuper_satellite for Angles'
            CALL OBS_AGRI%ObsThinning(XbRef, Y1, mpObs, Thin_w_bkg, .FALSE.)
            PRINT *, 'thinning for TBB finished'

          END IF

        END IF
        PRINT *, 'FINISH ObsSuper_satellite'
        CALL mpddGlob%barrier
        ! Concat thinning obs:
        ! CALL ObsConcat_s((/Y1, Y2/), Y)
        Y = Y1

        IF (sg%isActiveProc()) THEN
          CALL R%initialize(configFile, Y, X%sg) ! Initialize R
        END IF

        ! X%Fields(numvars)%data = 0.0D0
        CALL H%initialize(configFile, X, Y)

        ! ! Initialize J Function
        PRINT *, 'Initialize J Function'
        CALL JFunc%initialize(configFile, X, Y, H, B, R, sg)

        ! ! Run minimization
        IF (sg%isActiveProc()) THEN
          PRINT *, 'START minimization'
          CALL miniSolver%run(X, JFunc, sg)
          PRINT *, 'FINISH minimization'
        END IF

        ! Connect to coarser grid.
        IF (i .NE. mgEnd) THEN
          ! ===========================> prolongate to finer grid
          ASSOCIATE (sgFiner => geometry%mg%sg(i + 1))
            xANA = State_t(configFile, sgFiner)
            ! CALL prolongationMG(X, xANA, XbMG(i + 1), geometry%mg, IncFlag)
            CALL prolongationMGInc(X, xANA, XbMG(i), XbMG(i + 1), sg, sgFiner, geometry%mg)
          END ASSOCIATE
        ELSE
          xANA = X ! Return the state fields directly

        END IF

        PRINT *, 'WRITE OUT model space bak, ana'
        CALL Output_NC_State_AV(X, TRIM(outputDir), TRIM(test_name)//"_AGRI_ana", .TRUE., .TRUE.)
        CALL Output_NC_State_AV(XbMG(i), TRIM(outputDir), TRIM(test_name)//"_AGRI_bak", .TRUE., .TRUE.)

      END BLOCK

    END ASSOCIATE

  END DO

  ! CALL OBS_AGRI%ObsDestroy()

  ! Scaling back
  XbRef = XbMG(mgEnd)
  CALL Ctl2State%transFwdNonLinear(xANA)
  CALL Ctl2State%transFwdNonLinear(XbRef)

  CALL ioGrapes%m_write_Xm_into_bcg(xANA, geometry%mg%sg(mgEnd), XbMG(mgEnd))

  ! CALL mpddGlob%barrier
  CALL CPU_TIME(t2)
  PRINT *, 'Time cost:', t2 - t1, 's'

  ! Finalize
  DEALLOCATE (miniSolver)
  DEALLOCATE (XbMG)
  DEALLOCATE (HxMG)
  DEALLOCATE (varList)

  ! CALL mpddGlob%finalize

  IF (t2 - t1 > 1) THEN
    PRINT *, "Test passed!"
  ELSE
    PRINT *, "Test failed!"
  END IF

  CALL mpddGlob%finalize

END PROGRAM Test_wrfinput_App
