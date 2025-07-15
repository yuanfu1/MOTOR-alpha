!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-QC.FY4_AGRI_App
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Yali Wu
! VERSION           : V 0.1
! HISTORY           : 2023-12-22, created by Yali Wu, for testing Non-Linear features of the RTTOV forward operator.
!
!   Created by Yali Wu (wuyali@gbamwf.com), 2023/12/22, @GBA-MWF, Shenzhen
!
!!--------------------------------------------------------------------------------------------------

!> @brief
!!
PROGRAM NL_features_RTTOV
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
  USE YAMLRead_m
  USE parameters_m
  USE Obs2State_m
  USE RTTOV_utils_m, ONLY: RTTOV_utils_t
  USE Ctl2State_m, ONLY: Ctl2State_t
  
  IMPLICIT NONE
  TYPE(mpddGlob_t), TARGET :: mpddGlob
  TYPE(geometry_t), TARGET :: geometry
  TYPE(Ctl2State_t) :: Ctl2State
  TYPE(State_t) :: XRes
  REAL(r_kind) :: t1, t2
  INTEGER(i_kind) :: curr_sg, iv
  LOGICAL :: IncFlag = .True.
  TYPE(MPObs_t), TARGET :: mpObs

  TYPE(State_t) :: X
  TYPE(RTTOV_utils_t) :: RTTOV_utils
  CHARACTER(LEN=1024) :: configFile, ObsFile, ProfFile, StaticDir
  INTEGER(i_kind)     :: numChans, numVars

  ! Local: for tests only. Delete later.
  INTEGER(i_kind) :: i, j, k, ivar, numObs, numObsTotal
  INTEGER(i_kind) :: n_insts, i_inst, i_chan, ilevel, ifile
  INTEGER(i_kind) :: nchans

  ! Local variables
  REAL(r_kind) :: zenangle = 0.0, azangle = 0.0 , sunzenangle = 0.78539, sunazangle = 0.52359
  REAL(r_kind) :: temp_pert = 0.0D0, qv_pert = 0.0D0, qcloud_pert = 0.0D0, qice_pert = 0.0D0

  ! Get the configFile
  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", StaticDir)
  configFile = TRIM(StaticDir)//"/UnitTest/UnitTest_RTTOV.yaml"
  ProfFile = TRIM(StaticDir)//"/Satellite/prof.dat"
  PRINT *, 'NL_features_RTTOV Config: ', TRIM(configFile)
  PRINT *, 'numVars = ', numVars

  ! Initializer
  CALL mpddGlob%initialize()
  
  CALL geometry%initialize(configFile, mpddGlob)

  ifile = yaml_get_var(TRIM(configFile), 'FY4-AGRI', 'curr_sg', curr_sg)

  ! Initialize the Scaling implementation
  CALL Ctl2State%initialize(configFile)

  ASSOCIATE (sg => geometry%mg%sg(curr_sg))
    CALL X%initialize(configFile, sg)

    PRINT *, 'Read RTTOV profiles'
    ProfFile = TRIM(StaticDir)//"/Satellite/prof.dat"
    CALL RTTOV_utils%initialize(configFile, X)
    CALL RTTOV_utils%rttovProfile2Bkgd(TRIM(ProfFile), X)
    ! CALL Ctl2State%transBackward(X)

    ! temporarily test
    ! X%fields(X%getVarIdx("qcloud"))%data = 0.0D0
    ! X%fields(X%getVarIdx("qice"))%data = 0.0D0

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
    INTEGER(i_kind) :: iobs, it
    
    CALL mpObs%initializeMPObs(sg)    ! Initialize the observation parallel processing proc
    
    CALL Set_MP_ObsFields(configFile, Y1, mpObs, 220.0D0, zenangle, azangle, sunzenangle, sunazangle, sg, 'tbb', 7)
    CALL Y1%ObsFields(1)%Set_ObsType('fy4_1_agri')
    CALL Y1%ObsFields(2)%Set_ObsType('fy4_2_agri')
    CALL Y1%ObsFields(3)%Set_ObsType('fy4_3_agri')
    CALL Y1%ObsFields(4)%Set_ObsType('fy4_4_agri')
    CALL Y1%ObsFields(5)%Set_ObsType('fy4_5_agri')
    CALL Y1%ObsFields(6)%Set_ObsType('fy4_6_agri')
    CALL Y1%ObsFields(7)%Set_ObsType('fy4_7_agri')
    
    Y = Y1

    ! X%Fields(numvars)%data = ZERO
    CALL H%initialize(configFile, X, Y)

    ! BLOCK
    !   USE RTTOV_diag_out_m
    !   TYPE(State_t) :: RR
    !   INTEGER :: iv

    !   PRINT *, '-------------RTTOV_diag_out START-------------'
    !   CALL write_diag_vars(configFile, X, 'fy4_1', 'agri', Y, RR)
    !   CALL Output_NC_State_AV(RR, "/Users/yaliwu/Downloads/allsky/", "fy4_1-agri_diag_clear", .true. , .true.)
    !   PRINT *, '-------------RTTOV_diag_out OVER-------------'
    ! END BLOCK

    BLOCK
      TYPE(State_t) :: dX, Xnew
      REAL :: max_v = 0.1, min_v = -0.1, interval_v = 0.01
      ! REAL :: max_v = 0.01, min_v = 0.001, interval_v = 0.001
      ! REAL :: max_v = 1.5, min_v = 0, interval_v = 0.1
      REAL, ALLOCATABLE :: ratio(:), NLvalue(:), tout(:)
      REAL, ALLOCATABLE :: TLvalue(:)
      TYPE(ObsSet_t) :: NL, TL, NL0, NL1, NL2
      INTEGER :: nvalues, ivalue, nlevels, ivalue_qc, ivalue_qi
      INTEGER :: level_temp = 20, level_qv = 27, level_qc = 17, level_qi = 22
      ! 30-10km; 41-16km; 40-15.8km

      nvalues = int((max_v - min_v)/interval_v) + 1
      nlevels = SIZE(X%fields(X%getVarIdx("qcloud"))%data,1)
      PRINT *, 'nvalues = ', nvalues
      ALLOCATE(NLvalue(nvalues))
      ALLOCATE(TLvalue(nvalues))
      ALLOCATE(tout(nvalues))
      ALLOCATE(ratio(nvalues))

      DO ivalue = 1, nvalues
        ratio(ivalue) = min_v + interval_v * (ivalue - 1)
      END DO

      ! DO ivalue = 10,35
      !   Xnew = X
      !   level_qc = nlevels - ivalue + 1
      !   ! Xnew%fields(Xnew%getVarIdx("qcloud"))%data(level_qc, :, :) = 0.1 * MAXVAL(X%fields(X%getVarIdx("qcloud"))%data(:, 1, 1)) + &
      !   !   X%fields(X%getVarIdx("qcloud"))%data(level_qc, :, :)
      !   Xnew%fields(Xnew%getVarIdx("qcloud"))%data(level_qc, :, :) = (1+1.5) * X%fields(X%getVarIdx("qcloud"))%data(level_qc, :, :)
      !   NL = H%fwdNL_opr(Xnew)
      !   NL1 = H%fwdNL_opr(X)
      !   PRINT *, ivalue, level_qc, NL1%ObsFields(2)%values(1), NL%ObsFields(2)%values(1)
      ! END DO

      dX = X
      DO ivalue = 1, nvalues  
        Xnew = X
        level_qc = 18
        level_qc = nlevels - level_qc + 1
        ! print*, Xnew%fields(Xnew%getVarIdx("temp"))%data(:, 1, 1)
        Xnew%fields(Xnew%getVarIdx("qcloud"))%data(level_qc, :, :) = (1+ratio(ivalue)) * X%fields(X%getVarIdx("qcloud"))%data(level_qc, :, :)
        ! Xnew%fields(Xnew%getVarIdx("qcloud"))%data(level_qc, :, :) = 0.1*ratio(ivalue) * MAXVAL(X%fields(X%getVarIdx("qcloud"))%data(:, 1, 1)) + &
        !                             X%fields(X%getVarIdx("qcloud"))%data(level_qc, :, :)
        NL = H%fwdNL_opr(Xnew)
        print *, 'NLvalues: ', ratio(ivalue), Xnew%fields(Xnew%getVarIdx("qcloud"))%data(level_qc, 1, 1), NL%ObsFields(2)%values(1)
        CALL dX%setAllFieldData(0.0D0)
        dX = Xnew - X
        TL = H%fwdTL_opr(dX, X)
        NL0 = H%fwdNL_opr(X) 
        NLvalue(ivalue) = NL%ObsFields(2)%values(1)
        TLvalue(ivalue) = TL%ObsFields(2)%values(1) + NL0%ObsFields(2)%values(1)
        tout(ivalue) = Xnew%fields(Xnew%getVarIdx("qcloud"))%data(level_qc, 1, 1)
      END DO
      OPEN (12, file = "./NL_features_QC.txt", form="formatted")
      ! WRITE(12, *) "tout at 30th level, qvout at 30th level, RTTOV forward value"
      DO ivalue = 1, nvalues
        WRITE(12, *) tout(ivalue), NLvalue(ivalue),TLvalue(ivalue)
      END DO
      CLOSE(12)
      OPEN (12, file = "./NL_features_QC.bin", form="unformatted", access='stream')
      WRITE(12) tout
      WRITE(12) NLvalue
      WRITE(12) TLvalue
      CLOSE(12)
      DEALLOCATE(NLvalue, TLvalue, tout, ratio)

    END BLOCK

  END BLOCK
  END ASSOCIATE

! Finalize
CALL mpddGlob%finalize

END PROGRAM NL_features_RTTOV
