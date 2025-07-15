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
    !   CALL Output_NC_State_AV(RR, "/Users/yaliwu/Downloads/allsky/", "fy4_1-agri_diag", .true. , .true.)
    !   PRINT *, '-------------RTTOV_diag_out OVER-------------'
    ! END BLOCK

    BLOCK
      TYPE(State_t) :: dX, Xnew
      ! REAL :: max_v = 0.2, min_v = -0.1, interval_v = 0.02
      ! REAL :: max_v2 = 10.0, min_v2 = -0.8, interval_v2 = 0.1
      ! REAL :: max_v2 = 0.5, min_v2 = -0.5, interval_v2 = 0.1
      REAL :: max_v = 1.5, min_v = 0, interval_v = 0.1
      REAL :: max_v2 = 2.2, min_v2 = -0.8, interval_v2 = 0.2
      REAL, ALLOCATABLE :: ratio1(:,:), ratio2(:,:), NLvalue(:,:), tout(:,:), qvout(:,:), qcout(:,:), qiout(:,:)
      REAL, ALLOCATABLE :: TLvalue(:,:)
      TYPE(ObsSet_t) :: NL, TL, NL0
      INTEGER :: nvalues1, nvalues2, ivalue, nlevels, ivalue_qc, ivalue_qi
      ! INTEGER :: level_temp = 20, level_qv = 27, level_qc = 33, level_qi = 33
      INTEGER :: level_temp = 20, level_qv = 20, level_qc = 22, level_qi = 22
      ! 30-10km; 41-16km; 40-15.8km

      nvalues1 = int((max_v - min_v)/interval_v) + 1
      nvalues2 = int((max_v2 - min_v2)/interval_v2) + 1
      nlevels = SIZE(X%fields(X%getVarIdx("temp"))%data,1)
      PRINT *, 'nvalues = ', nvalues1, nvalues2
      ALLOCATE(NLvalue(nvalues1, nvalues2))
      ALLOCATE(TLvalue(nvalues1, nvalues2))
      ALLOCATE(tout(nvalues1, nvalues2))
      ALLOCATE(qvout(nvalues1, nvalues2))
      ALLOCATE(qcout(nvalues1, nvalues2))
      ALLOCATE(qiout(nvalues1, nvalues2))
      ALLOCATE(ratio1(nvalues1, nvalues2))
      ALLOCATE(ratio2(nvalues1, nvalues2))

      ! DO ivalue_qc = 1, nvalues
      !   DO ivalue_qi = 1, nvalues
      !     ratio(ivalue_qc, ivalue_qi) = min_v + interval_v * (ivalue_qc - 1) * (ivalue_qi-1)
      !   END DO
      ! END DO
      DO ivalue = 1, nvalues1
        ratio1(ivalue,:) = min_v + interval_v * (ivalue - 1)
      END DO
      DO ivalue = 1, nvalues2
        ratio2(:,ivalue) = min_v2 + interval_v2 * (ivalue - 1)
      END DO

      print *, 'nchans = ', SIZE(Y%ObsFields,1)
      PRINT *, 'pres = ', X%sg%FGPres(:,1,1)
      PRINT *, 'qcloud = ', X%fields(X%getVarIdx("qcloud"))%data(:, 1,1)
      dX = X
      DO ivalue_qc = 1, nvalues1
        DO ivalue_qi = 1, nvalues2
          Xnew = X
          level_qc = 18
          level_qc = nlevels - level_qc + 1
          ! Xnew%fields(Xnew%getVarIdx("qcloud"))%data(level_qc, :, :) = 0.1*ratio1(ivalue_qc, ivalue_qi) * MAXVAL(X%fields(X%getVarIdx("qcloud"))%data(:, 1, 1)) + &
          !                           X%fields(X%getVarIdx("qcloud"))%data(level_qc, :, :)
          Xnew%fields(Xnew%getVarIdx("qcloud"))%data(level_qc, :, :) = (1+ratio1(ivalue_qc, ivalue_qi)) * X%fields(X%getVarIdx("qcloud"))%data(level_qc, :, :)
          level_qi = 22
          level_qi = nlevels - level_qi + 1
          Xnew%fields(Xnew%getVarIdx("qice"))%data(level_qi, :, :) = (1+ratio2(ivalue_qc, ivalue_qi)) * X%fields(X%getVarIdx("qice"))%data(level_qi, :, :)
          ! PRINT *, Xnew%fields(Xnew%getVarIdx("qcloud"))%data(level_qi, 1, 1), Xnew%fields(Xnew%getVarIdx("qice"))%data(level_qi, 1, 1)
          NL = H%fwdNL_opr(Xnew)
          ! print *, 'NLvalues: ', NL%ObsFields(2)%values(1)
          CALL dX%setAllFieldData(0.0D0)
          dX = Xnew - X
          TL = H%fwdTL_opr(dX, X)
          NL0 = H%fwdNL_opr(X) 
          NLvalue(ivalue_qc, ivalue_qi) = NL%ObsFields(2)%values(1)
          TLvalue(ivalue_qc, ivalue_qi) = TL%ObsFields(2)%values(1) + NL0%ObsFields(2)%values(1)
          qcout(ivalue_qc, ivalue_qi) = Xnew%fields(Xnew%getVarIdx("qcloud"))%data(level_qc, 1, 1)
          qiout(ivalue_qc, ivalue_qi) = Xnew%fields(Xnew%getVarIdx("qice"))%data(level_qi, 1, 1)
        END DO
      END DO

      OPEN (12, file = "./NL_features_QCQI.txt", form="formatted")
      ! WRITE(12, *) "tout at 30th level, qvout at 30th level, RTTOV forward value"
      DO ivalue_qc = 1, nvalues1
        DO ivalue_qi = 1, nvalues2
          WRITE(12, *) qcout(ivalue_qc, ivalue_qi), qiout(ivalue_qc, ivalue_qi), NLvalue(ivalue_qc, ivalue_qi), TLvalue(ivalue_qc, ivalue_qi)
        END DO
      END DO
      CLOSE(12)
      OPEN (12, file = "./NL_features_QCQI.bin", form="unformatted", access='stream')
      WRITE(12) qcout
      WRITE(12) qiout
      WRITE(12) NLvalue
      WRITE(12) TLvalue
      CLOSE(12)
      DEALLOCATE(NLvalue, TLvalue, tout, qvout, qcout, qiout, ratio1, ratio2)

    END BLOCK

  END BLOCK
  END ASSOCIATE

! Finalize
CALL mpddGlob%finalize

END PROGRAM NL_features_RTTOV
