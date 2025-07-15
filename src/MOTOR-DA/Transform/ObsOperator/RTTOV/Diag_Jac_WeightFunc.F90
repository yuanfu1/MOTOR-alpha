!!--------------------------------------------------------------------------------------------------
! PROJECT           : Diag_Jac_WeightFunc
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Yali Wu
! VERSION           : V 0.1
! HISTORY           : 2024-1-11, created by Yali Wu, for testing Non-Linear features of the RTTOV forward operator.
!
!   Created by Yali Wu (wuyali@gbamwf.com), 2023/12/22, @GBA-MWF, Shenzhen
!
!!--------------------------------------------------------------------------------------------------

!> @brief
!!
PROGRAM Diag_Jac_WeightFunc
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
  USE Satellite_utils_m
  
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

  ! Get the configFile
  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", StaticDir)
  configFile = TRIM(StaticDir)//"/UnitTest/UnitTest_RTTOV.yaml"
  ProfFile = TRIM(StaticDir)//"/Satellite/prof.dat"
  PRINT *, 'Diag_Jac_WeightFunc Config: ', TRIM(configFile)
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
    ! CALL RTTOV_utils%rttovProfile2Bkgd(TRIM(ProfFile), X)
    CALL RTTOV_utils%USStandardProfile2Bkgd(X)
    ! CALL Ctl2State%transBackward(X)

    ! ! temporarily test
    ! X%fields(X%getVarIdx("qcloud"))%data = 0.0D0
    ! X%fields(X%getVarIdx("qice"))%data = 0.0D0

  END ASSOCIATE

  ASSOCIATE (sg => geometry%mg%sg(curr_sg))
  ! Run 3DVAR in each single grid.
  BLOCK
    USE RTTOV_diag_out_m
    TYPE(State_t) :: DD
    TYPE(MPObs_t), TARGET :: mpObs
    TYPE(ObsSet_t) :: Y1, YAngles
    TYPE(ObsSet_t) :: Y2
    TYPE(ObsSet_t) :: Y
    TYPE(C2O_t) :: H
    INTEGER(i_kind) :: iobs, it
    INTEGER(i_kind) :: nchans, ichan
    INTEGER(i_kind), ALLOCATABLE :: chan_lists(:), rttov_chan_lists(:)
    ! CHARACTER(len=10) :: platform = 'fy3_4', inst='mwts2', charvalue
    CHARACTER(len=10) :: platform = 'fy3_5', inst='mwts3', charvalue
    ! CHARACTER(len=10) :: platform = 'fy4_1', inst='agri', charvalue
    TYPE(State_t) :: RR
    INTEGER :: iv
    
    CALL mpObs%initializeMPObs(sg)    ! Initialize the observation parallel processing proc
    
    CALL Get_rttov_chan_info(trim(platform)//'-'//trim(inst), nchans)
    ALLOCATE (chan_lists(nchans), rttov_chan_lists(nchans))
    CALL Get_rttov_chan_info(trim(platform)//'-'//trim(inst), nchans, chan_lists, rttov_chan_lists)
    
    CALL Set_MP_ObsFields(configFile, Y1, mpObs, 220.0D0, zenangle, azangle, sunzenangle, sunazangle, sg, 'tbb', nchans)
    DO ichan = 1, nchans
      IF (ichan < 10) WRITE(charvalue, '(I1)') ichan
      IF (ichan > 10 .AND. ichan < 100) WRITE(charvalue, '(I2)') ichan
      print *, 'charvalue = ', charvalue
      CALL Y1%ObsFields(ichan)%Set_ObsType(TRIM(platform)//'_'//TRIM(charvalue)//'_'//TRIM(inst))
    END DO
    DEALLOCATE(chan_lists, rttov_chan_lists)   

    Y = Y1

    ! X%Fields(numvars)%data = ZERO
    CALL H%initialize(configFile, X, Y)

    PRINT *, '-------------RTTOV_diag_out START-------------'
    CALL write_diag_vars(configFile, X, TRIM(platform), TRIM(inst), Y, RR)
    CALL Output_NC_State_AV(RR, "./", TRIM(platform)//'-'//TRIM(inst)//"_diag_clear", .true. , .true.)
    PRINT *, '-------------RTTOV_diag_out OVER-------------'

  END BLOCK
  END ASSOCIATE

! Finalize
CALL mpddGlob%finalize

END PROGRAM Diag_Jac_WeightFunc
