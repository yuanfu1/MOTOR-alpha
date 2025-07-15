!!--------------------------------------------------------------------------------------------------
! PROJECT           : Application
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           : 2022-03-19   Created by Yuanfu Xie
!
!   Created by Yuanfu Xie (xieyf@gbamwf.com), 2022/03/19, @GBA-MWF, Shenzhen
!   Modified by Yuanfu Xie (yuanfu_xie@yahoo.com), 2022/06/03, for providing more detailed debugging
!     information on ideal test function. It can mark the analytic values at a given region for
!     checking the observations in the replace routine.
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module provides an initialization of an analysis, including namelist, grid setting etc.
! @warning abc
MODULE Applications_HybInc_m
  USE SolverLBFGS_m, ONLY: SolverLBFGS_t
  USE SolverFRCG_HybInc_m, ONLY: SolverFRCG_HybInc_t
  USE MiniSolver_m, ONLY: MiniSolver_t
  USE geometry_m, ONLY: geometry_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE MPObs_m, ONLY: MPObs_t
  USE State_m, ONLY: State_t
  USE ObsSet_m, ONLY: ObsSet_t
  USE ObsField_m, ONLY: ObsField_t
  USE C2O_m, ONLY: C2O_t
  USE BMatrix_m, ONLY: BMatrix_t
  USE RMatrix_m, ONLY: RMatrix_t
  USE JFunc_HybInc_m, ONLY: JFunc_HybInc_t
  USE State2NC_m
  USE Mock_m
  USE MGOpts_m
  USE parameters_m, ONLY: degree2radian
  USE YAMLRead_m
  USE UV2W_m, ONLY: UV2W_t
  USE Ctl2State_m, ONLY: Ctl2State_t

  USE IOWRF_m, ONLY: IOWRF_t
  USE IOGrapes_m, ONLY: IOGrapes_t
  USE IOGrapesModelvar_m, ONLY: IOGrapesModelvar_t
  USE IOGrapesPostvar_m, ONLY: IOGrapesPostvar_t
  USE IOModel_m, ONLY: IOModel_t

  USE EnLoc_m, ONLY: EnLoc_t

  USE ObsSurface_m, ONLY: ObsSurface_t
  USE ObsSound_m, ONLY: ObsSound_t
  USE ObsVwpw_m, ONLY: ObsVwpw_t
  USE ObsRadarRef_m, Only: ObsRadarRef_t
  USE ObsRadarVel_m, Only: ObsRadarVel_t
  USE ObsSatellite_m, ONLY: ObsSatellite_t
  USE Satob_m, ONLY: Satob_t
  USE ObsGWST_m, ONLY: ObsGWST_t
  USE ObsLBUOY_m, ONLY: ObsLBUOY_t
  USE ObsSING_m, ONLY: ObsSING_t

  USE ObsUtilities_m

  USE ObsBase_m
  USE Obs2State_m

  ! added bt TS@230131
  USE obsTools_m

  ! Define the Application type:
  TYPE Applications_HybInc_t
    CHARACTER(LEN=1024) :: configFile, ncOutputFile
    CHARACTER(LEN=20), ALLOCATABLE :: varList(:)
    CHARACTER(LEN=20) :: BECsolver
    CHARACTER(LEN=20) :: mode = 'Debug' !< mode: Debug, UnitTest, Alpha, Beta, Release
    CHARACTER(LEN=20) :: task
    LOGICAL           :: Hybrid, HybInc, UpdateBkgd ! Yuanfu Xie added UpdateBkgd parameter controlling the MG scheme
    ! INTEGER(i_kind) :: modeNO = 0
    INTEGER(i_kind)   :: numVar, ensNum
    INTEGER(i_kind)   :: mgStart, mgEnd
    INTEGER(i_kind), ALLOCATABLE :: mgGridSmooth(:, :)
    REAL(r_kind), ALLOCATABLE :: verifyLatlon(:)
    TYPE(mpddGlob_t) :: mpddGlob
    TYPE(geometry_t) :: geometry
    TYPE(State_t), ALLOCATABLE :: XbMG(:), XbRef(:)
    TYPE(State_t), ALLOCATABLE :: XbMGEns(:, :)
    TYPE(State_t), ALLOCATABLE :: Xb_enperts(:, :)
    CHARACTER(len=50), ALLOCATABLE :: platform_name(:), inst_name(:)
    LOGICAL, ALLOCATABLE :: turnOn(:)
    CHARACTER(len=20) :: framework = 'FullState'
    TYPE(State_t) :: analysis

    TYPE(ObsSurface_t)  :: sfc ! Surfaces
    TYPE(ObsSound_t)    :: snd ! Sounding
    TYPE(ObsVwpw_t)     :: pro ! Profiler
    TYPE(ObsRadarRef_t) :: ref ! Radar reflectivity
    TYPE(ObsRadarVel_t) :: vel ! Radar radial wind
    TYPE(ObsSatellite_t), ALLOCATABLE :: ObsSatellites(:)
    TYPE(Satob_t)       :: cdw ! cloud derived wind

    TYPE(ObsGWST_t)     :: shp ! ship data
    TYPE(ObsLBUOY_t)    :: boy ! buoy data
    TYPE(ObsSING_t)     :: air ! aircraft data

    TYPE(SolverFRCG_HybInc_t), ALLOCATABLE :: miniSolver
    CLASS(IOModel_t), ALLOCATABLE :: IOModel

    TYPE(ObsSet_t), ALLOCATABLE :: thinnedObs(:)  ! Thinned observation for each multigrid level
    TYPE(MPObs_t), ALLOCATABLE :: mpObs(:)

    TYPE(Ctl2State_t) :: Ctl2State

    ! Background field fill-in rampRanges:
    INTEGER(i_kind), ALLOCATABLE :: iRange_ramp(:, :)  ! space and time dimensions : variables

  CONTAINS
    PROCEDURE :: initial
    PROCEDURE :: backgrd
    PROCEDURE :: backgrdEns
    PROCEDURE :: calEnsB
    PROCEDURE :: readObs
    PROCEDURE :: analyss
    PROCEDURE :: analyze
    PROCEDURE :: verifyA
    PROCEDURE :: destroy
    PROCEDURE :: DASolve
    PROCEDURE :: varAnal
    PROCEDURE :: replace
    PROCEDURE :: dumpAna

    PROCEDURE, PRIVATE :: bkgdFillin
  END TYPE Applications_HybInc_t

CONTAINS

  !> @brief
  !================================================================
  !  This routine initializes an application:
  !  Input:
  !       namelist: character string of namelist filename, assuming
  !                 it is under static directory
  !       min:      integer option of minimization method,
  !                 1: LBFGS; 2: FRCG
  !  Output:
  !       istatus:  integer reporting if initialization successes
  !                 0: success; 1: undefined minimization method
  !
  SUBROUTINE initial(this, configFile, min, istatus)
    IMPLICIT NONE
    CLASS(Applications_HybInc_t) :: this
    CHARACTER(LEN=*), INTENT(IN) :: configFile
    INTEGER(i_kind), INTENT(IN) :: min
    INTEGER(i_kind), INTENT(OUT) :: istatus

    ! Local variables:
    INTEGER(i_kind) :: i
    REAL(r_kind) :: t1, t2
    INTEGER(i_kind) :: ifile, n_plats, n_insts, i_inst

    CALL CPU_TIME(t1)

    istatus = 0

    ! Get the configFile
    this%configFile = TRIM(configFile)

    ! Read in parameters:
    WRITE (*, 30)
30  FORMAT('Application - initial: reading yaml info...')

    istatus = yaml_get_var(this%configFile, 'IO', 'output_dir', this%ncOutputFile)
    IF (istatus .NE. 0) THEN
      WRITE (*, 157)
157   FORMAT('App initial: output directory is not provided, please check your yaml and rerun!')
      STOP
    END IF
    WRITE (*, 158) TRIM(this%ncOutputFile)
158 FORMAT('App initial: output file directory: ', A)
    this%ncOutputFile = TRIM(this%ncOutputFile)

    istatus = yaml_get_var(TRIM(this%configFile), 'modelState', 'varList', this%varList)
    IF (istatus .NE. 0) THEN
      this%numVar = 0
    ELSE
      this%numVar = UBOUND(this%varList, 1)
      WRITE (*, 11) this%numVar
11    FORMAT("Number of analysis vars: ", I2)
    END IF

    istatus = yaml_get_var(TRIM(this%configFile), 'geometry', 'mgStart', this%mgStart)
    IF (istatus .NE. 0) THEN
      WRITE (*, 51)
51    FORMAT('App initial: multigrid start level is required, please your yaml and provide it and rerun')
      STOP
    END IF
    istatus = yaml_get_var(TRIM(this%configFile), 'geometry', 'mgEnd', this%mgEnd)
    IF (istatus .NE. 0) THEN
      WRITE (*, 52)
52    FORMAT('App initial: multigrid end level is required, please your yaml and provide it and rerun')
      STOP
    END IF
    IF (this%mpddGlob%isBaseProc()) WRITE (*, 13) this%mgStart, this%mgEnd
13  FORMAT("mgStart and mgEnd: ", 2I2)

    istatus = yaml_get_var(TRIM(this%configFile), 'RunMode', 'Mode', this%mode)
    IF (istatus .NE. 0) THEN
      WRITE (*, 53)
53    FORMAT('App initial: run mode is required: please check your yaml setting RunMode/Mode and rerun')
      STOP
    END IF
    PRINT *, 'RunMode is: ', this%mode, istatus

    istatus = yaml_get_var(TRIM(this%configFile), 'RunMode', 'Framework', this%framework)
    IF (istatus .NE. 0) THEN
      WRITE (*, 29)
29    FORMAT('App initial: RunMode Framework is required: please check your yaml setting RunMode/Framework and rerun')
      STOP
    END IF
    PRINT *, 'RunMode Framework is: ', this%framework, istatus
    ! IF (TRIM(this%mode) == 'Debug') THEN
    !   this%modeNO = 0
    ! ELSE IF (TRIM(this%mode) == 'UnitTest') THEN
    !   this%modeNO = 0
    ! ELSE IF (TRIM(this%mode) == 'Alpha') THEN
    !   this%modeNO = 1
    ! ELSE IF (TRIM(this%mode) == 'Beta') THEN
    !   this%modeNO = 2
    ! ELSE IF (TRIM(this%mode) == 'Release') THEN
    !   this%modeNO = 3
    ! END IF

    IF (ALLOCATED(this%iRange_ramp)) DEALLOCATE (this%iRange_ramp)
    ALLOCATE (this%iRange_ramp(3, this%numVar))

    istatus = yaml_get_var(TRIM(this%configFile), 'RunMode', 'Task', this%Task)
    IF (istatus .NE. 0) THEN
      WRITE (*, 54)
54    FORMAT('App initial: Task is a required parameter. Check your yaml for RunMode/Task and rerun')
      STOP
    END IF
    PRINT *, 'Task is: ', this%Task, istatus

    ! Prepare for the loop of different satellite instruments
    ifile = yaml_get_var(TRIM(this%configFile), 'RTTOV', 'inst_name', this%inst_name)
    ifile = yaml_get_var(TRIM(this%configFile), 'RTTOV', 'platform_name', this%platform_name)
    ifile = yaml_get_var(TRIM(this%configFile), 'RTTOV', 'turnOn', this%turnOn)

    n_plats = size(this%platform_name, 1)
    n_insts = size(this%inst_name, 1)
    IF ( n_plats .NE. n_insts ) THEN
      PRINT *, 'STOP: Numbers of platform_name and inst_name are inconsistent'
      RETURN
    ELSE
      ALLOCATE(this%ObsSatellites(n_insts))
    END IF

    ! Auxtypes
    CALL this%mpddGlob%initialize()                                   ! Initialize the mpdd

    PRINT *, 'Finish initialization of mpddGlob: ', this%mpddGlob%myrank

!     ! Reading the verification area latlon: Yuanfu Xie added on 2022-06-11
!     IF (this%mpddGlob%isBaseProc()) THEN
!       WRITE (*, 19)
! 19    FORMAT('#================================================================================#')
!       WRITE (*, 18)
! 18    FORMAT('# If you have not set Verify parameters in your yaml file, this may crash!!! #', /, &
!              '# If so, please set it and rerun!                                                #', /, &
!              '# otherwise, you can ignore this message.                                        #')
!       WRITE (*, 19)
!     END IF
!     istatus = yaml_get_var(TRIM(configFile), 'Verify', 'VA_latlon', this%verifyLatlon)
!     IF (istatus .NE. 0) THEN
!       WRITE (*, 55)
! 55    FORMAT('App initial: verification area is not specified and then no verification is applied')
!       this%verifyLatlon = 0.0D0
!     END IF
!     istatus = yaml_get_var(TRIM(configFile), 'Verify', 'VA_time', this%verifyTime)
!     IF (istatus .NE. 0) THEN
!       WRITE (*, 56)
! 56    FORMAT('App initial: verification time is not specified and then no verification is applied')
!       this%verifyTime = 0.0D0
!     END IF
!     IF (istatus .EQ. 0 .AND. this%mpddGlob%isBaseProc()) WRITE (*, 15) this%verifyLatlon, this%verifyTime
! 15  FORMAT('Applications_HybInc/Initial - verifying area latlon: ', 4D12.4, ' time window: ', 2E12.4)

    ! Initialize geometry
    CALL this%geometry%initialize(this%configFile, this%mpddGlob)

    ! Initialize the geometry
    ALLOCATE (this%XbMG(this%geometry%mg%mg_coarsest:this%geometry%mg%mg_finest))
    ALLOCATE (this%XbRef(this%geometry%mg%mg_coarsest:this%geometry%mg%mg_finest))

    istatus = yaml_get_var(TRIM(this%configFile), 'BMat', 'ensNum', this%ensNum)
    IF (istatus .NE. 0) THEN
      this%ensNum = 0
    ELSE
      WRITE (*, 12) this%ensNum, istatus
12    FORMAT("Number of ensembles: ", I2, ' yaml status', I2)
    END IF

    istatus = yaml_get_var(TRIM(this%configFile), 'BMat', 'Hybrid', this%Hybrid)
    IF (istatus .NE. 0) THEN
      PRINT *, 'App initial: Hybrid is not specified. Use of a default .FALSE.'
      this%Hybrid = .FALSE.
    END IF

    istatus = yaml_get_var(TRIM(this%configFile), 'BMat', 'HybInc', this%HybInc)
    IF (istatus .NE. 0) THEN
      PRINT *, 'App initial: HybInc is not specified. Use of a default .FALSE.'
      this%HybInc = .FALSE.
    END IF

    IF (this%ensNum .GT. 1) THEN
      ALLOCATE (this%XbMGEns(this%geometry%mg%mg_coarsest:this%geometry%mg%mg_finest, this%ensNum))

      IF (this%HybInc) THEN
        ALLOCATE (this%Xb_enperts(this%geometry%mg%mg_coarsest:this%geometry%mg%mg_finest, this%ensNum))
      END IF
    END IF

    ! Multigrid scheme:
    istatus = yaml_get_var(TRIM(this%configFile), 'MultigridOptions', 'UpdateBkgd', this%UpdateBkgd)
    IF (istatus .NE. 0) THEN
      PRINT *, 'App initial: UpdateBkgd is not specified. Use of a default .FALSE.'
      this%UpdateBkgd = .FALSE.

      IF (this%HybInc) THEN
        this%UpdateBkgd = .TRUE.
        PRINT *, 'App initial: UpdateBkgd is setted as .TRUE. due to HybInc is .TRUE.'
      END IF
    END IF
    PRINT *, 'App initial: UpdateBkgd is setted as:', this%UpdateBkgd

   ALLOCATE(this%miniSolver)

        CALL this%miniSolver%initialize(configFile)

    istatus = yaml_get_var(TRIM(this%configFile), 'BMat', 'BECsolver', this%BECsolver)
    IF (istatus .NE. 0) THEN
      WRITE (*, 57)
57    FORMAT('App initial: BECsolvet is not specified. Use of a default Laplace')
      this%BECsolver = 'Laplace'
    END IF

    ! Initialize the Scaling implementation
    CALL this%Ctl2State%initialize(this%configFile)

    ! Get background fillin information:
    BLOCK
      INTEGER(i_kind), ALLOCATABLE :: irange(:)

      istatus = yaml_get_var(TRIM(this%configFile), 'modelState', 'Fillin_range', irange)
      IF (istatus .EQ. 0 .AND. MINVAL(irange) .GE. 0) THEN
        IF (UBOUND(irange, 1) .NE. 3 * this%numVar) THEN
          WRITE (*, 16) UBOUND(irange, 1), this%numVar
16        FORMAT('The background fillin parameters:', I2, ' do not match the number of model states: ', I2, /, &
                 'Check your yaml file under modelState named as Fillin-range and rerun!')
          STOP
        ELSE
          this%iRange_ramp = RESHAPE(irange, (/3, this%numVar/))
        END IF
      ELSE
        this%iRange_ramp(1, :) = this%geometry%mg%sg(this%geometry%mg%mg_finest)%vLevel + 100
        this%iRange_ramp(2, :) = this%geometry%mg%sg(this%geometry%mg%mg_finest)%num_cell + 100
        this%iRange_ramp(3, :) = this%geometry%mg%sg(this%geometry%mg%mg_finest)%tSlots + 100
      END IF

      IF (ALLOCATED(irange)) DEALLOCATE (irange)

      ! Read in multigrid smoothing info:
      istatus = yaml_get_var(TRIM(this%configFile), 'modelState', 'mgGridSmooth', irange)
      IF (istatus .EQ. 0 .AND. UBOUND(irange, 1) .EQ. (this%mgEnd - this%mgStart + 1) * this%numVar) THEN
        PRINT *, 'IRANGE: ', irange, this%numVar, this%mgEnd - this%mgStart + 1
        ALLOCATE (this%mgGridSmooth(this%numVar, this%mgStart:this%mgEnd))
        this%mgGridSmooth(:, this%mgStart:this%mgEnd) = RESHAPE(irange, (/this%numVar, this%mgEnd - this%mgStart + 1/))
        DO i = this%mgStart, this%mgEnd
          WRITE (*, 20) i, this%mgGridSmooth(:, i)
20        FORMAT('mgGridSmooth at G', I2, ' number of smooth: ', 20I2)
        END DO
      ELSE
        IF (istatus .EQ. 0) THEN
          WRITE (*, 21) UBOUND(irange, 1), (this%mgEnd - this%mgStart + 1) * this%numVar
21        FORMAT('No multigrid background smoothing as yaml parameters do not match with multigrid variables: ', 2I4)
        ELSE
          WRITE (*, 22)
22        FORMAT('No multigrid background smoothing parameters specified in your yaml...')
        END IF
      END IF

      IF (ALLOCATED(irange)) DEALLOCATE (irange)
    END BLOCK

    IF (this%mpddGlob%isBaseProc()) THEN
      DO i = 1, this%numVar
        WRITE (*, 17) this%iRange_ramp(:, i), this%varList(i)
17      FORMAT('Background Fillin_range: ', I4, I8, I4, ' for variable: ', A)
      END DO
    END IF

    istatus = 0 ! initialization successful
    CALL CPU_TIME(t2)
    WRITE (*, 1) t2 - t1
1   FORMAT('Time spent in application initial: ', D12.4)

  END SUBROUTINE initial

  !> @brief
  !================================================================
  !  This routine sets up background fields for an application:
  !  Input:
  !       mgStart:    integer of the start G level number
  !       mgEnd:      integer of the end G level number
  !
  SUBROUTINE backgrd(this, mgStart, mgEnd, unitTest)
    IMPLICIT NONE
    CLASS(Applications_HybInc_t) :: this
    INTEGER(i_kind), INTENT(IN) :: mgStart, mgEnd
    LOGICAL, INTENT(IN) :: unitTest

    ! Local variables:
    INTEGER(i_kind) :: i, ii
    REAL(r_kind) :: t1, t2

    TYPE(IOWRF_t), TARGET :: ioWRF
    TYPE(IOGrapes_t), TARGET :: ioGrapes
    TYPE(IOGrapesPostvar_t), TARGET :: ioGrapesPostvar

    ! Calculate the time used in background ingest:
    CALL CPU_TIME(t1)

    ! Check if a requested multigrid is valid:
    IF (mgStart .LT. this%geometry%mg%mg_coarsest .OR. &
        mgEnd .GT. this%geometry%mg%mg_finest) THEN
      WRITE (*, 2) mgStart, mgEnd, this%geometry%mg%mg_coarsest, this%geometry%mg%mg_finest
2     FORMAT('Your requested multigrid ', 2I2, ' exceed the geometry setting: (', I2, ', ', I2, /, &
             'Please resetting the geometry and rerun!')
      STOP
    END IF

    ! Get background fields:
    DO i = mgEnd, mgStart, -1
      IF (i .EQ. mgEnd) THEN
        ! Initialize a zeros background field at finest grid.
        ASSOCIATE (sg => this%geometry%mg%sg(mgEnd))
          CALL this%XbMG(mgEnd)%initialize(this%configFile, sg)
          CALL this%XbMG(mgEnd)%ShowInfo

          IF (.NOT. unitTest) THEN
            BLOCK
              CHARACTER(LEN=20) :: bkModel
              INTEGER(i_kind):: fileStat

              fileStat = yaml_get_var(this%configFile, 'IO', 'bk_model', bkModel)
              PRINT *, 'Background field is: ', bkModel
              IF (TRIM(bkModel) == 'WRF') THEN
                PRINT *, 'In Reading ', TRIM(bkModel)
                ALLOCATE (IOWRF_t:: this%IOModel)

                ASSOCIATE (IOModel => this%IOModel)
                  SELECT TYPE (IOModel)
                  TYPE IS (IOWRF_t)
                    CALL IOModel%initialize(this%configFile, this%geometry)
                  END SELECT
                END ASSOCIATE
                CALL this%IOModel%m_read_bcg_into_Xm_Ens(this%XbMG(mgEnd), sg)
                PRINT *, 'WRF reading is done!'

              ELSE IF (TRIM(bkModel) == 'GRAPES') THEN
                PRINT *, 'In Reading ', TRIM(bkModel)
                ALLOCATE (IOGrapes_t:: this%IOModel)

                ASSOCIATE (IOModel => this%IOModel)
                  SELECT TYPE (IOModel)
                  TYPE IS (IOGrapes_t)
                    CALL IOModel%initialize(this%configFile, this%geometry)
                  END SELECT
                END ASSOCIATE

                PRINT *, 'READING GRAPES...'
                CALL this%IOModel%m_read_bcg_into_Xm(this%XbMG(mgEnd), sg)

                ! this%XbMG(mgEnd)%fields(this%XbMG(mgEnd)%getVarIdx('uwnd'))%data = 0.0D0
                ! this%XbMG(mgEnd)%fields(this%XbMG(mgEnd)%getVarIdx('vwnd'))%data = 0.0D0

                PRINT *, 'GRAPES reading is done!'

              ELSE IF (TRIM(bkModel) == 'POSTVAR') THEN
                PRINT *, 'In Reading ', TRIM(bkModel)
                ALLOCATE (IOGrapesPostvar_t:: this%IOModel)

                ASSOCIATE (IOModel => this%IOModel)
                  SELECT TYPE (IOModel)
                  TYPE IS (IOGrapesPostvar_t)
                    CALL IOModel%initialize(this%configFile, this%geometry)
                  END SELECT
                END ASSOCIATE

                PRINT *, 'READING GrapesPostvar...'
                PRINT *, 'READING GrapesPostvar...'
                CALL this%IOModel%m_read_bcg_into_Xm(this%XbMG(mgEnd), sg)

                PRINT *, 'GrapesPostvar reading is done!'
              END IF
            END BLOCK

          END IF
        END ASSOCIATE
      ELSE
        ! Restrict to each coarser grid
        ASSOCIATE (sgFiner => this%geometry%mg%sg(i + 1), sgCoarser => this%geometry%mg%sg(i))
          CALL this%XbMG(i)%initialize(this%configFile, sgCoarser)
          PRINT *, 'Before restraction'
          CALL restrictionMG(this%XbMG(i), this%XbMG(i + 1), this%geometry%mg)
          PRINT *, 'After restraction'
          CALL this%geometry%mg%restrictionOfStatics(sgFiner, sgCoarser)
          CALL update_restrictionOfStatics(this%configFile, this%XbMG(i), sgCoarser)
        END ASSOCIATE
      END IF
      ! PRINT *, 'GRAPES reading is done!',maxval(this%XbMG(i)%Fields(this%XbMG(i)%getVarIdx('pres'))%data)

      IF (TRIM(this%mode) == 'Debug') THEN
        BLOCK
          ! TYPE(UV2W_t) :: UV2W

          ! UV2W = UV2W_t(this%configFile, this%XbMG(i))
          ! CALL UV2W%transFwdNonLinear(this%XbMG(i))
          PRINT *, "run motor in debug mode ----------*************"
          CALL Output_NC_State_AV(this%XbMG(i), this%ncOutputFile, &
                                  TRIM(this%task)//"_bak", .TRUE., .TRUE.)

          ! CALL this%XbMG(i)%rmVar('wwnd')
        END BLOCK
      END IF

      IF (unitTest) THEN
        DO ii = 1, this%numVar
          this%XbMG(i)%fields(ii)%DATA = 0.8D0
        END DO
      END IF
    END DO

    ! Add the Variable replacement
    DO i = mgEnd, mgStart, -1
      ! IF(this%geometry%mg%sg(i)%isActiveProc()) THEN
      !   BLOCK
      !     REAL(r_kind), ALLOCATABLE :: scaleQvapor(:)
      !     INTEGER :: j

      !     ALLOCATE(scaleQvapor(this%geometry%mg%sg(i)%vLevel))

      !     DO j = 1,this%geometry%mg%sg(i)%vLevel
      !       scaleQvapor(j) = SUM(this%XbMG(i)%Fields(this%XbMG(i)%getVarIdx('qvapor'))%data(j, 1:this%geometry%mg%sg(i)%num_icell, :))&
      !                        /this%geometry%mg%sg(i)%num_icell/this%geometry%mg%sg(i)%tSlots
      !       CALl this%XbMG(i)%sg%mpddInfo_sg%AllReduceSumReal(scaleQvapor(j), this%geometry%mg%sg(i)%s_qvapor(j))
      !       this%XbMG(i)%sg%s_qvapor(j) = this%XbMG(i)%sg%s_qvapor(j)/this%XbMG(i)%sg%mpddInfo_sg%nProc
      !     END DO

      !     IF(this%geometry%mg%sg(i)%isBaseProc())PRINT*,'scaleQvapor: ', scaleQvapor
      !       DEALLOCATE(scaleQvapor)

      !   END BLOCK
      ! ENDIF

      this%XbRef(i) = this%XbMG(i)

      IF ( TRIM(this%framework) .EQ. 'Incremental') THEN
        CALL this%Ctl2State%SubState(this%XbMG(i), this%XbRef(i))
      END IF
      CALL this%Ctl2State%transBackward(this%XbMG(i))
      
      BLOCK
        CHARACTER(LEN=2) :: iGrid
        WRITE (iGrid, "(I2.2)") this%XbMG(i)%sg%gLevel
        IF (SIZE(this%XbMG(i)%sg%SVapor, 2) .GT. 0) THEN
          OPEN (90, file=TRIM(this%ncOutputFile)//'/sv_g'//TRIM(iGrid)//'.txt', form='formatted')
          ! write(90,*) this%XbMG(i)%sg%SVapor1D
          WRITE (90, *) this%XbMG(i)%sg%SVapor(:, 1)
          CLOSE (90)
        END IF
      END BLOCK

      ! Smooth the background fields in the control space:
      IF (ALLOCATED(this%mgGridSmooth)) THEN
        CALL this%XbMG(i)%hStateFilter(this%mgGridSmooth(:, i))
        IF (TRIM(this%mode) == 'Debug') THEN
          CALL Output_NC_State_AV(this%XbMG(i), this%ncOutputFile, &
                                  TRIM(this%task)//"_smoothedBkgd", .TRUE., .TRUE.)
        END IF
      END IF

    END DO

    CALL CPU_TIME(t2)
    ! WRITE (*, 1) t2 - t1,UBOUND(this%XbMG(5)%fields,1), this%mpddGlob%myrank
1   FORMAT('backgrd: time spent ingesting background fields: ', D12.4, ' NUM-FIELDS: ', I2, ' at proc: ', I3)
  END SUBROUTINE backgrd

  SUBROUTINE dumpAna(this, unitTest)
    USE Export2SelDomain_m, ONLY: Export2SelDomain_t
    USE Export2HASCoordInSelDomain_m, ONLY: Export2HASCoordInSelDomain_t
    USE PostProcTools_m
    IMPLICIT NONE
    CLASS(Applications_HybInc_t) :: this
    LOGICAL, INTENT(IN) :: unitTest
    INTEGER(i_kind) :: i
    TYPE(State_t) :: XbRef

    IF (.NOT. unitTest) THEN
      ! Scaling back with control variable to physical space:
      XbRef = this%XbMG(this%mgEnd)
      CALL this%Ctl2State%transFwdNonLinear(this%analysis)
      CALL this%Ctl2State%transFwdNonLinear(XbRef)
      IF (TRIM(this%framework) .EQ. 'Incremental') THEN
        CALL this%Ctl2State%SumState(this%analysis, this%XbRef(this%mgEnd))
        CALL this%Ctl2State%SumState(XbRef, this%XbRef(this%mgEnd))
      END IF

      IF (TRIM(this%task) == '3DVar') THEN
        ASSOCIATE (ioModel => this%ioModel)
          SELECT TYPE (ioModel)
          TYPE IS (IOGrapes_t)
            PRINT *, 'Dumping Grapes background files.'
            CALL ioModel%m_write_Xm_into_bcg(this%analysis, this%geometry%mg%sg(this%mgEnd), XbRef)
          TYPE IS (IOWRF_t)
            ! Some wrf codes.
          END SELECT
        END ASSOCIATE
      ELSE IF (TRIM(this%task) == 'SfcAna') THEN
        BLOCK
          TYPE(Export2SelDomain_t) ::  Export2SelDomain

          ! CALL add_uv_to_state(this%analysis)
          Export2SelDomain = Export2SelDomain_t(this%configFile, this%analysis%sg)
          CALL Export2SelDomain%Export_State_AVLNST(this%analysis, this%ncOutputFile, &
                                                    TRIM(this%task), 1, .FALSE., .FALSE.)

          DO i = LBOUND(this%analysis%fields, 1), UBOUND(this%analysis%fields, 1)
            CALL Export2SelDomain%Export_State_SVLNST(this%analysis, this%ncOutputFile, &
                                                      TRIM(this%task), this%analysis%fields(i)%Get_Name(), &
                                                      1, .FALSE., .FALSE.)
          END DO
        END BLOCK
      ELSE IF (TRIM(this%task) == '3DAna') THEN
        BLOCK
          TYPE(Export2HASCoordInSelDomain_t) ::  Export2HASCoordInSelDomain
          ! CALL add_10min_uv_to_state(this%analysis)
          ! CALL add_wind_pressure(this%analysis)

          Export2HASCoordInSelDomain = Export2HASCoordInSelDomain_t(this%configFile, this%analysis%sg, this%analysis)
          CALL Export2HASCoordInSelDomain%Export_State_AVLNST(this%analysis, this%ncOutputFile, &
                                                              TRIM(this%task), 1, .FALSE., .FALSE.)

          DO i = LBOUND(this%analysis%fields, 1), UBOUND(this%analysis%fields, 1)
            CALL Export2HASCoordInSelDomain%Export_State_SVLNST(this%analysis, this%ncOutputFile, &
                                                                TRIM(this%task), this%analysis%fields(i)%Get_Name(), &
                                                                1, .FALSE., .FALSE.)
          END DO
        END BLOCK
      END IF

    END IF
  END SUBROUTINE dumpAna

  !> @brief
  !========================================================================
  !  This routine sets up ensemble background fields for an application:
  !  Up to now, only works for WRF and GrapesModelvar
  !  Input:
  !       mgStart:    integer of the start G level number
  !       mgEnd:      integer of the end G level number
  !
  SUBROUTINE backgrdEns(this)
    IMPLICIT NONE
    CLASS(Applications_HybInc_t) :: this
    ! INTEGER(i_kind), INTENT(IN) :: mgStart,mgEnd
    ! LOGICAL, INTENT(IN) :: unitTest

    ! Local variables:
    TYPE(IOWRF_t), ALLOCATABLE, TARGET :: ioWRF
    TYPE(IOGrapesModelvar_t), ALLOCATABLE, TARGET :: ioGrapesModelvar
    INTEGER(i_kind) :: i, ii, iens, mgStart, mgEnd, numCtrl
    REAL(r_kind) :: t1, t2
    CHARACTER(LEN=20) :: ensModel
    INTEGER(i_kind)   :: istatus
    CHARACTER(LEN=1)  :: ensID

    TYPE(State_t) :: Xbmean_tmp

    CALL CPU_TIME(t1)

    istatus = yaml_get_var(TRIM(this%configFile), 'BMat', 'ensModel', ensModel)

    mgStart = this%mgStart
    mgEnd = this%mgEnd
    ! Get background fields:
    DO i = mgEnd, mgStart, -1
      IF (i .EQ. mgEnd) THEN
        ! Initialize a zeros background field at finest grid.
        ASSOCIATE (sg => this%geometry%mg%sg(mgEnd))
          ! this%CALL XbMG(mgEnd) %initialize(this%configFile, sg)
          ! Initialize the ioWRF
          IF (TRIM(ensModel) == "WRF") THEN
            ALLOCATE (ioWRF)
            CALL ioWRF%initialize(this%configFile, this%geometry)
            CALL this%XbMGEns(i, 1)%initialize(this%configFile, sg)
            DO iens = 2, this%ensNum
              this%XbMGEns(i, iens) = this%XbMGEns(i, 1)
            END DO
            DO iens = 1, this%ensNum
              CALL ioWRF%m_read_bcg_into_Xm_Ens(this%XbMGEns(i, iens), sg, iens)   ! read wrfout into XbMG
            END DO
            DEALLOCATE (ioWRF)

          ELSE IF (TRIM(ensModel) == "GrapesModelvar") THEN
            ALLOCATE (ioGrapesModelvar)
            CALL ioGrapesModelvar%initialize(this%configFile, this%geometry)
            CALL this%XbMGEns(i, 1)%initialize(this%configFile, sg)
            DO iens = 2, this%ensNum
              this%XbMGEns(i, iens) = this%XbMGEns(i, 1)
            END DO
            DO iens = 1, this%ensNum
              CALL ioGrapesModelvar%m_read_bcg_into_Xm_Ens(this%XbMGEns(i, iens), sg, iens)   ! read wrfout into XbMG
            END DO
            DEALLOCATE (ioGrapesModelvar)

          END IF
        END ASSOCIATE

      ELSE
        ! Restrict to each coarser grid
        ASSOCIATE (sgFiner => this%geometry%mg%sg(i + 1), sgCoarser => this%geometry%mg%sg(i))
          CALL this%XbMGEns(i, 1)%initialize(this%configFile, sgCoarser)
          DO iens = 2, this%ensNum
            this%XbMGEns(i, iens) = this%XbMGEns(i, 1)
          END DO
          DO iens = 1, this%ensNum
            CALL restrictionMG(this%XbMGEns(i, iens), this%XbMGEns(i + 1, iens), this%geometry%mg)
          END DO
        END ASSOCIATE
      END IF
    END DO

    ! ! Add the Variable replacement
    ! DO i = mgEnd, mgStart, -1
    !   DO iens = 1, this%ensNum
    !     CALL this%Ctl2State%transBackward(this%XbMGEns(i,iens))
    !   END DO
    ! END DO

    ! Calculate the ensemble perturbations
    IF (this%HybInc) THEN
      DO i = mgEnd, mgStart, -1
        Xbmean_tmp = this%XbMGEns(i, 1)%zeroCopy()
        DO iens = 1, this%ensNum
          Xbmean_tmp = Xbmean_tmp + this%XbMGEns(i, iens) / (this%ensNum * 1.0D0)
        END DO
        ! CALL Output_NC_State_AV(Xbmean_tmp, this%ncOutputFile, TRIM(this%task)//"_EnsMean", .TRUE., .TRUE.)
        DO iens = 1, this%ensNum
          this%Xb_enperts(i, iens) = this%XbMGEns(i, iens)%zeroCopy()
          ! this%Xb_enperts(i, iens) = ( this%XbMGEns(i, iens) - Xbmean_tmp ) / SQRT( (this%ensNum - 1) * 1.0D0 )
          numCtrl = this%geometry%mg%sg(i)%vLevel * this%geometry%mg%sg(i)%num_icell_global * this%geometry%mg%sg(i)%tSlots
          this%Xb_enperts(i, iens) = (this%XbMGEns(i, iens) - Xbmean_tmp) / SQRT((this%ensNum) * 1.0D0)
          ! this%Xb_enperts(i, iens) = ( this%XbMGEns(i, iens) - Xbmean_tmp ) / SQRT( numCtrl * 1.0D0 )
          ! WRITE (ensID, '(I1)') iens
          ! CALL Output_NC_State_AV(this%Xb_enperts(i, iens), this%ncOutputFile, TRIM(this%task)//"_EnsPertMem"//TRIM(ensID), .TRUE., .TRUE.)
          ! CALL Output_NC_State_AV(this%XbMGEns(i, iens), this%ncOutputFile, TRIM(this%task)//"_EnsMem"//TRIM(ensID), .TRUE., .TRUE.)
        END DO
      END DO
    END IF

    CALL CPU_TIME(t2)
    WRITE (*, 1) t2 - t1
1   FORMAT('Time spent on background ensemble: ', D12.4)
  END SUBROUTINE backgrdEns

  !> @brief
  !! Calculate BEC through QR decomposition method:
  SUBROUTINE calEnsB(this)
    IMPLICIT NONE
    CLASS(Applications_HybInc_t) :: this

    ! Local variables:
    INTEGER(i_kind) :: iv, i
    REAL(r_kind) :: ana_max, ana_min, ana_max_global, ana_min_global, t1, t2
    REAL(r_kind), ALLOCATABLE :: ana_max_AV(:), ana_min_AV(:)
    REAL(r_kind), PARAMETER :: passthreshold_max = 6.0D0, passthreshold_min = 0.6D0
    TYPE(EnLoc_t) :: EnLoc

    CALL CPU_TIME(t1)

    CALL this%backgrdEns()

    IF (this%mpddGlob%isBaseProc()) THEN
      DO i = this%mgStart, this%mgEnd
        CALL EnLoc%b_getEnsdata(this%configFile, this%XbMGEns(i, :))
        PRINT *, 'Done EnLoc%b_getEnsdata in grid', i

        CALL EnLoc%b_qrDecomp(i)
        PRINT *, 'Done EnLoc%b_qrDecomp in grid', i

        CALL EnLoc%b_qrRinv()
        PRINT *, 'Done EnLoc%b_qrRinv in grid', i

        CALL EnLoc%b_ADJ(i)
        PRINT *, 'Done EnLoc%b_ADJ in grid', i

        CALL EnLoc%b_destroy()
        PRINT *, 'Done EnLoc%b_destroy in grid', i
      END DO
    END IF

    IF (ALLOCATED(this%XbMGEns)) DEALLOCATE (this%XbMGEns)

    CALL CPU_TIME(t2)
    WRITE (*, 1) t2 - t1
1   FORMAT('Time spent in calEnsB: ', D12.4)

  END SUBROUTINE calEnsB

  !> @brief
  !================================================================
  !  This routine reads all observation data for an application:
  !  Input:
  !       obsList:    character string array in length 8 containing
  !                   all observations:
  !                   surfaces (surface obs);
  !                   sounding (RAOB obs);
  !                   profiler (wind profiler obs)
  !                   ......
  !
  SUBROUTINE readObs(this, obsList, mgStart, mgEnd, unitTest)
    IMPLICIT NONE
    CLASS(Applications_HybInc_t) :: this
    CHARACTER(LEN=8), INTENT(IN) :: obsList(:)
    INTEGER(i_kind), INTENT(IN) :: mgStart, mgEnd
    LOGICAL, INTENT(IN) :: unitTest
    TYPE(State_t) :: XbRef

    ! Local variables:
    INTEGER(i_kind) :: i, i_inst, j
    REAL(r_kind) :: t1, t2

    ! XbRef = this%Ctl2State%fwdNL_opr(this%XbMG(mgEnd))
    XbRef = this%XbRef(mgEnd)

    ! Calculating CPU time used for obs ingest:
    CALL CPU_TIME(t1)

    DO i = 1, UBOUND(obsList, 1)
      SELECT CASE (obsList(i))
      CASE ('surfaces')
        CALL this%sfc%ObsInitial(this%configFile)
        CALL this%sfc%ObsIngest(XbRef)
        IF (unitTest) CALL this%replace(this%sfc%numobs, &
                                        this%sfc%numVars, this%sfc%olatlon, &
                                        this%sfc%obsTime, this%geometry%mg%sg(mgEnd), &
                                        this%sfc%obsData)
        PRINT *, 'readObs - sfc min/max: ', &
          MINVAL(this%sfc%obsData), MAXVAL(this%sfc%obsData)
      CASE ('sounding')
        CALL this%snd%ObsInitial(this%configFile)
        CALL this%snd%ObsIngest(XbRef)
        IF (unitTest) CALL this%replace(this%snd%numobs, &
                                        this%snd%numVars, this%snd%olatlon, &
                                        this%snd%obsTime, this%geometry%mg%sg(mgEnd), &
                                        this%snd%obsData)
        PRINT *, 'readObs - snd min/max: ', &
          MINVAL(this%snd%obsData), MAXVAL(this%snd%obsData)
      CASE ('profiler')
        CALL this%pro%ObsInitial(this%configFile)
        CALL this%pro%ObsIngest(XbRef)
        IF (unitTest) CALL this%replace(this%pro%numobs, &
                                        this%pro%numVars, this%pro%olatlon, &
                                        this%pro%obsTime, this%geometry%mg%sg(mgEnd), &
                                        this%pro%obsData)
        PRINT *, 'readObs - pro min/max: ', &
          MINVAL(this%pro%obsData), MAXVAL(this%pro%obsData)
      CASE ('radarRef')
        CALL this%ref%ObsInitial(this%configFile)
        CALL this%ref%ObsIngest(XbRef)
        PRINT *, 'readObs - ref min/max: ', &
          MINVAL(this%ref%obsData), MAXVAL(this%ref%obsData)
      CASE ('radarVel')
        CALL this%vel%ObsInitial(this%configFile)
        CALL this%vel%ObsIngest(XbRef)
        ! PRINT *,'readObs - vel min/max: ', &
        !    minval(this%vel%obsData),maxval(this%vel%obsData) ! Radar vel is not fill this variable
        ! Yuanfu Xie 2022-11-27 added Dr. Zhang Hua's cloud derived wind:
        ! Yuanfu Xie 2022-11-27 added Dr. Zhang Hua's cloud derived wind:
      CASE ('cloudWnd')
        CALL this%cdw%ObsInitial(this%configFile)
        CALL this%cdw%ObsIngest(XbRef)
      CASE ('shipRpts')
        CALL this%shp%ObsInitial(this%configFile)
        CALL this%shp%ObsIngest(XbRef)
      CASE ('bouyance')
        CALL this%boy%ObsInitial(this%configFile)
        CALL this%boy%ObsIngest(XbRef)
      CASE ('pilotRpt')
        CALL this%air%ObsInitial(this%configFile)
        CALL this%air%ObsIngest(XbRef)
      END SELECT
    END DO

    IF (ALL(.NOT. this%turnOn)) THEN
      PRINT *, " =========== NO satellite data will be assimilated ==========="
      RETURN
    ELSE
      DO i_inst = 1, SIZE(this%ObsSatellites)
        IF (this%turnOn(i_inst)) THEN
          PRINT *, 'readObs: ',TRIM(this%platform_name(i_inst))//'-'//TRIM(this%inst_name(i_inst)),' is being processed'
          CALL this%ObsSatellites(i_inst)%ObsInitial(this%configFile, i_inst)
          CALL this%ObsSatellites(i_inst)%ObsIngest(XbRef)
        END IF
      END DO
    END IF

    CALL cpu_time(t2)
    WRITE (*, 1) t2 - t1, this%mpddGlob%myrank
1   FORMAT('readObs: time spent ingesting observation data: ', D12.4, ' at proc: ', I3)

  END SUBROUTINE readObs

  SUBROUTINE varAnal(this, X, obsList, mg)
    CLASS(Applications_HybInc_t) :: this
    TYPE(State_t) :: X, XbCV
    INTEGER(i_kind) :: mg
    CHARACTER(LEN=8), INTENT(IN) :: obsList(:)

    TYPE(BMatrix_t) :: B
    TYPE(RMatrix_t) :: R
    TYPE(C2O_t) :: H
    TYPE(JFunc_HybInc_t) :: JFunc_HybInc

    INTEGER(i_kind) :: ii, j, ia
    REAL(r_kind) :: t1, t2

    ! IF(sg%isActiveProc()) THEN
    ! Initialize X, Y, B, H
    ! CALL B%initialize(this%configFile, X%sg)     ! Initialize the B matrix
    CALL B%initialize(this%configFile, X%sg, this%BECsolver)

    ! Thinning observations:
    CALL CPU_TIME(t1)       ! Timing observation thinning

    CALL H%initialize(this%configFile, X, this%thinnedObs(mg))

    CALL R%initialize(this%configFile, this%thinnedObs(mg), X%sg) ! Initialize R

    ! Initialize J Function
    SELECT CASE (TRIM(this%framework))
    CASE ('FullState')
      XbCV = this%XbMG(X%sg%gLevel)
    CASE ('Incremental')
      XbCV = this%XbRef(X%sg%gLevel)
      CALL this%Ctl2State%transBackward(XbCV)
      CALL this%Ctl2State%SumBothCtlState(XbCV, this%XbMG(X%sg%gLevel))
    END SELECT
    CALL JFunc_HybInc%initialize(this%configFile, X, this%thinnedObs(mg), H, B, R, X%sg, XbCV)

    ! Run minimization
    CALL CPU_TIME(t1)
    IF (X%sg%gLevel <= 5) THEN
      CALL this%miniSolver%run_HyyInc(X, JFunc_HybInc, X%sg, 50)
    ELSE
      CALL this%miniSolver%run_HyyInc(X, JFunc_HybInc, X%sg)
    END IF
    CALL CPU_TIME(t2)
    WRITE (*, 3) t2 - t1, mg, X%sg%mpddInfo_sg%myrank, X%sg%num_cell
3   FORMAT('Time spent in varAnal miniSolver: ', D10.2, ' MG: ', I2, ' at proc: ', I1, ' with Ncell: ', I8)

  END SUBROUTINE varAnal

  SUBROUTINE DASolve(this, X, obsList, sg, Nocoarest)
    CLASS(Applications_HybInc_t) :: this
    TYPE(State_t) :: X, XbRef, XbCV
    TYPE(SingleGrid_t) :: sg
    CHARACTER(LEN=*), INTENT(IN) :: obsList(:)
    INTEGER(i_kind) :: Nocoarest(:)
    INTEGER, PARAMETER :: NYmax = 15

    TYPE(BMatrix_t) :: B
    TYPE(BMatrix_t) :: B_e, B_einc
    TYPE(RMatrix_t) :: R
    TYPE(MPObs_t), TARGET :: mpObs
    TYPE(ObsSet_t) :: Yi(NYmax)
    TYPE(ObsSet_t) :: Y, Z
    TYPE(C2O_t) :: H, HY
    TYPE(JFunc_HybInc_t) :: JFunc_HybInc
    INTEGER(i_kind) :: ii, j, i, NObs_tmp, NY
    REAL(r_kind) :: t1, t2
    ! TYPE(UV2W_t) :: UV2W

    CALL mpObs%initializeMPObs(sg)    ! Initialize the observation parallel processing proc

    IF (.NOT. sg%isActiveProc()) RETURN


!     write(*,555) this%XbMG(sg%gLevel)%getVarIdx('pres_log'), &
!       minval(this%XbMG(sg%gLevel)%fields(this%XbMG(sg%gLevel)%getVarIdx('pres_log'))%data(1,:,:)), &
!       maxval(this%XbMG(sg%gLevel)%fields(this%XbMG(sg%gLevel)%getVarIdx('pres_log'))%data(1,:,:)), &
!       this%XbMG(sg%gLevel)%sg%gLevel,this%XbMG(sg%gLevel)%sg%mpddInfo_sg%myrank
! 555 format('DASolve check controls pres_log: ',I2,' min-max:',2D14.6,' G',I2,' pc',I2)
    ! XbRef = this%Ctl2State%fwdNL_opr(this%XbMG(sg%gLevel))
    XbRef = this%XbRef(sg%gLevel)
    ! UV2W = UV2W_t(this%configFile, XbRef)
    ! CALL UV2W%transFwdNonLinear(XbRef)

    ! Thinning observations:
    CALL CPU_TIME(t1)       ! Timing observation thinning

    DO j = 1, UBOUND(obsList, 1)
      SELECT CASE (TRIM(obsList(j)))
      CASE ('surfaces')
        CALL this%sfc%ObsThinning(XbRef, Yi(1), mpObs, .TRUE., .FALSE.)
        ! A temporary processing to remove analysis on precipitation on coarser grid
        IF (sg%gLevel < this%mgEnd - 3) THEN
          CALL Yi(1)%rmVar('SYNOP_pcpa')
          CALL Yi(1)%rmVar('SYNOP_pcpa5min')
        END IF

        IF (sg%gLevel < this%mgEnd - 2) THEN
          CALL Yi(1)%rmVar('SYNOP_temp')
          CALL Yi(1)%rmVar('SYNOP_qvapor')
          CALL Yi(1)%rmVar('SYNOP_pres')
          CALL Yi(1)%rmVar('SYNOP_psl')
        END IF
        ! A temporary processing to remove analysis on precipitation on coarser grid
      CASE ('sounding')
        CALL this%snd%ObsThinning(XbRef, Yi(2), mpObs, .TRUE., .FALSE.)

        IF (sg%gLevel < this%mgEnd - 1) THEN
          CALL Yi(2)%rmVar('SOUND_qvapor')
        END IF
        CALL Yi(2)%rmVar('SOUND_pres') ! remove the pressure in sounding

      CASE ('profiler')
        CALL this%pro%ObsThinning(XbRef, Yi(3), mpObs, .TRUE., .FALSE.)
      CASE ('radarRef')
        IF (sg%gLevel >= this%mgEnd) THEN
          CALL this%ref%ObsThinning(XbRef, Yi(4), mpObs, .TRUE., .FALSE.)
        END IF
      CASE ('radarVel')
        IF (sg%gLevel >= this%mgEnd - 2) THEN
          CALL this%vel%ObsPrepareForSg(XbRef)
          CALL this%vel%ObsThinning(XbRef, Yi(5), mpObs, .TRUE., .FALSE.)
        END IF
      CASE ('cloudWnd')
        CALL this%cdw%ObsThinning(XbRef, Yi(6), mpObs, .TRUE., .FALSE.)
      CASE ('shipRpts')
        CALL this%shp%ObsThinning(XbRef, Yi(7), mpObs, .TRUE., .FALSE.)
      CASE ('buoyance')
        CALL this%boy%ObsThinning(XbRef, Yi(8), mpObs, .TRUE., .FALSE.)
      CASE ('pilotRpt')
        CALL this%air%ObsThinning(XbRef, Yi(9), mpObs, .TRUE., .FALSE.)

      END SELECT
    END DO

    ! PRINT*, 'DASolve: ', TRIM(this%task), ' is being processed'
    ! STOP

    IF (ALL(.NOT. this%turnOn)) THEN
      PRINT *, " =========== NO satellite data will be assimilated ==========="
    ELSE
      NY = 9
      DO i_inst = 1, SIZE(this%ObsSatellites)
        IF (this%turnOn(i_inst)) THEN
          PRINT *, 'DASolve: ', TRIM(this%platform_name(i_inst))//'-'//TRIM(this%inst_name(i_inst)), ' is being processed'
          NY = NY + 1

          CALL this%ObsSatellites(i_inst)%ObsPrepareForSg(XbRef)
          CALL this%ObsSatellites(i_inst)%ObsThinning(XbRef, Yi(NY), mpObs, .TRUE., .False.)

          IF (sg%gLevel < this%mgEnd - 1) THEN
            DO i = LBOUND(Yi(NY)%obsFields, 1), UBOUND(Yi(NY)%ObsFields, 1)
              ! PRINT *, 'check Yi(NY) Get_Name: ', Yi(NY)%ObsFields(i)%Get_Name()
              ! PRINT *, 'check Yi(NY) Get_ObsType: ', Yi(NY)%ObsFields(i)%Get_ObsType()
              ! IF (Yi(NY)%ObsFields(i)%Get_Name() == 'tbb') THEN
              IF(INDEX(TRIM(Yi(NY)%ObsFields(i)%Get_ObsType()), 'fy4') .NE. 0) THEN
                ! PRINT *, 'ONLY assimialte fy4 on fine grids'
                CALL Yi(NY)%ObsFields(i)%Set_Name('tdbb')
              END IF
            END DO
          END IF

        END IF
      END DO
    END IF

    CALL cpu_time(t2)
    WRITE (*, 2) t2 - t1, sg%mpddInfo_sg%myrank, sg%num_cell
2   FORMAT('Time spent obsThinning: ', D14.7, ' at proc: ', I1, ' numCell: ', I5)

    ! Concat thinning obs:
    Y = ObsSet_t(this%configFile, mpObs)

    CALL ObsConcat_s((/Yi(1), Yi(2), Yi(3), Yi(4), Yi(5), Yi(6), Yi(7), Yi(8), Yi(9), Yi(10), Yi(11), Yi(12), Yi(13)/), Y)
    WRITE(*,12) UBOUND(Y%ObsFields,1), &
      (TRIM(Y%ObsFields(j)%get_name()),j=1,UBOUND(Y%obsFields,1)), &
      (TRIM(Y%ObsFields(j)%Get_Id_Name()),j=1,UBOUND(Y%obsFields,1))
12  FORMAT('Concatenated obs vars: ', I2, 20(1X,A))

    ! Model state variables to fill in data sparse areas
    Z = ObsSet_t(this%configFile, mpObs)
    CALL this%bkgdFillin(sg%gLevel, Y, Z)

    ! update FGPres here if you want
    ! CALL update_FGPres(XbRef, X)

    ! Initialize H operator for Z filled with background and Y no background fill in:
    IF (sg%gLevel .EQ. this%mgStart) THEN
      DO i = LBOUND(Y%ObsFields, 1), UBOUND(Y%ObsFields, 1)
        NObs_tmp = SIZE(Y%ObsFields(i)%values)
        CALL Y%mpObs%AllReduceSumInt(NObs_tmp, Nocoarest(i))
      END DO
    ELSE
      PRINT *, 'Check Nocoarest ', Nocoarest
    END IF

    ! Initialize X, Y, B, H
    ! CALL B%initialize(this%configFile, sg)     ! Initialize the B matrix
    ! CALL B%initialize(this%configFile, sg, this%BECsolver)     ! Initialize the B matrix
    IF (this%Hybrid) THEN
      CALL B%initialize(this%configFile, sg, 'Laplace')
      CALL B_e%initialize(this%configFile, sg, 'EnLoc')
    ELSE IF (this%HybInc) THEN
      CALL B%initialize(this%configFile, sg, 'Laplace')
      CALL B_einc%initialize(this%configFile, sg, 'HybInc')
    ELSE
      CALL B%initialize(this%configFile, sg, 'Laplace')
    END IF

    CALL H%initialize(this%configFile, X, Z, this%XbRef(sg%gLevel))
    WRITE (*, 13) X%sg%glevel, X%sg%mpddInfo_sg%myrank
13  FORMAT('H is initialized....at G', I2, ' pc ', I2)
    CALL HY%initialize(this%configFile, X, Y, this%XbRef(sg%gLevel))
    WRITE (*, 14) X%sg%glevel, X%sg%mpddInfo_sg%myrank
14  FORMAT('HY is initialized....at G', I2, ' pc ', I2)

    ! CALL B%initialize(this%configFile, sg)     ! Initialize the B matrix

    ! CALL B%initialize(this%configFile, sg, Y=Y)     ! Initialize the B matrix

    ! ----
    ! CALL this%Ctl2State%transFwdNonLinear(X)

    ! DO j = LBOUND(X%fields, 1), UBOUND(X%fields, 1)
    !   DO i = LBOUND(Y%ObsFields, 1), UBOUND(Y%ObsFields, 1)
    !     IF ((TRIM(Y%ObsFields(i)%Get_Name())) .EQ. (TRIM(X%fields(j)%Get_Name()))) THEN
    !       DO k = LBOUND(Y%ObsFields(i)%idx, 1), UBOUND(Y%ObsFields(i)%idx, 1)
    !         CALL X%fields(j)%Set_Value(Y%ObsFields(i)%idx(k), Y%ObsFields(i)%values(k))
    !       END DO
    !     END IF
    !   END DO
    ! END DO

    ! CALL this%Ctl2State%transBackward(X)
    ! ----

    IF (TRIM(this%mode) == 'Debug') THEN
      BLOCK
        TYPE(State_t) :: YY, DD, ZZ, BB

        YY = Obs2State_BaseTypeName(sg, Y)
        DD = Obs2State_BaseTypeName(sg, Y - HY%transFwdNonLinear_opr(this%XbMG(sg%gLevel)))
        WRITE (*, 15) X%sg%glevel, X%sg%mpddInfo_sg%myrank
15      FORMAT('DD is calculated...at G', I2, ' pc ', I2)
        BB = Obs2State_BaseTypeName(sg, HY%transFwdNonLinear_opr(this%XbMG(sg%gLevel)))

        CALL Output_NC_State_AV(YY, this%ncOutputFile, &
                                TRIM(this%task)//"_obsThinned", .TRUE., .TRUE.)
        WRITE (*, 18) X%sg%glevel, X%sg%mpddInfo_sg%myrank
18      FORMAT('Thinned is plotted...at G', I2, ' pc ', I2)
        CALL Output_NC_State_AV(DD, this%ncOutputFile, &
                                TRIM(this%task)//"_innovation", .TRUE., .TRUE.)
        WRITE (*, 28) X%sg%glevel, X%sg%mpddInfo_sg%myrank
28      FORMAT('Innovation at G', I2, ' pc ', I2)
        CALL Output_NC_State_AV(BB, this%ncOutputFile, &
                                TRIM(this%task)//"_bcgAtObs", .TRUE., .TRUE.)

        CALL Output_NC_State_AV(this%XbMG(sg%gLevel), this%ncOutputFile, &
                                TRIM(this%task)//"_bakInDA", .TRUE., .TRUE.)

        ! Plot fillin data:
        ZZ = Obs2State_BaseTypeName(sg, Z)
        CALL Output_NC_State_AV(ZZ, this%ncOutputFile, &
                                TRIM(this%task)//"_obsFilled", .TRUE., .TRUE.)
        WRITE (*, 38) X%sg%glevel, sg%gLevel, X%sg%mpddInfo_sg%myrank
38      FORMAT('Filled is plotted...at G', 2I2, ' pc ', I2)
      END BLOCK
    END IF

    WRITE (*, 16) X%sg%glevel, X%mpddGlob%myrank
16  FORMAT('Initializing R......at G', I2, ' pc ', I2)
    ! CALL R%initialize(this%configFile, Z, X%sg, Nocoarest) ! Initialize R
    CALL R%initialize(this%configFile, Z, X%sg) ! Initialize R

    CALL sg%mpddInfo_sg%barrier

    ! Initialize J Function
    WRITE (*, 17) X%sg%glevel, X%mpddGlob%myrank
17  FORMAT('Initializing JFunc_HybInc......at G', I2, ' pc ', I2)

    ! UpdateBkgd means the finer grid background is updated with its coarser grid analysis increment:
    IF (this%UpdateBkgd) THEN
      XbRef = X
    ELSE
      XbRef = this%XbMG(sg%gLevel)
    END IF

    SELECT CASE (TRIM(this%framework))
    CASE ('FullState')
      XbCV = this%XbMG(X%sg%gLevel)
    CASE ('Incremental')
      XbCV = this%XbRef(X%sg%gLevel)
      CALL this%Ctl2State%transBackward(XbCV)
      CALL this%Ctl2State%SumBothCtlState(XbCV, this%XbMG(X%sg%gLevel))
    END SELECT
    
    IF (this%Hybrid) THEN
      CALL JFunc_HybInc%initialize(this%configFile, X, Z, H, B, R, sg, Xb=XbCV, B_e=B_e)
    ELSE IF (this%HybInc) THEN
      CALL JFunc_HybInc%initialize(this%configFile, X, Z, H, B, R, sg, Xb=XbCV, B_e=B_einc, Xb_enperts=this%Xb_enperts(sg%gLevel, :))
    ELSE
      CALL JFunc_HybInc%initialize(this%configFile, X, Z, H, B, R, sg, Xb=XbCV)
    END IF

    CALL sg%mpddInfo_sg%barrier

    ! Run minimization
    CALL CPU_TIME(t1)
    IF (sg%gLevel <= 5) THEN
      CALL this%miniSolver%run_HyyInc(X, JFunc_HybInc, sg, 100)
      ! CALL this%miniSolver%run(X, JFunc_HybInc, sg)
    ELSE
      CALL this%miniSolver%run_HyyInc(X, JFunc_HybInc, sg)
    END IF
    CALL CPU_TIME(t2)
    
    WRITE (*, 3) t2 - t1, sg%mpddInfo_sg%myrank, sg%num_cell
3   FORMAT('Time spent in miniSolver: ', D10.2, ' at proc: ', I1, ' with Ncell: ', I8)

    IF (TRIM(this%mode) == 'Debug') THEN
      BLOCK
        TYPE(State_t) :: AA, ANA

        AA = Obs2State_BaseTypeName(sg, Y - H%transFwdNonLinear_opr(X))
        ANA = Obs2State_BaseTypeName(sg, H%transFwdNonLinear_opr(X))

        CALL Output_NC_State_AV(AA, this%ncOutputFile, &
                                TRIM(this%task)//"_AMO", .true., .TRUE.)
        CALL Output_NC_State_AV(ANA, this%ncOutputFile, &
                                TRIM(this%task)//"_AnaAtObs", .true., .TRUE.)

      END BLOCK

    !  BLOCK
    !    USE RTTOV_diag_out_m
    !    TYPE(State_t) :: RR
    !    INTEGER :: iv
 
    !    PRINT *, '-------------RTTOV_diag_out START-------------'
    !    CALL write_diag_vars(this%configFile, XbRef, 'fy4_1', 'agri', Y, RR)
    !    CALL Output_NC_State_AV(RR, this%ncOutputFile, "fy4_1-agri_diag", .true. , .true.)
    !    PRINT *, '-------------RTTOV_diag_out OVER-------------'
    !  END BLOCK
    END IF

  END SUBROUTINE DASolve

  !> @brief
  !================================================================
  !! This routine checks the obsSet Y to fill in background in area
  !! where no observation associated with the background field.
  SUBROUTINE bkgdFillin(this, glvl, Y, Z)
    IMPLICIT NONE
    CLASS(Applications_HybInc_t) :: this
    INTEGER(i_kind), INTENT(IN) :: glvl
    TYPE(ObsSet_t), INTENT(IN) :: Y
    TYPE(ObsSet_t), INTENT(INOUT) :: Z

    ! Local variables:
    INTEGER(i_kind) :: i, j, k, iu, iv, ix, nin, nv, nh, nt, imax, imin, jmax, jmin, kmax, kmin, nfill, idd
    REAL(i_kind) :: a
    TYPE(State_t) :: X, S

    INTEGER(i_kind), ALLOCATABLE :: idc(:, :)
    REAL(r_kind), ALLOCATABLE :: O(:)

    X = this%XbMG(glvl)%zeroCopy()
    CALL this%Ctl2State%transFwdNonLinear(X)
    IF (TRIM(this%framework) .EQ. 'Incremental') THEN
      CALL this%Ctl2State%SumState(X, this%XbRef(glvl))
    END IF

    ! Check the rampRange:
    WRITE (*, 1) X%sg%gLevel, X%sg%vLevel, this%iRange_ramp(1, :), X%sg%mpddInfo_sg%myrank
    WRITE (*, 2) X%sg%gLevel, X%sg%num_cell, this%iRange_ramp(2, :), X%sg%mpddInfo_sg%myrank
    WRITE (*, 3) X%sg%gLevel, X%sg%tSlots, this%iRange_ramp(3, :), X%sg%mpddInfo_sg%myrank
1   FORMAT('BkgcFillin - G', I2, ' Vertical  levels:', I6, ' ramp: ', 20I6)
2   FORMAT('BkgcFillin - G', I2, ' Horizon num_cell:', I6, ' ramp: ', 20I6)
3   FORMAT('BkgcFillin - G', I2, ' Time Framework slots:', I6, ' ramp: ', 20I6)

    IF (MINVAL(this%iRange_ramp(1, :)) .GT. X%sg%vLevel .AND. &
        MINVAL(this%iRange_ramp(2, :)) .GT. X%sg%num_cell .AND. &
        MINVAL(this%iRange_ramp(3, :)) .GT. X%sg%tSlots) THEN
      Z = Y
      WRITE (*, 4) X%sg%gLevel, X%sg%mpddInfo_sg%myrank
4     FORMAT('No background filled in Glevel:', I3, ' pc', I2)
      RETURN ! No fill in
    END IF

    ! distributing Y obs to X:
    DO i = 1, UBOUND(Y%ObsFields, 1)
      ! Associated obs:
      IF (TRIM(Y%ObsFields(i)%Get_Name()) .EQ. 'cdir' .OR. &
          TRIM(Y%ObsFields(i)%Get_Name()) .EQ. 'sdir' .OR. &
          TRIM(Y%ObsFields(i)%Get_Name()) .EQ. 'wspd' .OR. &
          TRIM(Y%ObsFields(i)%Get_Name()) .EQ. 'vel') THEN
        iu = X%getVarIdx('uwnd')
        iv = X%getVarIdx('vwnd')
        DO j = 1, UBOUND(Y%ObsFields(i)%idx, 1)
          X%fields(iu)%DATA(Y%ObsFields(i)%idx(j)%vIdx, &
                            Y%ObsFields(i)%idx(j)%hIdx, &
                            Y%ObsFields(i)%idx(j)%tIdx) = 1.0D0
          X%fields(iv)%DATA(Y%ObsFields(i)%idx(j)%vIdx, &
                            Y%ObsFields(i)%idx(j)%hIdx, &
                            Y%ObsFields(i)%idx(j)%tIdx) = 1.0D0
        END DO
      ELSE
        ix = X%getVarIdx(TRIM(Y%ObsFields(i)%Get_Name()))
        IF (ix .GT. 0) THEN ! Found model state var
          DO j = 1, UBOUND(Y%ObsFields(i)%idx, 1)
            X%fields(ix)%DATA(Y%ObsFields(i)%idx(j)%vIdx, &
                              Y%ObsFields(i)%idx(j)%hIdx, &
                              Y%ObsFields(i)%idx(j)%tIdx) = 1.0D0
          END DO
        END IF
      END IF
    END DO

    ! Exchange the halo values to make sure each process see all possible obs:
    CALL X%exHalo()

    ! Spread the markers to fill the ramp regions of the horizontal grid:
    DO iv = 1, UBOUND(X%fields, 1)
      DO j = 1, this%iRange_ramp(2, iv)
        ! These markers values may over-count because of the neighbors overlaps
        CALL X%fields(iv)%Sum_Neighbors()
        ! After each ramp spread of the obs, the halo values need updated:
        ! CALL X%exhalo()
        CALL X%sg%ExchangeMatOnHaloForFieldGrid(X%sg%tSlots, X%sg%vLevel, &
                                                X%Fields(iv)%DATA)
      END DO
    END DO

    ! Fill in:
    nin = UBOUND(Y%obsFields, 1)
    Z = Y   ! set Z to Y as default and add background to the existing states
    S = this%XbMG(glvl)
    CALL this%Ctl2State%transFwdNonLinear(S)
    IF (TRIM(this%framework) .EQ. 'Incremental') THEN
      CALL this%Ctl2State%SumState(S, this%XbRef(glvl))
    END IF

    ! Allocate O, and idc: setting them to the maximum possible fillins:
    ALLOCATE (O(UBOUND(X%fields(1)%DATA, 1) * UBOUND(X%fields(1)%DATA, 2) * UBOUND(X%fields(1)%DATA, 3)), &
              idc(3, UBOUND(X%fields(1)%DATA, 1) * UBOUND(X%fields(1)%DATA, 2) * UBOUND(X%fields(1)%DATA, 3)))

    idd = 0
    DO iv = 1, UBOUND(X%fields, 1)
      nfill = 0
      O = 0.0D0
      idc = 0
      nv = UBOUND(X%fields(iv)%DATA, 1)
      nh = X%sg%num_icell ! Only on the interior points to consider fillin
      nt = UBOUND(X%fields(iv)%DATA, 3)
      DO i = 1, 1 ! nv
        imin = 1 ! MAX(1, i - this%iRange_ramp(1, iv))
        imax = 1 !MIN(nv, i + this%iRange_ramp(1, iv))
        DO j = 1, nh
          jmin = MAX(1, j - this%iRange_ramp(2, iv))
          jmax = MIN(nh, j + this%iRange_ramp(2, iv))
          DO k = 1, nt
            kmin = MAX(1, k - this%iRange_ramp(3, iv))
            kmax = MIN(nt, k + this%iRange_ramp(3, iv))
            A = SUM(X%fields(iv)%DATA(imin:imax, j, kmin:kmax))

#ifdef TRACK_Fillins
            ! Debugging 2022-9-03:
            ! j values are at the gridpoint on different processor runs: e.g. 1 or 4 processors
            IF (iv .EQ. 1 .AND. j .EQ. 17 * 16 + 15 .OR. j .EQ. 34 * 16 + 32) THEN
              WRITE (*, 12) i, j, k, a, X%sg%cell_cntr(:, j), X%sg%mpddInfo_sg%myrank
12            FORMAT('Check the A values: ', 3I4, ' A: ', D10.2, ' LL', 2D15.8, ' pc:', I2)
            END IF
#endif

            IF (A .LT. 1.0D0) THEN
              ! No obs found in this ramp ranges:
              nfill = nfill + 1
              O(nfill) = S%fields(iv)%DATA(i, j, k)
              idc(1, nfill) = i
              idc(2, nfill) = j
              idc(3, nfill) = k
            END IF
          END DO
        END DO
      END DO

      ! Add or create Z ObsField:
      nv = 0
      iu = Y%getObsIdx(TRIM(X%fields(iv)%Get_Name()))
      IF (iu .GT. 0) THEN
        ! Y contains direct obs:
        nv = UBOUND(Y%ObsFields(iu)%values, 1)

        ! Release the data of Z and allocate addition memory for backround fillin:
        IF (ALLOCATED(Z%ObsFields(iu)%values)) DEALLOCATE (Z%ObsFields(iu)%values)
        IF (ALLOCATED(Z%ObsFields(iu)%idx)) DEALLOCATE (Z%ObsFields(iu)%idx)
        ALLOCATE (Z%ObsFields(iu)%values(nv + nfill), Z%ObsFields(iu)%idx(nv + nfill))
        Z%ObsFields(iu)%values(1:nv) = Y%ObsFields(iu)%values(1:nv)
        Z%ObsFields(iu)%values(nv + 1:nv + nfill) = O(1:nfill)
        Z%ObsFields(iu)%idx(1:nv) = Y%ObsFields(iu)%idx(1:nv)
        Z%ObsFields(iu)%idx(nv + 1:nv + nfill)%vIdx = idc(1, 1:nfill)
        Z%ObsFields(iu)%idx(nv + 1:nv + nfill)%hIdx = idc(2, 1:nfill)
        Z%ObsFields(iu)%idx(nv + 1:nv + nfill)%tIdx = idc(3, 1:nfill)

        CALL Z%ObsFields(iu)%Set_Name(TRIM(Y%ObsFields(iu)%Get_Name()))
        CALL Z%ObsFields(iu)%Set_ObsType(TRIM(Y%ObsFields(iu)%Get_ObsType()))
        Z%ObsFields(iu)%locObs = Y%ObsFields(iu)%locObs
      ELSE
        ! idd = idd + 1 ! One more var is additional to obs:
        ! ALLOCATE (Z%ObsFields(nin + idd)%values(nfill), Z%ObsFields(nin + idd)%idx(nfill))
        ! Z%ObsFields(nin + idd)%values(1:nfill) = O(1:nfill)
        ! Z%ObsFields(nin + idd)%idx(1:nfill)%vIdx = idc(1, 1:nfill)
        ! Z%ObsFields(nin + idd)%idx(1:nfill)%hIdx = idc(2, 1:nfill)
        ! Z%ObsFields(nin + idd)%idx(1:nfill)%tIdx = idc(3, 1:nfill)

        ! CALL Z%ObsFields(nin + idd)%Set_Name(TRIM(X%fields(iv)%Get_Name()))
        ! CALL Z%ObsFields(nin + idd)%Set_ObsType('bkgd')
        ! Z%ObsFields(nin + idd)%locObs(1:2) = X%sg%cell_cntr(1:2, 1)
        ! Z%ObsFields(nin + idd)%locObs(3) = 0.0D0 ! Temporarily set
        WRITE (*, 5) TRIM(X%fields(iv)%Get_Name()), X%mpddSub%myrank
5       FORMAT('No observation for this state variable: ', A, 1X, ' No need to fillin, pc: ', I4)
      END IF
    END DO

    ! Deallocate local variables:
    DEALLOCATE (O, idc)

  END SUBROUTINE bkgdFillin

  !> @brief
  !================================================================
   !!  This routine performs a MOTOR-DA analysis: using thinning from
   !!  the finest grid thinned obs to improve efficiency.
   !!
   !!  Input:
   !!       mgStart:    integer G level to start a multigrid analysis
   !!       mgEnd:      integer G level to end the multigrid analysis
  !
  SUBROUTINE analyze(this, obsList, mgStart, mgEnd)
    IMPLICIT NONE
    CLASS(Applications_HybInc_t) :: this
    CHARACTER(LEN=*), INTENT(IN) :: obsList(:)
    INTEGER(i_kind), INTENT(IN) :: mgStart, mgEnd

    ! Local variables:
    TYPE(State_t) :: XbRef
    CHARACTER(LEN=20) :: varName
    INTEGER(i_kind) :: i, j, iv, k
    LOGICAL :: IncFlag = .TRUE.
    REAL(r_kind) :: ana_max, ana_min, ana_max_global, ana_min_global
    REAL(r_kind), ALLOCATABLE :: ana_max_AV(:), ana_min_AV(:)
    TYPE(State_t) :: XTru
    TYPE(ObsSet_t) :: yObs(7) ! Hard coded for possible observation data types
    REAL(r_kind) :: tt1, tt2

    ! Calculate total time used in analysis:
    CALL CPU_TIME(tt1)

    ! General thinned observations for each analysis level:
    ALLOCATE (this%thinnedObs(mgStart:mgEnd), this%mpObs(mgStart:mgEnd))

    ! Initialize the observation parallel processing proc at the finest level grid:
    CALL this%mpObs(mgEnd)%initializeMPObs(this%geometry%mg%sg(mgEnd))

    ! XbRef = this%Ctl2State%fwdNL_opr(this%XbMG(mgEnd))
    XbRef = this%XbRef(mgEnd)

    ! Thinning the observations at the finest level:
    DO i = 1, UBOUND(obsList, 1)
      SELECT CASE (TRIM(obsList(i)))
      CASE ('surfaces')
        CALL this%sfc%ObsThinning(XbRef, yObs(1), this%mpObs(mgEnd), .TRUE., .FALSE.)
      CASE ('sounding')
        CALL this%snd%ObsThinning(XbRef, yObs(2), this%mpObs(mgEnd), .TRUE., .FALSE.)
      CASE ('profiler')
        CALL this%pro%ObsThinning(XbRef, yObs(3), this%mpObs(mgEnd), .TRUE., .FALSE.)
      CASE ('radarRef')
        CALL this%ref%ObsThinning(XbRef, yObs(4), this%mpObs(mgEnd), .TRUE., .FALSE.)
      CASE ('radarVel')
        CALL this%vel%ObsPrepareForSg(XbRef)
        CALL this%vel%ObsThinning(XbRef, yObs(5), this%mpObs(mgEnd), .TRUE., .FALSE.)
      END SELECT
    END DO
    ! Concat thinning obs:
    this%thinnedObs(mgEnd) = ObsSet_t(this%configfile, this%mpobs(mgEnd))
    CALL ObsConcat_s((/yObs(1), yObs(2), yObs(3), yObs(4), yObs(5), yObs(6), yObs(7)/), this%thinnedObs(mgEnd))

    ! General all other thinned observations from the finest:
    BLOCK
      USE Obs2State_m
      TYPE(State_t) :: DD
      INTEGER(i_kind) :: ia
      CHARACTER(len=20) :: varList

      ! Plot the thinned obs at the finest level if requested:
      IF (TRIM(this%mode) == 'Debug') THEN
        PRINT *, 'output thinned obs... ', this%XbMG(mgEnd)%sg%gLevel, this%XbMG(mgEnd)%sg%mpddInfo_sg%myrank
        DD = Obs2State_BaseTypeName(this%XbMG(mgEnd)%sg, this%thinnedObs(mgEnd))
        CALL Output_NC_State_AV(DD, this%ncOutputFile, &
                                TRIM(this%task)//"_obsThinned", .TRUE., .TRUE.)
        WRITE (*, 2) this%XbMG(mgEnd)%sg%gLevel, this%XbMG(mgEnd)%sg%mpddInfo_sg%myrank
2       FORMAT('Analyze - thinned obs output at G', I2, ' proc: ', I2)
      END IF

      DO i = mgEnd - 1, mgStart, -1 ! From the finest to the coarsest
        ! Use sfc function to handle all obs, this is not only for surface data:
        ASSOCIATE (sg => this%geometry%mg%sg(i))

          CALL this%mpObs(i)%initializeMPObs(this%geometry%mg%sg(i))
          ! modified by TS@230131, move griddedObsMap subroutine from obsBase_m to a new Module (obsTools_m)
          ! CALL this%sfc%griddedObsMap(this%geometry%mg%sg(i + 1), &
          !                            this%geometry%mg%sg(i), this%thinnedObs(i + 1), &
          !                            this%thinnedObs(i), this%mpObs(i))
          CALL griddedObsMap(this%sfc, this%geometry%mg%sg(i + 1), this%geometry%mg%sg(i), this%thinnedObs(i + 1), this%thinnedObs(i), this%mpObs(i))

          DO j = 1, UBOUND(this%thinnedObs(i)%ObsFields, 1)
            varList = TRIM(this%thinnedObs(i)%ObsFields(j)%Get_Id_Name())
            SELECT CASE (varList)
            CASE ('SYNOP_pcpa')
              IF (sg%gLevel < mgEnd - 3) this%thinnedObs(i)%ObsFields(j)%values = 0.0D0
            CASE ('SOUND_qvapor')
              IF (sg%gLevel < mgEnd - 1) CALL this%thinnedObs(i)%ObsFields(j)%Set_Name('qvapor_d')
            CASE ('fy4_1_agri_ch_2_tbb')
              IF (sg%gLevel < mgEnd - 1) CALL this%thinnedObs(i)%ObsFields(j)%Set_Name('tdbb')
            CASE ('fy4_1_agri_ch_3_tbb')
              IF (sg%gLevel < mgEnd - 1) CALL this%thinnedObs(i)%ObsFields(j)%Set_Name('tdbb')
            END SELECT
            ! PRINT *,'varList(i): ', this%thinnedObs(i)%ObsFields(j)%Get_Id_Name()
          END DO
        END ASSOCIATE
        IF (TRIM(this%mode) == 'Debug') THEN
          DD = Obs2State_BaseTypeName(this%XbMG(i)%sg, this%thinnedObs(i))
          CALL Output_NC_State_AV(DD, this%ncOutputFile, &
                                  TRIM(this%task)//"_obsThinned", .TRUE., .TRUE.)
          WRITE (*, 2) this%XbMG(i)%sg%gLevel, this%XbMG(i)%sg%mpddInfo_sg%myrank
        END IF
      END DO
    END BLOCK

    ! Give the initial value at the coarsest grid.
    this%analysis = this%XbMG(mgStart)

    DO k = 1, 1
      DO i = mgStart, mgEnd
        ! Run 3DVAR in each single grid.
        BLOCK
          TYPE(State_t) :: X
          ASSOCIATE (sg => this%geometry%mg%sg(i))
            X = this%analysis

            ! Run the minimization
            CALL this%varAnal(X, obsList, i)

            CALL this%mpddGlob%barrier

            IF (TRIM(this%mode) == 'Debug') THEN
              BLOCK
                ! USE UV2W_m, ONLY: UV2W_t
                ! TYPE(UV2W_t) :: UV2W
                TYPE(state_t) :: tempX

                ! UV2W = UV2W_t(this%configFile, this%XbMG(i))
                ! CALL UV2W%transFwdNonLinear(this%XbMG(i))
                ! CALL UV2W%transFwdNonLinear(X)

                tempX = X - this%XbMG(i)
                CALL this%Ctl2State%transFwdNonLinear(tempX)
                CALL Output_NC_State_AV(tempX, this%ncOutputFile, &
                                        TRIM(this%task)//"_diff", .TRUE., .TRUE.)

                tempX = X
                CALL this%Ctl2State%transFwdNonLinear(tempX)
                CALL Output_NC_State_AV(tempX, this%ncOutputFile, &
                                        TRIM(this%task)//"_ana", .TRUE., .TRUE.)

                CALL this%XbMG(i)%rmVar('wwnd')
                CALL X%rmVar('wwnd')
              END BLOCK
            END IF

            ! Connect to coarser grid.
            IF (i .NE. mgEnd) THEN
              ! ===========================> prolongate to finer grid
              ASSOCIATE (sgFiner => this%geometry%mg%sg(i + 1))
                ! CALL this%analysis%initialize(this%configFile, sgFiner)
                CALL prolongationMGInc(X, this%analysis, this%XbMG(i), this%XbMG(i + 1), sg, sgFiner, this%geometry%mg)
              END ASSOCIATE
            ELSE
              this%analysis = X ! Return the state fields directly

            END IF

            IF (sg%isActiveProc()) THEN
              IF (TRIM(this%mode) == 'Debug') THEN
                ALLOCATE (ana_max_AV(this%numVar), ana_min_AV(this%numVar))
                DO iv = 1, this%numVar
                  ana_max_AV(iv) = MAXVAL(this%analysis%fields(iv)%DATA)
                  ana_min_AV(iv) = MINVAL(this%analysis%fields(iv)%DATA)
                END DO
                ana_max = MAXVAL(ana_max_AV)
                ana_min = MINVAL(ana_min_AV)
                DEALLOCATE (ana_max_AV, ana_min_AV)

                CALL sg%mpddInfo_sg%AllReduceMaxReal(ana_max, ana_max_global)
                CALL sg%mpddInfo_sg%AllReduceMinReal(ana_min, ana_min_global)
                IF (this%mpddGlob%isBaseProc()) &
                  WRITE (*, 1) ana_max_global, ana_min_global, i, this%analysis%sg%vLevel
1               FORMAT('max/min of ana : ', 2D12.4, ' at Grid lvl: ', I2, ' with vLevels: ', I3)
              END IF
            END IF

          END ASSOCIATE
        END BLOCK
      END DO
    END DO

    ! Deallocate unused arrays:
    DEALLOCATE (this%thinnedObs)

    CALL CPU_TIME(tt2)
    WRITE (*, 4) tt2 - tt1, this%mpddGlob%myrank
4   FORMAT('analyze: time spent in analysis total: ', D12.4, ' at proc: ', I3)
  END SUBROUTINE analyze

  !> @brief
  !================================================================
  !  This routine performs a MOTOR-DA analysis:
  !  Input:
  !       mgStart:    integer G level to start a multigrid analysis
  !       mgEnd:      integer G level to end the multigrid analysis
  !
  SUBROUTINE analyss(this, obsList, mgStart, mgEnd)
    IMPLICIT NONE
    CLASS(Applications_HybInc_t) :: this
    CHARACTER(LEN=*), INTENT(IN) :: obsList(:)
    INTEGER(i_kind), INTENT(IN) :: mgStart, mgEnd

    ! Local variables:
    CHARACTER(LEN=20) :: varName
    INTEGER(i_kind) :: i, iv, k
    LOGICAL :: IncFlag = .TRUE.
    REAL(r_kind) :: ana_max, ana_min, ana_max_global, ana_min_global
    REAL(r_kind), ALLOCATABLE :: ana_max_AV(:), ana_min_AV(:)
    TYPE(State_t) :: XTru
    REAL(r_kind) :: tt1, tt2
    INTEGER(i_kind) :: Nocoarest(100) ! max is 100 types of obs

    ! Calculate total time used in analysis:
    CALL CPU_TIME(tt1)

    ! Give the initial value at the coarsest grid.
    this%analysis = this%XbMG(mgStart)

    DO k = 1, 1
      DO i = mgStart, mgEnd
        ! Run 3DVAR in each single grid.
        IF (i > mgStart) THEN
          CALL this%mpddGlob%bcast(Nocoarest)
        ELSE
          Nocoarest = 0
        END IF

        BLOCK
          TYPE(State_t) :: X
          TYPE(State_t) :: Xb
          ASSOCIATE (sg => this%geometry%mg%sg(i))
            X = this%analysis

            IF (X%getVarIdx('qcloud') .GT. 0) THEN
              WHERE(X%fields(X%getVarIdx('qcloud'))%data <= 1.0D-12) X%fields(X%getVarIdx('qcloud'))%data = 0.0D0
              WHERE(X%fields(X%getVarIdx('qice'))%data <= 1.0D-12) X%fields(X%getVarIdx('qice'))%data = 0.0D0
            END IF
            IF (X%getVarIdx('qcloud_ctl') .GT. 0) THEN
              WHERE(X%fields(X%getVarIdx('qcloud_ctl'))%data <= 1.0D-12) X%fields(X%getVarIdx('qcloud_ctl'))%data = 0.0D0
              WHERE(X%fields(X%getVarIdx('qice_ctl'))%data <= 1.0D-12) X%fields(X%getVarIdx('qice_ctl'))%data = 0.0D0
            END IF

            ! Run the minimization
            CALL this%DASolve(X, obsList, sg, Nocoarest)

            IF (TRIM(this%mode) == 'Debug') THEN
              BLOCK
                ! USE UV2W_m, ONLY: UV2W_t
                ! TYPE(UV2W_t) :: UV2W
                TYPE(state_t) :: tempX
                INTEGER :: j

                TYPE(state_t) :: DF

                ! UV2W = UV2W_t(this%configFile, this%XbMG(i))
                ! CALL UV2W%transFwdNonLinear(this%XbMG(i))
                ! CALL UV2W%transFwdNonLinear(X)

                DF = X - this%XbMG(i)
                CALL Output_NC_State_AV(X - this%XbMG(i), this%ncOutputFile, &
                                        TRIM(this%task)//"_diff", .TRUE., .TRUE.)

                tempX = X
                CALL this%Ctl2State%transFwdNonLinear(tempX) ! output in CV or BKG space
                IF (TRIM(this%framework) .EQ. 'Incremental') THEN
                  CALL this%Ctl2State%SumState(tempX, this%XbRef(i))
                END IF

                CALL Output_NC_State_AV(tempX, this%ncOutputFile, &
                                        TRIM(this%task)//"_ana", .TRUE., .TRUE.)

                ! Debugging the landmask: ! Yuanfu Xie added for plotting 2D horizonal statics on 2023/10/16
                ! tempX%fields(1)%DATA(1, :, 1) = sg%landmask
                ! tempX%fields(2)%DATA(1, :, 1) = sg%topo
                ! tempX%fields(3)%DATA(1, :, 1) = sg%horizonSimilarity  ! Yuanfu Xie added on 2023/10/18 for plotting
                ! CALL this%Ctl2State%transFwdNonLinear(tempX)
                ! CALL Output_NC_State_AV(tempX, this%ncOutputFile, &
                !                         TRIM(this%task)//"_mask", .TRUE., .TRUE.)

                ! CALL this%XbMG(i)%rmVar('wwnd')
                ! CALL X%rmVar('wwnd')
              END BLOCK
            END IF

            ! Connect to coarser grid.
            IF (i .NE. mgEnd) THEN
              ! ===========================> prolongate to finer grid
              ASSOCIATE (sgFiner => this%geometry%mg%sg(i + 1))
                this%analysis = this%XbMG(i + 1)%zeroCopy()

                ! CALL prolongationMG(X, this%analysis, this%XbMG(i + 1), this%geometry%mg, IncFlag)
                CALL prolongationMGInc(X, this%analysis, this%XbMG(i), this%XbMG(i + 1), sg, sgFiner, this%geometry%mg)

                BLOCK
                  INTEGER(i_kind) :: ii, jj, kk
                  DO ii = 1, SIZE(this%analysis%Fields)
                    DO jj = 1, sgFiner%num_cell
                      DO kk = 1, sgFiner%vLevel
                        IF (sgFiner%zHght(kk, jj) < sgFiner%topo(jj)) THEN
                          this%analysis%Fields(ii)%DATA(kk, jj, :) = 0.0D0
                        END IF
                      END DO
                    END DO
                  END DO
                END BLOCK

              END ASSOCIATE
            ELSE
              this%analysis = X ! Return the state fields directly
            END IF

            IF (TRIM(this%mode) == 'Debug') THEN
              IF (sg%isActiveProc()) THEN
                IF (TRIM(this%mode) == 'Debug') THEN
                  ALLOCATE (ana_max_AV(this%numVar), ana_min_AV(this%numVar))
                  DO iv = 1, this%numVar
                    ana_max_AV(iv) = MAXVAL(this%analysis%fields(iv)%DATA)
                    ana_min_AV(iv) = MINVAL(this%analysis%fields(iv)%DATA)
                  END DO
                  ana_max = MAXVAL(ana_max_AV)
                  ana_min = MINVAL(ana_min_AV)
                  DEALLOCATE (ana_max_AV, ana_min_AV)
                  ! PRINT *, 'ISSUE DEBUG --- MAX/MIN of ANA:', sg%mpddInfo_sg%myrank, ana_max, ana_min
                  CALL sg%mpddInfo_sg%AllReduceMaxReal(ana_max, ana_max_global)
                  CALL sg%mpddInfo_sg%AllReduceMinReal(ana_min, ana_min_global)
                  IF (this%mpddGlob%isBaseProc()) &
                    WRITE (*, 1) ana_max_global, ana_min_global, i, this%analysis%sg%vLevel
1                 FORMAT('max/min of ana : ', 2D12.4, ' at Grid lvl: ', I2, ' with vLevels: ', I3)
                END IF
              END IF

            END IF
          END ASSOCIATE
        END BLOCK

      END DO
    END DO

    CALL CPU_TIME(tt2)
    WRITE (*, 4) tt2 - tt1, this%mpddGlob%myrank
4   FORMAT('analysis: time spent in analysis total: ', D12.4, ' at proc: ', I3)
  END SUBROUTINE analyss

  !> @brief
  !! Validate the analysis:
  SUBROUTINE verifyA(this, obsList, mgEnd, unitTest)
    IMPLICIT NONE
    CLASS(Applications_HybInc_t) :: this
    CHARACTER(LEN=8), INTENT(IN) :: obsList(:)
    INTEGER(i_kind), INTENT(IN) :: mgEnd
    LOGICAL, INTENT(IN) :: unitTest

    ! Local variables:
    CHARACTER(LEN=20) :: varNames(20)
    INTEGER(i_kind) :: iv, j, idx_max(20), numVars
    REAL(r_kind) :: ana_max, ana_min, ana_max_global
    REAL(r_kind), ALLOCATABLE :: ana_max_AV(:), ana_min_AV(:), xmm(:, :)
    REAL(r_kind), PARAMETER :: passthreshold_max = 6.0D0, passthreshold_min = 0.6D0

    INTEGER(i_kind) :: numObsVars

    numObsVars = 50 ! Temporarily hard coded for handling model state variables are not identical to obs variables

    ALLOCATE (ana_max_AV(numObsVars), ana_min_AV(numObsVars), xmm(2, numObsVars))

    ! Verify the analysis over the verification area:
    ana_max_AV = 0.0D0
    DO j = 1, UBOUND(obsList, 1)
      SELECT CASE (obsList(j))
        ! Modified by TS@230131
      CASE ('surfaces')
        ! CALL this%sfc%ObsMinusState(this%analysis, this%verifyLatlon, numVars, ana_max_AV, idx_max, varNames, xmm)
        CALL ObsMinusState(this%sfc, this%analysis, this%verifyLatlon, numVars, ana_max_AV, idx_max, varNames, xmm)
      CASE ('sounding')
        ! CALL this%snd%ObsMinusState(this%analysis, this%verifyLatlon, numVars, ana_max_AV, idx_max, varNames, xmm)
        CALL ObsMinusState(this%snd, this%analysis, this%verifyLatlon, numVars, ana_max_AV, idx_max, varNames, xmm)
      CASE ('profiler')
        ! CALL this%pro%ObsMinusState(this%analysis, this%verifyLatlon, numVars, ana_max_AV, idx_max, varNames, xmm)
        CALL ObsMinusState(this%pro, this%analysis, this%verifyLatlon, numVars, ana_max_AV, idx_max, varNames, xmm)
      CASE ('radarRef')
        ! CALL this%ref%ObsMinusState(this%analysis, this%verifyLatlon, numVars, ana_max_AV, idx_max, varNames, xmm)
        CALL ObsMinusState(this%ref, this%analysis, this%verifyLatlon, numVars, ana_max_AV, idx_max, varNames, xmm)
      CASE ('radarVel')
        ! CALL this%vel%ObsMinusState(this%analysis, this%verifyLatlon, numVars, ana_max_AV, idx_max, varNames, xmm)
        CALL ObsMinusState(this%vel, this%analysis, this%verifyLatlon, numVars, ana_max_AV, idx_max, varNames, xmm)
      END SELECT
      IF (numVars .GT. numObsVars) THEN
        WRITE (*, 4) numVars, numObsVars
4       FORMAT('verifyA - the numVars: ', I2, ' in obsSet is greater than the hard coded numObsVars: ', I2, /, &
               ' Resetting the numObsVars in verifyA and rerun!')
        STOP
      END IF
      DO iv = 1, numVars
        IF (idx_max(iv) .GT. 0) THEN
!                WRITE(*,1) TRIM(obsList(j)),ana_max_AV(iv),idx_max(iv),TRIM(varNames(iv)),this%mpddGlob%myrank
! 1              FORMAT('Verify - Type: ',A,' Max ABS(O-A): ',D12.4,' at obs: ',I8,' varName: ',A,' proc: ',I2)
!                WRITE(*,2) TRIM(obsList(j)),xmm(1:2,iv),TRIM(varNames(iv)),this%mpddGlob%myrank
! 2              FORMAT('Verify - Type: ',A,' Max: ',D12.4,' Min Analysis: ',D12.4,' varName: ',A,' proc: ',I2)
        END IF
        ana_max = (ana_max_AV(iv))
        ana_min = (ana_max_AV(iv))

        CALL this%geometry%mg%sg(mgEnd)%mpddInfo_sg%AllReduceMaxReal(ana_max, ana_max_global)
        IF (this%mpddGlob%isBaseProc()) &
          WRITE (*, 6) TRIM(obsList(j)), TRIM(varNames(iv)), ana_max_global, iv
6       FORMAT('Max increment in Verify obsType:', A, ' varName: ', A, ' MaxIncrement: ', E12.4, ' VarId: ', I2)

      END DO
      CALL this%mpddGlob%barrier
    END DO
    ana_max = MAXVAL(ana_max_AV)

    CALL this%geometry%mg%sg(mgEnd)%mpddInfo_sg%AllReduceMaxReal(ana_max, ana_max_global)
    IF (this%mpddGlob%isBaseProc()) &
      PRINT *, 'GLOBAL max/min of ana:', ana_max_global, this%numVar

    DEALLOCATE (ana_max_AV, ana_min_AV, xmm)

    IF (unitTest) THEN
      IF (ana_min .GT. passthreshold_min .AND. ana_max .LT. passthreshold_max) THEN
        IF (this%mpddGlob%isBaseProc()) WRITE (*, 31)
31      FORMAT('Test passed')
      ELSE
        IF (this%mpddGlob%isBaseProc()) WRITE (*, 32)
32      FORMAT('Test failed')
      END IF
    END IF
  END SUBROUTINE verifyA

  !> @brief
  !! Destory the instance of geometry class.
  SUBROUTINE destroy(this)
    IMPLICIT NONE
    CLASS(Applications_HybInc_t) :: this

    CALL this%Geometry%destroy

    CALL this%mpddGlob%barrier

    DEALLOCATE (this%miniSolver)
    DEALLOCATE (this%XbMG)
    DEALLOCATE (this%XbRef)
    DEALLOCATE (this%platform_name, this%inst_name, this%turnOn)
    DEALLOCATE (this%ObsSatellites)
    ! Deallocate memory:
    IF (ALLOCATED(this%mgGridSmooth)) DEALLOCATE (this%mgGridSmooth)

    ! Finalize
    CALL this%mpddGlob%finalize

  END SUBROUTINE destroy

  !> @brief replace observation values by an analytic function
  !====================================================================
  !  This routine replaces observation values of a set of obs.
  !
  SUBROUTINE replace(this, numobs, numVar, latlon, obtime, sg, obdata)
    USE unitTest_sfc_m, ONLY: unitTest_sfc_t
    IMPLICIT NONE
    CLASS(Applications_HybInc_t) :: this

    INTEGER(i_kind), INTENT(IN) :: numobs, numVar, obtime(numobs)
    REAL(r_kind), INTENT(IN) :: latlon(2, numobs)
    TYPE(SingleGrid_t), INTENT(IN) :: sg
    REAL(r_kind), INTENT(OUT) :: obdata(numobs, numVar)

    ! Local variables:
    INTEGER(i_kind) :: i, j, k, ncount, mcount(sg%tSlots)
    REAL(r_kind) :: xyt(3, numobs), tmin, tmax, omin, omax
    REAL(r_kind) :: minGlobLat, maxGlobLat, minGlobLon, maxGlobLon, swap
    REAL(r_kind) :: minLat, maxLat, minLon, maxLon
    TYPE(unitTest_sfc_t) :: unitTest

    swap = MINVAL(sg%cell_cntr(1, 1:sg%num_icell))
    CALL sg%mpddInfo_sg%AllReduceMinReal(swap, minGlobLat)

    swap = MAXVAL(sg%cell_cntr(1, :))
    CALL sg%mpddInfo_sg%AllReduceMaxReal(swap, maxGlobLat)

    swap = MINVAL(sg%cell_cntr(2, :))
    CALL sg%mpddInfo_sg%AllReduceMinReal(swap, minGlobLon)

    swap = MAXVAL(sg%cell_cntr(2, :))
    CALL sg%mpddInfo_sg%AllReduceMaxReal(swap, maxGlobLon)

    WRITE (*, 51) minGlobLat / degree2radian, maxGlobLat / degree2radian, &
      minGlobLon / degree2radian, maxGlobLon / degree2radian
51  FORMAT('Analytic -- ranges of Global LL: ', 4D12.4)
    WRITE (*, 52) sg%mpddInfo_sg%myrank, &
      MINVAL(sg%cell_cntr(1, :)) / degree2radian, &
      MAXVAL(sg%cell_cntr(1, :)) / degree2radian, &
      MINVAL(sg%cell_cntr(2, :)) / degree2radian, &
      MAXVAL(sg%cell_cntr(2, :)) / degree2radian
52  FORMAT('Analytic -- At proc: ', I1, ' minmax lat: ', 2D12.4, ' lon: ', 2D12.4)
    ! Normalize the locations:
    xyt(1, :) = (latlon(1, 1:numObs) - minGlobLat) / &
                (maxGlobLat - minGlobLat)
    xyt(2, :) = (latlon(2, 1:numObs) - minGlobLon) / &
                (maxGlobLon - minGlobLon)
    xyt(3, :) = (obtime(1:numObs) - sg%tt(1)) / &
                (sg%tt(sg%tSlots) - sg%tt(1))

    ! Check obs time:
    tmin = 2.0D0
    tmax = -2.0D0
    ncount = 0
    DO i = 1, numobs
      IF (xyt(3, i) .GE. 0.0D0 .AND. tmin .GT. xyt(3, i)) tmin = xyt(3, i)
      IF (xyt(3, i) .LE. 1.0D0 .AND. tmax .LT. xyt(3, i)) tmax = xyt(3, i)
      IF (latlon(1, i) .GE. MINVAL(sg%cell_cntr(1, :)) .AND. latlon(1, i) .LE. MAXVAL(sg%cell_cntr(1, :)) .AND. &
          latlon(2, i) .GE. MINVAL(sg%cell_cntr(2, :)) .AND. latlon(2, i) .LE. MAXVAL(sg%cell_cntr(2, :)) .AND. &
          obtime(i) .GE. sg%tt(1) .AND. obtime(i) .LE. sg%tt(sg%tSlots)) THEN
        ncount = ncount + 1

        ! Output time localtion info for given processor:
        IF (sg%mpddInfo_sg%myrank .EQ. 1) WRITE (*, 99) sg%mpddInfo_sg%myrank, xyt(3, i), i
99      FORMAT('Analytic -- At proc: ', I2' scaled time: ', D12.4, ' with obs i: ', I6)
      END IF
    END DO
    WRITE (*, 58) ncount, sg%mpddInfo_sg%myrank
58  FORMAT('Analytic -- Total valid obs count = ', I5, ' at proc: ', I2)
    WRITE (*, 54) MINVAL(xyt(3, :)), MAXVAL(xyt(3, :)), tmin, tmax, sg%mpddInfo_sg%myrank, sg%tt(sg%tSlots) - sg%tt(1)
54  FORMAT('Analytic -- Obstime range: ', 2D10.3, ' validRange: ', 2D10.3, ' proc: ', I2, ' window: ', D10.3)
    DO i = 1, sg%tSlots
      ncount = 0
      DO j = 1, numobs
        IF (ABS(xyt(3, j) - (sg%tt(i) - sg%tt(1)) / (sg%tt(sg%tSlots) - sg%tt(1))) .LE. 0.5 / (sg%tt(sg%tSlots) - sg%tt(1)) .AND. &
            latlon(1, j) .GE. MINVAL(sg%cell_cntr(1, :)) .AND. latlon(1, j) .LE. MAXVAL(sg%cell_cntr(1, :)) .AND. &
            latlon(2, j) .GE. MINVAL(sg%cell_cntr(2, :)) .AND. latlon(2, j) .LE. MAXVAL(sg%cell_cntr(2, :))) THEN
          ncount = ncount + 1
        END IF
      END DO
      WRITE (*, 56) & ! ncount,i,sg%mpddInfo_sg%myrank,(sg%tt(i)-sg%tt(1))/(sg%tt(sg%tSlots)-sg%tt(1))
        COUNT(ABS(xyt(3, :) - (sg%tt(i) - sg%tt(1)) / (sg%tt(sg%tSlots) - sg%tt(1))) .LE. 0.5 / (sg%tt(sg%tSlots) - sg%tt(1)) .AND. &
              latlon(1, :) .GE. MINVAL(sg%cell_cntr(1, :)) .AND. latlon(1, :) .LE. MAXVAL(sg%cell_cntr(1, :)) .AND. &
              latlon(2, :) .GE. MINVAL(sg%cell_cntr(2, :)) .AND. latlon(2, :) .LE. MAXVAL(sg%cell_cntr(2, :))), &
        i, sg%mpddInfo_sg%myrank, (sg%tt(i) - sg%tt(1)) / (sg%tt(sg%tSlots) - sg%tt(1))
56    FORMAT('Analytic -- ObsCount: ', I5, ' at timeFrame: ', I2, ' at proc: ', I2, ' ttt: ', D12.4)
    END DO

    ! Replacing the obs values with the analytic function:
    !DO i = 1, numVar
    !   CALL unitTest%analytic(numObs, xyt, obData(:, i))
    !END DO

    ! Check a special region for observation validility:
    tmin = 2.0D0
    tmax = -2.0D0
    ncount = 0
    mcount = 0
    omax = -1.0D3
    omin = 1.0D8
    IF (sg%mpddInfo_sg%myrank .GE. 1) THEN
      DO i = 1, numobs
        IF (latlon(1, i) / degree2radian .GE. this%verifyLatlon(1) .AND. &
            latlon(1, i) / degree2radian .LE. this%verifyLatlon(2) .AND. &
            latlon(2, i) / degree2radian .GE. this%verifyLatlon(3) .AND. &
            latlon(2, i) / degree2radian .LE. this%verifyLatlon(4) .AND. xyt(3, i) .GE. 0.0D0 .AND. &
            latlon(1, i) .GE. MINVAL(sg%cell_cntr(1, :)) .AND. latlon(1, i) .LE. MAXVAL(sg%cell_cntr(1, :)) .AND. &
            latlon(2, i) .GE. MINVAL(sg%cell_cntr(2, :)) .AND. latlon(2, i) .LE. MAXVAL(sg%cell_cntr(2, :))) THEN !
          IF (obData(i, 5) .LT. 1.0D8) WRITE (*, 53) sg%mpddInfo_sg%myrank, i, obData(i, 5), xyt(:, i) ! latlon(1:2, i)/degree2radian
53        FORMAT('Analytic -- Obsvalues at proc: ', I1, ' order: ', I6, ' value: ', D10.3, ' xyt: ', 3D10.2)
          IF (tmin .GT. xyt(3, i)) tmin = xyt(3, i)
          IF (tmax .LT. xyt(3, i)) tmax = xyt(3, i)

          ! Determine time slot for this obs time:
          K = NINT(xyt(3, i) * (sg%tSlots - 1)) + 1
          mcount(K) = mcount(K) + 1

          ! Mark the region of obs with special values to identify: when needed to identify the region --
          IF (omax .LT. obdata(i, 5) .AND. obdata(i, 5) .LT. 1.0D8) omax = obdata(i, 5)
          IF (omin .GT. obdata(i, 5)) omin = obdata(i, 5)

          ! Replace the observation with an obslict value to mark on the analysis:
          obData(i, :) = 100.0D0

          ncount = ncount + 1
        END IF
      END DO
      WRITE (*, 55) tmin, tmax, ncount, sg%mpddInfo_sg%myrank, omax, omin
55    FORMAT('Analytic -- Special region obs time range: ', 2D12.4, ' numObs: ', I4, ' at proc: ', I2, ' max/min obs: ', 2D11.3)

    END IF
    DO i = 1, sg%tSlots
      IF (mcount(i) .GT. 0) WRITE (*, 88) i, mcount(i), sg%mpddInfo_sg%myrank
88    FORMAT('Analytic -- Obs in the special region -- At timeFramework: ', I2, ' there are obs: ', I5, ' at proc: ', I2)
    END DO
  END SUBROUTINE replace

END MODULE Applications_HybInc_m
