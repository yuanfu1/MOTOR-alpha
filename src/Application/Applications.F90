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
MODULE Applications_m
  USE SolverLBFGS_m, ONLY: SolverLBFGS_t
  USE SolverFRCG_m, ONLY: SolverFRCG_t
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
  USE JFunc_m, ONLY: JFunc_t
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
  USE IOGrapesModelvar_500m_m, ONLY: IOGrapesModelvar_500m_t
  USE IOGrapesPostvar_m, ONLY: IOGrapesPostvar_t
  USE IOECM_m, ONLY: IOECM_t
  USE IOERA5_m, ONLY: IOERA5_t
  USE IOModel_m, ONLY: IOModel_t

  USE EnLoc_m, ONLY: EnLoc_t

  USE ObsSurface_m, ONLY: ObsSurface_t
  USE ObsSound_m, ONLY: ObsSound_t
  USE ObsVwpw_m, ONLY: ObsVwpw_t
  USE ObsRadarRef_m, ONLY: ObsRadarRef_t
  USE ObsRadarVel_m, ONLY: ObsRadarVel_t
  USE ObsSatellite_m, ONLY: ObsSatellite_t
  USE Satob_m, ONLY: Satob_t
  USE ObsGWST_m, ONLY: ObsGWST_t
  USE ObsLBUOY_m, ONLY: ObsLBUOY_t
  USE ObsSING_m, ONLY: ObsSING_t
  USE ObsGNSS_m, ONLY: ObsGNSS_t

  USE ObsUtilities_m

  USE ObsBase_m
  USE Obs2State_m

  ! added bt TS@230131
  USE obsTools_m

  ! Define the Application type:
  TYPE Applications_t
    CHARACTER(LEN=1024) :: configFile, ncOutputFile
    CHARACTER(LEN=20), ALLOCATABLE :: varList(:)
    CHARACTER(LEN=20) :: BECsolver
    CHARACTER(LEN=20) :: mode = 'Debug' !< mode: Debug, UnitTest, Alpha, Beta, Release
    CHARACTER(LEN=20) :: task
    LOGICAL           :: Hybrid, UpdateBkgd ! Yuanfu Xie added UpdateBkgd parameter controlling the MG scheme
    ! INTEGER(i_kind) :: modeNO = 0
    INTEGER(i_kind)   :: numVar, ensNum
    INTEGER(i_kind)   :: mgStart, mgEnd
    INTEGER(i_kind), ALLOCATABLE :: mgGridSmooth(:, :)
    REAL(r_kind), ALLOCATABLE :: verifyLatlon(:)
    TYPE(mpddGlob_t) :: mpddGlob
    TYPE(geometry_t) :: geometry
    TYPE(State_t), ALLOCATABLE :: XbMG(:), XbRef(:), XbMGInitial(:)
    TYPE(State_t), ALLOCATABLE :: XbMGEns(:, :)
    TYPE(State_t), ALLOCATABLE :: Xb_enperts(:, :)
    CHARACTER(len=50), ALLOCATABLE :: platform_name(:), inst_name(:)
    LOGICAL, ALLOCATABLE :: turnOn(:)
    CHARACTER(len=20) :: framework = 'FullState'
    LOGICAL :: FixResolution = .FALSE.
    INTEGER :: OuterLoops = 1
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
    TYPE(ObsGNSS_t)     :: gnssro ! gnss ro data

    CLASS(MiniSolver_t), ALLOCATABLE :: miniSolver
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
    PROCEDURE :: destroy
    PROCEDURE :: DASolve
    PROCEDURE :: dumpAna

  END TYPE Applications_t

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
    CLASS(Applications_t) :: this
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

    IF (yaml_get_var(this%configFile, 'IO', 'output_dir', this%ncOutputFile) .NE. 0) THEN
      WRITE (*, 157)
157   FORMAT('App initial: output directory is not provided, please check your yaml and rerun!')
      STOP
    END IF
    WRITE (*, 158) TRIM(this%ncOutputFile)
158 FORMAT('App initial: output file directory: ', A)

    IF (yaml_get_var(TRIM(this%configFile), 'modelState', 'varList', this%varList) .NE. 0) THEN
      this%numVar = 0
    ELSE
      this%numVar = UBOUND(this%varList, 1)
      WRITE (*, 11) this%numVar
11    FORMAT("Number of analysis vars: ", I2)
    END IF

    IF (yaml_get_var(TRIM(this%configFile), 'geometry', 'mgStart', this%mgStart) .NE. 0) THEN
      WRITE (*, 51)
51    FORMAT('App initial: multigrid start level is required, please your yaml and provide it and rerun')
      STOP
    END IF

    IF (yaml_get_var(TRIM(this%configFile), 'geometry', 'mgEnd', this%mgEnd) .NE. 0) THEN
      WRITE (*, 52)
52    FORMAT('App initial: multigrid end level is required, please your yaml and provide it and rerun')
      STOP
    END IF
    IF (this%mpddGlob%isBaseProc()) WRITE (*, 13) this%mgStart, this%mgEnd
13  FORMAT("mgStart and mgEnd: ", 2I2)

    IF (yaml_get_var(TRIM(this%configFile), 'RunMode', 'Mode', this%mode) .NE. 0) THEN
      WRITE (*, 53)
53    FORMAT('App initial: run mode is required: please check your yaml setting RunMode/Mode and rerun')
      STOP
    END IF

    PRINT *, 'RunMode is: ', this%mode, istatus

    IF (yaml_get_var(TRIM(this%configFile), 'RunMode', 'Framework', this%framework) .NE. 0) THEN
      WRITE (*, 29)
29    FORMAT('App initial: RunMode Framework is required: please check your yaml setting RunMode/Framework and rerun')
      STOP
    END IF
    PRINT *, 'RunMode Framework is: ', this%framework, istatus

    IF (yaml_get_var(TRIM(this%configFile), 'RunMode', 'FixResolution', this%FixResolution) .NE. 0) this%FixResolution = .FALSE.
    IF (yaml_get_var(TRIM(this%configFile), 'RunMode', 'OuterLoops', this%OuterLoops) .NE. 0) this%OuterLoops = 1

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

    ALLOCATE (this%iRange_ramp(3, this%numVar))

    IF (yaml_get_var(TRIM(this%configFile), 'RunMode', 'Task', this%Task) .NE. 0) THEN
      WRITE (*, 54)
54    FORMAT('App initial: Task is a required parameter. Check your yaml for RunMode/Task and rerun')
      STOP
    END IF
    PRINT *, 'Task is: ', this%Task, istatus

    ! Prepare for the loop of different satellite instruments
    ifile = yaml_get_var(TRIM(this%configFile), 'RTTOV', 'inst_name', this%inst_name)
    ifile = yaml_get_var(TRIM(this%configFile), 'RTTOV', 'platform_name', this%platform_name)
    ifile = yaml_get_var(TRIM(this%configFile), 'RTTOV', 'turnOn', this%turnOn)

    n_plats = SIZE(this%platform_name, 1)
    n_insts = SIZE(this%inst_name, 1)
    IF (n_plats .NE. n_insts) THEN
      PRINT *, 'STOP: Numbers of platform_name and inst_name are inconsistent'
      RETURN
    ELSE
      ALLOCATE (this%ObsSatellites(n_insts))
    END IF

    ! Auxtypes
    CALL this%mpddGlob%initialize()                                   ! Initialize the mpdd

    PRINT *, 'Finish initialization of mpddGlob: ', this%mpddGlob%myrank, &
      ' with updateBkgd: ', this%UpdateBkgd

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
! 15  FORMAT('Applications/Initial - verifying area latlon: ', 4D12.4, ' time window: ', 2E12.4)

    ! Initialize geometry
    CALL this%geometry%initialize(this%configFile, this%mpddGlob)

    ! Initialize the geometry
    ALLOCATE (this%XbMG(this%geometry%mg%mg_coarsest:this%geometry%mg%mg_finest))
    ALLOCATE (this%XbRef(this%geometry%mg%mg_coarsest:this%geometry%mg%mg_finest))
    ALLOCATE (this%XbMGInitial(this%geometry%mg%mg_coarsest:this%geometry%mg%mg_finest))

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

    IF (this%ensNum .GT. 1) THEN
      ALLOCATE (this%XbMGEns(this%geometry%mg%mg_coarsest:this%geometry%mg%mg_finest, this%ensNum))
    END IF

    ! Multigrid scheme:
    istatus = yaml_get_var(TRIM(this%configFile), 'MultigridOptions', 'UpdateBkgd', this%UpdateBkgd)
    IF (istatus .NE. 0) THEN
      PRINT *, 'App initial: UpdateBkgd is not specified. Use of a default .FALSE.'
      this%UpdateBkgd = .FALSE.
    END IF
    PRINT *, 'App initial: UpdateBkgd is setted as:', this%UpdateBkgd

    ! Initialize solver
    SELECT CASE (min)
    CASE (1)
      ALLOCATE (SolverLBFGS_t::this%miniSolver)
      ! this%CALL miniSolver%initialize(this%configFile)
    CASE (2)
      ALLOCATE (SolverFRCG_t::this%miniSolver)
      ! this%CALL miniSolver%initialize(this%configFile)
    CASE DEFAULT
      ! Make sure to print warning message at last:
      CALL this%mpddGlob%barrier
      IF (this%mpddGlob%isBaseProc()) WRITE (*, 14) min
14    FORMAT('This minimization method is not implemented: ', I1, / &
             'Check your option in app%initial and rerun!')
      istatus = 1
    END SELECT

    ASSOCIATE (miniSolver => this%miniSolver)
      SELECT TYPE (miniSolver)
      TYPE IS (SolverFRCG_t)
        CALL miniSolver%initialize(configFile)
      TYPE IS (SolverLBFGS_t)
        CALL miniSolver%initialize(configFile)
      END SELECT
    END ASSOCIATE

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
    CLASS(Applications_t) :: this
    INTEGER(i_kind), INTENT(IN) :: mgStart, mgEnd
    LOGICAL, INTENT(IN) :: unitTest

    ! Local variables:
    INTEGER(i_kind) :: i, ii
    REAL(r_kind) :: t1, t2

    TYPE(IOWRF_t), TARGET :: ioWRF
    TYPE(IOECM_t), TARGET :: ioECM
    TYPE(IOGrapes_t), TARGET :: ioGrapes
    TYPE(IOGrapesPostvar_t), TARGET :: ioGrapesPostvar
    TYPE(IOGrapesModelvar_500m_t), TARGET :: IOGrapesModelvar_500m

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

              SELECT CASE (TRIM(bkModel))
              CASE ('WRF')
                ALLOCATE (IOWRF_t:: this%IOModel)
              CASE ('GRAPES')
                ALLOCATE (IOGrapes_t:: this%IOModel)
              CASE ('POSTVAR')
                ALLOCATE (IOGrapesPostvar_t:: this%IOModel)
              CASE ('ECM')
                ALLOCATE (IOECM_t:: this%IOModel)
              CASE ('MODELVAR')
                ALLOCATE (IOGrapesModelvar_500m_t:: this%IOModel)
              END SELECT

              ASSOCIATE (IOModel => this%IOModel)
                SELECT TYPE (IOModel)
                TYPE IS (IOWRF_t)
                  CALL IOModel%initialize(this%configFile, this%geometry)
                TYPE IS (IOGrapes_t)
                  CALL IOModel%initialize(this%configFile, this%geometry)
                TYPE IS (IOGrapesPostvar_t)
                  CALL IOModel%initialize(this%configFile, this%geometry)
                TYPE IS (IOECM_t)
                  CALL IOModel%initialize(this%configFile, this%geometry)
                TYPE IS (IOGrapesModelvar_500m_t)
                  CALL IOModel%initialize(this%configFile, this%geometry)
                END SELECT
              END ASSOCIATE

              PRINT *, "In reading ", TRIM(bkModel)
              CALL this%IOModel%m_read_bcg_into_Xm(this%XbMG(mgEnd), sg)
              PRINT *, "Done reading ", TRIM(bkModel)
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

      IF (TRIM(this%mode) == 'Debug') THEN
        PRINT *, "Output the background field. ----------*************"
        CALL Output_NC_State_AV(this%XbMG(i), this%ncOutputFile, &
                                TRIM(this%task)//"_bak", .TRUE., .TRUE.)
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

      CALL this%Ctl2State%transBackward(this%XbMG(i))
      IF (TRIM(this%framework) .EQ. 'Incremental') CALL this%XbMG(i)%setAllFieldData(0.0D0)

      BLOCK
        CHARACTER(LEN=2) :: iGrid
        WRITE (iGrid, "(I2.2)") this%XbMG(i)%sg%gLevel
        IF (SIZE(this%XbMG(i)%sg%SVapor, 2) .GT. 0) THEN
          OPEN (90, file=TRIM(this%ncOutputFile)//'/sv_g'//TRIM(iGrid)//'.txt', form='formatted')
          ! write(90,*) this%XbMG(i)%sg%SVapor1D
          WRITE (90, *) 'SVapor: '
          WRITE (90, *) this%XbMG(i)%sg%SVapor(:, 1)
          WRITE (90, *) 'SIce: '
          WRITE (90, *) this%XbMG(i)%sg%SIce(:, 1)
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

      this%XbMGInitial(i) = this%XbMG(i)

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
    CLASS(Applications_t) :: this
    LOGICAL, INTENT(IN) :: unitTest
    INTEGER(i_kind) :: i
    TYPE(State_t) :: XbRef

    IF (.NOT. unitTest) THEN
      ! Scaling back with control variable to physical space:
      XbRef = this%XbMGInitial(this%mgEnd)
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

            BLOCK
              LOGICAL :: testRepUV = .FALSE., testRepClouds = .FALSE.
              TYPE(IOERA5_t) :: ioERA5
              TYPE(state_t) :: XbERA5

              IF (yaml_get_var(TRIM(this%configFile), 'Verify', 'testRepUV', testRepUV) == 0) THEN
                IF (testRepUV) THEN
                  CALL ioERA5%initialize(this%configFile, this%geometry)
                  XbERA5 = XbRef
                  CALL ioERA5%m_read_bcg_into_Xm(XbERA5, this%geometry%mg%sg(this%mgEnd))
                  this%Analysis%Fields(this%Analysis%getVarIdx('uwnd'))%DATA = XbERA5%Fields(XbERA5%getVarIdx('uwnd'))%DATA
                  this%Analysis%Fields(this%Analysis%getVarIdx('vwnd'))%DATA = XbERA5%Fields(XbERA5%getVarIdx('vwnd'))%DATA
                  this%Analysis%Fields(this%Analysis%getVarIdx('qvapor'))%DATA = XbERA5%Fields(XbERA5%getVarIdx('qvapor'))%DATA
                  this%Analysis%Fields(this%Analysis%getVarIdx('temp'))%DATA = XbERA5%Fields(XbERA5%getVarIdx('temp'))%DATA
                  CALL Output_NC_State_AVST(this%Analysis, '/Users/yaliwu/Desktop/MOTOR/MOTOR/build/', 'testIOEra5', this%geometry%mg%sg(this%mgEnd)%tSlots, .TRUE., .TRUE.)
                END IF
              END IF
              IF (yaml_get_var(TRIM(this%configFile), 'Verify', 'testRepClouds', testRepClouds) == 0) THEN
                IF (testRepClouds) THEN
                  CALL ioERA5%initialize(this%configFile, this%geometry)
                  XbERA5 = XbRef
                  CALL ioERA5%m_read_bcg_into_Xm(XbERA5, this%geometry%mg%sg(this%mgEnd))
                  this%Analysis%Fields(this%Analysis%getVarIdx('temp'))%DATA = XbERA5%Fields(XbERA5%getVarIdx('temp'))%DATA
                  this%Analysis%Fields(this%Analysis%getVarIdx('qvapor'))%DATA = XbERA5%Fields(XbERA5%getVarIdx('qvapor'))%DATA
                  this%Analysis%Fields(this%Analysis%getVarIdx('qcloud'))%DATA = XbERA5%Fields(XbERA5%getVarIdx('qcloud'))%DATA
                  this%Analysis%Fields(this%Analysis%getVarIdx('qice'))%DATA = XbERA5%Fields(XbERA5%getVarIdx('qice'))%DATA
                  CALL Output_NC_State_AVST(this%Analysis, '/Users/yaliwu/Desktop/MOTOR/MOTOR/build/', 'testIOEra5', this%geometry%mg%sg(this%mgEnd)%tSlots, .TRUE., .TRUE.)
                END IF
              END IF

            END BLOCK

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
    CLASS(Applications_t) :: this
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

    CALL CPU_TIME(t2)
    WRITE (*, 1) t2 - t1
1   FORMAT('Time spent on background ensemble: ', D12.4)
  END SUBROUTINE backgrdEns

  !> @brief
  !! Calculate BEC through QR decomposition method:
  SUBROUTINE calEnsB(this)
    IMPLICIT NONE
    CLASS(Applications_t) :: this

    ! Local variables:
    INTEGER(i_kind) :: iv, i
    REAL(r_kind) :: ana_max, ana_min, ana_max_global, ana_min_global, t1, t2
    REAL(r_kind), ALLOCATABLE :: ana_max_AV(:), ana_min_AV(:)
    REAL(r_kind), PARAMETER :: passthreshold_max = 6.0D0, passthreshold_min = 0.6D0
    TYPE(EnLoc_t) :: EnLoc
    INTEGER(i_kind) :: count_start, count_end, count_rate
    REAL(r_kind) :: ts

    CALL CPU_TIME(t1)

    CALL SYSTEM_CLOCK(count_start, count_rate)
    CALL this%backgrdEns()
    CALL SYSTEM_CLOCK(count_end)
    ts = REAL(count_end - count_start) / REAL(count_rate)
    PRINT *, 'Done backgrdEns, cost ', ts

    IF (this%mpddGlob%isBaseProc()) THEN
      DO i = this%mgStart, this%mgEnd
        PRINT *, 'begin EnLoc at Gid:', i
        CALL SYSTEM_CLOCK(count_start, count_rate)
        CALL EnLoc%b_getEnsData(this%configFile, this%XbMGEns(i, :))
        CALL SYSTEM_CLOCK(count_end)
        ts = REAL(count_end - count_start) / REAL(count_rate)
        PRINT *, 'Done EnLoc%b_getEnsdata in grid ', i, ', cost ', ts

        CALL SYSTEM_CLOCK(count_start, count_rate)
        CALL EnLoc%b_qrDecomp(this%configFile, this%XbMGEns(i, :), i)
        CALL SYSTEM_CLOCK(count_end)
        ts = REAL(count_end - count_start) / REAL(count_rate)
        PRINT *, 'Done EnLoc%b_qrDecomp in grid', i, ', cost ', ts

        CALL SYSTEM_CLOCK(count_start, count_rate)
        CALL EnLoc%b_qrRinv(this%XbMGEns(i, 1))
        CALL SYSTEM_CLOCK(count_end)
        ts = REAL(count_end - count_start) / REAL(count_rate)
        PRINT *, 'Done EnLoc%b_qrRinv in grid', i, ', cost ', ts

        CALL SYSTEM_CLOCK(count_start, count_rate)
        CALL EnLoc%b_ADJ(this%configFile, this%XbMGEns(i, 1), i)
        CALL SYSTEM_CLOCK(count_end)
        ts = REAL(count_end - count_start) / REAL(count_rate)
        PRINT *, 'Done EnLoc%b_ADJ in grid', i, ', cost ', ts

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
    CLASS(Applications_t) :: this
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
        PRINT *, 'readObs - sfc min/max: ', &
          MINVAL(this%sfc%obsData), MAXVAL(this%sfc%obsData)
      CASE ('sounding')
        CALL this%snd%ObsInitial(this%configFile)
        CALL this%snd%ObsIngest(XbRef)
        PRINT *, 'readObs - snd min/max: ', &
          MINVAL(this%snd%obsData), MAXVAL(this%snd%obsData)
      CASE ('profiler')
        CALL this%pro%ObsInitial(this%configFile)
        CALL this%pro%ObsIngest(XbRef)
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
      CASE ('gnssrefr')
        CALL this%gnssro%ObsInitial(this%configFile)
        CALL this%gnssro%ObsIngest(XbRef)
      END SELECT
    END DO

    IF (ALL(.NOT. this%turnOn)) THEN
      PRINT *, " =========== NO satellite data will be assimilated ==========="
      RETURN
    ELSE
      DO i_inst = 1, SIZE(this%ObsSatellites)
        IF (this%turnOn(i_inst)) THEN
          PRINT *, 'readObs: ', TRIM(this%platform_name(i_inst))//'-'//TRIM(this%inst_name(i_inst)), ' is being processed'
          CALL this%ObsSatellites(i_inst)%ObsInitial(this%configFile, i_inst)
          CALL this%ObsSatellites(i_inst)%ObsIngest(XbRef)
        END IF
      END DO
    END IF

    CALL CPU_TIME(t2)
    WRITE (*, 1) t2 - t1, this%mpddGlob%myrank
1   FORMAT('readObs: time spent ingesting observation data: ', D12.4, ' at proc: ', I3)

  END SUBROUTINE readObs

  SUBROUTINE DASolve(this, X, obsList, sg, Nocoarest, k)
    CLASS(Applications_t) :: this
    TYPE(State_t) :: X, XbRef, XbCV
    TYPE(SingleGrid_t) :: sg
    CHARACTER(LEN=*), INTENT(IN) :: obsList(:)
    INTEGER(i_kind) :: Nocoarest(:)
    INTEGER, PARAMETER :: NYmax = 15
    INTEGER(i_kind) :: k

    TYPE(BMatrix_t) :: B
    TYPE(BMatrix_t) :: B_e, B_einc
    TYPE(RMatrix_t) :: R
    TYPE(MPObs_t), TARGET :: mpObs
    TYPE(ObsSet_t) :: Yi(NYmax)
    TYPE(ObsSet_t) :: Y, Z
    TYPE(C2O_t) :: H, HY
    TYPE(JFunc_t) :: JFunc
    INTEGER(i_kind) :: ii, j, i, NObs_tmp, NY
    REAL(r_kind) :: t1, t2

    CALL mpObs%initializeMPObs(sg)    ! Initialize the observation parallel processing proc

    IF (.NOT. sg%isActiveProc()) RETURN

    ! XbRef = this%Ctl2State%fwdNL_opr(this%XbMG(sg%gLevel))
    XbRef = this%XbRef(sg%gLevel)

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

        IF (TRIM(this%task) == '3DVar') THEN
          CALL Yi(1)%rmVar('SYNOP_qvapor')

          IF (sg%gLevel < this%mgEnd - 2) THEN
            CALL Yi(1)%rmVar('SYNOP_temp')
          END IF

          IF (sg%gLevel < this%mgEnd - 2) THEN
            CALL Yi(1)%rmVar('SYNOP_pres')
            CALL Yi(1)%rmVar('SYNOP_psl')
          END IF
        END IF

        ! A temporary processing to remove analysis on precipitation on coarser grid
      CASE ('sounding')
        CALL this%snd%ObsThinning(XbRef, Yi(2), mpObs, .TRUE., .FALSE.)
        CALL Yi(2)%rmVar('SOUND_qvapor')
        CALL Yi(2)%rmVar('SOUND_pres') ! remove the pressure in sounding
      CASE ('profiler')
        ! IF (sg%gLevel >= this%mgEnd - 1) THEN
        CALL this%pro%ObsThinning(XbRef, Yi(3), mpObs, .TRUE., .FALSE.)
        ! END IF
      CASE ('radarRef')
        IF (sg%gLevel >= this%mgEnd) THEN
          CALL this%ref%ObsThinning(XbRef, Yi(4), mpObs, .TRUE., .FALSE.)
        END IF
      CASE ('radarVel')
        IF (sg%gLevel >= this%mgEnd) THEN
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
      CASE ('gnssrefr')
        IF (sg%gLevel >= this%mgEnd - 3) THEN
          CALL this%gnssro%ObsThinning(XbRef, Yi(14), mpObs, .TRUE., .FALSE.)
        END IF
      END SELECT
    END DO

    IF (ALL(.NOT. this%turnOn)) THEN
      PRINT *, " =========== NO satellite data will be assimilated ==========="
    ELSE
      NY = 9
      DO i_inst = 1, SIZE(this%ObsSatellites)
        IF (this%turnOn(i_inst)) THEN
          PRINT *, 'DASolve: ', TRIM(this%platform_name(i_inst))//'-'//TRIM(this%inst_name(i_inst)), ' is being processed'
          ! IF (sg%gLevel .GE. this%mgEnd - 1) THEN
          NY = NY + 1
          CALL this%ObsSatellites(i_inst)%ObsPrepareForSg(XbRef)
          CALL this%ObsSatellites(i_inst)%ObsThinning(XbRef, Yi(NY), mpObs, .TRUE., .FALSE.)
          ! END IF
        END IF
      END DO
    END IF

    CALL CPU_TIME(t2)
    WRITE (*, 2) t2 - t1, sg%mpddInfo_sg%myrank, sg%num_cell
2   FORMAT('Time spent obsThinning: ', D14.7, ' at proc: ', I1, ' numCell: ', I5)

    ! Concat thinning obs:
    Y = ObsSet_t(this%configFile, mpObs)

    CALL ObsConcat_s((/Yi(1), Yi(2), Yi(3), Yi(4), Yi(5), Yi(6), Yi(7), Yi(8), Yi(9), Yi(10), Yi(11), Yi(12), Yi(13), Yi(14)/), Y)
    WRITE (*, 12) UBOUND(Y%ObsFields, 1), &
      (TRIM(Y%ObsFields(j)%get_name()), j=1, UBOUND(Y%obsFields, 1)), &
      (TRIM(Y%ObsFields(j)%Get_Id_Name()), j=1, UBOUND(Y%obsFields, 1))
12  FORMAT('Concatenated obs vars: ', I2, 20(1X, A))

    ! Model state variables to fill in data sparse areas

    BLOCK
      TYPE(State_t) :: X_Ctl
      X_Ctl = this%XbMG(sg%gLevel)

      Z = ObsSet_t(this%configFile, mpObs)
      CALL this%Ctl2State%transFwdNonLinear(X_Ctl)
      IF (TRIM(this%framework) .EQ. 'Incremental') THEN
        CALL this%Ctl2State%SumState(X_Ctl, this%XbRef(sg%gLevel))
      END IF

      CALL bkgdFillToObs(X_Ctl, Y, Z, this%iRange_ramp)
    END BLOCK

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
    ELSE
      CALL B%initialize(this%configFile, sg, 'Laplace')
    END IF

    CALL H%initialize(this%configFile, X, Z, this%XbRef(sg%gLevel))
    WRITE (*, 13) X%sg%glevel, X%sg%mpddInfo_sg%myrank
13  FORMAT('H is initialized....at G', I2, ' pc ', I2)
    CALL HY%initialize(this%configFile, X, Y, this%XbRef(sg%gLevel))
    WRITE (*, 14) X%sg%glevel, X%sg%mpddInfo_sg%myrank
14  FORMAT('HY is initialized....at G', I2, ' pc ', I2)

    IF (TRIM(this%mode) == 'Debug') THEN
      BLOCK
        TYPE(State_t) :: YY, DD, ZZ, BB
        CHARACTER(len=2) :: kstr

        YY = Obs2State_BaseTypeName(sg, Y)
        DD = Obs2State_BaseTypeName(sg, Y - HY%transFwdNonLinear_opr(this%XbMG(sg%gLevel)))

        IF (k < 10) WRITE (kstr, '(I1)') k
        IF (k < 100 .AND. k >= 10) WRITE (kstr, '(I2)') k

!         WRITE (*, 15) X%sg%glevel, X%sg%mpddInfo_sg%myrank
! 15      FORMAT('DD is calculated...at G', I2, ' pc ', I2)
        BB = Obs2State_BaseTypeName(sg, HY%transFwdNonLinear_opr(this%XbMG(sg%gLevel)))

        CALL Output_NC_State_AV(YY, this%ncOutputFile, &
                                TRIM(this%task)//"_obsThinned_"//TRIM(kstr), .TRUE., .TRUE.)
!         WRITE (*, 18) X%sg%glevel, X%sg%mpddInfo_sg%myrank
! 18      FORMAT('Thinned is plotted...at G', I2, ' pc ', I2)
        CALL Output_NC_State_AV(DD, this%ncOutputFile, &
                                TRIM(this%task)//"_innovation_"//TRIM(kstr), .TRUE., .TRUE.)
!         WRITE (*, 28) X%sg%glevel, X%sg%mpddInfo_sg%myrank
! 28      FORMAT('Innovation at G', I2, ' pc ', I2)
        CALL Output_NC_State_AV(BB, this%ncOutputFile, &
                                TRIM(this%task)//"_bcgAtObs_"//TRIM(kstr), .TRUE., .TRUE.)

!         CALL Output_NC_State_AV(this%XbMG(sg%gLevel), this%ncOutputFile, &
!                                 TRIM(this%task)//"_bakInDA_"//TRIM(kstr), .TRUE., .TRUE.)

!         ! Plot fillin data:
!         ZZ = Obs2State_BaseTypeName(sg, Z)
!         CALL Output_NC_State_AV(ZZ, this%ncOutputFile, &
!                                 TRIM(this%task)//"_obsFilled_"//TRIM(kstr), .TRUE., .TRUE.)
!         WRITE (*, 38) X%sg%glevel, sg%gLevel, X%sg%mpddInfo_sg%myrank
! 38      FORMAT('Filled is plotted...at G', 2I2, ' pc ', I2)
      END BLOCK
    END IF

    WRITE (*, 16) X%sg%glevel, X%mpddGlob%myrank
16  FORMAT('Initializing R......at G', I2, ' pc ', I2)
    ! CALL R%initialize(this%configFile, Z, X%sg, Nocoarest) ! Initialize R
    CALL R%initialize(this%configFile, Z, X%sg) ! Initialize R

    CALL sg%mpddInfo_sg%barrier

    ! Initialize J Function
    WRITE (*, 17) X%sg%glevel, X%mpddGlob%myrank
17  FORMAT('Initializing JFunc......at G', I2, ' pc ', I2)

    SELECT CASE (TRIM(this%framework))
    CASE ('FullState')
      XbCV = this%XbMG(X%sg%gLevel)
    CASE ('Incremental')
      XbCV = this%XbRef(X%sg%gLevel)
      CALL this%Ctl2State%transBackward(XbCV)
      CALL this%Ctl2State%SumBothCtlState(XbCV, this%XbMG(X%sg%gLevel))
    END SELECT

    IF (this%Hybrid) THEN
      CALL JFunc%initialize(this%configFile, X, Z, H, B, R, sg, Xb=XbCV, B_e=B_e)
    ELSE
      CALL JFunc%initialize(this%configFile, X, Z, H, B, R, sg, Xb=XbCV)
    END IF

    CALL sg%mpddInfo_sg%barrier

    ! Run minimization
    CALL CPU_TIME(t1)
    IF (sg%gLevel <= 5) THEN
      ! CALL this%miniSolver%run(X, JFunc, sg, 100)
      CALL this%miniSolver%run(X, JFunc, sg)
    ELSE
      CALL this%miniSolver%run(X, JFunc, sg)
    END IF
    CALL CPU_TIME(t2)

    WRITE (*, 3) t2 - t1, sg%mpddInfo_sg%myrank, sg%num_cell
3   FORMAT('Time spent in miniSolver: ', D10.2, ' at proc: ', I1, ' with Ncell: ', I8)

    IF (TRIM(this%mode) == 'Debug') THEN
      BLOCK
        TYPE(State_t) :: AA, ANA, tempX
        CHARACTER(len=2) :: kstr

        AA = Obs2State_BaseTypeName(sg, Y - H%transFwdNonLinear_opr(X))
        ANA = Obs2State_BaseTypeName(sg, H%transFwdNonLinear_opr(X))

        IF (k < 10) WRITE (kstr, '(I1)') k
        IF (k < 100 .AND. k >= 10) WRITE (kstr, '(I2)') k

        CALL Output_NC_State_AV(AA, this%ncOutputFile, &
                                TRIM(this%task)//"_AMO_"//TRIM(kstr), .TRUE., .TRUE.)
        CALL Output_NC_State_AV(ANA, this%ncOutputFile, &
                                TRIM(this%task)//"_AnaAtObs_"//TRIM(kstr), .TRUE., .TRUE.)

        tempX = X
        CALL this%Ctl2State%transFwdNonLinear(tempX) ! output in CV or BKG space
        IF (TRIM(this%framework) .EQ. 'Incremental') THEN
          CALL this%Ctl2State%SumState(tempX, this%XbRef(i))
        END IF
        CALL tempX%fillPresWithHydrostatic()

        CALL Output_NC_State_AV(tempX, this%ncOutputFile, &
                                TRIM(this%task)//"_ana_"//TRIM(kstr), .TRUE., .TRUE.)
      END BLOCK
    END IF

  END SUBROUTINE DASolve

  !> @brief
  !================================================================
  !  This routine performs a MOTOR-DA analysis:
  !  Input:
  !       mgStart:    integer G level to start a multigrid analysis
  !       mgEnd:      integer G level to end the multigrid analysis
  !
  SUBROUTINE analyss(this, obsList, mgStart, mgEnd)
    IMPLICIT NONE
    CLASS(Applications_t) :: this
    CHARACTER(LEN=*), INTENT(IN) :: obsList(:)
    INTEGER(i_kind), INTENT(IN) :: mgStart, mgEnd

    ! Local variables:
    CHARACTER(LEN=20) :: varName
    INTEGER(i_kind) :: i, iv, k, loops, idx = 0
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

    IF (this%FixResolution) THEN
      loops = this%OuterLoops
    ELSE
      loops = 1
    END IF
    DO k = 1, loops
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
              WHERE (X%fields(X%getVarIdx('qcloud'))%DATA <= 1.0D-12) X%fields(X%getVarIdx('qcloud'))%DATA = 0.0D0
              WHERE (X%fields(X%getVarIdx('qice'))%DATA <= 1.0D-12) X%fields(X%getVarIdx('qice'))%DATA = 0.0D0
            END IF
            IF (X%getVarIdx('qcloud_ctl') .GT. 0) THEN
              WHERE (X%fields(X%getVarIdx('qcloud_ctl'))%DATA <= 1.0D-12) X%fields(X%getVarIdx('qcloud_ctl'))%DATA = 0.0D0
              WHERE (X%fields(X%getVarIdx('qice_ctl'))%DATA <= 1.0D-12) X%fields(X%getVarIdx('qice_ctl'))%DATA = 0.0D0
            END IF

            ! Run the minimization
            CALL this%DASolve(X, obsList, sg, Nocoarest, k)

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
                ! print*, 'this%XbMG(i): ', maxval(this%XbMG(i)%fields(7)%data), minval(this%XbMG(i)%fields(7)%data)
                CALL Output_NC_State_AV(X - this%XbMG(i), this%ncOutputFile, &
                                        TRIM(this%task)//"_diff", .TRUE., .TRUE.)
                CALL Output_NC_State_AV(this%XbMG(i), this%ncOutputFile, &
                                        TRIM(this%task)//"_diff_before", .TRUE., .TRUE.)
                tempX = X
                CALL this%Ctl2State%transFwdNonLinear(tempX) ! output in CV or BKG space
                IF (TRIM(this%framework) .EQ. 'Incremental') THEN
                  CALL this%Ctl2State%SumState(tempX, this%XbRef(i))
                END IF
                CALL tempX%fillPresWithHydrostatic()

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
                ! CALL prolongationMGInc(X, this%analysis, this%XbMG(i), this%XbMG(i + 1), sg, sgFiner, this%geometry%mg)
                CALL prolongationMGInc(X, this%analysis, this%XbMGInitial(i), this%XbMG(i + 1), sg, sgFiner, this%geometry%mg)

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
              idx = i + 1
            ELSE
              this%analysis = X ! Return the state fields directly
              idx = i
            END IF

            ! UpdateBkgd means the finer grid background is updated with its coarser grid analysis increment:
            IF (this%UpdateBkgd) THEN
              this%XbMG(idx) = this%analysis
              ! print*, 'this%XbMG(idx): ', maxval(this%XbMG(idx)%fields(7)%data), minval(this%XbMG(idx)%fields(7)%data)

              BLOCK
                TYPE(State_t) :: XCV, XbCV

                SELECT CASE (TRIM(this%framework))
                CASE ('FullState')
                  this%XbRef(idx) = this%Ctl2State%transFwdNonLinear_opr(this%analysis)
                CASE ('Incremental')
                  XCV = this%analysis
                  XbCV = this%XbRef(idx)
                  CALL this%Ctl2State%transBackward(XbCV)
                  CALL this%Ctl2State%SumBothCtlState(XCV, XbCV)
                  this%XbRef(idx) = this%Ctl2State%transFwdNonLinear_opr(XCV)
                END SELECT
              END BLOCK
              IF (this%XbRef(idx)%getVarIdx('qvapor') .GT. 0) THEN
                WHERE (this%XbRef(idx)%fields(this%XbRef(idx)%getVarIdx('qvapor'))%DATA <= 1.0D-12) &
                  this%XbRef(idx)%fields(this%XbRef(idx)%getVarIdx('qvapor'))%DATA = 1.0D-12
              END IF
              IF (this%XbRef(idx)%getVarIdx('qcloud') .GT. 0) THEN
                WHERE (this%XbRef(idx)%fields(this%XbRef(idx)%getVarIdx('qcloud'))%DATA <= 1.0D-12) &
                  this%XbRef(idx)%fields(this%XbRef(idx)%getVarIdx('qcloud'))%DATA = 0.0D0
              END IF
              IF (this%XbRef(idx)%getVarIdx('qice') .GT. 0) THEN
                WHERE (this%XbRef(idx)%fields(this%XbRef(idx)%getVarIdx('qice'))%DATA <= 1.0D-12) &
                  this%XbRef(idx)%fields(this%XbRef(idx)%getVarIdx('qice'))%DATA = 0.0D0
              END IF

              IF (TRIM(this%framework) .EQ. 'Incremental') THEN
                CALL this%XbMG(idx)%setAllFieldData(0.0D0)
              END IF

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
  !! Destory the instance of geometry class.
  SUBROUTINE destroy(this)
    IMPLICIT NONE
    CLASS(Applications_t) :: this

    CALL this%Geometry%destroy

    CALL this%mpddGlob%barrier

    DEALLOCATE (this%miniSolver)
    DEALLOCATE (this%XbMG)
    DEALLOCATE (this%XbRef)
    DEALLOCATE (this%XbMGInitial)
    DEALLOCATE (this%platform_name, this%inst_name, this%turnOn)
    DEALLOCATE (this%ObsSatellites)
    ! Deallocate memory:
    IF (ALLOCATED(this%mgGridSmooth)) DEALLOCATE (this%mgGridSmooth)

    ! Finalize
    CALL this%mpddGlob%finalize

  END SUBROUTINE destroy

END MODULE Applications_m
