!!--------------------------------------------------------------------------------------------------
! PROJECT           : Application of MOTOR-DA 4DVAR
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : V 1.0
! HISTORY           : 2024-01-17, created by Yuanfu Xie.
!
!   Created by Yuanfu Xie (yuanfu_xie@yahoo.com), 2024/01/17, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------
!!
!> @brief
!! This is a MOTOR-DA 4DVAR with analysis
MODULE App_4DVAR_m
  USE SolverLBFGS_m, ONLY: SolverLBFGS_t
  USE SolverFRCG_m, ONLY: SolverFRCG_t
  USE MiniSolver_m, ONLY: MiniSolver_t
  USE FRCGSolver_m, ONLY: FRCGSolver_t
  USE variationalSolvers_m, ONLY: variationalSolvers_t
  USE poissonSolver_m, ONLY: poissonSolver_t
  USE geometry_m, ONLY: geometry_t
  USE GenContainers_m, ONLY: GenContainers_t    !< Jilong's wrapper initializing a geometry
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE State_m, ONLY: State_t

  ! Control to model mapping:
  USE C2MBase_m, ONLY: C2MBase_t
  USE c4DVar_m, ONLY: c4DVar_t
  USE cDefault_m, ONLY: cDefault_t

  USE ObsSet_m, ONLY: ObsSet_t
  USE ObsField_m, ONLY: ObsField_t
  USE BMatrix_m, ONLY: BMatrix_t
  USE RMatrix_m, ONLY: RMatrix_t
  USE JFunc_m, ONLY: JFunc_t
  USE State2NC_m
  USE parameters_m, ONLY: degree2radian
  USE YAMLRead_m

  USE obsMG_m, ONLY: obsMG_t
  USE bkgMG_m, ONLY: bkgMG_t
  USE ObsUtilities_m
  USE dyCoresBase_m, ONLY: dyCoresBase_t
  USE dyCoreGZM_m, ONLY: dyCoreGZM_t
  USE rhsBase_m, ONLY: rhsBase_t
  USE rhsShallowWater_m, ONLY: rhsShallowWater_t
  USE TimeIntegrationBase_m, ONLY: TimeIntegrationBase_t
  USE TimeIntegrationAB3_m, ONLY: TimeIntegrationAB3_t

  !USE Ctl2State_m, ONLY: Ctl2State_t ! Control to state converter, turned off for new design

  ! Yaml list:
  TYPE yamlVars_t
    CHARACTER(LEN=20), ALLOCATABLE :: varList(:), c2m(:)
    CHARACTER(LEN=8) :: minimization, bkgdModel
    CHARACTER(LEN=1024) :: yamlFile, ncOutputFile
    CHARACTER(LEN=20) :: mode = 'Debug'
    CHARACTER(LEN=20) :: task
    LOGICAL :: updateBkgd = .FALSE.   ! Default as no background update
    REAL(r_kind), ALLOCATABLE :: weights(:)
  END TYPE yamlVars_t

  ! Yuanfu Xie added a new type for the 4DVar application:
  TYPE choice4dvar_t
    CLASS(C2MBase_t), POINTER :: ctl         ! Control variables, C2MBase_t is an abstract type
    CLASS(dyCoresBase_t), POINTER :: dyc   ! Dycore for the 4DVar, can be a dyCoreGZM_t or other dycores
    CLASS(rhsBase_t), POINTER :: rhs
    CLASS(TimeIntegrationBase_t), POINTER :: tim  ! Time integration for the dycore
  END TYPE choice4dvar_t

  TYPE App_4DVAR_t

    INTEGER(i_kind) ::numVar
    INTEGER(i_kind), ALLOCATABLE :: mgGridSmooth(:, :)

    ! Yaml list:
    TYPE(yamlVars_t) :: yamlVars
    ! Analysis: Yuanfu Xie turned off the analysis variable in the new design 2025-06-26
    ! TYPE(State_t) :: analysis ! In the new design, the analysis is done in the control variable
    ! Multi-processors handler:
    TYPE(mpddGlob_t) :: mpddGlob
    ! Geometry:
    TYPE(geometry_t) :: geometry
    ! Controls for a multigrid:
    ! Minimization solver:
    ! CLASS(MiniSolver_t), ALLOCATABLE :: miniSolver
    ! CLASS(variationalSolvers_t), ALLOCATABLE :: varSolver
    TYPE(FRCGSolver_t) :: varSolver  ! Variational solver, FRCG or LBFGS

    ! 4DVar choices
    CLASS(choice4dvar_t), POINTER :: Choice4dv(:)  ! Control variables, dycore, rhs and time integration

    TYPE(PoissonSolver_t) :: poissonSolver  ! Poisson solver

    !TYPE(Ctl2State_t) :: Ctl2State   ! Control to state converter, turned off for new design
    ! Background field fill-in rampRanges:
    INTEGER(i_kind), ALLOCATABLE :: iRange_ramp(:, :)  ! space and time dimensions : variables

  CONTAINS
    PROCEDURE, PRIVATE :: readInYaml_s
    PROCEDURE :: initialize_s
    PROCEDURE :: destroy_s

    PROCEDURE :: analysis_s
    PROCEDURE :: singleLvlAnal_s
    PROCEDURE :: bkgdFillin_s
  END TYPE App_4DVAR_t

CONTAINS
  !> @brief
  !！This routine reads all yaml information, no other places to do so.
  SUBROUTINE readInYaml_s(this)
    IMPLICIT NONE
    CLASS(App_4DVAR_t) :: this

    ! Local variables:
    INTEGER(i_kind) :: i, istatus

    istatus = yaml_get_var(TRIM(this%yamlVars%yamlFile), 'IO', 'bk_model', this%yamlVars%bkgdModel)
    istatus = yaml_get_var(TRIM(this%yamlVars%yamlFile), 'Minimization', 'method', this%yamlVars%minimization)
    istatus = yaml_get_var(TRIM(this%yamlVars%yamlFile), 'RunMode', 'Mode', this%yamlVars%mode)
    istatus = yaml_get_var(TRIM(this%yamlVars%yamlFile), 'RunMode', 'Task', this%yamlVars%Task)
    ! Multigrid scheme:
    istatus = yaml_get_var(TRIM(this%yamlVars%yamlFile), 'MultigridOptions', 'UpdateBkgd', this%yamlVars%UpdateBkgd)
    istatus = yaml_get_var(TRIM(this%yamlVars%yamlFile), 'IO', 'output_dir', this%yamlVars%ncOutputFile)
    istatus = yaml_get_var(TRIM(this%yamlVars%yamlFile), 'BMat', 'bkgdWeights', this%yamlVars%weights)

    ! Yuanfu Xie added a yaml variable for c2m mapping: 2025-06-26
    istatus = yaml_get_var(TRIM(this%yamlVars%yamlFile), 'modelState', 'control2Model', this%yamlVars%c2m)
    IF (istatus .NE. 0) THEN
      WRITE (*, 10)
10    FORMAT('No c2m mapping specified in your yaml file, default will be used!')
    ELSE
      WRITE(*, 12) UBOUND(this%yamlVars%c2m,1), (TRIM(this%yamlVars%c2m(i)),i=1,UBOUND(this%yamlVars%c2m,1))
12    FORMAT('c2m mapping: ', I2, ' c2m: ', 100(A,' '))
    END IF

    istatus = yaml_get_var(TRIM(this%yamlVars%yamlFile), 'modelState', 'varList', this%yamlVars%varList)
    IF (istatus .NE. 0) THEN
      this%numVar = 0
    ELSE
      this%numVar = UBOUND(this%yamlVars%varList, 1)
      WRITE (*, 11) this%numVar
11    FORMAT("Number of analysis vars: ", I2)
    END IF
  END SUBROUTINE readInYaml_s

  !> @brief
  !！This routine initializes analysis modules, geometry, mpdd etc.

!#define jilongGeometry
  SUBROUTINE initialize_s(this, istatus)
    IMPLICIT NONE
    CLASS(App_4DVAR_t) :: this
    INTEGER(i_kind), INTENT(OUT) :: istatus

    ! Local variables:
    INTEGER(i_kind) :: i
    TYPE(GenContainers_t) :: GenContainers  ! Temporarily use of Jilong's GenContainers to initialize a geometry

#ifndef jilongGeometry
    CALL this%mpddGlob%initialize()
#endif

    ! Get the configuration file
    CALL getarg(1, this%yamlVars%yamlFile)

    ! Initialize geometry
#ifdef jilongGeometry
    GenContainers = GenContainers_t(TRIM(this%yamlVars%yamlFile))
    CALL GenContainers%GenGeometry(this%geometry)
    this%mpddGlob = GenContainers%mpddGlob
#endif

#ifndef jilongGeometry
    CALL this%geometry%initialize(this%yamlVars%yamlFile, this%mpddGlob)
#endif

    IF (this%mpddGlob%isBaseProc() .AND. TRIM(this%yamlVars%yamlFile) .EQ. '') THEN
    ! IF (GenContainers%mpddGlob%isBaseProc() .AND. TRIM(this%yamlVars%yamlFile) .EQ. '') THEN
      PRINT *, TRIM(this%yamlVars%yamlFile)

      CALL getarg(0, this%yamlVars%yamlFile)
      PRINT *, TRIM(this%yamlVars%yamlFile)
      WRITE (*, 1)
1     FORMAT('Usage of this driver: mpirun -n <n> Debug/App_CMA_GD.exe yamlFile', /, &
             ' Check your configure file and rerun!')
    ELSE
      WRITE (*, 2) TRIM(this%yamlVars%yamlFile)
2     FORMAT('The yamlFile is: ', A)
    END IF

    CALL this%mpddGlob%barrier
    ! CALL GenContainers%mpddGlob%barrier

    ! Get yaml parameters:
    CALL this%readInYaml_s

    ! Initialize the geometry
    ! ALLOCATE (this%XbMG(this%geometry%mg%mg_coarsest:this%geometry%mg%mg_finest))

    ! Initialize solver
    SELECT CASE (this%yamlVars%minimization)
    ! CASE ('LBFGS')
    !   ALLOCATE (SolverLBFGS_t::this%miniSolver)
    CASE ('FRCG')
      ! ALLOCATE (SolverFRCG_t::this%miniSolver)
    CASE DEFAULT
      ! Make sure to print warning message at last:
      CALL this%mpddGlob%barrier
      IF (this%mpddGlob%isBaseProc()) WRITE (*, 14) TRIM(this%yamlVars%minimization)
14    FORMAT('This minimization method is not implemented: ', A8, / &
             'Check your option in app%initial and rerun!')
      istatus = 1
    END SELECT

    ! CALL this%miniSolver%initialize(this%yamlVars%yamlFile)
    CALL this%varSolver%initialize_s(this%yamlVars%yamlFile)

    ! Get background fillin information:
    BLOCK
      INTEGER(i_kind) :: i
      INTEGER(i_kind), ALLOCATABLE :: irange(:)

      istatus = yaml_get_var(TRIM(this%yamlVars%yamlFile), 'modelState', 'Fillin_range', irange)
      IF (ALLOCATED(this%iRange_ramp)) DEALLOCATE (this%iRange_ramp)
      ALLOCATE (this%iRange_ramp(3, this%numVar))
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
      istatus = yaml_get_var(TRIM(this%yamlVars%yamlFile), 'modelState', 'mgGridSmooth', irange)
      IF (istatus .EQ. 0 .AND. UBOUND(irange, 1) .EQ. (this%geometry%mg%mg_finest - this%geometry%mg%mg_coarsest + 1) * this%numVar) THEN
        ALLOCATE (this%mgGridSmooth(this%numVar, this%geometry%mg%mg_coarsest:this%geometry%mg%mg_finest))
        this%mgGridSmooth(:, this%geometry%mg%mg_coarsest:this%geometry%mg%mg_finest) = RESHAPE(irange, (/this%numVar, this%geometry%mg%mg_finest - this%geometry%mg%mg_coarsest + 1/))
        DO i = this%geometry%mg%mg_coarsest, this%geometry%mg%mg_finest
          WRITE (*, 20) i, this%mgGridSmooth(:, i)
20        FORMAT('mgGridSmooth at G', I2, ' number of smooth: ', 20I2)
        END DO
      ELSE
        IF (istatus .EQ. 0) THEN
          WRITE (*, 21) UBOUND(irange, 1), (this%geometry%mg%mg_finest - this%geometry%mg%mg_coarsest + 1) * this%numVar
21        FORMAT('No multigrid background smoothing as yaml parameters do not match with multigrid variables: ', 2I4)
        ELSE
          WRITE (*, 22)
22        FORMAT('No multigrid background smoothing parameters specified in your yaml...')
        END IF
      END IF

      IF (ALLOCATED(irange)) DEALLOCATE (irange)

      ! Initialize the Ctl2State_t: Turned off for new design:
      !CALL this%Ctl2State%initialize(this%yamlVars%yamlFile)
    END BLOCK

  END SUBROUTINE initialize_s

  !> @brief
  !! Destory the instance of geometry class.
  SUBROUTINE destroy_s(this)
    IMPLICIT NONE
    CLASS(App_4DVAR_t) :: this

    ! Local variables:
    INTEGER(i_kind) :: i

    IF (ASSOCIATED(this%Choice4dv)) THEN
      DO i= LBOUND(this%Choice4dv, 1), UBOUND(this%Choice4dv, 1)
        SELECT CASE (TRIM(this%yamlVars%c2m(i)))
        CASE ('zsw4dvar')
          IF (ASSOCIATED(this%Choice4dv(i)%rhs)) DEALLOCATE (this%Choice4dv(i)%rhs)
          IF (ASSOCIATED(this%Choice4dv(i)%tim)) DEALLOCATE (this%Choice4dv(i)%tim)
          IF (ASSOCIATED(this%Choice4dv(i)%dyc)) DEALLOCATE (this%Choice4dv(i)%dyc)
        END SELECT
        IF (ASSOCIATED(this%Choice4dv(i)%ctl)) DEALLOCATE (this%Choice4dv(i)%ctl)
      END DO
    END IF

    CALL this%Geometry%destroy

    CALL this%mpddGlob%barrier

    ! Finalize
    CALL this%mpddGlob%finalize
  END SUBROUTINE destroy_s

  ! The analysis using a multigrid from mgStart to mgEnd:
  SUBROUTINE analysis_s(this, obs, bkg, mgStart, mgEnd, istatus, dt)

    USE MGOpts_m

    IMPLICIT NONE

    CLASS(App_4DVAR_t) :: this
    ! Note: An abstract type, so CLASS is used to declare
    ! CLASS(C2MBase_t), INTENT(INOUT) :: ctl          ! Control variable
    TYPE(ObsMG_t), INTENT(IN) :: obs                  ! Observations
    CLASS(bkgMG_t), INTENT(INOUT) :: bkg              ! Background
    INTEGER(i_kind), INTENT(IN) :: mgStart, mgEnd     ! Grid level to analyze
    INTEGER(i_kind), INTENT(OUT) :: istatus
    REAL(r_kind), INTENT(IN), OPTIONAL :: dt  ! Time step for the analysis optional

    ! Local variables:
    INTEGER(i_kind) :: img

    ! Check the multigrid parameters with the grid geometry:
    IF (mgStart .LT. this%geometry%mg%mg_coarsest .OR. &
        mgEnd .GT. this%geometry%mg%mg_finest) THEN
      WRITE (*, 1) mgStart, mgEnd, this%geometry%mg%mg_coarsest, &
        this%geometry%mg%mg_finest
1     FORMAT('Multigrid analysis levels do not match the geometry setting, Check and rerun!', &
             ' Requested: ', 2I2, ' Geometry setting: ', 2I2)
      STOP
    END IF

    ! Initialize the analysis control variable for all grid level:
    IF (UBOUND(this%yamlVars%c2m, 1) .LT. mgEnd) THEN
      WRITE (*, 2) UBOUND(this%yamlVars%c2m, 1), mgEnd
2     FORMAT('c2m mapping does not match the multigrid levels, Check and rerun!', &
             ' Number of c2m mappings: ', I2, ' mgEnd: ', I2,' all c2m set to default!!!')
      ! Set all c2m to default:
      IF (ALLOCATED(this%yamlVars%c2m)) DEALLOCATE (this%yamlVars%c2m)
      ALLOCATE (this%yamlVars%c2m(mgEnd))
      this%yamlVars%c2m = 'default'
    END IF

    ! Choices of 4DVar:
    IF (ASSOCIATED(this%Choice4dv)) DEALLOCATE (this%Choice4dv)
    ALLOCATE (this%Choice4dv(mgStart:mgEnd))
    DO img = mgStart, mgEnd
      ! Initialize the analysis control variable at the current grid level based on the user
      ! options set in yaml file:
      this%choice4dv(img)%ctl => NULL()
      this%choice4dv(img)%rhs => NULL() 
      this%choice4dv(img)%tim => NULL()
      this%choice4dv(img)%dyc => NULL()

      WRITE(*, 5) img, TRIM(this%yamlVars%c2m(img))
5     FORMAT('Initializing analysis_s at level: ', I2, ' with c2m mapping: ', A)
      SELECT CASE (TRIM(this%yamlVars%c2m(img)))
      CASE ('zsw4dvar')
        ASSOCIATE(sg => this%geometry%mg%sg(img))
          ! Allocate the control variable:
          ALLOCATE(c4DVar_t :: this%choice4dv(img)%ctl)
          ALLOCATE(rhsShallowWater_t :: this%choice4dv(img)%rhs)
          ALLOCATE(dyCoreGZM_t :: this%choice4dv(img)%dyc)
          ALLOCATE(TimeIntegrationAB3_t :: this%choice4dv(img)%tim)

          CALL this%choice4dv(img)%rhs%initialize_s(sg, this%poissonSolver)
          CALL this%choice4dv(img)%tim%initialization_s(this%choice4dv(img)%rhs, bkg%XbMG(img), this%yamlVars%yamlFile)
          CALL this%choice4dv(img)%dyc%initialize_s(this%yamlVars%yamlFile, this%choice4dv(img)%tim, this%choice4dv(img)%rhs, sg, dt)
          ! Initialize the control variable:
          CALL this%choice4dv(img)%ctl%initialize_s(sg, bkg%XbMG(img), this%choice4dv(img)%dyc)
        END ASSOCIATE
      CASE ('vortDivg')
        ALLOCATE(cDefault_t :: this%choice4dv(img)%ctl)
        PRINT*, 'Using vortDivg for analysis_s at level: ', img
        CALL this%choice4dv(img)%ctl%initialize_s(this%geometry%mg%sg(img), bkg%XbMG(img))
        PRINT*, 'vortDivg_t initialized with default temporarily at level: ', img
      CASE DEFAULT
        ! If the c2m mapping is not recognized, use the default cDefault_t:
        WRITE (*, 3) TRIM(this%yamlVars%c2m(img))
3       FORMAT('c2m mapping not recognized: ', A, ' using default cDefault_t')
        ALLOCATE(cDefault_t :: this%choice4dv(img)%ctl)
        CALL this%choice4dv(img)%ctl%initialize_s(this%geometry%mg%sg(img), bkg%XbMG(img))
      END SELECT
    END DO

    ! Start the analysis:
    WRITE(*, 4) mgStart, mgEnd
4   FORMAT('Starting analysis_s from level: ', I2, ' to level: ', I2)
    DO img = mgStart, mgEnd
      ! Solve the analysis at level img:
      PRINT *, 'Starting singleLvlAnal_s... at level ', img
      CALL this%singleLvlAnal_s(this%choice4dv(img)%ctl, obs, bkg%XbMG(img), this%geometry%mg%sg(img), istatus)
      PRINT*,'Finised single grid analysis at level: ',img

      ! Output NC files for debugging: it calls the model forecast for a 4DVar:
      ! IF (TRIM(this%yamlVars%mode) == 'Debug') THEN
      !   PRINT *, 'Outputting diff in NC files at level: ', img
      !   CALL Output_NC_State_AV(ctl%forwardOpr_s(ctl%analysis) - bkg%XbMG(img), &
      !                           TRIM(this%yamlVars%ncOutputFile), &
      !                           TRIM(this%yamlVars%task)//"_diff", .TRUE., .TRUE.)
      !   PRINT *, 'Outputting ana in NC files at level: ', img
      !   CALL Output_NC_State_AV(ctl%forwardOpr_s(ctl%analysis), &
      !                           TRIM(this%yamlVars%ncOutputFile), &
      !                           TRIM(this%yamlVars%task)//"_ana", .TRUE., .TRUE.)
      ! END IF

      ! Multigrid prolongation if needed:
      IF (img .LT. mgEnd) THEN
        PRINT*,'Prolongation at level: ',img,mgEnd
        CALL prolongationMGInc(this%choice4dv(img)%ctl%analysis, this%choice4dv(img + 1)%ctl%analysis, &
                               bkg%XbMG(img), bkg%XbMG(img + 1), &
                               this%geometry%mg%sg(img), this%geometry%mg%sg(img + 1), this%geometry%mg)

        ! Check if model backgroud needs to be updated: If updateBkgd is set to true, recalculate the background
        ! from the controls interpolated from coarser grid analysis 
        IF (this%yamlVars%updateBkgd) THEN
          PRINT*,'Analysis_s - forwardOpr_s at level: ',img+1,UBOUND(this%choice4dv(img + 1)%ctl%analysis%fields,1)
          bkg%XbMG(img + 1) = this%choice4dv(img+1)%ctl%forwardOpr_s(this%choice4dv(img+1)%ctl%analysis)
        END IF
      END IF
    END DO
  END SUBROUTINE analysis_s

  ! A single level solver:
  SUBROUTINE singleLvlAnal_s(this, ctl, obs, bkg, sg, istatus)

    USE MPObs_m, ONLY: MPObs_t
    USE C2O_m, ONLY: C2O_t
    USE costFuncGrad_m, ONLY: costFuncGrad_t

    IMPLICIT NONE

    CLASS(App_4DVAR_t) :: this
    CLASS(C2MBase_t), INTENT(INOUT) :: ctl   ! Control variable
    TYPE(ObsMG_t), INTENT(IN) :: obs      ! Observations
    TYPE(state_t), INTENT(INOUT) :: bkg   ! Background
    TYPE(SingleGrid_t), INTENT(IN) :: sg  ! Single grid information
    INTEGER(i_kind), INTENT(OUT) :: istatus

    ! Local variables:
    TYPE(BMatrix_t) :: B      ! Static background covariance
    TYPE(BMatrix_t) :: B_e    ! Ensemble background covariance
    !TYPE(RMatrix_t) :: R      ! Observation error covariance
    TYPE(ObsSet_t) :: Y, Z
    TYPE(C2O_t) :: H
    TYPE(state_t) :: ctlBkgd
    TYPE(costFuncGrad_t) cost
    REAL(r_kind) :: t1, t2

    INTEGER(i_kind) :: i

    istatus = 0 ! Normal flag

    ! Pass mask information to model background:
    DO i = 1, UBOUND(bkg%fields, 1)
      bkg%fields(i)%maskHorizontal = ctl%analysis%fields(i)%maskHorizontal
      bkg%fields(i)%maskVertical = ctl%analysis%fields(i)%maskVertical
      bkg%fields(i)%maskTemporal = ctl%analysis%fields(i)%maskTemporal
    END DO

    ! Note: model background needs to convert to control variable bkgd

    ctlBkgd = ctl%bckwardOpr_s(bkg)

    IF (.NOT. sg%isActiveProc()) RETURN

    CALL obs%thinningSG_s(bkg)

    ! Model state variables to fill in data sparse areas
    Z = ObsSet_t(this%yamlVars%yamlFile, obs%mpObs(sg%gLevel))
    ! Thinned obs is saved in obs thinnedObs now:
    CALL this%bkgdFillin_s(sg%gLevel, bkg, obs, Z)

    ! Background covariance:
    CALL B%initialize(this%yamlVars%yamlFile, sg, 'Laplace')

    ! Need to decide if H is needed ... temporarily kept here
    CALL H%initialize(this%yamlVars%yamlFile, ctl%analysis, Z)

    ! Observation covariance:
    ! Use the obs R matrix loader:
    CALL obs%R_Matrices_s(sg)

    ! Define the cost function:
    CALL cost%initialize_s(this%yamlVars%yamlFile, ctlBkgd, &
                            ctl, obs, H, B, B, sg, ctlBkgd)

    ! Synchronize
    CALL sg%mpddInfo_sg%barrier

    ! Solve the variational analysis with increment as control
    IF (this%yamlVars%updateBkgd) THEN
      ctl%analysis = ctl%analysis%zeroCopy()
    ELSE
      ctl%analysis = ctl%analysis - ctlBkgd
    END IF
    CALL CPU_TIME(t1)
    WRITE(*, 2) sg%mpddInfo_sg%myrank, sg%num_cell
2     FORMAT('Starting variationalSolver at proc: ', I1, ' with Ncell: ', I8)
    CALL this%varSolver%solve_s(ctlBkgd, this%yamlVars%weights, ctl%analysis, cost, sg)
    CALL CPU_TIME(t2)
    WRITE (*, 3) t2 - t1, sg%mpddInfo_sg%myrank, sg%num_cell
3     FORMAT('Time spent in variationalSolver: ', D10.2, ' at proc: ', I1, ' with Ncell: ', I8)
    ! Return to full state variable:
    ctl%analysis = ctl%analysis + ctlBkgd

  END SUBROUTINE singleLvlAnal_s

  SUBROUTINE bkgdFillin_s(this, glvl, bkg, obs, Z)
    IMPLICIT NONE
    CLASS(App_4DVAR_t) :: this
    INTEGER(i_kind), INTENT(IN) :: glvl
    TYPE(State_t), INTENT(IN) :: bkg    ! Background
    TYPE(ObsMG_t), INTENT(IN) :: obs      ! Observations
    TYPE(ObsSet_t), INTENT(INOUT) :: Z

    ! Local variables:
    INTEGER(i_kind) :: i, j, k, iu, iv, ix, nin, nv, nh, nt, imax, imin, jmax, jmin, kmax, kmin, nfill, idd
    REAL(i_kind) :: a
    TYPE(State_t) :: X
    TYPE(ObsSet_t) :: Y   ! Save the concatenated obs following the old coding

    INTEGER(i_kind), ALLOCATABLE :: idc(:, :)
    REAL(r_kind), ALLOCATABLE :: O(:)

    X = bkg%zeroCopy()

    ! Check the rampRange:
    WRITE (*, 1) X%sg%gLevel, X%sg%vLevel, this%iRange_ramp(1, :), X%sg%mpddInfo_sg%myrank
    WRITE (*, 2) X%sg%gLevel, X%sg%num_cell, this%iRange_ramp(2, :), X%sg%mpddInfo_sg%myrank
    WRITE (*, 3) X%sg%gLevel, X%sg%tSlots, this%iRange_ramp(3, :), X%sg%mpddInfo_sg%myrank
1   FORMAT('bkgdFillin_s - G', I2, ' Vertical  levels:', I6, ' ramp: ', 20I6)
2   FORMAT('bkgdFillin_s - G', I2, ' Horizon num_cell:', I6, ' ramp: ', 20I6)
3   FORMAT('bkgdFillin_s - G', I2, ' Time frame slots:', I6, ' ramp: ', 20I6)

    IF (MINVAL(this%iRange_ramp(1, :)) .GT. X%sg%vLevel .AND. &
        MINVAL(this%iRange_ramp(2, :)) .GT. X%sg%num_cell .AND. &
        MINVAL(this%iRange_ramp(3, :)) .GT. X%sg%tSlots) THEN
      !Z = Y
      ! Concate the thinned obs for assimilation:
      CALL ObsConcat_s(obs%thinnedObs%obsArray, Z)
      WRITE (*, 4) X%sg%gLevel, X%sg%mpddInfo_sg%myrank
4     FORMAT('No background filled in Glevel:', I3, ' pc', I2)
      RETURN ! No fill in
    ELSE ! Filling is done on the concatenated obs:
      CALL ObsConcat_s(obs%thinnedObs%obsArray, Y)
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
              O(nfill) = bkg%fields(iv)%DATA(i, j, k)
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
        WRITE (*, 5) TRIM(X%fields(iv)%Get_Name()), X%mpddSub%myrank
5       FORMAT('No observation for this state variable: ', A, 1X, ' No need to fillin, pc: ', I4)
      END IF
    END DO

    ! Deallocate local variables:
    DEALLOCATE (O, idc)

  END SUBROUTINE bkgdFillin_s
END MODULE App_4DVAR_m
