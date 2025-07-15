!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2021/1/26, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!!
PROGRAM Test_VelUVW2RWndVar
  USE kinds_m, ONLY: i_kind, r_kind, r_double
  USE RadarRAW_m, ONLY: RadarRAW_t
  USE ObsRadarVel_m, ONLY: ObsRadarVel_t
  USE geometry_m, ONLY: geometry_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE MPObs_m, ONLY: MPObs_t
  USE State_m, ONLY: State_t
  USE ObsSet_m, ONLY: ObsSet_t
  USE C2O_m, ONLY: C2O_t
  USE MiniSolver_m, ONLY: MiniSolver_t
  USE SolverLBFGS_m, ONLY: SolverLBFGS_t
  USE SolverFRCG_m, ONLY: SolverFRCG_t
  USE JFunc_m, ONLY: JFunc_t
  USE BMatrix_m, ONLY: BMatrix_t
  USE RMatrix_m, ONLY: RMatrix_t
  USE State2NC_m
  USE MGOpts_m
  USE IOGrapes_m, ONLY: IOGrapes_t
  USE YAMLRead_m
  USE Obs2State_m

  IMPLICIT NONE

  ! Define types
  TYPE(mpddGlob_t), TARGET :: mpddGlob
  TYPE(geometry_t), TARGET :: geometry
  TYPE(JFunc_t) :: JFunc
  TYPE(IOGrapes_t), TARGET :: ioGrapes

  CHARACTER(LEN=1024) :: inputDir
  CHARACTER(LEN=1024) :: outputDir
  CHARACTER(len=1024) :: filenameInput
  CHARACTER(len=1024) :: configFile
  CHARACTER(LEN=1024) :: staticDir
  TYPE(ObsRadarVel_t) :: ObsRadarVel
  CLASS(MiniSolver_t), ALLOCATABLE :: miniSolver
  TYPE(BMatrix_t) :: B
  INTEGER(i_kind) :: i, j, k, mgStart, mgEnd
  TYPE(State_t), ALLOCATABLE :: XbMG(:)
  TYPE(State_t) :: XRes
  LOGICAL :: IncFlag = .TRUE.  !.True., .False.

  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", staticDir)
  configFile = TRIM(staticDir)//"/UnitTest"//"/Test_VelUV2RWndVar.yaml"

  IF (yaml_get_var(configFile, 'IO', 'output_dir', outputDir) /= 0) STOP

  ! Initializer
  ! Auxtypes
  CALL mpddGlob%initialize()                                   ! Initialize the mpdd

  ! Initialize geometry
  CALL geometry%initialize(configFile, mpddGlob)               ! Initialize the geometry
  ALLOCATE (XbMG(geometry%mg%mg_coarsest:geometry%mg%mg_finest))

  ! Initialize solver
  ALLOCATE (SolverFRCG_t::miniSolver)
  SELECT TYPE (miniSolver)
  TYPE IS (SolverFRCG_t)
    CALL miniSolver%initialize(configFile)
  TYPE IS (SolverLBFGS_t)
    CALL miniSolver%initialize(configFile)
  END SELECT

  mgStart = 4      ! start index of multigrid
  mgEnd = 7        ! end index of multigrid

  DO i = mgEnd, mgStart, -1
    PRINT *, 'i', i
    IF (i .EQ. mgEnd) THEN
      ! Initialize a zeros background field at finest grid.
      ASSOCIATE (sg => geometry%mg%sg(mgEnd))

        CALL XbMG(mgEnd)%initialize(configFile, sg)
        CALL ioGrapes%initialize(configFile, geometry)

        ! ! Read the grapesinput into the Xm
        CALL ioGrapes%m_read_bcg_into_Xm(XbMG(mgEnd), sg)

      END ASSOCIATE
    ELSE
      ! Restrict to each coarser grid
      ASSOCIATE (sgFiner => geometry%mg%sg(i + 1), sgCoarser => geometry%mg%sg(i))
        CALL XbMG(i)%initialize(configFile, sgCoarser)
        CALL restrictionMG(XbMG(i), XbMG(i + 1), geometry%mg)
      END ASSOCIATE
    END IF
  END DO

  ! Give the initial value at the coarsest grid.
  XRes = XbMG(mgStart)

  ! Initial
  CALL ObsRadarVel%ObsInitial(configFile)
  CALL ObsRadarVel%ObsIngest(XRes)

  ! MultiGrid
  DO i = mgStart, mgEnd
    ASSOCIATE (sg => geometry%mg%sg(i))
      BLOCK

        TYPE(MPObs_t), TARGET :: mpObs
        TYPE(State_t) :: X, DD
        TYPE(ObsSet_t) :: Y
        TYPE(C2O_t) :: H
        TYPE(RMatrix_t) :: R

        CALL mpObs%initializeMPObs(sg)    ! Initialize the observation parallel processing proc
        X = XRes                ! Initialize the state

        Y = ObsSet_t(configFile, mpObs)
        CALL B%initialize(configFile, sg)     ! Initialize the B matrix

        CALL ObsRadarVel%ObsPrepareForSg(X)
        ! CALL ObsRadarVel%ObsSuper(X, Y, mpObs)
        CALL ObsRadarVel%ObsThinning(X, Y, mpObs, .TRUE., .FALSE.)

        !       DD= Obs2State_BaseTypeName(sg, Y)
        ! CALL Output_NC_State_AV(DD, outputDir, "TestRadarRef"//TRIM(DD%Fields(1)%Get_Name()))

        CALL R%initialize(configFile, Y, X%sg) ! Initialize R
        CALL H%initialize(configFile, X, Y)

        ! Initialize J Function
        CALL JFunc%initialize(configFile, X, Y, H, B, R, sg)

        ! Run minimization
        CALL miniSolver%run(X, JFunc, sg)

        ! Sync
        CALL mpddGlob%barrier

        ! Connect to coarser grid.
        IF (i .NE. mgEnd) THEN
          ! ===========================> prolongate to finer grid
          ASSOCIATE (sgFiner => geometry%mg%sg(i + 1))
            CALL XRes%initialize(configFile, sgFiner)
            CALL prolongationMG(X, XRes, XbMG(i + 1), geometry%mg, IncFlag)
          END ASSOCIATE
        ELSE
          XRes = X ! Return the state fields directly
          CALL H%UV2W%fwdNL(XRes)
        END IF

        ! PRINT *, 'U:', XRes%fields(1)%data(10, 6435, 1)
        ! PRINT *, 'V:', XRes%fields(2)%data(10, 6435, 1)
        ! PRINT *, 'W:', XRes%fields(3)%data(10, 6435, 1)
        ! PRINT *, ''

        !   DO k = 1, SIZE(Y%ObsFields)
        !     DO j = 1, SIZE(Y%ObsFields(k)%values)
        !       IF (Y%ObsFields(k)%idx(j)%vIdx == 10 .AND. Y%ObsFields(k)%idx(j)%hIdx == 6435 .AND. &
        !           Y%ObsFields(k)%idx(j)%tIdx == 1) THEN
        !         PRINT *, 'Radar', k, 'radial wind at this point: '
        !         PRINT *, 'Thining value is:', Y%ObsFields(k)%values(j)
        !         PRINT *, 'UVW projection coefficients:', H%OprRadarVel%attrRadars(k)%factor(:, j)
        !  PRINT *, 'Radial wind cal from U, V, W: ', sum(H%OprRadarVel%attrRadars(k)%factor(:, j)*(/XRes%fields(1)%data(10, 6435, 1), &
        !                                                          XRes%fields(2)%data(10, 6435, 1), XRes%fields(3)%data(10, 6435, 1)/))
        !         PRINT *, ''
        !       END IF
        !     END DO
        !   END DO

      END BLOCK
    END ASSOCIATE
  END DO

  CALL Output_NC_State_AV(XRes, outputDir, "Test_VelUV2RWndVar")

  ! Destory
  CALL ObsRadarVel%ObsDestroy()
  PRINT *, 'It is done'

  CALL mpddGlob%finalize
END PROGRAM Test_VelUVW2RWndVar
