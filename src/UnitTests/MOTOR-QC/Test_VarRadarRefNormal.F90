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
PROGRAM Test_RadarRefVar
  USE kinds_m, ONLY: i_kind, r_kind, r_double
  USE RadarRAW_m, ONLY: RadarRAW_t
  USE ObsRadarRef_m, ONLY: ObsRadarRef_t
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
  USE RhoRCtl2RhoR_m, ONLY: RhoRCtl2RhoR_t
  USE OprRadarRef_m, ONLY: OprRadarRef_t

  IMPLICIT NONE

  ! Define types
  TYPE(mpddGlob_t), TARGET :: mpddGlob
  TYPE(geometry_t), TARGET :: geometry
  TYPE(JFunc_t) :: JFunc
  TYPE(IOGrapes_t), TARGET :: ioGrapes

  CHARACTER(LEN=1024) :: inputDir
  CHARACTER(LEN=1024) :: outputDir
  CHARACTER(len=1024) :: filenameInput, filenameOutput
  CHARACTER(len=1024) :: configFile
  CHARACTER(LEN=1024) :: staticDir
  TYPE(ObsRadarRef_t) :: ObsRadarRef
  CLASS(MiniSolver_t), ALLOCATABLE :: miniSolver
  TYPE(BMatrix_t) :: B
  INTEGER(i_kind) :: i, j, k, mgStart, mgEnd
  TYPE(State_t), ALLOCATABLE :: XbMG(:)
  TYPE(State_t) :: XRes
  LOGICAL :: IncFlag = .TRUE.  !.True., .False.

  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", staticDir)
  configFile = TRIM(staticDir)//"/UnitTest"//"/testRadarRAWNormal.yaml"

  IF (yaml_get_var(configFile, 'IO', 'output_dir', outputDir) /= 0) STOP
  filenameOutput = TRIM(outputDir)//"/testRadar.nc"

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

  mgStart = 7      ! start index of multigrid
  mgEnd = 7        ! end index of multigrid

  DO i = mgEnd, mgStart, -1
    PRINT *, 'i', i
    IF (i .EQ. mgEnd) THEN
      ! Initialize a zeros background field at finest grid.
      ASSOCIATE (sg => geometry%mg%sg(mgEnd))

        CALL XbMG(mgEnd)%initialize(configFile, sg)
        CALL ioGrapes%initialize(configFile, geometry)

        ! Read the grapesinput into the Xm
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
  CALL ObsRadarRef%ObsInitial(configFile)
  CALL ObsRadarRef%ObsIngest(XbMG(mgEnd))

  ! MultiGrid
  DO i = mgStart, mgEnd
    ASSOCIATE (sg => geometry%mg%sg(i))
      BLOCK

        TYPE(MPObs_t), TARGET :: mpObs
        TYPE(State_t) :: X, DD
        TYPE(ObsSet_t) :: Y
        TYPE(C2O_t) :: H
        TYPE(RMatrix_t) :: R
        TYPE(RhoRCtl2RhoR_t) :: RhoRCtl2RhoR
        TYPE(OprRadarRef_t) :: OprRadarRef

        CALL mpObs%initializeMPObs(sg)    ! Initialize the observation parallel processing proc
        X = XRes                ! Initialize the state

        Y = ObsSet_t(configFile, mpObs)
        CALL B%initialize(configFile, sg)     ! Initialize the B matrix

        CALL H%initialize(configFile, X, Y)

        CALL ObsRadarRef%ObsPrepareForSg(X)
        CALL ObsRadarRef%ObsSuper(X, Y, mpObs)
        !       DD= Obs2State_BaseTypeName(sg, Y)
        ! CALL Output_NC_State_AV(DD, outputDir, "TestRadarRef"//TRIM(DD%Fields(1)%Get_Name()))

        ! Initialize the qrain
        IF (i .EQ. mgStart) THEN

          RhoRCtl2RhoR = RhoRCtl2RhoR_t(configFile, X)
          OprRadarRef = OprRadarRef_t(configFile, X)

          CALL RhoRCtl2RhoR%transFwdNonLinear(X)

          CALL OprRadarRef%calRhoRFromRefRho(X, Y)
          PRINT *, 'MAXVAL: ', MAXVAL(X%Fields(1)%DATA)
          X%Fields(1)%DATA = 0.00001

          CALL RhoRCtl2RhoR%transBackward(X)
        END IF

        CALL R%initialize(configFile, Y, X%sg) ! Initialize R

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
        END IF
      END BLOCK
    END ASSOCIATE
  END DO

  BLOCK
    USE RhoRCtl2RhoR_m, ONLY: RhoRCtl2RhoR_t
    TYPE(RhoRCtl2RhoR_t) :: RhoRCtl2RhoR

    RhoRCtl2RhoR = RhoRCtl2RhoR_t(configFile, XRes)

    CALL RhoRCtl2RhoR%transFwdNonLinear(XRes)
  END BLOCK

  CALL Output_NC_State_AV(XRes, outputDir, "TestRadarVarNormal_G9-9"//TRIM(XRes%Fields(1)%Get_Name()))

  ! Destory
  CALL ObsRadarRef%ObsDestroy()
  PRINT *, 'It is done'

  CALL mpddGlob%finalize
END PROGRAM Test_RadarRefVar
