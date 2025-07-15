!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-QC.ObsConvention
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Jiongming Pang
! VERSION           : V 0.0
! HISTORY           : 2022-01-20, created by Jiongming Pang for test of thinning Obsercations.
!
!   Created by Jiongming Pang for Vwpw obs (pang.j.m@hotmail.com), 2022/1/20, @GBA-MWF, Shenzhen
!
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This program implements a test for reading observation.
!
PROGRAM Test_ObsConcat

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
  USE JFunc_m, ONLY: JFunc_t
  USE State2NC_m
  USE Mock_m
  USE MGOpts_m
  USE parameters_m, ONLY: degree2radian
  !USE NMLRead_m
  USE YAMLRead_m

  USE ObsBase_m, ONLY: ObsBase_t
  USE ObsSurface_m, ONLY: ObsSurface_t
  USE ObsSound_m, ONLY: ObsSound_t
  USE ObsVwpw_m, ONLY: ObsVwpw_t

  USE ObsUtilities_m
  USE Obs2State_m

  IMPLICIT NONE

  ! Define types
  TYPE(mpddGlob_t), TARGET :: mpddGlob
  TYPE(geometry_t), TARGET :: geometry
  TYPE(MPObs_t), TARGET :: mpObs
  TYPE(State_t) :: XRes, XTru
  TYPE(State_t), ALLOCATABLE :: XbMG(:)
  TYPE(C2O_t) :: H
  REAL(r_kind) :: t1, t2
  INTEGER(i_kind) :: mgStart, mgEnd

  CLASS(MiniSolver_t), ALLOCATABLE :: miniSolver

  CHARACTER(LEN=1024) :: configFile, ncOutputFile !< Config and output file name.

  INTEGER(i_kind) :: i, iii, imx, iv
  REAL(r_kind) :: amx

  TYPE(ObsSurface_t) :: OBS1 !OBS_surface
  TYPE(ObsSound_t) :: OBS2 !OBS_sound
  TYPE(ObsVwpw_t) :: OBS3 !OBS_vwpw

  LOGICAL :: IncFlag = .TRUE.
  REAL(r_kind), ALLOCATABLE :: xyt(:, :)

  CHARACTER(LEN=20), ALLOCATABLE :: varList(:)
  CHARACTER(LEN=20) :: varName
  INTEGER(i_kind)   :: varNum
  INTEGER(i_kind)   :: istatus

  ! Get the configFile
  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", configFile)
  configFile = TRIM(configFile)//"/UnitTest/testMiniSolver_newObs.yaml"

  istatus = yaml_get_var(configFile, 'IO', 'output_dir', ncOutputFile)
  ncOutputFile = TRIM(ncOutputFile)

  !CALL namelist_read(TRIM(configFile), "varList", varList)
  istatus = yaml_get_var(TRIM(configFile), 'modelState', 'varList', varList)
  varNum = UBOUND(varList, 1)
  PRINT *, "num of analysis vars: ", varNum

  ! Initializer
  ! Auxtypes
  CALL mpddGlob%initialize()                                   ! Initialize the mpdd

  CALL mpddGlob%barrier
  CALL CPU_TIME(t1)

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
!  CALL miniSolver%initialize(configFile)

  !CALL namelist_read(TRIM(configFile), "mgStart", mgStart)
  !CALL namelist_read(TRIM(configFile), "mgEnd", mgEnd)
  istatus = yaml_get_var(TRIM(configFile), 'geometry', 'mgStart', mgStart)
  istatus = yaml_get_var(TRIM(configFile), 'geometry', 'mgEnd', mgEnd)

  DO i = mgEnd, mgStart, -1
    PRINT *, 'i', i
    IF (i .EQ. mgEnd) THEN
      ! Initialize a zeros background field at finest grid.
      ASSOCIATE (sg => geometry%mg%sg(mgEnd))
        CALL XbMG(mgEnd)%initialize(configFile, sg)

      END ASSOCIATE
    ELSE
      ! Restrict to each coarser grid
      ASSOCIATE (sgFiner => geometry%mg%sg(i + 1), sgCoarser => geometry%mg%sg(i))
        CALL XbMG(i)%initialize(configFile, sgCoarser)
        CALL restrictionMG(XbMG(i), XbMG(i + 1), geometry%mg)
      END ASSOCIATE
    END IF
    ! XbMG(i)%fields(1)%data = 0.8D0
    DO iv = 1, varNum
      XbMG(i)%Fields(iv)%DATA = 0.0D0
    END DO
  END DO

  ! Give the initial value at the coarsest grid.
  XRes = XbMG(mgStart)

  CALL OBS1%ObsInitial(configFile)
  CALL OBS1%ObsIngest(XbMG(mgEnd))

  CALL OBS2%ObsInitial(configFile)
  CALL OBS2%ObsIngest(XbMG(mgEnd))

  CALL OBS3%ObsInitial(configFile)
  CALL OBS3%ObsIngest(XbMG(mgEnd))

  ! MultiGrid
  DO i = mgStart, mgEnd
    ! Run 3DVAR in each single grid.
    BLOCK
      TYPE(State_t) :: X
      TYPE(BMatrix_t) :: B
      TYPE(MPObs_t), TARGET :: mpObs
      TYPE(ObsSet_t) :: Y1
      TYPE(ObsSet_t) :: Y2
      TYPE(ObsSet_t) :: Y3
      TYPE(ObsSet_t) :: Y
      TYPE(C2O_t) :: H
      TYPE(JFunc_t) :: JFunc

      ASSOCIATE (sg => geometry%mg%sg(i))
        CALL mpObs%initializeMPObs(sg)    ! Initialize the observation parallel processing proc

        ! Initialize X, Y, B, H
        X = Xres     ! Initialize the state
        CALL B%initialize(configFile, sg, 'Laplace')     ! Initialize the B matrix

        ! Thinning observations:
        CALL OBS1%ObsThinning(XbMG(i), Y1, mpObs, .TRUE., .FALSE.)
        CALL OBS2%ObsThinning(XbMG(i), Y2, mpObs, .TRUE., .FALSE.)
        CALL OBS3%ObsThinning(XbMG(i), Y3, mpObs, .TRUE., .FALSE.)

        ! Concat thinning obs:
        CALL ObsConcat_s((/Y1, Y2, Y3/), Y)

        CALL mpddGlob%barrier

        PRINT *, 'Max/min super obs values: ', sg%mpddInfo_sg%myRank, MAXVAL(Y%ObsFields(1)%values(:)), MINVAL(Y%ObsFields(1)%values(:))

        PRINT *, 'Number of thinned obs: ', UBOUND(Y%ObsFields(1)%values, 1)

        CALL H%initialize(configFile, X, Y)

        BLOCK
          TYPE(State_t) :: DD
          CHARACTER(LEN=128) :: temp

          WRITE (temp, "(I2.2)") i
          DD = Obs2State_BaseTypeName(sg, Y)

          DO iv = 1, SIZE(DD%Fields)
            varName = DD%Fields(iv)%Get_Name()
            CALL Output_NC_State_SV(DD, ncOutputFile, "MGTest_ConcatOBS_G"//TRIM(temp)//"_"//TRIM(varName), TRIM(varName))
          END DO
        END BLOCK

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

      END ASSOCIATE
    END BLOCK
  END DO

  CALL mpddGlob%barrier

  CALL CPU_TIME(t2)
  IF (mpddGlob%isBaseProc()) PRINT *, 'Time cost:', t2 - t1, 's'

  DEALLOCATE (miniSolver)

  ! Finalize
  CALL mpddGlob%finalize

END PROGRAM Test_ObsConcat
