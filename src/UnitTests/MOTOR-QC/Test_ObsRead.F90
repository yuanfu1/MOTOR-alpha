!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-QC.ObsConvention
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Jiongming Pang
! VERSION           : V 0.0
! HISTORY           : 2022-01-20, created by Jiongming Pang for test of reading Obsercations.
!
!   Created by Jiongming Pang for Vwpw obs (pang.j.m@hotmail.com), 2022/1/20, @GBA-MWF, Shenzhen
!
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This program implements a test for reading observation.
!
PROGRAM Test_ObsRead

  USE kinds_m, ONLY: i_kind, r_kind, r_double
  USE ObsBase_m, ONLY: ObsBase_t
  USE State_m, ONLY: State_t
  USE geometry_m, ONLY: geometry_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE MPObs_m, ONLY: MPObs_t
  USE ObsSet_m, ONLY: ObsSet_t
  USE C2O_m, ONLY: C2O_t
  !USE NMLRead_m
  USE YAMLRead_m
  USE State2NC_m
  USE ObsSurface_m, ONLY: ObsSurface_t
  USE ObsSound_m, ONLY: ObsSound_t
  USE ObsVwpw_m, ONLY: ObsVwpw_t
  ! Yuanfu Xie added the following for creating a multigrid background:
  USE MGOpts_m
  USE YAMLRead_m

  IMPLICIT NONE

  CHARACTER(LEN=1024) :: configFile, outputDir
  CHARACTER(LEN=1024), ALLOCATABLE :: filename(:)
  CHARACTER(LEN=20), ALLOCATABLE :: varname(:)
  CHARACTER(LEN=128) :: temp
  INTEGER(i_kind) :: numVars, i, obsFilesNum
  TYPE(State_t) :: state
  TYPE(State_t) :: state2
  TYPE(mpddGlob_t), TARGET :: mpddGlob
  TYPE(geometry_t), TARGET :: geometry
  TYPE(MPObs_t), TARGET :: mpObs
  TYPE(ObsSet_t) :: obsThinned
  TYPE(ObsSet_t) :: obsThinned2
  ! TYPE(ObsSurface_t) :: OBS
  ! TYPE(ObsSound_t) :: OBS
  TYPE(ObsVwpw_t) :: OBS
  TYPE(C2O_t) :: H
  TYPE(State_t), ALLOCATABLE :: XbMG(:)
  INTEGER(i_kind) :: fileStat

  ! Yuanfu Xie modified the test drive to include the XbMG for the modified obs ingest:
  CHARACTER(LEN=20), ALLOCATABLE :: varList(:)
  INTEGER(i_kind) :: iv, mgStart, mgEnd, varNum, ifile
  REAL(r_kind) :: t1, t2

  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", configFile)
  ! configFile = TRIM(configFile)//"/UnitTest/testMiniSolver_newObs.yaml"
  configFile = TRIM(configFile)//"/UnitTest/testMiniSolver_newObs_RVWPW.yaml"
  PRINT *, 'Config: ', TRIM(configFile)

  fileStat = yaml_get_var(configFile, 'IO', 'output_dir', outputDir)
  PRINT *, "outputDir: ", outputDir

  !CALL namelist_read(TRIM(configFile), "varList", varList)
  ifile = yaml_get_var(TRIM(configFile), 'modelState', 'varList', varList)
  varNum = UBOUND(varList, 1)

  ! Initializer
  ! Auxtypes
  CALL mpddGlob%initialize()                                   ! Initialize the mpdd

  CALL mpddGlob%barrier
  CALL CPU_TIME(t1)

  ! Initialize geometry
  CALL geometry%initialize(configFile, mpddGlob)               ! Initialize the geometry
  ALLOCATE (XbMG(geometry%mg%mg_coarsest:geometry%mg%mg_finest))
  ifile = yaml_get_var(TRIM(configFile), 'geometry', 'mgStart', mgStart)
  ifile = yaml_get_var(TRIM(configFile), 'geometry', 'mgEnd', mgEnd)

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
    DO iv = 1, varNum
      XbMG(i)%Fields(iv)%DATA = 0.0D0
    END DO
  END DO

  CALL OBS%ObsInitial(configFile)
  CALL OBS%ObsIngest(XbMG(mgEnd))

END PROGRAM Test_ObsRead
