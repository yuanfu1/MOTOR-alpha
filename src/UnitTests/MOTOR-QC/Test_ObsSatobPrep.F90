!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-QC.ObsConvention
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Hua Zhang
! VERSION           : V 0.0
! HISTORY           : 2022-10-05, created by Hua Zhang for test of reading Satob Obsercations.
!
!   Created by Hua Zhang for Satob obs (zhangh@cma.gov.cn), 2022/10/05, @GBA-MWF, Beijing
!
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This program implements a test for reading observation.
!
PROGRAM Test_ObsSatobRead

  USE kinds_m, ONLY: i_kind, r_kind, r_double
  USE ObsBase_m, ONLY: ObsBase_t
  USE State_m, ONLY: State_t
  USE geometry_m, ONLY: geometry_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE MPObs_m, ONLY: MPObs_t
  USE ObsSet_m, ONLY: ObsSet_t
  USE C2O_m, ONLY: C2O_t
  USE NMLRead_m
  USE State2NC_m
  USE ObsSurface_m, ONLY: ObsSurface_t
  USE ObsSound_m, ONLY: ObsSound_t
  USE ObsVwpw_m, ONLY: ObsVwpw_t
  USE Meteoro_constants
  USE Meteoro_type_define
  USE extreme_value_m
  USE Satob_m, ONLY: Satob_t
  USE GtsQC_m
  USE interp_m
  USE YAMLRead_m
  USE IOModel_m, ONLY: IOModel_t
  USE IOGrapes_m, ONLY: IOGrapes_t

  ! Yuanfu Xie added the following for creating a multigrid background:
  USE MGOpts_m

  IMPLICIT NONE

  CHARACTER(LEN=1024) :: configFile, outputDir, inputDir
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
  TYPE(ObsSurface_t) :: OBS
  TYPE(Satob_t) :: ObsSatob
  ! TYPE(ObsSound_t) :: OBS
  ! TYPE(ObsVwpw_t) :: OBS
  TYPE(C2O_t) :: H
  TYPE(State_t), ALLOCATABLE :: XbMG(:)

  TYPE(IOGrapes_t) :: IOModel

  ! Yuanfu Xie modified the test drive to include the XbMG for the modified obs ingest:
  CHARACTER(LEN=20), ALLOCATABLE :: varList(:)
  INTEGER(i_kind) :: iv, mgStart, mgEnd, varNum, istatus
  REAL(r_kind) :: t1, t2

!   CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", inputDir)
!   !inputDir='/sources/MOTOR/static/UnitTest/'
!   configFile = TRIM(inputDir)//"/UnitTest/testSatob.nl"
!   PRINT *,'Config: ',TRIM(configFile)
  CALL getarg(1, configFile)
  IF (TRIM(configFile) .EQ. '') THEN
    PRINT *, TRIM(configFile)

    CALL getarg(0, configFile)
    PRINT *, TRIM(configFile)
    WRITE (*, 1)
1   FORMAT('Usage of this driver: mpirun -n <n> Debug/Test_ObsSatobPrep.exe configFile', /, &
           ' Check your configure file and rerun!')
    configFile = './App_3DVarVerification.yaml'
  ELSE
    WRITE (*, 2) TRIM(configFile)
2   FORMAT('ConfigFile is: ', A)
  END IF

  CALL GET_ENVIRONMENT_VARIABLE("OUTPUT_DIR", outputDir)
  PRINT *, "outputDir: ", TRIM(outputDir)

  !CALL namelist_read(TRIM(configFile), "varList", varList)
  !varNum = UBOUND(varList, 1)
  istatus = yaml_get_var(TRIM(configFile), 'modelState', 'varList', varList)
  IF (istatus .NE. 0) THEN
    varNum = 0
  ELSE
    varNum = UBOUND(varList, 1)
    WRITE (*, 11) varNum
11  FORMAT("Number of analysis vars: ", I2)
  END IF

  ! Initializer
  ! Auxtypes
  ! mpddGlob = mpddGlob_t()                                   ! Initialize the mpdd
  CALL mpddGlob%initialize()

  CALL mpddGlob%barrier
  CALL CPU_TIME(t1)

  ! Initialize geometry
  ! geometry = geometry_t(configFile, mpddGlob)               ! Initialize the geometry
  CALL geometry%initialize(configFile, mpddGlob)
  ALLOCATE (XbMG(geometry%mg%mg_coarsest:geometry%mg%mg_finest))
  CALL XbMG(geometry%mg%mg_finest)%initialize(configFile, geometry%mg%sg(geometry%mg%mg_finest))
  CALL XbMG(geometry%mg%mg_finest)%ShowInfo

  !ALLOCATE (IOModel)
  CALL IOModel%initialize(configFile, geometry)
  CALL IOModel%m_read_bcg_into_Xm(XbMG(geometry%mg%mg_finest), geometry%mg%sg(geometry%mg%mg_finest))

  PRINT *, 'SG-latlon: ', XbMG(geometry%mg%mg_finest)%sg%minLatGlob / degree2radian, &
    XbMG(geometry%mg%mg_finest)%sg%maxLatGlob / degree2radian, &
    XbMG(geometry%mg%mg_finest)%sg%minLonGlob / degree2radian, &
    XbMG(geometry%mg%mg_finest)%sg%maxLonGlob / degree2radian
  PRINT *, 'Start satob_setup'
  CALL ObsSatob%ObsInitial(configFile)
  PRINT *, 'Start satob_read'
  CALL ObsSatob%ObsIngest(XbMG(geometry%mg%mg_finest))

  IF (ObsSatob%numObs .GT. 0) THEN
    PRINT *, 'ObsSatob successfully processed!'
  ELSE
    WRITE (*, 3) XbMG(geometry%mg%mg_finest)%mpddGlob%myrank
3   FORMAT('No cloud drift wind processed at proc: ', I3)
  END IF

  ! Finalize
  CALL mpddGlob%finalize

END PROGRAM Test_ObsSatobRead
