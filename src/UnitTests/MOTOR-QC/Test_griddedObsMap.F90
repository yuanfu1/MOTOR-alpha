!!--------------------------------------------------------------------------------------------------
! PROJECT           : Gridded obs Test
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           : 2022-06-21, created by Yuanfu Xie.
!
!   Created by Yuanfu Xie (yuanfu_xie@yahoo.com), 2022/06/21, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------
!!
!> @brief
!! This is a test of MOTOR-DA gridded obs mapping
!
PROGRAM Test_griddedObsMap

  USE geometry_m, ONLY: geometry_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE ObsBase_m
  USE ObsSurface_m, ONLY: ObsSurface_t
  USE ObsSet_m, ONLY: ObsSet_t
  USE YAMLRead_m
  USE Obs2State_m
  USE State2NC_m

  ! Added by TS@230131
  USE obsTools_m, ONLY: griddedObsMap

  IMPLICIT NONE

  INCLUDE "mpif.h"

  CHARACTER(LEN=1024) :: configFile, ncOutputFile
  INTEGER(i_kind) :: mgStart = 1, mgEnd = 5
  INTEGER(i_kind) :: istatus, i
  CLASS(SingleGrid_t), POINTER :: sgFine, sgCoarse
  TYPE(ObsSet_t) :: yFine
  TYPE(ObsSet_t) :: yCoarse
  CLASS(ObsBase_t), POINTER :: obs
  TYPE(ObsSurface_t) :: sfc
  TYPE(mpddGlob_t), TARGET :: mpddGlob
  TYPE(geometry_t), TARGET :: geometry
  INTEGER(i_kind)   :: ifile, iv
  CHARACTER(LEN=20), ALLOCATABLE :: varList(:)
  CHARACTER(LEN=20) :: varName
  INTEGER(i_kind)   :: varNum
  INTEGER(i_kind), ALLOCATABLE :: NRanks(:)
  TYPE(State_t) :: X
  TYPE(MPObs_t), TARGET, ALLOCATABLE :: mpObs(:)
  TYPE(MPObs_t), POINTER :: mpObsFine, mpObsCoarse
  TYPE(State_t) :: dCoarse, dFine

  ! Get the configFile
  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", configFile)
  configFile = TRIM(configFile)//"/UnitTest/testMiniSolver_newObs.yaml"

  ifile = yaml_get_var(configFile, 'IO', 'output_dir', ncOutputFile)
  ncOutputFile = TRIM(ncOutputFile)

  ifile = yaml_get_var(TRIM(configFile), 'modelState', 'varList', varList)
  varNum = UBOUND(varList, 1)
  PRINT *, "num of analysis vars: ", varNum

  ! Initializer
  ! Auxtypes
  CALL mpddGlob%initialize()
  PRINT *, 'Total procs: ', mpddGlob%nProc
  ALLOCATE (NRanks(mpddGlob%nProc))

  istatus = yaml_get_var(TRIM(configFile), 'geometry', 'mgStart', mgStart)
  istatus = yaml_get_var(TRIM(configFile), 'geometry', 'mgEnd', mgEnd)
  IF (mpddGlob%isBaseProc()) WRITE (*, 3) mgStart, mgEnd
3 FORMAT("mgStart and mgEnd: ", 2I2)

  ! Allocate mpObs:
  ALLOCATE (mpObs(mgStart:mgEnd))

  ! Initialize geometry
  CALL geometry%initialize(configFile, mpddGlob)

  DO i = mgStart, mgEnd
    CALL mpObs(i)%initializeMPObs(geometry%mg%sg(i))
  END DO

  sgFine => geometry%mg%sg(mgEnd - 1)
  sgCoarse => geometry%mg%sg(mgEnd - 2)
  mpObsFine => mpObs(mgEnd - 1)
  mpObsCoarse => mpObs(mgEnd - 2)
  DO i = 1, UBOUND(sgCoarse%c_tc_idx, 2)
    WRITE (*, 1) i, sgCoarse%c_tc_idx(:, i), sgCoarse%mpddInfo_sg%myrank
1   FORMAT('MAP Coarse check c_tc_idx at ', I3, ' idx: ', 4I6, ' proc: ', I2)
  END DO
  WRITE (*, 5) sgCoarse%mpddInfo_sg%myrank, sgCoarse%num_cell, sgCoarse%num_icell
5 FORMAT('Num_cell Coarse: ', I2, 2I4)
  WRITE (*, 6) sgFine%mpddInfo_sg%myrank, sgFine%num_cell, sgFine%num_icell, UBOUND(sgFine%c_tc_idx, 2)
6 FORMAT('Num_cell fine: ', I2, 2I4, ' Size of c_tc_idx: ', I6)
  WRITE (*, 7) geometry%mg%sg(3)%mpddInfo_sg%myrank, geometry%mg%sg(3)%num_cell, geometry%mg%sg(3)%num_icell
7 FORMAT('Num_cell 3: ', I2, 2I4)
  DO i = 1, UBOUND(sgFine%c_tc_idx, 2)
    WRITE (*, 2) i, sgFine%c_tc_idx(:, i), sgFine%mpddInfo_sg%myrank
2   FORMAT('MAP Fine check c_tc_idx at ', I3, ' idx: ', 4I6, ' proc: ', I2)
  END DO
  ! No use of these arrays:
  PRINT *, 'Allocated p_tc_idx? : ', ALLOCATED(sgFine%p_tc_idx), ALLOCATED(sgFine%sp_t_g_idx_toFiner)

  WRITE (*, 8) sgFine%mpddInfo_sg%myrank, UBOUND(sgFine%sp_t_g_idx_toFiner, 1), sgFine%sp_t_g_idx_toFiner(1:20)
8 FORMAT('sp_t_g_idx_toFiner at proc: ', 2I2, ' idx: ', 20I4)
  WRITE (*, 9) sgFine%mpddInfo_sg%myrank, UBOUND(sgFine%g_t_sp_idx_toFiner, 1), sgFine%g_t_sp_idx_toFiner
9 FORMAT('g_t_sp_idx_toFiner at proc: ', 2I4, ' idx: ', 36I3)
  WRITE (*, 10) sgFine%mpddInfo_sg%myrank, UBOUND(sgFine%sp_t_g_idx_toCoarser, 1), sgFine%sp_t_g_idx_toCoarser
10 FORMAT('sp_t_g_idx_toCoarser at proc: ', I2, I4' idx: ', 20I4)
  WRITE (*, 11) sgFine%mpddInfo_sg%myrank, UBOUND(sgFine%g_t_sp_idx_toCoarser, 1), sgFine%g_t_sp_idx_toCoarser
11 FORMAT('g_t_sp_idx_toCoarser at proc: ', I2, I4, ' idx: ', 36I3)
  WRITE (*, 12) sgFine%mpddInfo_sg%myrank, UBOUND(sgFine%sp_t_g_idx, 1), sgFine%sp_t_g_idx
12 FORMAT('sp_t_g_idx at           proc: ', I2, I6, ' idx: ', 20I4)
  WRITE (*, 13) sgFine%mpddInfo_sg%myrank, UBOUND(sgFine%g_t_sp_idx, 1), sgFine%g_t_sp_idx(1)
13 FORMAT('g_t_sp_idx at           proc: ', I2, I6, ' idx: ', 36I3)

  CALL MPI_ALLGATHER(sgFine%mpddInfo_sg%myrank, 1, MPI_INTEGER, Nranks, 1, MPI_INTEGER, sgFine%mpddInfo_sg%comm, istatus)

  PRINT *, 'NRanks: ', NRanks, sgCoarse%mpddInfo_sg%comm, sgFine%mpddInfo_sg%comm

  CALL X%initialize(configfile, sgFine)

  CALL sfc%ObsInitial(configFile)
  CALL sfc%ObsIngest(X)
  CALL sfc%ObsThinning(X, yFine, mpObsFine, .TRUE., .FALSE.)

  dFine = Obs2State_BaseName(sgFine, yFine)
  DO iv = 1, SIZE(dFine%Fields)
    WRITE (*, 23) iv, sgFine%mpddInfo_sg%myrank, sgFine%num_cell, MAXVAL(ABS(dFine%fields(iv)%DATA(:, :, :)))
23  FORMAT('MaxABS of combined thinned obs of var: ', I1, ' at proc: ', I1, I8, ' with value: ', D12.4)
    CALL Output_NC_State_SV(dFine, ncOutputFile, &
                            "testMO_obsThinnedMap_"//TRIM(dFine%Fields(iv)%Get_Name()), TRIM(dFine%Fields(iv)%Get_Name()), .TRUE.)
  END DO

  ! yCoarse = ObsSet_t(configfile, mpobsCoarse)
  PRINT *, 'ObsFields Allocated? ', ALLOCATED(yCoarse%ObsFields)
  ! Modified by TS@230131
  ! CALL sfc%griddedObsMap(sgFine,sgCoarse,yFine,yCoarse,mpobsCoarse)
  CALL griddedObsMap(sfc, sgFine, sgCoarse, yFine, yCoarse, mpobsCoarse)
  dCoarse = Obs2State_BaseName(sgCoarse, yCoarse)
  DO iv = 1, SIZE(dCoarse%Fields)
    WRITE (*, 33) iv, sgCoarse%mpddInfo_sg%myrank, sgCoarse%num_cell, MAXVAL(ABS(dCoarse%fields(iv)%DATA(:, :, :)))
33  FORMAT('MaxABS of combined thinned obs of var: ', I1, ' at proc: ', I1, I8, ' with value: ', D12.4)
    CALL Output_NC_State_SV(dCoarse, ncOutputFile, &
                            "testMO_obsThinnedMap_"//TRIM(dCoarse%Fields(iv)%Get_Name()), TRIM(dCoarse%Fields(iv)%Get_Name()), .TRUE.)
  END DO

  CALL mpddGlob%barrier

  ! Test passed if come to this point:
  PRINT *, 'Test passed'

  DEALLOCATE (NRanks, mpObs)

  ! Finalize
  CALL mpddGlob%finalize

END PROGRAM Test_griddedObsMap
