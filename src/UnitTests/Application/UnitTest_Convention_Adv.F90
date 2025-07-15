!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Jiongming Pang, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           : 2022-01-19, created by Jiongming Pang.
!                     2022-01-20, unified the modules and types of OBS by Jiongming Pang
!
!   Created by Jiongming Pang (pang.j.m@hotmail.com), 2022/01/19, @GBA-MWF, Shenzhen
!   Modified by Jiongming Pang (pang.j.m@hotmail.com), 2022/01/20, @GBA-MWF, Shenzhen
!   Modified by Zilong Qin (zilong.qin@gmail.com), 2022/3/13, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------
!!
!> @brief
!! This is a unit test of MOTOR-DA analysis by using conventional observations.
!!   Uses a module of a multiscale analytic function to test the multiscale capabity of MOTOR-DA
!! This is a unit test using real conventional observation stations replacing their observation values with
!!   the analytic test function values.
!
PROGRAM unitTest_Convention_Adv
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
  USE RMatrix_m, ONLY: RMatrix_t
  USE JFunc_m, ONLY: JFunc_t
  USE State2NC_m
  USE Mock_m
  USE MGOpts_m
  USE parameters_m, ONLY: degree2radian
  USE YAMLRead_m

  ! Unified the modules and types of OBS. -Jiongming Pang 2022-01-20
  USE ObsSurface_m, ONLY: ObsSurface_t
  USE ObsSound_m, ONLY: ObsSound_t

  USE unitTest_sfc_m, ONLY: unitTest_sfc_t
  USE Obs2State_m
  USE Obs2State_m

  IMPLICIT NONE

  ! Define types
  TYPE(mpddGlob_t), TARGET :: mpddGlob
  TYPE(geometry_t), TARGET :: geometry
  TYPE(MPObs_t), TARGET :: mpObs
  TYPE(State_t) :: XRes, XTru
  TYPE(State_t), ALLOCATABLE :: XbMG(:)
  TYPE(ObsSet_t) :: Y
  TYPE(C2O_t) :: H
  !TYPE(BMatrix_t) :: B
  TYPE(JFunc_t) :: JFunc
  REAL(r_kind) :: t1, t2
  INTEGER(i_kind) :: mgStart, mgEnd, kk

  CLASS(MiniSolver_t), ALLOCATABLE :: miniSolver

  CHARACTER(LEN=1024) :: configFile, ncOutputFile !< Config and output file name.

  ! Yuanfu Xie test:
  CHARACTER(LEN=1024) :: filename
  CHARACTER(LEN=20), ALLOCATABLE :: varname(:)
  INTEGER(i_kind) :: numVars, i, iii, imx, ifile
  REAL(r_kind) :: amx
  TYPE(unitTest_sfc_t) :: unitTest

  ! TYPE(ObsSurface_t) :: OBS
  TYPE(ObsSound_t) :: OBS

  LOGICAL :: IncFlag = .TRUE.
  REAL(r_kind), ALLOCATABLE :: xyt(:, :)

  ! Get the configFile
  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", configFile)
  configFile = TRIM(configFile)//"/UnitTest/testMiniSolver_newObs.yaml"

  IF (yaml_get_var(configFile, 'IO', 'output_dir', ncOutputFile) /= 0) STOP
  ncOutputFile = TRIM(ncOutputFile)

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

  !CALL namelist_read(TRIM(configFile), "mgStart", mgStart)
  !CALL namelist_read(TRIM(configFile), "mgEnd", mgEnd)
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
    XbMG(i)%fields(1)%DATA = 0.8D0
  END DO

  ! ! Give the initial value at the coarsest grid.
  ! XRes = XbMG(mgStart)

  ! Give the initial value at the finest grid.
  XRes = XbMG(mgEnd)

  CALL OBS%ObsInitial(configFile)
  CALL OBS%ObsIngest(XbMG(mgEnd))

  DO kk = 1, 10
    ! MultiGrid
    DO i = mgEnd, mgStart, -1
      ! Run 3DVAR in each single grid.
      BLOCK
        TYPE(State_t) :: X
        TYPE(BMatrix_t) :: B
        TYPE(RMatrix_t) :: R
        TYPE(MPObs_t), TARGET :: mpObs
        TYPE(ObsSet_t) :: Y
        TYPE(C2O_t) :: H
        TYPE(JFunc_t) :: JFunc
        ! Yuanfu Xie test:
        CHARACTER(LEN=1024) :: filename(1)
        CHARACTER(LEN=20), ALLOCATABLE :: varname(:)
        INTEGER(i_kind) :: numVars

        ASSOCIATE (sg => geometry%mg%sg(i))
          CALL mpObs%initializeMPObs(sg)    ! Initialize the observation parallel processing proc

          ! Initialize X, Y, B, H
          X = Xres     ! Initialize the state
          CALL B%initialize(configFile, sg)     ! Initialize the B matrix
          Y = ObsSet_t(configFile, mpObs)

          ! Mock observation data:
          filename = '../data/20210906_1210'
          numVars = 1
          ALLOCATE (varname(numVars))
          varname(1) = 'temperature'

          ! Use of unit test analytic observation data:
          ALLOCATE (xyt(3, OBS%numObs))

          BLOCK
            REAL(r_kind) :: minGlobLat, maxGlobLat, minGlobLon, maxGlobLon, swap

            swap = MINVAL(sg%cell_cntr(1, 1:sg%num_icell))
            CALL sg%mpddInfo_sg%AllReduceMinReal(swap, minGlobLat)

            swap = MAXVAL(sg%cell_cntr(1, :))
            CALL sg%mpddInfo_sg%AllReduceMaxReal(swap, maxGlobLat)

            swap = MINVAL(sg%cell_cntr(2, :))
            CALL sg%mpddInfo_sg%AllReduceMinReal(swap, minGlobLon)

            swap = MAXVAL(sg%cell_cntr(2, :))
            CALL sg%mpddInfo_sg%AllReduceMaxReal(swap, maxGlobLon)

            WRITE (*, 52) minGlobLat / degree2radian, maxGlobLat / degree2radian, &
              minGlobLon / degree2radian, maxGlobLon / degree2radian
52          FORMAT('ranges of Global LL: ', 4D12.4)

            xyt(1, :) = (OBS%olatlon(1, 1:OBS%numObs) * degree2radian - minGlobLat) / &
                        (maxGlobLat - minGlobLat)
            xyt(2, :) = (OBS%olatlon(2, 1:OBS%numObs) * degree2radian - minGlobLon) / &
                        (maxGlobLon - minGlobLon)
            xyt(3, :) = (OBS%obsTime(1:OBS%numObs) - sg%tt(1)) / &
                        (sg%tt(sg%tSlots) - sg%tt(1))
            CALL unitTest%analytic(OBS%numObs, xyt, OBS%obsData)
          END BLOCK
          PRINT *, 'Max/min analytic obs values: ', sg%mpddInfo_sg%myRank, MAXVAL(OBS%obsData(:, :)), MINVAL(OBS%obsData(:, :))

          ! Thinning observations:
          CALL OBS%ObsSuper(XbMG(i), Y, mpObs)

          PRINT *, 'Max/min super obs values: ', sg%mpddInfo_sg%myRank, MAXVAL(Y%ObsFields(1)%values(:)), MINVAL(Y%ObsFields(1)%values(:))

          PRINT *, 'Number of thinned obs: ', UBOUND(Y%ObsFields(1)%values, 1)
          amx = 0.0D0
          DO iii = 1, UBOUND(Y%ObsFields(1)%values, 1)
            IF (sg%cell_cntr(1, Y%ObsFields(1)%idx(iii)%hIdx) / degree2radian .GT. 24.0D0 .AND. &
                sg%cell_cntr(2, Y%ObsFields(1)%idx(iii)%hIdx) / degree2radian .LT. 112.0D0 .AND. &
                Y%ObsFields(1)%values(iii) .LT. 1.5D0) &
              WRITE (*, 112) sg%mpddInfo_sg%myRank, i, iii, Y%ObsFields(1)%values(iii), &
              sg%cell_cntr(1, Y%ObsFields(1)%idx(iii)%hIdx) / degree2radian, &
              sg%cell_cntr(2, Y%ObsFields(1)%idx(iii)%hIdx) / degree2radian
112         FORMAT('SuperObs too small at procs: ', I1, ' MG: ', I2, I6, D12.4, ' at: ', 2D12.4)
            IF (amx .LT. Y%ObsFields(1)%values(iii)) THEN
              amx = Y%ObsFields(1)%values(iii)
              imx = iii
            END IF
          END DO
          WRITE (*, 212) sg%mpddInfo_sg%myRank, imx, amx
212       FORMAT('Thinned obs max at procs: ', I1, I6, D12.4)

          ! Debug:
          IF (i .EQ. 117) THEN
            CALL mpddGlob%barrier
            STOP
          END IF

          CALL H%initialize(configFile, X, Y)

          BLOCK
            CHARACTER(len=10) :: temp

            WRITE (temp, "(I2.2)") i
            CALL Output_NC_State_SV(Obs2State_BaseTypeName(sg, Y), ncOutputFile, "MGTest_Yt_newOBSAdv_G"//TRIM(temp), "SYNOP_temp")
          END BLOCK

          CALL R%initialize(configFile, Y, X%sg) ! Initialize R

          ! Initialize J Function
          CALL JFunc%initialize(configFile, X, Y, H, B, R, sg)                              !

          ! Run minimization
          CALL miniSolver%run(X, JFunc, sg)

          ! Sync
          CALL mpddGlob%barrier

          ! Connect to coarser grid.
          IF (i .NE. mgStart) THEN
            ! ===========================> prolongate to finer grid
            ASSOCIATE (sgCoarser => geometry%mg%sg(i - 1))
              CALL XRes%initialize(configFile, sgCoarser)
              CALL restrictionMG(XRes, X, geometry%mg)
            END ASSOCIATE
          ELSE
            XRes = X ! Return the state fields directly
          END IF

          PRINT *, 'Max/min analysis values: ', sg%mpddInfo_sg%myrank, MAXVAL(XRes%fields(1)%DATA), MINVAL(XRes%fields(1)%DATA)

          ! Deallocate xyt:
          DEALLOCATE (xyt)

        END ASSOCIATE
      END BLOCK
    END DO

    ! MultiGrid
    DO i = mgStart, mgEnd
      ! Run 3DVAR in each single grid.
      BLOCK
        TYPE(State_t) :: X
        TYPE(BMatrix_t) :: B
        TYPE(RMatrix_t) :: R
        TYPE(MPObs_t), TARGET :: mpObs
        TYPE(ObsSet_t) :: Y
        TYPE(C2O_t) :: H
        TYPE(JFunc_t) :: JFunc
        ! Yuanfu Xie test:
        CHARACTER(LEN=1024) :: filename(1)
        CHARACTER(LEN=20), ALLOCATABLE :: varname(:)
        INTEGER(i_kind) :: numVars

        ASSOCIATE (sg => geometry%mg%sg(i))
          CALL mpObs%initializeMPObs(sg)    ! Initialize the observation parallel processing proc

          ! Initialize X, Y, B, H
          X = Xres     ! Initialize the state
          CALL B%initialize(configFile, sg)     ! Initialize the B matrix
          Y = ObsSet_t(configFile, mpObs)

          ! Mock observation data:
          filename = '../data/20210906_1210'
          numVars = 1
          ALLOCATE (varname(numVars))
          varname(1) = 'temperature'

          ! Use of unit test analytic observation data:
          ALLOCATE (xyt(3, OBS%numObs))

          BLOCK
            REAL(r_kind) :: minGlobLat, maxGlobLat, minGlobLon, maxGlobLon, swap

            swap = MINVAL(sg%cell_cntr(1, 1:sg%num_icell))
            CALL sg%mpddInfo_sg%AllReduceMinReal(swap, minGlobLat)

            swap = MAXVAL(sg%cell_cntr(1, :))
            CALL sg%mpddInfo_sg%AllReduceMaxReal(swap, maxGlobLat)

            swap = MINVAL(sg%cell_cntr(2, :))
            CALL sg%mpddInfo_sg%AllReduceMinReal(swap, minGlobLon)

            swap = MAXVAL(sg%cell_cntr(2, :))
            CALL sg%mpddInfo_sg%AllReduceMaxReal(swap, maxGlobLon)

            WRITE (*, 5) minGlobLat / degree2radian, maxGlobLat / degree2radian, &
              minGlobLon / degree2radian, maxGlobLon / degree2radian
5           FORMAT('ranges of Global LL: ', 4D12.4)

            xyt(1, :) = (OBS%olatlon(1, 1:OBS%numObs) * degree2radian - minGlobLat) / &
                        (maxGlobLat - minGlobLat)
            xyt(2, :) = (OBS%olatlon(2, 1:OBS%numObs) * degree2radian - minGlobLon) / &
                        (maxGlobLon - minGlobLon)
            xyt(3, :) = (OBS%obsTime(1:OBS%numObs) - MINVAL(OBS%obsTime(1:OBS%numObs))) / &
                        (MAXVAL(OBS%obsTime(1:OBS%numObs)) - MINVAL(OBS%obsTime(1:OBS%numObs)))
            CALL unitTest%analytic(OBS%numObs, xyt, OBS%obsData)
          END BLOCK
          PRINT *, 'Max/min analytic obs values: ', sg%mpddInfo_sg%myRank, MAXVAL(OBS%obsData(:, :)), MINVAL(OBS%obsData(:, :))

          ! Thinning observations:
          CALL OBS%ObsSuper(XbMG(i), Y, mpObs)

          PRINT *, 'Max/min super obs values: ', sg%mpddInfo_sg%myRank, MAXVAL(Y%ObsFields(1)%values(:)), MINVAL(Y%ObsFields(1)%values(:))

          PRINT *, 'Number of thinned obs: ', UBOUND(Y%ObsFields(1)%values, 1)
          amx = 0.0D0
          DO iii = 1, UBOUND(Y%ObsFields(1)%values, 1)
            IF (sg%cell_cntr(1, Y%ObsFields(1)%idx(iii)%hIdx) / degree2radian .GT. 24.0D0 .AND. &
                sg%cell_cntr(2, Y%ObsFields(1)%idx(iii)%hIdx) / degree2radian .LT. 112.0D0 .AND. &
                Y%ObsFields(1)%values(iii) .LT. 1.5D0) &
              WRITE (*, 11) sg%mpddInfo_sg%myRank, i, iii, Y%ObsFields(1)%values(iii), &
              sg%cell_cntr(1, Y%ObsFields(1)%idx(iii)%hIdx) / degree2radian, &
              sg%cell_cntr(2, Y%ObsFields(1)%idx(iii)%hIdx) / degree2radian
11          FORMAT('SuperObs too small at procs: ', I1, ' MG: ', I2, I6, D12.4, ' at: ', 2D12.4)
            IF (amx .LT. Y%ObsFields(1)%values(iii)) THEN
              amx = Y%ObsFields(1)%values(iii)
              imx = iii
            END IF
          END DO
          WRITE (*, 21) sg%mpddInfo_sg%myRank, imx, amx
21        FORMAT('Thinned obs max at procs: ', I1, I6, D12.4)

          ! Debug:
          IF (i .EQ. 117) THEN
            CALL mpddGlob%barrier
            STOP
          END IF

          CALL H%initialize(configFile, X, Y)

          ! PRINT *, 'temp: ', Y%ObsFields(1)%values(1:5)
          ! CALL Y%ObsFields(1)%Set_Name('temp')
          ! Output Observations
          BLOCK
            CHARACTER(len=10) :: temp
            WRITE (temp, "(I2.2)") i
            CALL Output_NC_State_SV(Obs2State_BaseTypeName(sg, Y), ncOutputFile, "MGTest_Yt_newOBSAdv_G"//TRIM(temp), "SYNOP_temp")
          END BLOCK

          CALL R%initialize(configFile, Y, X%sg) ! Initialize R

          ! Initialize J Function
          CALL JFunc%initialize(configFile, X, Y, H, B, R, sg)                           ! Initialize H
          ! goto 1

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

          PRINT *, 'Max/min analysis values: ', sg%mpddInfo_sg%myrank, MAXVAL(XRes%fields(1)%DATA), MINVAL(XRes%fields(1)%DATA)

          ! Deallocate xyt:
          DEALLOCATE (xyt)

        END ASSOCIATE
      END BLOCK
    END DO

  END DO ! end kk loop

  CALL mpddGlob%barrier

  ! Analytic solution:
  ASSOCIATE (sg => geometry%mg%sg(mgEnd))
    CALL XTru%initialize(configFile, sg)
    CALL unitTest%analytic(sg%num_icell, sg%cell_cntr, Xtru%fields(1)%DATA)
    ALLOCATE (xyt(3, sg%num_icell))
    !CALL unitTest%checking(sg%num_icell,xyt,Xtru%fields(1)%data(1,:,1),XRes%fields(1)%data(1,:,1))
    DEALLOCATE (xyt)
  END ASSOCIATE

  ! Output state file to NC for view.
  CALL Output_NC_State_SV(XRes, ncOutputFile, 'test_newOBSAdv', "temp")
  CALL CPU_TIME(t2)
  IF (mpddGlob%isBaseProc()) PRINT *, 'Time cost:', t2 - t1, 's'

  DEALLOCATE (miniSolver)

  ! Finalize
  CALL mpddGlob%finalize

END PROGRAM unitTest_Convention_Adv
