!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           : 2021-10-27, created by Yuanfu Xie
!                     2022-01-20, unified the modules and types of OBS by Jiongming Pang
!
!   Created by Yuanfu Xie (yuanfu_xie@yahoo.com.com), 2021/10/27, @GBA-MWF, Shenzhen
!   Modified by Jiongming Pang (pang.j.m@hotmail.com), 2022/1/20, @GBA-MWF, Shenzhen
!   Modified by Zilong Qin (zilong.qin@gmail.com), 2022/3/13, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This is a unit test of MOTOR-DA surface analysis
!!   Uses a module of a multiscale surface analytic function to test the multiscale capabity of MOTOR-DA
!! This is an ideal test using randomly generated observation locations. This is suitable for any
!!   standalone test.
!
PROGRAM UniTest_sfc_randoms
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
  USE MGOpts_m
  USE parameters_m, ONLY: degree2radian

  USE ObsSurface_m, ONLY: ObsSurface_t
  USE unitTest_sfc_m, ONLY: unitTest_sfc_t
  USE YAMLRead_m
  USE Obs2State_m

  IMPLICIT NONE

  ! Define types
  TYPE(mpddGlob_t), TARGET :: mpddGlob
  TYPE(geometry_t), TARGET :: geometry
  TYPE(State_t) :: XRes, XTru
  TYPE(State_t), ALLOCATABLE :: XbMG(:)
  REAL(r_kind) :: t1, t2, ran
  INTEGER(i_kind) :: mgStart, mgEnd    ! Multigrid levels

  CLASS(MiniSolver_t), ALLOCATABLE :: miniSolver

  CHARACTER(LEN=1024) :: configFile, ncOutputFile !< Config and output file name.

  TYPE(unitTest_sfc_t) :: unitTest

  LOGICAL :: IncFlag = .TRUE., passed = .FALSE.
  INTEGER(i_kind) :: numVars, i, imx, iii, numRandomPoints, istatus
  REAL(r_kind) :: amx, passthreshold
  REAL(r_kind), ALLOCATABLE :: xyt(:, :)

  TYPE(ObsSurface_t) :: obs_in

  ! Get the configFile
  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", configFile)
  configFile = TRIM(configFile)//"/UnitTest/testSurfacRandom.yaml"

  istatus = yaml_get_var(configFile, 'IO', 'output_dir', ncOutputFile)
  ncOutputFile = TRIM(ncOutputFile)

  ! Initializer - Auxtypes
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

  mgStart = 3      ! start index of multigrid
  mgEnd = 8       ! end index of multigrid

  ! Initialize obs:
  numRandomPoints = 1000
  CALL obs_in%ObsInitial(configFile)
  CALL unitTest%randomObs(obs_in, numRandomPoints, geometry%mg%sg(mgEnd))

  DO i = mgEnd, mgStart, -1
    PRINT *, 'Multigrid level: ', i
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
    XbMG(i)%fields(1)%DATA = 0
  END DO

  ! Give the initial value at the coarsest grid.
  XRes = XbMG(mgStart)

  ! MultiGrid
  DO i = mgStart, mgEnd
    ! Run 3DVAR in each single grid.
    BLOCK
      TYPE(State_t) :: X, Xb
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

        ! Thinning observations:
        CALL obs_in%ObsThinning(XbMG(i), Y, mpObs, .FALSE., .FALSE.)

! #ifdef DEBUG
!           print*,'Max/min super obs values: ',sg%mpddInfo_sg%myRank, &
!             maxval(Y%ObsFields(1)%values(:)),minval(Y%ObsFields(1)%values(:))
!           print*,'Number of thinned obs: ',UBOUND(Y%ObsFields(1)%values,1)
!           amx = 0.0D0
!           do iii=1,UBOUND(Y%ObsFields(1)%values,1)
!             if (sg%cell_cntr(1,Y%ObsFields(1)%idx(iii)%hIdx)/degree2radian .GT. 24.0D0 .and. &
!                 sg%cell_cntr(2,Y%ObsFields(1)%idx(iii)%hIdx)/degree2radian .LT. 112.0D0 .and. &
!                 Y%ObsFields(1)%values(iii) .LT. 1.5D0) &
!                 write(*,11) sg%mpddInfo_sg%myRank,i,iii,Y%ObsFields(1)%values(iii), &
!                   sg%cell_cntr(1,Y%ObsFields(1)%idx(iii)%hIdx)/degree2radian, &
!                   sg%cell_cntr(2,Y%ObsFields(1)%idx(iii)%hIdx)/degree2radian
!   11            format('SuperObs too small at procs: ',I1,' MG: ',I2,I6,D12.4,' at: ',2D12.4)
!             if (amx .LT. Y%ObsFields(1)%values(iii)) THEN
!               amx = Y%ObsFields(1)%values(iii)
!               imx = iii
!             end if
!           end do
!           WRITE(*,21) sg%mpddInfo_sg%myRank,imx,amx
!   21      format('Thinned obs max at procs: ',I1,I6,D12.4)

!           ! Debug:
!           IF (i .EQ. 15) THEN
!             CALL mpddGlob%barrier
!             stop
!           END IF
! #endif

        CALL H%initialize(configFile, X, Y)

        ! CALL Y%ObsFields(1)%Set_Name('temp')

        BLOCK
          CHARACTER(len=10) :: temp
          WRITE (temp, "(I2.2)") i
          CALL Output_NC_State_SV(Obs2State_BaseTypeName(sg, Y), ncOutputFile, "MGTest_Yt_G"//TRIM(temp), "SYNOP_temp")
        END BLOCK

        CALL R%initialize(configFile, Y, X%sg) ! Initialize R

        ! Initialize J Function
        CALL JFunc%initialize(configFile, X, Y, H, B, R, sg, XRes)
        ! CALL JFunc%initialize(configFile, X, Y, H, B, R, sg, XbMG(sg%gLevel))

        ! Run minimization
        CALL miniSolver%run(X, JFunc, sg, 100)

        ! Sync
        CALL mpddGlob%barrier

        ! Connect to coarser grid.
        IF (i .NE. mgEnd) THEN
          ! ===========================> prolongate to finer grid
          ASSOCIATE (sgFiner => geometry%mg%sg(i + 1))
            ! CALL XRes%initialize(configFile, sgFiner)
            XRes = XbMG(i + 1)%zeroCopy()
            CALL prolongationMG(X, XRes, XbMG(i + 1), geometry%mg, IncFlag)
          END ASSOCIATE
        ELSE
          XRes = X ! Return the state fields directly
        END IF

        PRINT *, 'Max/min analysis values: ', sg%mpddInfo_sg%myrank, MAXVAL(XRes%fields(1)%DATA), MINVAL(XRes%fields(1)%DATA)

      END ASSOCIATE
    END BLOCK
  END DO
  PRINT*,'Dim of Xres: ',ubound(XRes%fields,1)

  CALL mpddGlob%barrier

  ! Analytic solution:
  ASSOCIATE (sg => geometry%mg%sg(mgEnd))
    CALL XTru%initialize(configFile, sg)
    ALLOCATE (xyt(3, sg%num_cell))
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
5     FORMAT('ranges of Global LL: ', 4D12.4)

      xyt(1, :) = (sg%cell_cntr(1, 1:sg%num_cell) - minGlobLat) / &
                  (maxGlobLat - minGlobLat)
      xyt(2, :) = (sg%cell_cntr(2, 1:sg%num_cell) - minGlobLon) / &
                  (maxGlobLon - minGlobLon)
    END BLOCK
    xyt(3, :) = 0.0D0
    CALL unitTest%analytic(sg%num_cell, xyt, Xtru%fields(1)%DATA(1, :, 1))
    !CALL unitTest%checking(sg%num_icell,xyt,Xtru%fields(1)%data(1,:,1),XRes%fields(1)%data(1,:,1))
    DEALLOCATE (xyt)
    amx = 0.0D0
    imx = 0
    DO iii = 1, sg%num_cell
      IF (amx .LT. ABS(XRes%fields(1)%DATA(1, iii, 1) - Xtru%fields(1)%DATA(1, iii, 1))) THEN
        amx = ABS(XRes%fields(1)%DATA(1, iii, 1) - Xtru%fields(1)%DATA(1, iii, 1))
        imx = iii
      END IF
    END DO
    WRITE (*, 1) amx, imx, sg%mpddInfo_sg%myRank, XRes%fields(1)%DATA(1, imx, 1), Xtru%fields(1)%DATA(1, imx, 1)
1   FORMAT('Maximum analysis error: ', D12.4, ' at grid: ', I8, ' at procs: ', I1, /, 2D20.10)

    ! Xtru%fields(1)%DATA = Xres%fields(1)%DATA-Xtru%fields(1)%DATA

  END ASSOCIATE

  ! Output state file to NC for view.
  CALL Output_NC_State_SV(XRes, ncOutputFile, 'unitest', "temp")
  XRes%fields(2)%DATA = XRes%fields(1)%DATA - Xtru%fields(1)%DATA
  CALL Output_NC_State_SV(XRes, ncOutputFile, 'truth', "diff")
  CALL Output_NC_State_SV(Xtru, ncOutputFile, 'truth', "temp")
  CALL CPU_TIME(t2)
  IF (mpddGlob%isBaseProc()) PRINT *, 'Time cost:', t2 - t1, 's'

  !CALL mpddGlob%barrier

  ! NOTE: as a careful calibration of the analysis and the truth, Yuanfu Xie temporarily requires
  ! the analysis error is less than half of the maximum value of the truth. on 2025-04-17
  ! Careful calibration requires
  ! 1. carefully adjust the analysis background error covariance matrix
  ! 2. carefully adjust the analytic function's spatial scales, the current analytic function
  !    has too fine spatial scales, which requres extremely large amount of observations to
  !    calibrate the analysis error.
  passThreshold = 0.5D0*MAXVAL(ABS(Xtru%fields(1)%DATA(1,:,1)))
  PRINT *, 'AMX XIE: ', amx, passThreshold
  IF (amx .LT. passThreshold) THEN
    IF (mpddGlob%isBaseProc()) WRITE (*, 2)
2   FORMAT('Test passed')
  ELSE
    IF (mpddGlob%isBaseProc()) WRITE (*, 3)
3   FORMAT('Test failed',E14.6)
  END IF

  DEALLOCATE (miniSolver)

  ! Finalize
  CALL mpddGlob%barrier
  CALL mpddGlob%finalize

END PROGRAM UniTest_sfc_randoms
