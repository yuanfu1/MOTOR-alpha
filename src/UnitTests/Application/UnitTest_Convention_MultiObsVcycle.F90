!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Jiongming Pang, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           : 2022-01-25, created by Jiongming Pang for mulitsource DA unit test.
!
!   Created by Jiongming Pang (pang.j.m@hotmail.com), 2022/01/19, @GBA-MWF, Shenzhen
!
!!--------------------------------------------------------------------------------------------------
!!
!> @brief
!! This is a unit test of MOTOR-DA analysis by using conventional observations.
!!   Uses a module of a multiscale analytic function to test the multiscale capabity of MOTOR-DA
!! This is a unit test using real conventional observation stations replacing their observation values with
!!   the analytic test function values.
!
PROGRAM UnitTest_Convention_MultiObsVcycle
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

  USE ObsSurface_m, ONLY: ObsSurface_t
  USE ObsSound_m, ONLY: ObsSound_t
  USE ObsVwpw_m, ONLY: ObsVwpw_t

  USE ObsUtilities_m

  USE unitTest_sfc_m, ONLY: unitTest_sfc_t
  USE Obs2State_m
  USE Obs2State_m

  IMPLICIT NONE

  ! Define types
  TYPE(mpddGlob_t), TARGET :: mpddGlob
  TYPE(geometry_t), TARGET :: geometry
  TYPE(State_t) :: XRes, XTru
  TYPE(State_t), ALLOCATABLE :: XbMG(:)
  REAL(r_kind) :: t1, t2
  INTEGER(i_kind) :: mgStart, mgEnd

  CLASS(MiniSolver_t), ALLOCATABLE :: miniSolver

  CHARACTER(LEN=1024) :: configFile, ncOutputFile !< Config and output file name.

  INTEGER(i_kind) :: i, iii, imx, iv, ii, kk
  REAL(r_kind) :: amx
  TYPE(unitTest_sfc_t) :: unitTest

  TYPE(ObsSurface_t) :: OBS1 !OBS_surface
  TYPE(ObsSound_t)   :: OBS2 !_sound
  TYPE(ObsVwpw_t)    :: OBS3 !_vwpw

  LOGICAL :: IncFlag = .TRUE.
  REAL(r_kind), ALLOCATABLE :: xyt(:, :)
  REAL(r_kind) :: passthreshold_max, passthreshold_min
  REAL(r_kind) :: ana_max, ana_min, ana_max_global, ana_min_global
  REAL(r_kind), ALLOCATABLE :: ana_max_AV(:), ana_min_AV(:)

  CHARACTER(LEN=20), ALLOCATABLE :: varList(:)
  CHARACTER(LEN=20) :: varName
  INTEGER(i_kind)   :: varNum
  INTEGER(i_kind)   :: ifile

  LOGICAL :: ObsFlag_Surface = .TRUE.
  LOGICAL :: ObsFlag_Sound = .TRUE.
  LOGICAL :: ObsFlag_Vwpw = .TRUE.

  ! Get the configFile
  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", configFile)
  configFile = TRIM(configFile)//"/UnitTest/testConvention_MultiObsVCycle.yaml"

  IF (yaml_get_var(configFile, 'IO', 'output_dir', ncOutputFile) /= 0) STOP
  ncOutputFile = TRIM(ncOutputFile)

  !CALL namelist_read(TRIM(configFile), "varList", varList)
  ifile = yaml_get_var(TRIM(configFile), 'modelState', 'varList', varList)
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

  ! CALL namelist_read(TRIM(configFile), "mgStart", mgStart)
  ! CALL namelist_read(TRIM(configFile), "mgEnd", mgEnd)
  !ifile =  yaml_get_var(TRIM(configFile), 'geometry', 'mgStart', mgStart)
  !ifile =  yaml_get_var(TRIM(configFile), 'geometry', 'mgEnd', mgEnd)
  mgStart = 3
  mgEnd = 6

  DO i = mgEnd, mgStart, -1
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
    DO ii = 1, varNum
      XbMG(i)%fields(ii)%DATA = 0.8D0
    END DO
  END DO

  ! Give the initial value at the coarsest grid.
  XRes = XbMG(mgStart)

  passthreshold_min = 0.6
  passthreshold_max = 6.0

  IF (ObsFlag_Surface) THEN

    CALL OBS1%ObsInitial(configFile)
    CALL OBS1%ObsIngest(XbMG(mgEnd))

    ! Use of unit test analytic observation data:
    ASSOCIATE (sg => geometry%mg%sg(mgEnd))
      ALLOCATE (xyt(3, OBS1%numObs))

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

        WRITE (*, 51) minGlobLat / degree2radian, maxGlobLat / degree2radian, &
          minGlobLon / degree2radian, maxGlobLon / degree2radian
51      FORMAT('ranges of Global LL: ', 4D12.4)

        xyt(1, :) = (OBS1%olatlon(1, 1:OBS1%numObs) - minGlobLat) / &
                    (maxGlobLat - minGlobLat)
        xyt(2, :) = (OBS1%olatlon(2, 1:OBS1%numObs) - minGlobLon) / &
                    (maxGlobLon - minGlobLon)
        xyt(3, :) = (OBS1%obsTime(1:OBS1%numObs) - sg%tt(1)) / &
                    (sg%tt(sg%tSlots) - sg%tt(1))
        DO ii = 1, UBOUND(OBS1%obsData, 2)
          CALL unitTest%analytic(OBS1%numObs, xyt, OBS1%obsData(:, ii))
        END DO

      END BLOCK

      PRINT *, 'MAX/MIN analytic obs_Surface values: ', sg%mpddInfo_sg%myRank, MAXVAL(OBS1%obsData), MINVAL(OBS1%obsData)

      DEALLOCATE (xyt)

    END ASSOCIATE

  END IF

  IF (ObsFlag_Sound) THEN

    CALL OBS2%ObsInitial(configFile)
    CALL OBS2%ObsIngest(XbMG(mgEnd))

    ! Use of unit test analytic observation data:
    ASSOCIATE (sg => geometry%mg%sg(mgEnd))
      ALLOCATE (xyt(3, OBS2%numObs))

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
52      FORMAT('ranges of Global LL: ', 4D12.4)

        xyt(1, :) = (OBS2%olatlon(1, 1:OBS2%numObs) - minGlobLat) / &
                    (maxGlobLat - minGlobLat)
        xyt(2, :) = (OBS2%olatlon(2, 1:OBS2%numObs) - minGlobLon) / &
                    (maxGlobLon - minGlobLon)
        xyt(3, :) = (OBS2%obsTime(1:OBS2%numObs) - MINVAL(OBS2%obsTime(1:OBS2%numObs))) / &
                    (MAXVAL(OBS2%obsTime(1:OBS2%numObs)) - MINVAL(OBS2%obsTime(1:OBS2%numObs)))
        DO ii = 1, UBOUND(OBS2%obsData, 2)
          CALL unitTest%analytic(OBS2%numObs, xyt, OBS2%obsData(:, ii))
        END DO

      END BLOCK

      PRINT *, 'MAX/MIN analytic obs_Sound values: ', sg%mpddInfo_sg%myRank, MAXVAL(OBS2%obsData), MINVAL(OBS2%obsData)

      DEALLOCATE (xyt)
    END ASSOCIATE
  END IF

  IF (ObsFlag_Vwpw) THEN

    CALL OBS3%ObsInitial(configFile)
    CALL OBS3%ObsIngest(XbMG(mgEnd))

    ! Use of unit test analytic observation data:
    ASSOCIATE (sg => geometry%mg%sg(mgEnd))
      ALLOCATE (xyt(3, OBS3%numObs))

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

        WRITE (*, 53) minGlobLat / degree2radian, maxGlobLat / degree2radian, &
          minGlobLon / degree2radian, maxGlobLon / degree2radian
53      FORMAT('ranges of Global LL: ', 4D12.4)

        xyt(1, :) = (OBS3%olatlon(1, 1:OBS3%numObs) - minGlobLat) / &
                    (maxGlobLat - minGlobLat)
        xyt(2, :) = (OBS3%olatlon(2, 1:OBS3%numObs) - minGlobLon) / &
                    (maxGlobLon - minGlobLon)
        xyt(3, :) = (OBS3%obsTime(1:OBS3%numObs) - MINVAL(OBS3%obsTime(1:OBS3%numObs))) / &
                    (MAXVAL(OBS3%obsTime(1:OBS3%numObs)) - MINVAL(OBS3%obsTime(1:OBS3%numObs)))
        DO ii = 1, UBOUND(OBS3%obsData, 2)
          CALL unitTest%analytic(OBS3%numObs, xyt, OBS3%obsData(:, ii))
        END DO

      END BLOCK

      PRINT *, 'MAX/MIN analytic obs_Vwpw values: ', sg%mpddInfo_sg%myRank, MAXVAL(OBS3%obsData), MINVAL(OBS3%obsData)

      DEALLOCATE (xyt)

    END ASSOCIATE

  END IF

  XRes = XbMG(mgEnd)
  DO kk = 1, 3

    ! MultiGrid
    DO i = mgEnd, mgStart, -1
      ! Run 3DVAR in each single grid.
      BLOCK
        TYPE(State_t) :: X
        TYPE(BMatrix_t) :: B
        TYPE(RMatrix_t) :: R
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
          CALL B%initialize(configFile, sg)     ! Initialize the B matrix

          ! Thinning observations:
          ! CALL OBS%ObsSuper(XbMG(i), Y, mpObs)
          IF (ObsFlag_Surface) CALL OBS1%ObsSuper(XbMG(i), Y1, mpObs)
          IF (ObsFlag_Sound) CALL OBS2%ObsSuper(XbMG(i), Y2, mpObs)
          IF (ObsFlag_Vwpw) CALL OBS3%ObsSuper(XbMG(i), Y3, mpObs)

          ! Concat thinning obs:
          CALL ObsConcat_s((/Y1, Y2, Y3/), Y)

! #ifdef DEBUG
!         amx = 0.0D0
!         do iii=1,UBOUND(Y%ObsFields(1)%values,1)
!           if (sg%cell_cntr(1,Y%ObsFields(1)%idx(iii)%hIdx)/degree2radian .GT. 24.0D0 .and. &
!               sg%cell_cntr(2,Y%ObsFields(1)%idx(iii)%hIdx)/degree2radian .LT. 112.0D0 .and. &
!               Y%ObsFields(1)%values(iii) .LT. 1.5D0) &
!               write(*,11) sg%mpddInfo_sg%myRank,i,iii,Y%ObsFields(1)%values(iii), &
!                 sg%cell_cntr(1,Y%ObsFields(1)%idx(iii)%hIdx)/degree2radian, &
!                 sg%cell_cntr(2,Y%ObsFields(1)%idx(iii)%hIdx)/degree2radian
! 11            format('SuperObs too small at procs: ',I1,' MG: ',I2,I6,D12.4,' at: ',2D12.4)
!           if (amx .LT. Y%ObsFields(1)%values(iii)) THEN
!             amx = Y%ObsFields(1)%values(iii)
!             imx = iii
!           end if
!         end do
!         WRITE(*,21) sg%mpddInfo_sg%myRank,imx,amx
! 21      format('Thinned obs max at procs: ',I1,I6,D12.4)

!         ! Debug:
!         IF (i .EQ. 117) THEN
!           CALL mpddGlob%barrier
!           stop
!         END IF
! #endif

          CALL H%initialize(configFile, X, Y)

          BLOCK
            TYPE(State_t) :: DD
            CHARACTER(len=10) :: temp

            WRITE (temp, "(I2.2)") i
            DD = Obs2State_BaseName(sg, Y)

            DO iv = 1, SIZE(DD%Fields)
              varName = DD%Fields(iv)%Get_Name()
              CALL Output_NC_State_SV(DD, ncOutputFile, "testMO_obsThinned_G"//TRIM(temp)//"_"//TRIM(varName), TRIM(varName))
            END DO
          END BLOCK

          CALL R%initialize(configFile, Y, X%sg) ! Initialize R

          ! Initialize J Function
          CALL JFunc%initialize(configFile, X, Y, H, B, R, sg)

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

          ! print*,'Max/min analysis values: ',sg%mpddInfo_sg%myrank,maxval(XRes%fields(1)%data),minval(XRes%fields(1)%data)

          ALLOCATE (ana_max_AV(varNum), ana_min_AV(varNum))
          DO iv = 1, varNum
            ana_max_AV(iv) = MAXVAL(XRes%fields(iv)%DATA)
            ana_min_AV(iv) = MINVAL(XRes%fields(iv)%DATA)
          END DO
          ana_max = MAXVAL(ana_max_AV)
          ana_min = MINVAL(ana_min_AV)
          DEALLOCATE (ana_max_AV, ana_min_AV)
          ! PRINT *, 'ISSUE DEBUG --- MAX/MIN of ANA:', sg%mpddInfo_sg%myrank, ana_max, ana_min
          CALL sg%mpddInfo_sg%AllReduceMaxReal(ana_max, ana_max_global)
          CALL sg%mpddInfo_sg%AllReduceMinReal(ana_min, ana_min_global)
          IF (mpddGlob%isBaseProc()) &
            PRINT *, 'max/min of ana at Grid:', ana_max_global, ana_min_global, i

          BLOCK
            CHARACTER(LEN=128) :: igrid

            WRITE (igrid, "(I2.2)") i
            DO iv = 1, varNum
              varName = varList(iv)
              CALL Output_NC_State_SV(XRes, ncOutputFile, "testMO_ana_G"//TRIM(igrid)//"_"//TRIM(varName), TRIM(varName))
            END DO
          END BLOCK

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
          CALL B%initialize(configFile, sg)     ! Initialize the B matrix

          ! Thinning observations:
          ! CALL OBS%ObsSuper(XbMG(i), Y, mpObs)
          IF (ObsFlag_Surface) CALL OBS1%ObsSuper(XbMG(i), Y1, mpObs)
          IF (ObsFlag_Sound) CALL OBS2%ObsSuper(XbMG(i), Y2, mpObs)
          IF (ObsFlag_Vwpw) CALL OBS3%ObsSuper(XbMG(i), Y3, mpObs)

          ! Concat thinning obs:
          CALL ObsConcat_s((/Y1, Y2, Y3/), Y)

! #ifdef DEBUG
!         amx = 0.0D0
!         do iii=1,UBOUND(Y%ObsFields(1)%values,1)
!           if (sg%cell_cntr(1,Y%ObsFields(1)%idx(iii)%hIdx)/degree2radian .GT. 24.0D0 .and. &
!               sg%cell_cntr(2,Y%ObsFields(1)%idx(iii)%hIdx)/degree2radian .LT. 112.0D0 .and. &
!               Y%ObsFields(1)%values(iii) .LT. 1.5D0) &
!               write(*,11) sg%mpddInfo_sg%myRank,i,iii,Y%ObsFields(1)%values(iii), &
!                 sg%cell_cntr(1,Y%ObsFields(1)%idx(iii)%hIdx)/degree2radian, &
!                 sg%cell_cntr(2,Y%ObsFields(1)%idx(iii)%hIdx)/degree2radian
! 11            format('SuperObs too small at procs: ',I1,' MG: ',I2,I6,D12.4,' at: ',2D12.4)
!           if (amx .LT. Y%ObsFields(1)%values(iii)) THEN
!             amx = Y%ObsFields(1)%values(iii)
!             imx = iii
!           end if
!         end do
!         WRITE(*,21) sg%mpddInfo_sg%myRank,imx,amx
! 21      format('Thinned obs max at procs: ',I1,I6,D12.4)

!         ! Debug:
!         IF (i .EQ. 117) THEN
!           CALL mpddGlob%barrier
!           stop
!         END IF
! #endif

          CALL H%initialize(configFile, X, Y)

          BLOCK
            TYPE(State_t) :: DD
            CHARACTER(len=10) :: temp

            WRITE (temp, "(I2.2)") i
            DD = Obs2State_BaseName(sg, Y)

            DO iv = 1, SIZE(DD%Fields)
              varName = DD%Fields(iv)%Get_Name()
              CALL Output_NC_State_SV(DD, ncOutputFile, "testMO_obsThinned_G"//TRIM(temp)//"_"//TRIM(varName), TRIM(varName))
            END DO
          END BLOCK

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

          ! print*,'Max/min analysis values: ',sg%mpddInfo_sg%myrank,maxval(XRes%fields(1)%data),minval(XRes%fields(1)%data)

          ALLOCATE (ana_max_AV(varNum), ana_min_AV(varNum))
          DO iv = 1, varNum
            ana_max_AV(iv) = MAXVAL(XRes%fields(iv)%DATA)
            ana_min_AV(iv) = MINVAL(XRes%fields(iv)%DATA)
          END DO
          ana_max = MAXVAL(ana_max_AV)
          ana_min = MINVAL(ana_min_AV)
          DEALLOCATE (ana_max_AV, ana_min_AV)
          ! PRINT *, 'ISSUE DEBUG --- MAX/MIN of ANA:', sg%mpddInfo_sg%myrank, ana_max, ana_min
          CALL sg%mpddInfo_sg%AllReduceMaxReal(ana_max, ana_max_global)
          CALL sg%mpddInfo_sg%AllReduceMinReal(ana_min, ana_min_global)
          IF (mpddGlob%isBaseProc()) &
            PRINT *, 'max/min of ana at Grid:', ana_max_global, ana_min_global, i

          BLOCK
            CHARACTER(LEN=128) :: igrid

            WRITE (igrid, "(I2.2)") i
            DO iv = 1, varNum
              varName = varList(iv)
              CALL Output_NC_State_SV(XRes, ncOutputFile, "testMO_ana_G"//TRIM(igrid)//"_"//TRIM(varName), TRIM(varName))
            END DO
          END BLOCK

        END ASSOCIATE
      END BLOCK
    END DO

  END DO

  CALL mpddGlob%barrier

! #ifdef DEBUG
!   ! Analytic solution:
!   ASSOCIATE (sg => geometry%mg%sg(mgEnd))
!     CALL XTru%initialize(configFile, sg)
!     CALL unitTest%analytic(sg%num_icell,sg%cell_cntr,Xtru%fields(1)%data)
!     ALLOCATE(xyt(3,sg%num_icell))
!     !CALL unitTest%checking(sg%num_icell,xyt,Xtru%fields(1)%data(1,:,1),XRes%fields(1)%data(1,:,1))
!     DEALLOCATE(xyt)
!   END ASSOCIATE
! #endif

  ! ! Output state file to NC for view.
  ! DO iv =1,varNum
  !   varName = varList(iv)
  !   CALL Output_NC_State_SV(XRes, ncOutputFile, "test_newOBS_"//TRIM(varName), TRIM(varName))
  ! END DO
  CALL CPU_TIME(t2)
  CALL mpddGlob%barrier
  IF (mpddGlob%isBaseProc()) PRINT *, 'Time cost:', t2 - t1, 's'

  ALLOCATE (ana_max_AV(varNum), ana_min_AV(varNum))
  DO iv = 1, varNum
    ana_max_AV(iv) = MAXVAL(XRes%fields(iv)%DATA(:, :, :))
    ana_min_AV(iv) = MINVAL(XRes%fields(iv)%DATA(:, :, :))
  END DO
  ana_max = MAXVAL(ana_max_AV)
  ana_min = MINVAL(ana_min_AV)
  DEALLOCATE (ana_max_AV, ana_min_AV)

  CALL geometry%mg%sg(mgEnd)%mpddInfo_sg%AllReduceMaxReal(ana_max, ana_max_global)
  CALL geometry%mg%sg(mgEnd)%mpddInfo_sg%AllReduceMinReal(ana_min, ana_min_global)
  IF (mpddGlob%isBaseProc()) &
    PRINT *, 'GLOBAL max/min of ana:', ana_max_global, ana_min_global

  IF (ana_min .GT. passthreshold_min .AND. ana_max .LT. passthreshold_max) THEN
    IF (mpddGlob%isBaseProc()) WRITE (*, 2)
2   FORMAT('Test passed')
  ELSE
    IF (mpddGlob%isBaseProc()) WRITE (*, 3)
3   FORMAT('Test failed')
  END IF

  DEALLOCATE (miniSolver)

  ! Finalize
  CALL mpddGlob%finalize

END PROGRAM UnitTest_Convention_MultiObsVcycle
