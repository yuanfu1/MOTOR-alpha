!!--------------------------------------------------------------------------------------------------
! PROJECT           : Test of Geostrophic balance
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Yuanfu Xie and Jia Wang
! VERSION           : V 0.0
! HISTORY           : 2022-09-14, modified by Yuanfu Xie based on Wang's Test_GeosBalTLAD.F90.
!
!   Created by Yuanfu Xie (yuanfu_xie@yahoo.com), 2022/09/14, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------
!!
!> @brief
!! This is a test of MOTOR-DA with a Geostrophic balance
!
PROGRAM Test_GeosBalTLAD
  USE Applications_m
  USE GeosBal_m, ONLY: GeosBal_t
  USE geometry_m, ONLY: geometry_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE MPObs_m, ONLY: MPObs_t
  USE State_m, ONLY: State_t
  ! USE ObsSet_m, ONLY: ObsSet_t
  USE IOGrapes_m, ONLY: IOGrapes_t
  USE State2NC_m
  USE GeosBal_m, ONLY: GeosBal_t
  USE YAMLRead_m

!#define TOPO

  IMPLICIT NONE

  ! Define local variables:
  CHARACTER(LEN=1024) :: configFile
  INTEGER(i_kind) :: istatus
  LOGICAL :: unitTest

  TYPE(mpddGlob_t) :: mpddGlob
  REAL(r_kind) :: rs1, rs2, f1, f2
  TYPE(IOGrapes_t), TARGET :: ioGrapes

  TYPE(Applications_t) :: app

  CHARACTER(LEN=1024) :: staticDir
  CHARACTER(LEN=1024) :: outputDir
  CHARACTER(LEN=20) :: task

  REAL(r_kind) :: t1, t2

  CALL GETARG(1, configFile)
  IF (TRIM(configFile) .EQ. '') THEN
    WRITE (*, 1)
1   FORMAT('Usage of this driver: mpirun -n <n> Debug/Test_NewGeosBal.exe configFile', /, &
           ' Check your configure file and rerun!')
    STOP
  ELSE
    WRITE (*, 2) TRIM(configFile)
2   FORMAT('ConfigFile is: ', A)
  END IF

  unitTest = .FALSE.

  ! Initialize app:
  CALL app%initial(configFile, 2, istatus)
  CALL app%Ctl2State%initialize(configFile)

  ! Read in the background for testing the terrain-following coordinate: Unit Test turned this off
  ! CALL app%backgrd(app%mgStart, app%mgEnd, unitTest)
  ! WRITE(*,3)
3 FORMAT('Background fields have been processed successively')

  BLOCK
    !TYPE(geometry_t) :: geometry
    ! Initialize geometry
    !geometry = geometry_t(configFile, mpddGlob)               ! Initialize the geometry

    ASSOCIATE (sg => app%geometry%mg%sg(3))
      BLOCK

        !TYPE(MPObs_t), TARGET :: mpObs
        TYPE(State_t) :: X, Y, tempX, tempY
        TYPE(GeosBal_t) :: H

        INTEGER(i_kind) :: i, j, k, ic, iv, imx
        REAL(r_kind) :: amx

        ! Insert topography:
        CALL GetTopo(sg)
        CALL sg%update_zHght_from_topo_and_sigma

        !mpObs = MPObs_t(sg)
        CALL X%initialize(configFile, sg)
        CALL tempX%initialize(configFile, sg)

        ! Get a test state:
        PRINT *, 'Getting a test State...', X%mpddSub%myrank
        CALL testState(X)
        CALL app%Ctl2State%transBackward(X)

        CALL CPU_TIME(t1)
        !H = GeosBal_t(configFile, X)
        CALL H%initial(configFile, X)
        CALL CPU_TIME(t2)
        WRITE (*, 42) t2 - t1, app%mpddGlob%myrank
42      FORMAT('Time spent in GeosBal construction:', D12.4, ' pc', I2)

        CALL CPU_TIME(t1)

        ! Y is an arbitrary vector to check adjoint operation: <Y*X, Z> = <X, Y^T Z>
        CALL Y%initialize(configFile, sg)
        CALL tempY%initialize(configFile, sg)
        CALL Y%setAllFieldData(0.0D0)
        CALL tempY%setAllFieldData(0.0D0)
        IF (X%mpddGlob%myrank .EQ. 1) THEN
          CALL Y%setAllFieldData(0.0D0)
          CALL tempY%setAllFieldData(0.0D0)
        END IF

        ic = 7; iv = 2 !X%sg%vLevel

        IF (X%mpddGlob%myrank .EQ. 11) THEN
          WRITE (*, 55) X%sg%num_icell
55        FORMAT('Which grid cell to check: enter an integer between: 1', I3)
          READ (*, *) ic
          WRITE (*, 56) ic
56        FORMAT('Checking ', I2, ' cell')
        END IF
        CALL X%mpddGlob%barrier()

        ! Single gridpoint check:
        IF (X%mpddGlob%myrank .EQ. 1) THEN
          Y%fields(H%iu)%DATA(:, :, :) = 1.0D0
          tempY%fields(H%iu)%DATA(:, :, :) = 1.0D0
          Y%fields(H%iv)%DATA(:, :, :) = 1.0D0
          tempY%fields(H%iv)%DATA(:, :, :) = 1.0D0
          Y%fields(H%it)%DATA(:, :, :) = 1.0D0
          tempY%fields(H%it)%DATA(:, :, :) = 1.0D0
          Y%fields(Y%getVarIdx('qvapor'))%DATA(:, :, :) = 1.0D0
          tempY%fields(Y%getVarIdx('qvapor'))%DATA(:, :, :) = 1.0D0
        END IF
        ! Full grid check:
        ! CALL testState(Y)
        ! CALL testState(tempY)
        CALL app%Ctl2State%transBackward(Y)
        CALL app%Ctl2State%transBackward(tempY)
        BLOCK
          PRINT *, "Calling Forward of Tangent Linear operator: ", X%mpddGlob%myrank
          tempX = X             ! X: pres/rho

          CALL H%fwdTL(tempX, X)  ! tempX: temp/rho/u/v

          ! Debugging check the Tangent linear result: designed for G3 testing
          WRITE (*, 111) H%iv, &
            ! maxval(tempX%fields(H%iv)%data),&
            ! minval(tempX%fields(H%iv)%data),&
            tempX%fields(H%iv)%DATA(iv, ic, 2), &
            tempX%fields(H%iv)%DATA(iv, ic + 5, 2), &
            tempX%fields(H%iv)%DATA(iv, ic + 10, 2), &
            tempX%fields(H%iv)%DATA(iv, ic + 15, 2), &
            !X%fields(X%getVarIdx('temp'))%data(2,ic,2), &
            !X%fields(X%getVarIdx('temp'))%data(2,12,2), &
            X%sg%cell_cntr(2, ic), X%sg%cell_cntr(2, ic + 5), &
            X%mpddGlob%myrank
111       FORMAT('Driver: ', I2, ' Hx', 4D12.4, ' lats: ', 2D14.6, ' pc', I2)

          ! Test topo plot using temp field:
#ifdef TOPO
! This section is used to plot topography and their derivatives only;
! Undefined TOPO at the beginning of the driver if a real test is run!
          DO i = 1, X%sg%vLevel
            tempX%fields(H%it)%DATA(i, :, 1) = sg%sigmax(:)
            tempX%fields(H%iu)%DATA(i, :, 1) = sg%sigmay(:)
            tempX%fields(H%iv)%DATA(i, :, 1) = sg%topo(:)
          END DO
          WRITE (*, 222) MINVAL(sg%sigmax), MAXVAL(sg%sigmax), &
            MINVAL(sg%sigmay), MAXVAL(sg%sigmay), sg%mpddGlob%myrank
222       FORMAT('sigmax: ', 2D12.4, ' sigmay: ', 2D12.4, ' pc', I2)
#endif

          CALL Output_NC_State_AV(tempX, TRIM(app%ncOutputFile), &
                                  TRIM(app%task)//"_HX", .TRUE., .TRUE.)

          ! (Y)^T (H X)
          rs1 = Y.DOT.tempX

          ! The adjoint operator: tempX contains pres/rho
          CALL H%adjMul(tempY, X)

          ! Debugging: 2022-09-11:
          IF (X%mpddGlob%myrank .EQ. 1) THEN
            imx = 0
            amx = 0.0D0
            DO i = 1, UBOUND(tempY%fields(Y%getVarIdx(TRIM('uwnd')))%DATA, 2)
              IF (ABS(tempX%fields(Y%getVarIdx(TRIM('uwnd')))%DATA(iv, i, 1)) .GT. amx) THEN
                imx = i
                amx = ABS(tempX%fields(Y%getVarIdx(TRIM('uwnd')))%DATA(iv, i, 1))
              END IF
            END DO
            WRITE (*, 11) amx, imx, tempX%fields(X%getVarIdx(TRIM('pres_ctl')))%DATA(iv, ic, 1)
11          FORMAT('Maximum abs value of HX:', D12.4, ' at cell:', I6, ' pres_ctl:', D12.4)
            DO i = iv, iv !iv-1,iv+1
              WRITE (*, 5) (tempY%fields(H%it)%DATA(i, j, 1), j=21, 25), i, sg%mpddInfo_sg%myrank
              WRITE (*, 5) (tempY%fields(H%it)%DATA(i, j, 1), j=16, 20), i, sg%mpddInfo_sg%myrank
              WRITE (*, 5) (tempY%fields(H%it)%DATA(i, j, 1), j=11, 15), i, sg%mpddInfo_sg%myrank
              WRITE (*, 5) (tempY%fields(H%it)%DATA(i, j, 1), j=6, 10), i, sg%mpddInfo_sg%myrank
              WRITE (*, 5) (tempY%fields(H%it)%DATA(i, j, 1), j=1, 5), i, sg%mpddInfo_sg%myrank
5             FORMAT('Adjoint temp: ', 5D12.4, ' vertical', I2, ' pc', I2)
            END DO
            DO i = iv, iv !iv-1,iv+1
              WRITE (*, 7) (tempY%fields(Y%getVarIdx(TRIM('pres_ctl')))%DATA(i, j, 1), j=21, 25), i, sg%mpddInfo_sg%myrank
              WRITE (*, 7) (tempY%fields(Y%getVarIdx(TRIM('pres_ctl')))%DATA(i, j, 1), j=16, 20), i, sg%mpddInfo_sg%myrank
              WRITE (*, 7) (tempY%fields(Y%getVarIdx(TRIM('pres_ctl')))%DATA(i, j, 1), j=11, 15), i, sg%mpddInfo_sg%myrank
              WRITE (*, 7) (tempY%fields(Y%getVarIdx(TRIM('pres_ctl')))%DATA(i, j, 1), j=6, 10), i, sg%mpddInfo_sg%myrank
              WRITE (*, 7) (tempY%fields(Y%getVarIdx(TRIM('pres_ctl')))%DATA(i, j, 1), j=1, 5), i, sg%mpddInfo_sg%myrank
7             FORMAT('Adjoint prro: ', 5D12.4, ' vertical', I2, ' pc', I2)
            END DO
          END IF

          rs2 = X.DOT.tempY   ! X^T H^T Y

          ! Output the <Hu, v> and <u, H^Tv> values:
          IF (X%mpddGlob%myrank .EQ. 1) THEN
            WRITE (*, 6) rs1, rs2
6           FORMAT('<HX,  Y> =', D20.10, '     <X,H^TY> =', D20.10)
          END IF

        END BLOCK
        CALL CPU_TIME(t2)
        WRITE (*, 22) t2 - t1, app%mpddGlob%myrank
22      FORMAT('Time spent in GeosBal plotting:', D12.4, ' pc', I2)

      END BLOCK
    END ASSOCIATE
  END BLOCK

  ! Check if test passes:
  IF (app%mpddGlob%isBaseProc()) THEN
    IF (ABS(rs1 - rs2) .LE. 1.0D-8) THEN
      PRINT *, 'Test passed'
    ELSE
      WRITE (*, 12) rs1, rs2, rs1 - rs2, app%mpddGlob%myrank
12    FORMAT('Test failed!', 2D20.12, ' diff: ', D20.12, I2)
    END IF
  END IF

  CALL app%mpddGlob%finalize

END PROGRAM Test_GeosBalTLAD

SUBROUTINE testState(X)
  USE kinds_m, ONLY: i_kind, r_kind
  USE State_m, ONLY: State_t
  IMPLICIT NONE
  TYPE(State_t), INTENT(INOUT) :: X

  ! Local variables:
  INTEGER(i_kind) :: i, j, ip, ir, it, ll
  TYPE(State_t) :: O

  !X = X%zeroCopy()
  O = X - X
  DO j = 1, X%sg%vLevel
    DO i = 1, X%sg%num_cell
      IF (ISNAN(X%sg%sigmax(j, i)) .OR. ISNAN(X%sg%sigmay(j, i))) THEN
        WRITE (*, 1) j, i, X%sg%sigmax(j, i), X%sg%sigmay(j, i), X%mpddGlob%myrank
1       FORMAT('testState - NaN sigma: ', I3, I8, ' sigmaX-Y: ', 2D12.4, ' proc: ', I2)
      END IF
    END DO
  END DO

  !PRINT*,'in side of testState... '
  !CALL flush
  !call X%mpddSub%barrier

  ir = X%getVarIdx(TRIM('rho'))
  ip = X%getVarIdx(TRIM('pres'))
  it = X%getVarIdx(TRIM('temp'))
  X = X%zeroCopy()
  ASSOCIATE (tem => X%Fields(it)%DATA)
    ll = 1
    DO i = 1, X%sg%num_cell
      DO j = 1, X%sg%vLevel
        tem(j, i, :) = X%sg%cell_cntr(ll, i) ! DBLE(X%sg%vLevel)-DBLE(j)+1 ! 1.0D0
        IF (ir .GT. 0) X%Fields(ir)%DATA(j, i, :) = X%sg%cell_cntr(ll, i) !1.0D0 DBLE(X%sg%vLevel)-DBLE(j)+1 !2.0D0
        IF (ip .GT. 0) X%Fields(ip)%DATA(j, i, :) = X%sg%cell_cntr(ll, i) !1.0D0 DBLE(X%sg%vLevel)-DBLE(j)+1 !2.0D0
      END DO
    END DO
  END ASSOCIATE
END SUBROUTINE testState

SUBROUTINE GetTopo(sg)
  USE kinds_m, ONLY: i_kind, r_kind
  USE SingleGrid_m, ONLY: SingleGrid_t
  IMPLICIT NONE
  TYPE(SingleGrid_t), INTENT(INOUT) :: sg

  ! Local variables:
  INTEGER(i_kind) :: i
  REAL(r_kind) :: topoHeight, topoRadius, x, y
  REAL(r_kind) :: swap, minGlobLat, maxGlobLat, minGlobLon, maxGlobLon

  ! Find the global latlon:
  swap = MINVAL(sg%cell_cntr(1, 1:sg%num_icell))
  CALL sg%mpddInfo_sg%AllReduceMinReal(swap, minGlobLat)

  swap = MAXVAL(sg%cell_cntr(1, :))
  CALL sg%mpddInfo_sg%AllReduceMaxReal(swap, maxGlobLat)

  swap = MINVAL(sg%cell_cntr(2, :))
  CALL sg%mpddInfo_sg%AllReduceMinReal(swap, minGlobLon)

  swap = MAXVAL(sg%cell_cntr(2, :))
  CALL sg%mpddInfo_sg%AllReduceMaxReal(swap, maxGlobLon)

  ! Assigning a Gaussian bell as a topography:
  topoHeight = 1.0D3 ! 1km
  topoRadius = 4.0D0 ! domain ratio, 0.5 is half the size of domain
  DO i = 1, sg%num_cell
    x = -1.0D0 + 2.0D0 * (sg%cell_cntr(2, i) - minGlobLon) / (maxGlobLon - minGlobLon)
    y = -1.0D0 + 2.0D0 * (sg%cell_cntr(1, i) - minGlobLat) / (maxGlobLat - minGlobLat)
    sg%topo(i) = topoHeight * EXP(-topoRadius * (x**2 + y**2))
  END DO

END SUBROUTINE GetTopo
