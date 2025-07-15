PROGRAM Test_GeosBalTLAD
  USE kinds_m, ONLY: i_kind, r_kind, r_double
  USE GeosBal_m, ONLY: GeosBal_t
  USE geometry_m, ONLY: geometry_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE MPObs_m, ONLY: MPObs_t
  USE State_m, ONLY: State_t
  USE ObsSet_m, ONLY: ObsSet_t
  USE IOGrapes_m, ONLY: IOGrapes_t
  USE State2NC_m
  USE GeosBal_m, ONLY: GeosBal_t
  USE YAMLRead_m

  IMPLICIT NONE

  ! Define types
  !TYPE(mpddGlob_t), TARGET :: mpddGlob
  TYPE(mpddGlob_t) :: mpddGlob
  REAL(r_kind) :: rs1, rs2, f1, f2
  TYPE(IOGrapes_t), TARGET :: ioGrapes

  CHARACTER(len=1024) :: configFile
  CHARACTER(LEN=1024) :: staticDir
  CHARACTER(LEN=1024) :: outputDir
  CHARACTER(LEN=20) :: task
  INTEGER(i_kind) :: istatus

  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", staticDir)
  configFile = TRIM(staticDir)//"/UnitTest"//"/Test_NewGeosBalTLAD.yaml"
  istatus = yaml_get_var(configFile, 'IO', 'output_dir', outputDir)
  istatus = yaml_get_var(TRIM(configFile), 'RunMode', 'Task', task)
  PRINT *, 'Task is: ', task, istatus

  ! Initializer
  ! Auxtypes
  CALL mpddGlob%initialize()                                   ! Initialize the mpdd

  BLOCK
    TYPE(geometry_t) :: geometry
    ! Initialize geometry
    CALL geometry%initialize(configFile, mpddGlob)               ! Initialize the geometry

    ASSOCIATE (sg => geometry%mg%sg(3))
      BLOCK

        TYPE(MPObs_t), TARGET :: mpObs
        TYPE(State_t) :: X, Y, tempX, tempY, XB
        TYPE(GeosBal_t) :: H

        INTEGER(i_kind) :: i, j, k, ic, iv
        REAL(r_kind) :: xsquare

        CALL mpObs%initializeMPObs(sg)
        CALL X%initialize(configFile, sg)
        CALL tempX%initialize(configFile, sg)

        ! Initialize state data:
        PRINT *, 'Assigning values to field data.... ', &
          ALLOCATED(sg%cell_cntr), ALLOCATED(sg%edgeNorm2), X%mpddGlob%myrank

        IF (mpddGlob%isBaseProc()) THEN
          CALL X%setAllFieldData(7.0D0) ! base processor uses different values
        ELSE
          CALL X%setAllFieldData(9.0D0)
        END IF

        CALL XB%initialize(configFile, sg)
        CALL XB%setAllFieldData(3.0D0)
        CALL tempX%setAllFieldData(5.0D0)

        ! Get a test state:
        PRINT *, 'Getting a test State...'
        CALL testState(X)

        PRINT *, 'Initialize the Geostrophic balance model... ', sg%cell_type(1), X%mpddGlob%myrank
        !H = GeosBal_t(configFile, X)   ! Replaced by its initial routine as follows
        CALL H%initial(configFile, X)
        PRINT *, 'Geostrophic balance has been initialized at proc: ', X%mpddGlob%myrank

        ! Y is an arbitrary vector to check adjoint operation: <Y*X, Z> = <X, Y^T Z>
        CALL Y%initialize(configFile, sg)
        CALL tempY%initialize(configFile, sg)
        CALL Y%setAllFieldData(0.0D0)
        CALL tempY%setAllFieldData(0.0D0)

        ic = 20; iv = 2

        IF (X%mpddGlob%myrank .EQ. 1) THEN
          Y%fields(Y%getVarIdx(TRIM('vwnd')))%DATA(iv, ic, 1) = 1.0D0
          tempY%fields(Y%getVarIdx(TRIM('vwnd')))%DATA(iv, ic, 1) = 1.0D0
        END IF
        CALL testState(Y)
        CALL testState(tempY)

        BLOCK
          PRINT *, "Calling Forward of Tangent Linear operator: ", X%mpddGlob%myrank
          tempX = X             ! X: pres/rho
          CALL H%fwdTL(tempX, X)  ! tempX: temp/rho/u/v
          PRINT *, "Forward tangent linear operator has called at proc: ", &
            tempX%fields(Y%getVarIdx(TRIM('uwnd')))%DATA(1, ic, 1), X%mpddGlob%myrank

          CALL Output_NC_State_AV(tempX, outputDir, &
                                  TRIM(task)//"_HX", .TRUE., .TRUE.)

          rs1 = Y.DOT.tempX ! (Y)^T (H X)

          CALL H%adjMul(tempY, X)  ! tempX: pres/rho

          ! Debugging: 2022-09-11:
          DO i = iv, iv !iv-1,iv+1
            WRITE (*, 5) (tempY%fields(Y%getVarIdx(TRIM('temp')))%DATA(i, j, 1), j=21, 25), i, sg%mpddInfo_sg%myrank
            WRITE (*, 5) (tempY%fields(Y%getVarIdx(TRIM('temp')))%DATA(i, j, 1), j=16, 20), i, sg%mpddInfo_sg%myrank
            WRITE (*, 5) (tempY%fields(Y%getVarIdx(TRIM('temp')))%DATA(i, j, 1), j=11, 15), i, sg%mpddInfo_sg%myrank
            WRITE (*, 5) (tempY%fields(Y%getVarIdx(TRIM('temp')))%DATA(i, j, 1), j=6, 10), i, sg%mpddInfo_sg%myrank
            WRITE (*, 5) (tempY%fields(Y%getVarIdx(TRIM('temp')))%DATA(i, j, 1), j=1, 5), i, sg%mpddInfo_sg%myrank
5           FORMAT('Adjoint tempY: ', 5D12.4, ' vertical', I2, ' pc', I2)
          END DO

          ! Test if tempY%rho is cancelled out: This test passes!
          ! CALL X%setAllFieldData(0.0D0)
          ! X%fields(2)%data = 1.0D0
          PRINT *, 'Value of rho: ', MINVAL(X%fields(2)%DATA), MAXVAL(X%fields(2)%DATA), sg%mpddInfo_sg%myrank

          rs2 = X.DOT.tempY   ! X^T H^T Y

          PRINT *, 'Y x (H(X)):', rs1, mpddGlob%myrank
          PRINT *, 'X x (H^TY):', rs2, mpddGlob%myrank

        END BLOCK

      END BLOCK
    END ASSOCIATE
  END BLOCK
  WRITE (*, 1) mpddGlob%myrank
1 FORMAT('Test successully completes... ', I2)

  CALL mpddGlob%finalize

!   IF (mpddGlob%isBaseProc()) THEN
!     IF (ABS(f1 - 1.0) < 1e-7 .AND. ABS(f2 - 1.0) < 1e-7 .AND. ABS(rs1 - rs2) < 1e-7) THEN
!       PRINT *, 'Test passed!'
!     ELSE
!       PRINT *, 'Test failed!'
!     END IF
!   END IF

END PROGRAM Test_GeosBalTLAD

SUBROUTINE testState(X)
  USE kinds_m, ONLY: i_kind, r_kind
  USE State_m, ONLY: State_t
  IMPLICIT NONE
  TYPE(State_t), INTENT(INOUT) :: X

  ! Local variables:
  INTEGER(i_kind) :: i, j
  TYPE(State_t) :: O

  !X = X%zeroCopy()
  O = X - X
  DO j = 1, UBOUND(X%fields, 1)
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

  ASSOCIATE (tem => X%Fields(X%getVarIdx(TRIM('temp')))%DATA, &
             rho => X%Fields(X%getVarIdx(TRIM('rho')))%DATA)
    DO i = 1, X%sg%num_cell
      tem(:, i, :) = X%sg%cell_cntr(2, i)
      rho(:, i, :) = X%sg%cell_cntr(1, i) !1.0D0
    END DO
  END ASSOCIATE
END SUBROUTINE testState
