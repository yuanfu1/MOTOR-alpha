!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather
! Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Jilong CHEN, Zilong Qin
! VERSION           : V 0.0
! HISTORY           :
!   Created by Jilong CHEN (jchen@link.cuhk.edu.hk), 2021/1/26, @GBA-MWF,
!   Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!!

PROGRAM Test_Cum
  !USE NMLRead_m
  USE YAMLRead_m
  USE PhyCon_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE State_m, ONLY: State_t
  USE geometry_m, ONLY: geometry_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE State2NC_m
  USE UnitTestData_m, ONLY: UnitTestData

  CHARACTER(LEN=1024) :: configFile, outputDir, &
                         DataPreFlag, CumFwdFlag
  TYPE(mpddGlob_t), TARGET :: mpddGlob
  TYPE(geometry_t), TARGET :: geometry
  ! TYPE(Cumulus_t) :: cumulus
  TYPE(State_t) :: X, dX
  TYPE(PhyCon_t) :: PhyCon

  INTEGER :: i, j, k

  ! Get the configFile
  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", configFile)
  configFile = TRIM(configFile)//"/UnitTest/testCum.yaml"

  IF (yaml_get_var(configFile, 'IO', 'output_dir', outputDir) /= 0) STOP

  ! Initialize
  CALL mpddGlob%initialize()
  CALL geometry%initialize(configFile, mpddGlob)          ! Initialize the geometry

  ASSOCIATE (sg => geometry%mg%sg(1))
    CALL X%initialize(configFile, sg)
    dCALL X%initialize(configFile, sg)

    CALL UnitTestData(X)

    X%Fields(1)%DATA(:, :, 2) = 1.01D0 * X%Fields(1)%DATA(:, :, 1)
    X%Fields(2)%DATA(:, :, 2) = 1.01D0 * X%Fields(2)%DATA(:, :, 1)
    X%Fields(3)%DATA(:, :, 2) = 1.0D0 * X%Fields(3)%DATA(:, :, 1)
    X%Fields(4)%DATA(:, :, 2) = 1.0D0 * X%Fields(4)%DATA(:, :, 1)
    X%Fields(5)%DATA(:, :, 2) = 1.0D0 * X%Fields(5)%DATA(:, :, 1)
    X%Fields(6)%DATA(:, :, 2) = 1.0D0 * X%Fields(6)%DATA(:, :, 1)

    dX%Fields(1)%DATA = 0.0D0
    dX%Fields(2)%DATA = 0.0D0
    dX%Fields(3)%DATA = 0.0D0
    dX%Fields(4)%DATA = 0.0D0
    dX%Fields(5)%DATA = 0.0D0
    dX%Fields(6)%DATA = 0.0D0

    dX%Fields(4)%DATA(:, :, 2) = 0.00D0
    dX%Fields(5)%DATA(20, 5, 2) = 0.0000000D0

    PRINT *, "num_cell, num_icell are", X%sg%num_cell, X%sg%num_icell

    CALL PhyCon%FwdDataPre(X)

    IF (dabs((PhyCon%uwnd_out(3, 5, 1) - PhyCon%uwnd_out(3, 5, 7)) / &
             (PhyCon%uwnd_out(3, 5, 1) - PhyCon%uwnd_out(3, 5, 13)) &
             - 0.5D0) < 1.0D-12) THEN

      DataPreFlag = "passed"
      PRINT *, "Subrotine FwdDataPre in PhyCon is tested ", DataPreFlag
    ELSE
      DataPreFlag = "failed"
      PRINT *, "Subrotine FwdDataPre in PhyCon is tested ", DataPreFlag
    END IF

    CALL PhyCon%RK4(X)

    IF (dabs(PhyCon%theta_out(20, 5, 13) - 312.71451627799695D0) < 1.0D-15 .AND. &
        dabs(PhyCon%qvapor_out(20, 5, 13) - 1.0740591796961579D-2) < 1.0D-15) THEN
      CumFwdFlag = "passed"
      PRINT *, "Subrotine CumFwd in PhyCon is tested ", CumFwdFlag
    ELSE
      CumFwdFlag = "failed"
      PRINT *, "Subrotine CumFwd in PhyCon is tested ", CumFwdFlag
    END IF

    IF (DataPreFlag .EQ. "passed" .AND. CumFwdFlag .EQ. "passed") THEN
      PRINT *, "tested passed"
    ELSE
      PRINT *, "tested failed"
    END IF

    CALL mpddGlob%finalize
  END ASSOCIATE

END PROGRAM Test_Cum
