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

PROGRAM Test_Cum_interface
  !USE NMLRead_m
  USE YAMLRead_m
  USE CumInterface_m, ONLY: CumInterface_t
  USE kinds_m, ONLY: i_kind, r_kind
  USE State_m, ONLY: State_t
  USE geometry_m, ONLY: geometry_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE State2NC_m
  USE UnitTestData_m, ONLY: UnitTestData
  USE Pres2Rho_m, ONLY: Pres2Rho_t
  USE Theta2Temp_m, ONLY: Theta2Temp_t
  USE Qvapor2Rhov_m, ONLY: Qvapor2Rhov_t

  IMPLICIT NONE
  CHARACTER(LEN=1024) :: configFile, outputDir, &
                         InitFlag, FwdFlag, AdjFlag
  TYPE(mpddGlob_t), TARGET :: mpddGlob
  TYPE(geometry_t), TARGET :: geometry
  ! TYPE(Cumulus_t) :: cumulus
  TYPE(State_t) :: X, dX
  TYPE(CumInterface_t) :: CumInterface
  TYPE(Pres2Rho_t) :: Pres2Rho
  TYPE(Theta2Temp_t) :: Theta2Temp
  TYPE(Qvapor2Rhov_t) :: Qvapor2Rhov
  INTEGER(i_kind) :: i, j, k
  REAL(r_kind) :: theta_ratio

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

    dX%Fields(4)%DATA(:, :, 2) = 0.005D0
    ! dX%Fields(5)%data(20, 5, 2) = 0.0000000D0

    PRINT *, X%sg%cell_type

    CALL CumInterface%CumForward(X)

    IF (dabs((CumInterface%uwnd_out(3, 5, 1) - CumInterface%uwnd_out(3, 5, 7)) / &
             (CumInterface%uwnd_out(3, 5, 1) - CumInterface%uwnd_out(3, 5, 13)) &
             - 0.5D0) < 1.0D-12) THEN
      InitFlag = "passed"
      PRINT *, "Subrotine Fwdinitial in the CumInterface%CumForward is tested ", InitFlag
    ELSE
      InitFlag = "failed"
      PRINT *, "Subrotine Initialize in the CumInterface is tested ", InitFlag
    END IF

    IF (dabs(X%Fields(4)%DATA(20, 5, 2) - 312.71451627799695D0) < 1.0D-14 .AND. &
        dabs(X%Fields(5)%DATA(20, 5, 2) - 1.0740591796961579D-2) < 1.0D-14) THEN
      FwdFlag = "passed"
      PRINT *, "Subrotine CumForward in the CumInterface is tested ", FwdFlag
    ELSE
      FwdFlag = "failed"
      PRINT *, "Subrotine CumForward in the CumInterface is tested ", FwdFlag
    END IF

    CALL CumInterface%CumAdjoint(X, dX)

    theta_ratio = dX%fields(4)%DATA(20, 5, 1) / CumInterface%theta_outb(20, 5, 1)

    IF (dabs(theta_ratio - 1.0D0) < 1.0D-15) THEN
      AdjFlag = "passed"
      PRINT *, "Subrotine Adjoint in the CumInterface is tested ", AdjFlag
    ELSE
      AdjFlag = "failed"
      PRINT *, "Subrotine Adjoint in the CumInterface is tested ", AdjFlag
    END IF

    IF (InitFlag .EQ. "passed" .AND. FwdFlag .EQ. "passed" .AND. AdjFlag .EQ. "passed") THEN
      PRINT *, "tested passed"
    ELSE
      PRINT *, "tested failed"
    END IF

  END ASSOCIATE

  PRINT *, "theta is ", X%fields(x%getVarIdx("theta"))%DATA(1, 3, 1)
  PRINT *, "qvapor is ", X%fields(x%getVarIdx("qvapor"))%DATA(1, 3, 1)
  PRINT *, "pres is ", X%fields(x%getVarIdx("pres"))%DATA(1, 3, 1)

  CALL X%addVar("rho")
  CALL X%addVar("temp")
  CALL X%addVar("rhov")

  CALL Theta2Temp%transBackward(X)
  CALL Pres2Rho%transBackward(X)
  CALL Qvapor2Rhov%transBackward(X)

  PRINT *, "temp is", X%fields(x%getVarIdx("temp"))%DATA(1, 3, 1)
  PRINT *, "rhov is", X%fields(x%getVarIdx("rhov"))%DATA(1, 3, 1)
  PRINT *, "rho is", X%fields(x%getVarIdx("rho"))%DATA(1, 3, 1)

  CALL Theta2Temp%fwdNL(X)
  CALL Pres2Rho%fwdNL(X)
  CALL Qvapor2Rhov%fwdNL(X)

  PRINT *, "theta after converting is", X%fields(x%getVarIdx("theta"))%DATA(1, 3, 1)
  PRINT *, "pres after converting is", X%fields(x%getVarIdx("pres"))%DATA(1, 3, 1)
  PRINT *, "qvapor after converting is", X%fields(x%getVarIdx("qvapor"))%DATA(1, 3, 1)

  PRINT *, "theta in dX is ", dX%fields(dX%getVarIdx("theta"))%DATA(2, 5, 1)
  PRINT *, "qvapor in dX is ", dX%fields(dX%getVarIdx("qvapor"))%DATA(2, 5, 1)
  PRINT *, "pres in dX is ", dX%fields(dX%getVarIdx("pres"))%DATA(2, 5, 1)

  CALL Theta2Temp%adjMul(dX, X)
  CALL Pres2Rho%adjMul(dX, X)
  CALL Qvapor2Rhov%adjMul(dX, X)

  PRINT *, "temp in dX is", dX%fields(dX%getVarIdx("temp"))%DATA(2, 5, 1)
  PRINT *, "rhov in dX is", dX%fields(dX%getVarIdx("rhov"))%DATA(2, 5, 1)
  PRINT *, "rho in dX is", dX%fields(dX%getVarIdx("rho"))%DATA(2, 5, 1)

END PROGRAM Test_Cum_interface
