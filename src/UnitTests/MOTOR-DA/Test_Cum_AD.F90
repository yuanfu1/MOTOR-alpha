!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather
! Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Jilong CHEN
! VERSION           : V 0.0
! HISTORY           :
!   Created by Jilong CHEN (jchen@link.cuhk.edu.hk), 2021/1/26, @GBA-MWF,
!   Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!!

PROGRAM Test_Cum_AD
  !USE NMLRead_m
  USE YAMLRead_m
  USE PhyCon_m, ONLY: PhyCon_t
  USE RK4TL_m, ONLY: CumTL_t
  USE kinds_m, ONLY: i_kind, r_kind
  USE State_m, ONLY: State_t
  USE geometry_m, ONLY: geometry_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE State2NC_m
  USE UnitTestData_m, ONLY: UnitTestData

  CHARACTER(LEN=1024) :: configFile, ncOutputFile
  TYPE(mpddGlob_t), TARGET :: mpddGlob
  TYPE(geometry_t), TARGET :: geometry
  TYPE(State_t) :: X, dX
  TYPE(PhyCon_t) :: PhyCon
  TYPE(CumTL_t) :: CumTL
  INTEGER(i_kind) :: i, j, k, size_3
  REAL(r_kind) :: Y_dot_Lx_q, LY_dot_x_q, &
                  Y_dot_Lx_theta, LY_dot_x_theta, &
                  theta_ratio, q_ratio

  ! Get the configFile
  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", configFile)
  configFile = TRIM(configFile)//"/UnitTest/testCum.yaml"

  CALL GET_ENVIRONMENT_VARIABLE("GRID_DIR", ncOutputFile)
  ncOutputFile = TRIM(ncOutputFile)//"/testCum.nc"

  ! Initialize
  CALL mpddGlob%initialize()
  CALL geometry%initialize(configFile, mpddGlob)          ! Initialize the geometry

  ASSOCIATE (sg => geometry%mg%sg(1))
    CALL X%initialize(configFile, sg)
    CALL dX%initialize(configFile, sg)

    ! put values
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

    CALL PhyCon%FwdDataPre(X)

    CALL PhyCon%RK4(X)

    CALL PhyCon%AdjDataPre(X, dX)

    CALL CumTL%TLDataPre(PhyCon)

    CumTL%qvapor_d(20, 5, 1) = 0.0000005D0

    CALL CumTL%RK4_TL(PhyCon, X)

    size_3 = SIZE(PhyCon%qvapor_out, dim=3)

    PhyCon%uwnd_outb(:, :, size_3) = CumTL%uwnd_d(:, :, size_3)
    PhyCon%vwnd_outb(:, :, size_3) = CumTL%vwnd_d(:, :, size_3)
    PhyCon%wwnd_outb(:, :, size_3) = CumTL%wwnd_d(:, :, size_3)
    PhyCon%theta_outb(:, :, size_3) = CumTL%theta_d(:, :, size_3)
    PhyCon%qvapor_outb(:, :, size_3) = CumTL%qvapor_d(:, :, size_3)
    PhyCon%pres_outb(:, :, size_3) = CumTL%pres_d(:, :, size_3)
    PhyCon%precipb(1, :, size_3) = CumTL%precipd(1, :, size_3)

    CALL PhyCon%RK4_AD(X)

    Y_dot_Lx_q = SUM(CumTL%theta_d(:, :, size_3) * CumTL%theta_d(:, :, size_3)) &
                 + SUM(CumTL%qvapor_d(:, :, size_3) * CumTL%qvapor_d(:, :, size_3)) &
                 + SUM(CumTL%precipd(:, :, size_3) * CumTL%precipd(:, :, size_3))
    LY_dot_x_q = SUM(PhyCon%qvapor_outb(:, 5, 1) * CumTL%qvapor_d(:, 5, 1))

    q_ratio = Y_dot_Lx_q / LY_dot_x_q

    PRINT *, "qvapor_ratio is", q_ratio

    CumTL%uwnd_d = 0.0D0
    CumTL%vwnd_d = 0.0D0
    CumTL%wwnd_d = 0.0D0
    CumTL%theta_d = 0.0D0
    CumTL%qvapor_d = 0.0D0
    CumTL%pres_d = 0.0D0
    CumTL%precipd = 0.0D0

    PhyCon%uwnd_outb = 0.0D0
    PhyCon%vwnd_outb = 0.0D0
    PhyCon%wwnd_outb = 0.0D0
    PhyCon%theta_outb = 0.0D0
    PhyCon%qvapor_outb = 0.0D0
    PhyCon%pres_outb = 0.0D0
    PhyCon%precipb = 0.0D0

    CumTL%theta_d(20, 5, 1) = 0.05D0

    CALL CumTL%RK4_TL(PhyCon, X)

    PhyCon%uwnd_outb(:, :, size_3) = CumTL%uwnd_d(:, :, size_3)
    PhyCon%vwnd_outb(:, :, size_3) = CumTL%vwnd_d(:, :, size_3)
    PhyCon%wwnd_outb(:, :, size_3) = CumTL%wwnd_d(:, :, size_3)
    PhyCon%theta_outb(:, :, size_3) = CumTL%theta_d(:, :, size_3)
    PhyCon%qvapor_outb(:, :, size_3) = CumTL%qvapor_d(:, :, size_3)
    PhyCon%pres_outb(:, :, size_3) = CumTL%pres_d(:, :, size_3)
    PhyCon%precipb(1, :, size_3) = CumTL%precipd(1, :, size_3)

    CALL PhyCon%RK4_AD(X)

    Y_dot_Lx_theta = SUM(CumTL%theta_d(:, :, size_3) * CumTL%theta_d(:, :, size_3)) &
                     + SUM(CumTL%qvapor_d(:, :, size_3) * CumTL%qvapor_d(:, :, size_3)) &
                     + SUM(CumTL%precipd(:, :, size_3) * CumTL%precipd(:, :, size_3))
    LY_dot_x_theta = SUM(PhyCon%theta_outb(:, 5, 1) * &
                         CumTL%theta_d(:, 5, 1))

    theta_ratio = Y_dot_Lx_theta / LY_dot_x_theta

    PRINT *, "theta_ratio is", theta_ratio

    IF (dabs(theta_ratio - 1.0D0) < 1.0D-11 .AND. dabs(q_ratio - 1.0D0) < 1.0D-11) THEN
      PRINT *, "tested passed"
    ELSE
      PRINT *, "tested failed"
    END IF

  END ASSOCIATE

  CALL mpddGlob%finalize

END PROGRAM Test_Cum_AD

