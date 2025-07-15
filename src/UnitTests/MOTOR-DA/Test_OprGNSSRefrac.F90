PROGRAM Test_OprGNSSRefrac
  USE kinds_m, ONLY: i_kind, r_kind
  USE geometry_m, ONLY: geometry_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE MPObs_m, ONLY: MPObs_t
  USE State_m, ONLY: State_t
  USE ObsSet_m, ONLY: ObsSet_t
  USE IOGrapes_m, ONLY: IOGrapes_t
  USE State2NC_m
  USE OprGNSSRefrac_m, ONLY: OprGNSSRefrac_t
  USE ObsGNSS_m, ONLY: ObsGNSS_t
  USE C2O_m, ONLY: C2O_t

  IMPLICIT NONE

  TYPE(ObsSet_t) :: Y_obs, Z_obs
  TYPE(ObsSet_t) :: tl_01, tl_02
  TYPE(State_t) :: X, dX
  TYPE(OprGNSSRefrac_t) :: GNSSRefrac
  TYPE(mpddGlob_t) :: mpddGlob
  TYPE(geometry_t) :: geometry
  TYPE(MPObs_t), TARGET :: mpObs
  TYPE(IOGrapes_t), TARGET :: ioGrapes
  TYPE(ObsGNSS_t) :: OBS_GNSS
  INTEGER(i_kind) :: k, i
  REAL(r_kind) :: rs1, rs2
  CHARACTER(LEN=1024) :: configFile
  TYPE(C2O_t) :: H
  REAL(r_kind) :: rhs, lhs, ad_check, tl_check_min, tl_check_max, tl_check_avg

  configFile = "App_3DVarVerification_gnss.yaml"

  PRINT *, '------ RUN TEST Opr GNSS CASE ------ with config file:', configFile

  ! Initialize the mpdd
  CALL mpddGlob%initialize()

  ! Initialize geometry:
  CALL geometry%initialize(configFile, mpddGlob)

  ASSOCIATE (sg => geometry%mg%sg(5))

    !! init X and read obs
    CALL X%initialize(configFile, sg)
    CALL mpObs%initializeMPObs(sg)
    CALL ioGrapes%initialize(configFile, geometry)
    CALL ioGrapes%m_read_bcg_into_Xm(X, sg)

    CALL OBS_GNSS%ObsInitial(configFile)
    CALL OBS_GNSS%ObsIngest(X)

    CALL OBS_GNSS%ObsThinning(X, Y_obs, mpObs, .TRUE., .FALSE.)

    PRINT *, 'after thinning obs count: ', UBOUND(Y_obs%obsFields(1)%values, 1)

    CALL H%initialize(configFile, X, Y_obs)
    GNSSRefrac = OprGNSSRefrac_t(configFile, X, Y_obs)

    !! Test forward
    Z_obs = Y_obs%zeroCopy()
    Z_obs = H%fwdNL_opr(X)

    DO i = LBOUND(Z_obs%obsFields(1)%values, 1), UBOUND(Z_obs%obsFields(1)%values, 1)
      PRINT *, Y_obs%obsFields(1)%values(i), Z_obs%obsFields(1)%values(i), Y_obs%obsFields(1)%values(i) - Z_obs%obsFields(1)%values(i)
    END DO

    PRINT *, ' *********** '
    PRINT *, ' **** T **** '
    PRINT *, ' **** E **** '
    PRINT *, ' **** S **** '
    PRINT *, ' **** T **** '
    PRINT *, ' ****   **** '
    PRINT *, ' **** T **** '
    PRINT *, ' **** L **** '
    PRINT *, ' *********** '

    !! TEST TL
    dX = X%zeroCopy()
    dX%fields(dX%getVarIdx('temp'))%DATA = X%fields(dX%getVarIdx('temp'))%DATA * 0.001D0
    ! dX%fields(dX%getVarIdx('pres'))%DATA = X%fields(dX%getVarIdx('pres'))%DATA * 0.001D0
    dX%fields(dX%getVarIdx('qvapor'))%DATA = X%fields(dX%getVarIdx('qvapor'))%DATA * 0.001D0
    tl_01 = ((H%fwdNL_opr(X + dX / 2.0D0)) - (H%fwdNL_opr(X - dX / 2.0D0)))
    tl_02 = H%fwdTL_opr(dX, X)

    DO i = LBOUND(tl_01%obsFields(1)%values, 1), UBOUND(tl_01%obsFields(1)%values, 1)
      PRINT *, tl_01%obsFields(1)%values(i), tl_02%obsFields(1)%values(i), tl_01%obsFields(1)%values(i) - tl_02%obsFields(1)%values(i)
    END DO

    PRINT *, ' *********** '
    PRINT *, ' **** T **** '
    PRINT *, ' **** E **** '
    PRINT *, ' **** S **** '
    PRINT *, ' **** T **** '
    PRINT *, ' ****   **** '
    PRINT *, ' **** A **** '
    PRINT *, ' **** D **** '
    PRINT *, ' *********** '

    lhs = Y_obs.DOT. (H%fwdTL_opr(dX, X))
    rhs = (H%adjMul_opr(Y_obs, X)) .DOT.dX
    ad_check = lhs - rhs
    PRINT *, "ad_check: ", lhs, rhs, ad_check

    ! !! TEST FORWARD END

    ! !! TEST ADJ START
    ! dX = X%zeroCopy()
    ! CALL dX%setAllFieldData(0.0D0)
    ! dX%fields(dX%getVarIdx('temp'))%DATA = X%fields(dX%getVarIdx('temp'))%DATA * 0.0001D0
    ! ! dX%fields(dX%getVarIdx('pres'))%DATA = X%fields(dX%getVarIdx('pres'))%DATA * 0.01D0
    ! ! dX%fields(dX%getVarIdx('qvapor'))%DATA = X%fields(dX%getVarIdx('qvapor'))%DATA * 0.01

    ! tl_01 = ((H%fwdNL_opr(X + dX / 2.0D0)) - (H%fwdNL_opr(X - dX/2.0D0)))
    ! tl_02 = (H%fwdTL_opr(dX, X))
    ! ! PRINT*, 'tl_01 :', tl_01%ObsFields(1)%values
    ! ! PRINT*, 'tl_02 :', tl_02%ObsFields(1)%values
    ! Y_tl_check =  tl_01 / tl_02

    ! ! tl_01 = GNSSRefrac%fwdNL_opr(X +  dX / 2.0D0) - GNSSRefrac%fwdNL_opr(X -  dX / 2.0D0)

    ! ! PRINT*, dX%fields(dX%getVarIdx('temp'))%DATA
    ! ! PRINT*, X%fields(dX%getVarIdx('temp'))%DATA
    ! ! tl_02 = GNSSRefrac%fwdTL_opr(dX, X)

    ! ! Y_tl_check = tl_01 / tl_02

    ! ! PRINT*, 'Y_tl_check:', Y_tl_check%ObsFields(1)%values

    ! tl_check_max = MAXVAL(Y_tl_check%ObsFields(1)%values)
    ! tl_check_min = MINVAL(Y_tl_check%ObsFields(1)%values)
    ! tl_check_avg = SUM(Y_tl_check%ObsFields(1)%values) / REAL(SIZE(Y_tl_check%ObsFields(1)%values))

    ! lhs = Y_obs.DOT. (H%fwdTL_opr(dX, X))
    ! X1 = H%adjMul_opr(Y_obs, X)
    ! ! X2 =
    ! rhs = X1%fields(X%getVarIdx('temp')) .DOT.dX%fields(dX%getVarIdx('temp'))
    ! ad_check = lhs - rhs
    ! ! PRINT*, dX_adj%fields(dX%gx/etVarIdx('temp'))%DATA - dX%fields(dX%getVarIdx('temp'))%DATA
    ! PRINT *, "ad_check: ", lhs, rhs,  ad_check

    ! IF (mpddGlob%isBaseProc()) THEN
    !   ! IF (ABS(tl_check_avg - 1.0D0) .LE. 1E-7 .AND. ABS(ad_check) .LE. 1E-7) THEN
    !     PRINT *, 'Test passed!'
    !     PRINT *, 'TL check accuracy: ', ABS(tl_check_avg - 1.0D0)
    !     PRINT *, 'AD check accuracy: ', ABS(ad_check)
    !   ! ELSE
    !   !   PRINT *, 'Test failed!'
    !   ! END IF
    ! END IF
  END ASSOCIATE

  CALL mpddGlob%barrier

  ! Finalize
  CALL mpddGlob%finalize

END PROGRAM Test_OprGNSSRefrac
