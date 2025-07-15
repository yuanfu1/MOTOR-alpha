PROGRAM Test_Halo
  USE kinds_m, ONLY: i_kind, r_kind
  USE geometry_m, ONLY: geometry_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE IOGrapes_m, ONLY: IOGrapes_t
  USE State_m, ONLY: State_t
  USE State2NC_m

  IMPLICIT NONE

  TYPE(State_t) :: X, dX
  TYPE(mpddGlob_t) :: mpddGlob
  TYPE(geometry_t) :: geometry
  TYPE(IOGrapes_t), TARGET :: ioGrapes
  INTEGER(i_kind) :: k, i
  REAL(r_kind) :: rs1, rs2
  CHARACTER(LEN=1024) :: configFile, outputDir

  configFile = "App_3DVarVerification_TestHalo.yaml"
  outputDir = "/home/ubt/workspace/code/MOTOR_feature_halo_width/build/Debug/"
  PRINT *, '------ RUN TEST Halo CASE ------ with config file:', configFile

  ! Initialize the mpdd
  CALL mpddGlob%initialize()

  ! Initialize geometry:
  CALL geometry%initialize(configFile, mpddGlob)

  ASSOCIATE (sg => geometry%mg%sg(7))

    BLOCK

      CALL X%initialize(configFile, sg)
      CALL ioGrapes%initialize(configFile, geometry)
      CALL ioGrapes%m_read_bcg_into_Xm(X, sg)

      CALL mpddGlob%barrier

      CALL Output_NC_State_AV(X, TRIM(outputDir), 'HaloTest_bcg', .TRUE., .FALSE.)

      CALL mpddGlob%barrier

      DO i = LBOUND(X%fields, 1), UBOUND(X%fields, 1)
        X%fields(i)%DATA = 1
      END DO
      CALL mpddGlob%barrier

      CALL Output_NC_State_AV(X, TRIM(outputDir), 'HaloTest_SET1', .TRUE., .FALSE.)

      CALL mpddGlob%barrier

      CALL X%exHaloRevSum()

      CALL mpddGlob%barrier

      CALL Output_NC_State_AV(X, TRIM(outputDir), 'HaloTest_RevSum', .TRUE., .FALSE.)

    END BLOCK
  END ASSOCIATE

  CALL mpddGlob%barrier

  ! Finalize
  CALL mpddGlob%finalize

END PROGRAM Test_Halo
