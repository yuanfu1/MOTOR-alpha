!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
! Created by Zilong Qin (zilong.qin@gmail.com), 2022/3/23, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!!
PROGRAM Test_UV2WParaTest
  USE kinds_m, ONLY: i_kind, r_kind, r_double
  USE UV2W_m, ONLY: UV2W_t
  USE geometry_m, ONLY: geometry_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE MPObs_m, ONLY: MPObs_t
  USE State_m, ONLY: State_t
  USE ObsSet_m, ONLY: ObsSet_t
  USE IOGrapes_m, ONLY: IOGrapes_t
  USE YAMLRead_m
  USE State2NC_m
  USE UV2W_m, ONLY: UV2W_t

  IMPLICIT NONE

  ! Define types
  TYPE(mpddGlob_t), TARGET :: mpddGlob
  TYPE(geometry_t), TARGET :: geometry
  REAL(r_kind) :: rs1, rs2, f1, f2
  TYPE(IOGrapes_t), TARGET :: ioGrapes

  CHARACTER(len=1024) :: configFile
  CHARACTER(LEN=1024) :: staticDir
  CHARACTER(LEN=1024) :: outputDir

  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", staticDir)
  configFile = TRIM(staticDir)//"/UnitTest"//"/Test_UV2WParaTest.yaml"

  IF (yaml_get_var(configFile, 'IO', 'output_dir', outputDir) /= 0) STOP

  ! Initializer
  ! Auxtypes
  CALL mpddGlob%initialize()                                   ! Initialize the mpdd

  ! Initialize geometry
  CALL geometry%initialize(configFile, mpddGlob)               ! Initialize the geometry
  ! Initial
  CALL ioGrapes%initialize(configFile, geometry)

  CALL mpddGlob%barrier
  ASSOCIATE (sg => geometry%mg%sg(8))
    BLOCK

      TYPE(MPObs_t), TARGET :: mpObs
      TYPE(State_t) :: X, Y, tempX
      TYPE(UV2W_t) :: H

      INTEGER(i_kind) :: i, j, k

      CALL mpObs%initializeMPObs(sg)
      CALL X%initialize(configFile, sg)
      ! CALL X%setAllFieldData(1.0D0)
      CALL ioGrapes%m_read_bcg_into_Xm(X, sg)

      PRINT *, 'Initialize'
      H = UV2W_t(configFile, X)
      PRINT *, 'End initialize'

      CALL H%fwdNL(X)

      X%Fields(3)%DATA(:, :, 1) = sg%zHght

      IF (mpddGlob%nProc == 1) THEN
        CALL Output_NC_State_AV(X, outputDir, "Test_UV2WParaTest_serial")
      ELSE
        CALL Output_NC_State_AV(X, outputDir, "Test_UV2WParaTest_parallel")
      END IF
    END BLOCK
  END ASSOCIATE

!   IF (mpddGlob%isBaseProc()) THEN
!     IF (ABS(f1 - 1.0) < 1e-7 .AND. ABS(f2 - 1.0) < 1e-7 .AND. ABS(rs1 - rs2) < 1e-7) THEN
!       PRINT *, 'Test passed!'
!     ELSE
!       PRINT *, 'Test failed!'
!     END IF
!   END IF

  CALL mpddGlob%finalize
END PROGRAM Test_UV2WParaTest
