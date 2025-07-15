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
PROGRAM Test_ObsRadarVel
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
  USE UV2Divg_m, ONLY: UV2Divg_t

  IMPLICIT NONE

  ! Define types
  TYPE(mpddGlob_t), TARGET :: mpddGlob
  TYPE(geometry_t), TARGET :: geometry
  REAL(r_kind) :: rs1, rs2, f1, f2
  TYPE(IOGrapes_t), TARGET :: ioGrapes

  CHARACTER(len=1024) :: configFile
  CHARACTER(LEN=1024) :: staticDir
  CHARACTER(LEN=1024) :: outputDir
  INTEGER(i_kind):: fileStat

  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", staticDir)
  configFile = TRIM(staticDir)//"/UnitTest"//"/Test_UV2WTLAD.yaml"
  configFile = "/Users/qzl/sources/MOTOR/input/220527_0000_3km/App_3DVarVerification.yaml"

  fileStat = yaml_get_var(configFile, 'IO', 'output_dir', outputDir)
  ! Initializer
  ! Auxtypes
  CALL mpddGlob%initialize()                                   ! Initialize the mpdd

  ! Initialize geometry
  CALL geometry%initialize(configFile, mpddGlob)               ! Initialize the geometry
  ! Initial

  CALL mpddGlob%barrier
  ASSOCIATE (sg => geometry%mg%sg(7))
    BLOCK
      TYPE(MPObs_t), TARGET :: mpObs
      TYPE(State_t) :: X, Y, tempX
      TYPE(UV2Divg_t) :: H

      INTEGER(i_kind) :: i, j, k

      CALL mpObs%initializeMPObs(sg)
      CALL X%initialize(configFile, sg)
      ! CALL X%setAllFieldData(1.0D0)
      CALL Y%initialize(configFile, sg)

      BLOCK
        PRINT *, 'Initialize'
        CALL ioGrapes%initialize(configFile, geometry)
        CALL ioGrapes%m_read_bcg_into_Xm(X, sg)
        PRINT *, 'End initialize'
      END BLOCK

      PRINT *, 'Initialize'
      H = UV2Divg_t(configFile, X)
      PRINT *, 'End initialize'

      ! CALL X%addVar('wwnd')
      ! CALL Y%addVar('wwnd')

      BLOCK
        ! Y = H%fwdNL_opr(X);
        CALL Y%addVar('divg')
        Y%Fields(Y%getVarIdx('divg'))%DATA = 1.0D0
        rs1 = Y .DOT. (H%fwdNL_opr(X))
        rs2 = (H%transAdjMultiply_opr(Y, X)) .DOT.X

        IF (mpddGlob%isBaseProc()) THEN
          PRINT *, 'Y*(H^-1/2*X):', rs1
          PRINT *, 'H^-1/2^T*X*Y:', rs2
        END IF
      END BLOCK

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
END PROGRAM Test_ObsRadarVel
