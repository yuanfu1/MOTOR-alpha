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
  USE UV2W_Poisson_m, ONLY: UV2W_Poisson_t

  IMPLICIT NONE

  ! Define types
  TYPE(mpddGlob_t), TARGET :: mpddGlob
  TYPE(geometry_t), TARGET :: geometry
  REAL(r_kind) :: rs1, rs2, f1, f2
  TYPE(IOGrapes_t), TARGET :: ioGrapes

  CHARACTER(len=1024) :: configFile
  CHARACTER(LEN=1024) :: staticDir
  CHARACTER(LEN=1024) :: outputDir

  INTEGER(i_kind) :: fileStat

  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", staticDir)
  configFile = TRIM(staticDir)//"/UnitTest"//"/Test_UV2WTLAD.yaml"
  configFile = "/Users/qzl/sources/MOTOR/input/220527_0000_3km/App_UV2W.yaml"

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
      CLASS(UV2W_t), ALLOCATABLE :: H

      INTEGER(i_kind) :: i, j, k

      CALL mpObs%initializeMPObs(sg)
      CALL X%initialize(configFile, sg)
      ! CALL X%setAllFieldData(1.0D0)

      BLOCK
        CALL ioGrapes%initialize(configFile, geometry)
        CALL ioGrapes%m_read_bcg_into_Xm(X, sg)
      END BLOCK

      PRINT *, 'Initialize'
      H = UV2W_Poisson_t(configFile, X)
      PRINT *, 'End initialize'

      CALL Y%initialize(configFile, sg)

      CALL H%fwdNL(X)
      CALL Output_NC_State_AV(X, outputDir, "Test_UV2WTLAD")

      CALL Y%addVar('wwnd')
      CALL Y%setAllFieldData(1.0D0)
      ! Y%getVar('wwnd')

      ! CALL Output_NC_State_AV(Y, outputDir, "Test_UV2WTLAD")

      ! Test 1:  TL test
      !   BLOCK
      !     TYPE(State_t) :: dX

      !     dCALL X%initialize(configFile, sg)
      !     CALL dX%setAllFieldData(0.1D0)

      !     fY = ((H.TL. (X + dX)) - (H.TL. (X - dX)))/((H.TL.dX)*2.0D0)

      !     f1 = mpddGlob%allRedMaxReal(MAXVAL(fY%ObsFields(1)%values))
      !     IF (mpddGlob%isBaseProc()) PRINT *, 'MAXVAL((H(X+dX)-H(X-dX))/2*H(dX)): ', f1

      !     f2 = mpddGlob%allRedMinReal(MINVAL(fY%ObsFields(1)%values))
      !     IF (mpddGlob%isBaseProc()) PRINT *, 'MINVAL((H(X+dX)-H(X-dX))/2*H(dX)): ', f2

      !   END BLOCK

      BLOCK

        rs1 = Y.DOT. (H%fwdNL_opr(X))
        rs2 = (H%transAdjMultiply_opr(Y, X)) .DOT.X

        IF (mpddGlob%isBaseProc()) THEN
          PRINT *, 'Y*(H(X)):', rs1
          PRINT *, 'H^T*Y*X:', rs2
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
