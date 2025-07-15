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
PROGRAM Test_UV2WForward
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
  USE UV2W_SurfInte_m, ONLY: UV2W_SurfInte_t
  USE UVW2Divg_m, ONLY: UVW2Divg_t

  IMPLICIT NONE

  ! Define types
  TYPE(mpddGlob_t), TARGET :: mpddGlob
  TYPE(geometry_t), TARGET :: geometry
  REAL(r_kind) :: rs1, rs2, f1, f2
  TYPE(IOGrapes_t), TARGET :: ioGrapes

  CHARACTER(len=1024) :: configFile
  CHARACTER(LEN=1024) :: staticDir
  CHARACTER(LEN=1024) :: outputDir

  INTEGER(i_kind):: fileStat, i

  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", staticDir)
  configFile = TRIM(staticDir)//"/UnitTest"//"/Test_UV2WTLAD.yaml"
  configFile = "/public/home/simi/optest/Case0908/2022090809/App_3DVarVerificationForTest.yaml"
  configFile = "/Users/qzl/sources/MOTOR/input/220527_0000_3km/App_UV2W.yaml"

  PRINT *, 'configFile: ', configFile
  fileStat = yaml_get_var(configFile, 'IO', 'output_dir', outputDir)

  ! Initializer
  ! Auxtypes
  CALL mpddGlob%initialize()                                   ! Initialize the mpdd

  ! Initialize geometry
  CALL geometry%initialize(configFile, mpddGlob)               ! Initialize the geometry
  ! Initial

  CALL ioGrapes%initialize(configFile, geometry)

  CALL mpddGlob%barrier
  ASSOCIATE (sg => geometry%mg%sg(7))
    BLOCK
      TYPE(MPObs_t), TARGET :: mpObs
      TYPE(State_t) :: X, Y, tempX
      CLASS(UV2W_t), ALLOCATABLE :: H
      TYPE(UVW2Divg_t) :: F
      LOGICAL, DIMENSION(sg%num_cell) :: mask
      INTEGER(i_kind) :: i, j, k
      INTEGER(i_kind), ALLOCATABLE :: idxSelected(:)

      CALL X%initialize(configFile, sg)
      CALL ioGrapes%m_read_bcg_into_Xm(X, sg)

      H = UV2W_SurfInte_t(configFile, X)
      F = UVW2Divg_t(configFile, X)
      PRINT *, 'End initialize'

      CALL H%fwdNL(X)
      CALL X%addVar('wwnd')
      CALL F%fwdNL(X)

      mask = (sg%cell_type == 0)

      idxSelected = PACK((/(i, i=1, sg%num_cell)/), mask)

      PRINT *, 'Output: ', SUM(ABS(X%fields(X%getVarIdx('divg'))%DATA(2:X%sg%vLevel - 1, idxSelected, 1))) &
        / ((X%sg%vLevel - 2) * X%sg%num_icell), MAXVAL(X%fields(X%getVarIdx('divg'))%DATA(2:X%sg%vLevel - 1, idxSelected, 1)), &
        MINVAL(X%fields(X%getVarIdx('divg'))%DATA(2:X%sg%vLevel - 1, idxSelected, 1)), mpddGlob%myrank

      BLOCK
        REAL(r_kind) :: sumEach, sumAll

        sumEach = SUM(X%fields(X%getVarIdx('wwnd'))%DATA(2:X%sg%vLevel - 1, idxSelected, 1))
        CALL mpddGlob%AllReduceSumReal(sumEach, sumAll)
        PRINT *, 'All sum is: ', sumAll

      END BLOCK

      CALL Output_NC_State_AV(X, outputDir, "Test_UV2WForward_SurfInte")

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
END PROGRAM Test_UV2WForward
