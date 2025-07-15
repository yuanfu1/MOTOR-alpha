!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2022/3/18, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!!
PROGRAM Test_ObsRadarVel
  USE kinds_m, ONLY: i_kind, r_kind, r_double
  USE RadarRAW_m, ONLY: RadarRAW_t
  USE ObsRadarVel_m, ONLY: ObsRadarVel_t
  USE geometry_m, ONLY: geometry_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE MPObs_m, ONLY: MPObs_t
  USE State_m, ONLY: State_t
  USE ObsSet_m, ONLY: ObsSet_t
  USE C2O_m, ONLY: C2O_t
  USE YAMLRead_m
  USE State2NC_m
  USE Obs2State_m
  IMPLICIT NONE

  ! Define types
  TYPE(mpddGlob_t), TARGET :: mpddGlob
  TYPE(geometry_t), TARGET :: geometry
  REAL(r_kind) :: rs1, rs2, f1, f2

  CHARACTER(LEN=1024) :: inputDir
  CHARACTER(LEN=1024) :: outputDir
  CHARACTER(len=1024) :: filenameInput, filenameOutput
  CHARACTER(len=1024) :: configFile
  CHARACTER(LEN=1024) :: staticDir
  TYPE(ObsRadarVel_t) ::  ObsRadarVel

  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", staticDir)
  configFile = TRIM(staticDir)//"/UnitTest"//"/Test_VelUVW2RwndTLAD.yaml"

  IF (yaml_get_var(configFile, 'IO', 'output_dir', outputDir) /= 0) STOP

  ! Initializer
  ! Auxtypes
  CALL mpddGlob%initialize()                                   ! Initialize the mpdd

  ! Initialize geometry
  CALL geometry%initialize(configFile, mpddGlob)               ! Initialize the geometry

  ASSOCIATE (sg => geometry%mg%sg(9))
    BLOCK

      TYPE(MPObs_t), TARGET :: mpObs
      TYPE(State_t) :: X
      TYPE(ObsSet_t) :: Y
      TYPE(C2O_t) :: H
      INTEGER(i_kind) :: i, j, k

      CALL mpObs%initializeMPObs(sg)
      CALL X%initialize(configFile, sg)
      CALL X%setAllFieldData(1.0D0)

      Y = ObsSet_t(configFile, mpObs)

      ! Thinning the radial wind
      ! Initial
      CALL ObsRadarVel%ObsInitial(configFile)
      CALL ObsRadarVel%ObsIngest(X)
      CALL ObsRadarVel%ObsPrepareForSg(X)
      CALL ObsRadarVel%ObsSuper(X, Y, mpObs)

      CALL H%initialize(configFile, X, Y)

      ! Test 1:  TL test
      BLOCK
        TYPE(State_t) :: dX
        TYPE(ObsSet_t) :: fY

        CALL X%initialize(configFile, sg)
        CALL dX%setAllFieldData(0.1D0)

        fY = ((H%fwdNL_opr(X + dX)) - (H%fwdNL_opr(X - dX))) / ((H%fwdNL_opr(dX)) * 2.0D0)

        f1 = mpddGlob%allRedMaxReal(MAXVAL(fY%ObsFields(1)%values))
        IF (mpddGlob%isBaseProc()) PRINT *, 'MAXVAL((H(X+dX)-H(X-dX))/2*H(dX)): ', f1

        f2 = mpddGlob%allRedMinReal(MINVAL(fY%ObsFields(1)%values))
        IF (mpddGlob%isBaseProc()) PRINT *, 'MINVAL((H(X+dX)-H(X-dX))/2*H(dX)): ', f2

      END BLOCK

      BLOCK
        ! Test 2:  AD test
        rs1 = Y.DOT. (H%fwdNL_opr(X))
        rs2 = (H%adjMul_opr(Y, X)) .DOT.X

        IF (mpddGlob%isBaseProc()) THEN
          PRINT *, 'Y*(H(X)):', rs1
          PRINT *, 'H^T*Y*X:', rs2
        END IF
      END BLOCK

    END BLOCK
  END ASSOCIATE

  IF (mpddGlob%isBaseProc()) THEN
    IF (ABS(f1 - 1.0) < 1E-7 .AND. ABS(f2 - 1.0) < 1E-7 .AND. ABS(rs1 - rs2) < 1E-7) THEN
      PRINT *, 'Test passed!'
    ELSE
      PRINT *, 'Test failed!'
    END IF
  END IF

  ! Destory
  CALL ObsRadarVel%ObsDestroy()
  ! PRINT *, 'It is done'

  CALL mpddGlob%finalize
END PROGRAM Test_ObsRadarVel
