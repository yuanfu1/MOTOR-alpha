!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA.ObsSet
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2021/1/26, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!!
! @note
! @warning
! @attention
PROGRAM Test_ObsSet
  USE geometry_m, ONLY: geometry_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE MPObs_m, ONLY: MPObs_t
  USE kinds_m, ONLY: i_kind, r_kind
  USE ObsSet_m, ONLY: ObsSet_t
  USE Mock_m

  TYPE(mpddGlob_t), TARGET :: mpddGlob
  TYPE(geometry_t), TARGET :: geometry
  TYPE(MPObs_t), TARGET :: mpObs
  TYPE(ObsSet_t) :: Y1, Y2

  CHARACTER(LEN=1024) :: configFile
  REAL(r_kind) :: J1, J2, J3

  ! Get the configFile
  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", configFile)
  configFile = TRIM(configFile)//"/UnitTest/testObsSet.yaml"

  ! Initialize
  CALL mpddGlob%initialize()
  CALL geometry%initialize(configFile, mpddGlob)          ! Initialize the geometry
  CALL mpObs%initializeMPObs(geometry%mg%sg(geometry%mg%mg_finest)) ! Initialize the observation parallel processing proc

  ! Mock observation data:
  CALL Set_Mock_Single_Pt_For_ObsSet(configFile, Y1, mpObs, 1.0D0, geometry%mg%sg(geometry%mg%mg_finest))

  J1 = Y1.dot.Y1
  PRINT *, 'J1 is: ', J1

! Mock observation data:
  CALL Set_Mock_Single_Pt_For_ObsSet(configFile, Y2, mpObs, 4.0D0, geometry%mg%sg(geometry%mg%mg_finest))

  J2 = Y2.dot.Y2
  PRINT *, 'J2 is: ', J2

  Y2 = Y2 - Y1
  J3 = Y2.dot.Y2
  PRINT *, 'J3 is: '

  IF ((J1 .EQ. 4) .AND. (J2 .EQ. 64) .AND. (J3 .EQ. 36)) THEN
    PRINT *, 'Test passed!'
  ELSE
    PRINT *, 'Test failed!'
  END IF

  ! Y%Fields(1)%para=ObsParaSAT_t(configFile)
  ! Y%Fields(1)%para=ObsParaSFC_t(configFile)

  ! Destroy the structures
  CALL mpddGlob%finalize
END PROGRAM Test_ObsSet
