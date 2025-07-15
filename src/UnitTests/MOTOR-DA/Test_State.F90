!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA.State_Xm
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2020/12/31, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!!
!! @copyright (C) 2020 GBA-MWF, All rights reserved.
PROGRAM Test_State
  USE State_m, ONLY: State_t
  USE geometry_m, ONLY: geometry_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE kinds_m, ONLY: i_kind, r_kind

  TYPE(mpddGlob_t), TARGET :: mpddGlob
  TYPE(geometry_t), TARGET :: geometry
  TYPE(State_t), TARGET :: Xm

  CHARACTER(LEN=1024) :: configFile
  REAL(r_kind) :: J

  ! Get the configFile
  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", configFile)
  configFile = TRIM(configFile)//"/UnitTest/testState.yaml"

  ! Initialize
  CALL mpddGlob%initialize()
  CALL geometry%initialize(configFile, mpddGlob)          ! Initialize the geometry
  CALL Xm%initialize(configFile, geometry%mg%sg(geometry%mg%mg_finest))

  Xm%fields(1)%DATA = 1.0D0
  Xm%fields(2)%DATA = 0.0D0
  Xm%fields(3)%DATA = 0.0D0
  Xm%fields(4)%DATA = 0.0D0
  Xm%fields(5)%DATA = 0.0D0

  ! Xm = Xm+Xm
  ! J = XX*XX
  J = (Xm.DOT.Xm)
  PRINT *, 'Programe is done! J is ', J

  ! For ctest
  IF (mpddGlob%isBaseProc()) THEN
    IF (J .EQ. 69476616) THEN
      PRINT *, 'Test passed!'
    ELSE
      PRINT *, 'Test failed!'
    END IF
  END IF

  Xm = Xm * 2.0D0 - Xm
  PRINT *, 'MAX(Xm*2.0D0 - Xm)', MAXVAL(Xm%Fields(1)%DATA)

  Xm = Xm / 2.0D0
  PRINT *, 'MAX(Xm/2.0D0)', MAXVAL(Xm%Fields(1)%DATA)

  Xm = Xm / (Xm / 2.0D0)
  PRINT *, 'Xm/(Xm/2.0d0)', MAXVAL(Xm%Fields(1)%DATA)

  ! Destroy the structures
  CALL mpddGlob%finalize

END PROGRAM Test_State
