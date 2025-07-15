!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2021/3/3, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!!
!! @author Zilong Qin
!! @copyright (C) 2020 GBA-MWF, All rights reserved.
PROGRAM testIOERA5
  USE geometry_m, ONLY: geometry_t
  USE IOGrapes_m, ONLY: IOGrapes_t
  USE IOERA5_m, ONLY: IOERA5_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE State_m, ONLY: State_t
  USE State2NC_m
  USE parameters_m, ONLY: spec_heat_const_pres, dry_air_gas_const
  USE YAMLRead_m

  ! Define types
  TYPE(mpddGlob_t), TARGET :: mpddGlob
  TYPE(geometry_t), TARGET :: geometry
  TYPE(IOGrapes_t), TARGET :: ioGrapes
  TYPE(IOERA5_t), TARGET :: ioERA5
  TYPE(State_t) :: Xm, XEra5
  REAL(r_kind) :: t1, t2

  CHARACTER(LEN=1024) :: configFile, outputDir

  ! Get the configFile
  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", configFile)
  configFile = "/Users/qzl/sources/MOTOR/input/2024110300_T10p1/App_3DVarVerification.yaml"

  !   IF(yaml_get_var(configFile, 'IO', 'output_dir', outputDir) /= 0) STOP
  outputDir = "/Users/qzl/sources/MOTOR/build"

  ! Initialize
  CALL mpddGlob%initialize()  ! Initialize the mpdd
  CALL geometry%initialize(configFile, mpddGlob)          ! Initialize the geometry

  ! Initialize the ioGrapes
  CALL ioGrapes%initialize(configFile, geometry)
  CALL ioERA5%initialize(configFile, geometry)

  ASSOCIATE (sg => geometry%mg%sg(geometry%mg%mg_finest))

    ! Initialize the Xm
    CALL Xm%initialize(configFile, sg)

    CALL mpddGlob%barrier
    CALL CPU_TIME(t1)

    ! Read the grapesinput into the Xm
    CALL ioGrapes%m_read_bcg_into_Xm(Xm, sg)

    XEra5 = Xm
    CALL ioERA5%m_read_bcg_into_Xm(XEra5, sg)

    CALL mpddGlob%barrier
    CALL CPU_TIME(t2)
    IF (mpddGlob%isBaseProc()) PRINT *, 'Time cost:', t2 - t1, 's'

    ! Xm%Fields(Xm%getVarIdx('pip'))%data = (Xm%Fields(Xm%getVarIdx('rho'))%data* &
    !                                        Xm%Fields(Xm%getVarIdx('temp'))%data*dry_air_gas_const/100000.0D0)** &
    !                                       (dry_air_gas_const/spec_heat_const_pres) - Xm%Fields(Xm%getVarIdx('pi'))%data

    ! Output Xm to NC files
    !  CALL Output_NC_State_AV(Xm, outputDir, 'testIOEra5')
    CALL Output_NC_State_AVST(Xm, outputDir, 'testIOEra5', sg%tSlots, .TRUE., .TRUE.)
    ! CALL Output_NC_State_AVST(XEra5 - Xm, outputDir, 'testDiff', sg%tSlots, .TRUE., .TRUE.)

  END ASSOCIATE

  ! Destroy the structures
  CALL mpddGlob%finalize

END PROGRAM
