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
PROGRAM testIOGrapes
  USE geometry_m, ONLY: geometry_t
  USE IOGrapes_m, ONLY: IOGrapes_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE State_m, ONLY: State_t
  USE State2NC_m
  USE parameters_m, ONLY: spec_heat_const_pres, dry_air_gas_const
  USE YAMLRead_m

  ! Define types
  TYPE(mpddGlob_t), TARGET :: mpddGlob
  TYPE(geometry_t), TARGET :: geometry
  TYPE(IOGrapes_t), TARGET :: ioGrapes
  TYPE(State_t) :: Xm
  REAL(r_kind) :: t1, t2

  CHARACTER(LEN=1024) :: configFile, outputDir

  ! Get the configFile
  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", configFile)
  ! configFile = TRIM(configFile)//"/UnitTest/Test_IOGrapes.yaml"
  ! configFile = TRIM(configFile)//"/UnitTest/Test_IOGrapes.yaml"
  configFile = TRIM(configFile)//"/UnitTest/Test_IOError.yaml"
  configFile = "/public/home/simi/optest/3DVarVerification/srf/220527/220527_0000_3km/App_3DVarVerification.yaml"

  !   IF(yaml_get_var(configFile, 'IO', 'output_dir', outputDir) /= 0) STOP

  ! Initialize
  CALL mpddGlob%initialize()  ! Initialize the mpdd

  CALL geometry%initialize(configFile, mpddGlob)          ! Initialize the geometry

  ! Initialize the ioGrapes
  CALL ioGrapes%initialize(configFile, geometry)

  ASSOCIATE (sg => geometry%mg%sg(7))

    ! Initialize the Xm
    CALL Xm%initialize(configFile, sg)

    CALL mpddGlob%barrier
    CALL CPU_TIME(t1)

    ! Read the grapesinput into the Xm
    CALL ioGrapes%m_read_bcg_into_Xm(Xm, sg)
    ! Xm%Fields(Xm%getVarIdx('uwnd')) %data = 0.0D0
    ! Xm%Fields(Xm%getVarIdx('vwnd')) %data = 0.0D0
    ! Xm%Fields(Xm%getVarIdx('temp')) %data = 0.0D0
    ! Xm%Fields(Xm%getVarIdx('rhov')) %data = 0.0D0
    ! Xm%Fields(Xm%getVarIdx('rho')) %data = 0.0D0

    CALL ioGrapes%m_write_Xm_into_bcg(Xm, sg, Xm)

    CALL mpddGlob%barrier
    CALL CPU_TIME(t2)
    IF (mpddGlob%isBaseProc()) PRINT *, 'Time cost:', t2 - t1, 's'

    ! Xm%Fields(Xm%getVarIdx('pip'))%data = (Xm%Fields(Xm%getVarIdx('rho'))%data* &
    !                                        Xm%Fields(Xm%getVarIdx('temp'))%data*dry_air_gas_const/100000.0D0)** &
    !                                       (dry_air_gas_const/spec_heat_const_pres) - Xm%Fields(Xm%getVarIdx('pi'))%data

    ! Output Xm to NC files
    CALL Output_NC_State_AV(Xm, outputDir, 'testIOGrapes')

  END ASSOCIATE

  ! Destroy the structures
  CALL mpddGlob%finalize

END PROGRAM testIOGrapes
