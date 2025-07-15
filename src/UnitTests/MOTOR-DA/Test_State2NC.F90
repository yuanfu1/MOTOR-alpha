!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA.State2NC
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
PROGRAM Test_State2NC
  USE geometry_m, ONLY: geometry_t
  ! USE IOGrapes_m, ONLY: IOGrapes_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE State_m, ONLY: State_t
  USE State2NC_m
  USE YAMLRead_m

  ! Define types
  TYPE(mpddGlob_t), TARGET :: mpddGlob
  TYPE(geometry_t), TARGET :: geometry
  ! TYPE(IOGrapes_t), TARGET :: ioGrapes
  TYPE(State_t) :: Xm
  INTEGER(i_kind) :: fileStat

  CHARACTER(LEN=1024) :: configFile, outputDir

  ! Get the configFile
  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", configFile)
  configFile = TRIM(configFile)//"/UnitTest/testState2NC.yaml"

  fileStat = yaml_get_var(configFile, 'IO', 'output_dir', outputDir)

  ! Initialize
  CALL mpddGlob%initialize()  ! Initialize the mpdd
  CALL geometry%initialize(configFile, mpddGlob)          ! Initialize the geometry

  ASSOCIATE (sg => geometry%mg%sg(8))
    ! Initialize the Xm
    CALL Xm%initialize(configFile, sg)
    Xm%fields(1)%DATA(1, :, 1) = 1.0D0
    Xm%fields(1)%DATA(2, :, 1) = 2.0D0

    ! Output Xm to NC files
    CALL Output_NC_State_AV(Xm, outputDir, 'testState2NC')
    CALL Output_NC_State_SVLNST(Xm, outputDir, 'testState2NC_SVLNST_uwnd_', 'uwnd', 4)
    CALL Output_NC_State_AVLNST(Xm, outputDir, 'testState2NC_AVLNST', 4)

  END ASSOCIATE

  ! Destroy the structures
  CALL mpddGlob%finalize

END PROGRAM Test_State2NC
