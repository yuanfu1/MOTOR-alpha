!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2021/1/26, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!!
PROGRAM Test_C2O
  USE geometry_m, ONLY: geometry_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE MPObs_m, ONLY: MPObs_t
  USE State_m, ONLY: State_t
  USE ObsSet_m, ONLY: ObsSet_t
  USE C2O_m, ONLY: C2O_t
  USE State2NC_m
  USE Mock_m
  USE RossbyHaurwitzSphere_m, ONLY: RossbyHaurwitzSphere_t ! rossbyhaurwitz function from utility
  USE YAMLRead_m

  ! Define types
  TYPE(mpddGlob_t), TARGET :: mpddGlob
  TYPE(geometry_t), TARGET :: geometry
  TYPE(MPObs_t), TARGET :: mpObs
  TYPE(State_t) :: Xm
  TYPE(ObsSet_t) :: Y
  TYPE(C2O_t) :: H
  TYPE(RossbyHaurwitzSphere_t) :: rossby
  CHARACTER(LEN=1024) :: configFile, outputDir
  INTEGER(i_kind) :: fileStat
  REAL(r_kind), ALLOCATABLE :: rhs(:,:)

  ! Get the configFile
  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", configFile)
  configFile = TRIM(configFile)//"/UnitTest/testC2O.yaml"

  fileStat = yaml_get_var(configFile, 'IO', 'output_dir', outputDir)

  ! Initializer
  CALL mpddGlob%initialize()                              ! Initialize the mpdd
  CALL geometry%initialize(configFile, mpddGlob)          ! Initialize the geometry
  CALL mpObs%initializeMPObs(geometry%mg%sg(geometry%mg%mg_finest)) ! Initialize the observation parallel processing proc
  CALL Xm%initialize(configFile, geometry%mg%sg(geometry%mg%mg_finest))

  ! Main code
  CALL rossby%initializ(Xm%sg%num_icell) ! Beaware that this initialization determine the cache size for generating the corresponding values

  ! Xm%fields(1)%data(1, 125000) = 100D0
  ALLOCATE(rhs(Xm%sg%num_icell,3))
  CALL rossby%UpdateVar(Xm%sg%num_icell, Xm%sg%cell_cntr, -600.0D0, Xm%fields(1)%DATA(1, 2, :), &
    Xm%fields(1)%DATA(1, 2, :),Xm%fields(3)%DATA(1,2,:),rhs)
  DEALLOCATE(rhs) ! Unused rhs
  Xm%fields(1)%DATA(1, 1, :) = Xm%fields(1)%DATA(1, 2, :)

  ! Mock observation data:
  CALL Set_Mock_Single_Pt(configFile, Y, mpObs, 0.0D0, geometry%mg%sg(geometry%mg%mg_finest), 'ua')
  CALL H%initialize(configFile, Xm, Y)

  ! Output state file to NC for view.
  CALL Output_NC_State_SV(Xm, outputDir, 'test', "ua")

  ! Calculate H*X
  Y = H%fwdNL_opr(Xm)

  ! Finalize
  CALL mpddGlob%finalize

END PROGRAM Test_C2O
