!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DP.Process_mwhs
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Yali Wu
! VERSION           : V 0.1
! HISTORY           : 2023-07-5, created by Yali Wu.
!
!   Created by Yali Wu (wuyali@gbamwf.com), 2023/07/5, @GBA-MWF, Shenzhen
!
!!--------------------------------------------------------------------------------------------------

!> @brief
!!
PROGRAM Process_mwhs
  USE geometry_m, ONLY: geometry_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE MPObs_m, ONLY: MPObs_t
  USE State_m, ONLY: State_t
  USE ObsSet_m, ONLY: ObsSet_t
  USE ObsField_m, ONLY: ObsField_t
  USE MGOpts_m
  USE State2NC_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE YAMLRead_m
  USE parameters_m
  USE Obs2State_m
  USE ReadWriteMWHS_m, ONLY: ReadWriteMWHS_t
  USE IOGrapes_m, ONLY: IOGrapes_t
  USE Ctl2State_m, ONLY: Ctl2State_t

  IMPLICIT NONE
  TYPE(mpddGlob_t), TARGET :: mpddGlob
  TYPE(geometry_t), TARGET :: geometry
  TYPE(State_t), ALLOCATABLE  :: XbMG(:)
  REAL(r_kind) :: t1, t2
  TYPE(MPObs_t), TARGET :: mpObs
  TYPE(IOGrapes_t), TARGET :: ioGrapes

  TYPE(State_t) :: XbRef
  TYPE(Ctl2State_t) :: Ctl2State
  CHARACTER(LEN=1024)   :: configFile
  TYPE(ReadWriteMWHS_t) :: ReadWriteMWHS
  CHARACTER(len=50), ALLOCATABLE :: platform_name(:), inst_name(:)
  INTEGER(i_kind) :: mgEnd, istatus, n_plats, n_insts, i_inst

  ! Get the configFile
  CALL getarg(1, configFile)
  IF (TRIM(configFile) .EQ. '') THEN
    PRINT *, TRIM(configFile)

    CALL getarg(0, configFile)
    PRINT *, TRIM(configFile)
    WRITE (*, 1)
1   FORMAT('Usage of this driver: mpirun -n <n> Debug/Process_mwhs.exe configFile', /, &
           ' Check your configure file and rerun!')
    configFile = './App_3DVarVerification.yaml'
  ELSE
    WRITE (*, 2) TRIM(configFile)
2   FORMAT('ConfigFile is: ', A)
  END IF

  istatus = yaml_get_var(TRIM(configFile), 'geometry', 'mgEnd', mgEnd)
  IF (yaml_get_var(TRIM(configFile), 'FY3-MWHS', 'satellite', platform_name) .NE. 0) STOP
  IF (yaml_get_var(TRIM(configFile), 'FY3-MWHS', 'instrument', inst_name) .NE. 0) STOP

  n_plats = SIZE(platform_name, 1)
  n_insts = SIZE(inst_name, 1)
  IF (n_plats .NE. n_insts) THEN
    PRINT *, 'STOP: Numbers of platform_name and inst_name are inconsistent'
    STOP
  END IF

  ! Initializer
  CALL mpddGlob%initialize()

  CALL geometry%initialize(configFile, mpddGlob)
  ALLOCATE (XbMG(geometry%mg%mg_coarsest:geometry%mg%mg_finest))

  ! Initialize the Scaling implementation
  CALL Ctl2State%initialize(configFile)

  ! Initialize a zeros background field at finest grid.
  ASSOCIATE (sg => geometry%mg%sg(mgEnd))

    CALL XbMG(mgEnd)%initialize(configFile, sg)
    CALL ioGrapes%initialize(configFile, geometry)
    CALL ioGrapes%m_read_bcg_into_Xm(XbMG(mgEnd), sg)

  END ASSOCIATE

  DO i_inst = 1, n_insts
    ! ! Read observations at the finest grid
    CALL ReadWriteMWHS%ObsInitial(configFile, i_inst)
    CALL ReadWriteMWHS%ObsIngest(XbMG(mgEnd))

    PRINT *, 'Num of OBS = ', ReadWriteMWHS%numObs
    CALL mpddGlob%barrier

    CALL Ctl2State%transBackward(XbMG(mgEnd))

    ASSOCIATE (sg => geometry%mg%sg(mgEnd))

      CALL mpObs%initializeMPObs(sg)    ! Initialize the observation parallel processing proc

      XbRef = Ctl2State%fwdNL_opr(XbMG(mgEnd))

      IF (XbRef%sg%isActiveProc()) THEN

        CALL ReadWriteMWHS%ObsPrepareForSg(XbRef)
        CALL ReadWriteMWHS%OutputForThinning(XbRef, mpObs)
        CALL ReadWriteMWHS%destructor()
        PRINT *, 'ReadWriteMWHS%OutputForThinning OVER'
        PRINT *, 'here ', ALLOCATED(ReadWriteMWHS%rttov_chan_lists), ALLOCATED(ReadWriteMWHS%ObsData)

      END IF

    END ASSOCIATE
  END DO

  ! CALL mpddGlob%barrier

  ! Finalize
  DEALLOCATE (XbMG)
  DEALLOCATE (platform_name, inst_name)
  CALL mpddGlob%finalize

END PROGRAM Process_mwhs
