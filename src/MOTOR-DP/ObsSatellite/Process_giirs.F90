!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DP.Process_giirs
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
PROGRAM Process_giirs
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
  USE ReadWriteGIIRS_m,  ONLY: ReadWriteGIIRS_t
  USE IOGrapes_m, ONLY: IOGrapes_t
  USE Ctl2State_m, ONLY: Ctl2State_t

  IMPLICIT NONE
  TYPE(mpddGlob_t), TARGET :: mpddGlob
  TYPE(geometry_t), TARGET :: geometry
  TYPE(State_t), ALLOCATABLE  :: XbMG(:)
  REAL(r_kind) :: t1, t2
  TYPE(MPObs_t), TARGET :: mpObs
  TYPE(IOGrapes_t), target :: ioGrapes

  TYPE(State_t) :: XbRef
  TYPE(Ctl2State_t) :: Ctl2State
  CHARACTER(LEN=1024)   :: configFile
  TYPE(ReadWriteGIIRS_t) :: ReadWriteGIIRS
  INTEGER(i_kind) :: mgEnd, istatus

  CALL cpu_time(t1)

  ! Get the configFile
  CALL getarg(1, configFile)
  IF (TRIM(configFile) .EQ. '') THEN
    PRINT*, TRIM(configFile)

    CALL getarg(0, configFile)
    PRINT*, TRIM(configFile)
    WRITE (*, 1)
1   FORMAT('Usage of this driver: mpirun -n <n> Debug/Process_giirs.exe configFile', /, &
           ' Check your configure file and rerun!')
    configFile = './App_3DVarVerification.yaml'
  ELSE
    WRITE (*, 2) TRIM(configFile)
2   FORMAT('ConfigFile is: ', A)
  END IF

  istatus = yaml_get_var(TRIM(configFile), 'geometry', 'mgEnd', mgEnd)

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

  ! ! Read observations at the finest grid
  CALL ReadWriteGIIRS%ObsInitial(configFile)
  CALL ReadWriteGIIRS%ObsIngest(XbMG(mgEnd))

  PRINT *, 'Num of OBS = ', ReadWriteGIIRS%numObs
  CALL mpddGlob%barrier

  CALL Ctl2State%transBackward(XbMG(mgEnd))

  ASSOCIATE (sg => geometry%mg%sg(mgEnd))
    
    CALL mpObs%initializeMPObs(sg)    ! Initialize the observation parallel processing proc   

    XbRef = Ctl2State%fwdNL_opr(XbMG(mgEnd))

    IF (XbRef%sg%isActiveProc()) THEN
      
      CALL ReadWriteGIIRS%ObsPrepareForSg(XbRef)  
      CALL ReadWriteGIIRS%OutputForThinning(XbRef, mpObs)
      CALL ReadWriteGIIRS%destructor()
      PRINT *, 'ReadWriteGIIRS%OutputForThinning OVER'
      PRINT *, 'here ', ALLOCATED(ReadWriteGIIRS%rttov_chan_lists), ALLOCATED(ReadWriteGIIRS%ObsData)
    
    END IF

  END ASSOCIATE
  
  ! CALL mpddGlob%barrier
  CALL cpu_time(t2)
  PRINT *, 'Time cost:', t2 - t1, 's'

  ! Finalize
  DEALLOCATE (XbMG)
  CALL mpddGlob%finalize

END PROGRAM Process_giirs
