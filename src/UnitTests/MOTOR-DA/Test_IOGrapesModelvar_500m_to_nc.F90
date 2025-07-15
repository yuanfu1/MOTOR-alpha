!----------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Centers
!                     for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Jiongming Pang
! VERSION           : V 0.0
! HISTORY           : 2023-07-01, created by Jiongming Pang for testing IOGrapesModelvar.
!
!   Created by Jiongming Pang (pang.j.m@hotmail.com), 2023-07-01, @GBA-MWF, Shenzhen
!
!----------------------------------------------------------------------------------------

PROGRAM Test_IOGrapesModelvar_500m_to_nc
    USE kinds_m
    USE State2NC_m
    USE YAMLRead_m
    USE InterpolateTime_m
    USE geometry_m, ONLY: geometry_t
    USE mpddGlob_m, ONLY: mpddGlob_t
    USE State_m, ONLY: State_t
    USE Export2SelDomain_m, ONLY: Export2SelDomain_t
    USE IOGrapesModelvar_500m_m, ONLY: IOGrapesModelvar_500m_t
    USE Export2HASCoordInSelDomain_m, ONLY: Export2HASCoordInSelDomain_t

    IMPLICIT NONE
    ! Define types
    TYPE(State_t)  :: Xm
    TYPE(mpddGlob_t), TARGET :: mpddGlob
    TYPE(geometry_t), TARGET :: geometry
    TYPE(Export2SelDomain_t) ::  Export2SelDomain
    TYPE(IOGrapesModelvar_500m_t), TARGET :: ioGrapesModelvar_500m
    TYPE(Export2HASCoordInSelDomain_t) ::  Export2HASCoordInSelDomain
  
    CHARACTER(LEN=1024) :: configFile, bakFilesDir, outputDir, arg_buffer, ncOutputFile
    CHARACTER(LEN=20), ALLOCATABLE :: varList(:)
    CHARACTER(LEN=20) :: varName, task
    CHARACTER(LEN=3)  :: memid
    INTEGER(i_kind)   :: varNum, iv, ifile,  arg_length, arg_status, i
  
    ! Get the configFile
    CALL GET_COMMAND_ARGUMENT(1, arg_buffer, arg_length, arg_status)
    IF (arg_status == 0 .AND. arg_length > 0) THEN
    configFile = TRIM(arg_buffer)
    ELSE
    WRITE(*,*) "error: no configFile specified"
    STOP 1
    END IF
  
    ifile = yaml_get_var(TRIM(configFile), 'IO', 'output_dir', outputDir)
    ! ifile = yaml_get_var(TRIM(configFile), 'RunMode', 'Task', task)
    ifile = yaml_get_var(TRIM(configFile), 'modelState', 'varList', varList)
    varNum = UBOUND(varList, 1)
    PRINT *, "num of analysis vars: ", varNum
  
    ! Initialize
    CALL mpddGlob%initialize()   ! Initialize the mpdd
    CALL geometry%initialize(configFile, mpddGlob)   ! Initialize the geometry
  
    ! Initialize the ioGrapesModelvar_500m
    CALL ioGrapesModelvar_500m%initialize(configFile, geometry)
  
    ASSOCIATE (sg => geometry%mg%sg(geometry%mg%mg_finest))
  
      ! Initialize the Xm
      CALL Xm%initialize(configFile, sg)
  
      ! Read the modelvar into the Xm
      CALL ioGrapesModelvar_500m%m_read_bcg_into_Xm(Xm, sg)
  
      ! Output Xm to NC files
      IF (Xm%sg%vLevel == 1) THEN
        ncOutputFile = outputDir // "Modelvar_500m_surfaces.nc"
        task = "Modelvar_500m_surfaces"
       
        Export2SelDomain = Export2SelDomain_t(configFile, sg)
        CALL Export2SelDomain%Export_State_AVST(Xm, TRIM(outputDir),  &
                                                  TRIM(task), 1, .FALSE., .FALSE.)

        ! Export2HASCoordInSelDomain = Export2HASCoordInSelDomain_t(configFile, sg, Xm)
        ! CALL Export2HASCoordInSelDomain%Export_State_AVLNST(Xm, ncOutputFile, &
        !                                                         TRIM(task), 1, .FALSE., .FALSE.)
        ! DO iv = 1, varNum
        !   varName = varList(iv)
        !   CALL Output_NC_State_AV(Xm, TRIM(outputDir), "Modelvar_500m_surfaces")
        ! END DO
      ELSE
        ncOutputFile = outputDir // "Modelvar_500m_levels.nc"
        task = "Modelvar_500m_levels"
        Export2HASCoordInSelDomain = Export2HASCoordInSelDomain_t(configFile, sg, Xm)
        CALL Export2HASCoordInSelDomain%Export_State_AVLNST(Xm, ncOutputFile, &
                                                                TRIM(task), 1, .FALSE., .FALSE.)
      END IF
     
      
    END ASSOCIATE
  
    CALL mpddGlob%finalize
  
  END PROGRAM Test_IOGrapesModelvar_500m_to_nc
  