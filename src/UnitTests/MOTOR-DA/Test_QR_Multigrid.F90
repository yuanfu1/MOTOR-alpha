!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA.EnLoc.Test_EnLoc
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Jiongming Pang
! VERSION           : V 0.0   Calculating QR in multigrid and save as temporary files.
! HISTORY           :
!   Created by Jiongming Pang (pang.j.m@hotmail.com), 2022/4/5, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!!
PROGRAM Test_QR_Multigrid

  USE kinds_m
  USE geometry_m, ONLY: geometry_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE State_m, ONLY: State_t
  USE IOWRF_m, ONLY: IOWRF_t
  USE EnLoc_m, ONLY: EnLoc_t
  USE State2NC_m
  !USE NMLRead_m
  USE YAMLRead_m
  USE MGOpts_m

  IMPLICIT NONE

  ! Define types
  TYPE(mpddGlob_t), TARGET :: mpddGlob
  TYPE(geometry_t), TARGET :: geometry

  TYPE(IOWRF_t), ALLOCATABLE, TARGET :: ioWRF
  ! TYPE(IOWRF_t),  target :: ioWRF
  TYPE(State_t), ALLOCATABLE         :: Xm(:)
  TYPE(State_t), ALLOCATABLE         :: XbMG(:, :)
  ! TYPE(EnLoc_t), ALLOCATABLE         :: EnLoc
  TYPE(EnLoc_t)       :: EnLoc
  CHARACTER(LEN=1024) :: configFile, bakFilesDir, outputDir
  CHARACTER(LEN=20), ALLOCATABLE :: varList(:)
  CHARACTER(LEN=20)   :: varName
  CHARACTER(LEN=3)    :: memid, gid
  INTEGER(i_kind)     :: i, iv, varNum, ensNum, iens, mgStart, mgEnd, ifile

  ! Get the configFile
  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", configFile)
  configFile = TRIM(configFile)//"/UnitTest/testQR_MulitiGrid.yaml"

  IF (yaml_get_var(configFile, 'IO', 'output_dir', outputDir) /= 0) STOP

  !CALL namelist_read(TRIM(configFile), "varList", varList)
  ifile = yaml_get_var(TRIM(configFile), 'modelState', 'varList', varList)
  varNum = UBOUND(varList, 1)
  PRINT *, "num of analysis vars: ", varNum

  !CALL namelist_read(TRIM(configFile), "ensNum", ensNum)
  ifile = yaml_get_var(TRIM(configFile), 'BMat', 'ensNum', ensNum)
  PRINT *, "num of ensembles: ", ensNum

  ! Initialize
  CALL mpddGlob%initialize()   ! Initialize the mpdd
  CALL geometry%initialize(configFile, mpddGlob)   ! Initialize the geometry

  ! ALLOCATE(Xm(ensNum))
  ALLOCATE (XbMG(geometry%mg%mg_coarsest:geometry%mg%mg_finest, ensNum))

  !CALL namelist_read(TRIM(configFile), "mgStart", mgStart)
  !CALL namelist_read(TRIM(configFile), "mgEnd", mgEnd)
  ifile = yaml_get_var(TRIM(configFile), 'geometry', 'mgStart', mgStart)
  ifile = yaml_get_var(TRIM(configFile), 'geometry', 'mgEnd', mgEnd)
  PRINT *, "mgStart and mgEnd:", mgStart, mgEnd

  DO i = mgEnd, mgStart, -1
    PRINT *, 'i', i
    IF (i .EQ. mgEnd) THEN
      ! Initialize a zeros background field at finest grid.
      ASSOCIATE (sg => geometry%mg%sg(mgEnd))
        ! Initialize the ioWRF
        ALLOCATE (ioWRF)
        CALL ioWRF%initialize(configFile, geometry)

        CALL XbMG(mgEnd, 1)%Initialize(configFile, sg)
        IF (ensNum > 1) THEN
          DO iens = 2, ensNum
            XbMG(mgEnd, iens) = XbMG(mgEnd, 1)
          END DO
        END IF
        DO iens = 1, ensNum
          CALL ioWRF%m_read_bcg_into_Xm_Ens(XbMG(mgEnd, iens), sg, iens)   ! read wrfout into XbMG
        END DO
        DEALLOCATE (ioWRF)
      END ASSOCIATE
    ELSE
      ! Restrict to each coarser grid
      ASSOCIATE (sgFiner => geometry%mg%sg(i + 1), sgCoarser => geometry%mg%sg(i))
        CALL XbMG(i, 1)%Initialize(configFile, sgCoarser)
        IF (ensNum > 1) THEN
          DO iens = 2, ensNum
            XbMG(i, iens) = XbMG(i, 1)
          END DO
        END IF
        DO iens = 1, ensNum
          CALL restrictionMG(XbMG(i, iens), XbMG(i + 1, iens), geometry%mg)
        END DO
      END ASSOCIATE
    END IF

    ! ! Output Xm to NC files
    ! write(gid, '(A1,I2.2)') "G", i
    ! DO iens=1,ensNum
    !   write(memid, '(A1,I2.2)') "m", iens
    !   DO iv =1,varNum
    !     varName = varList(iv)
    !     CALL Output_NC_State_SV(XbMG(i,iens), TRIM(outputDir), &
    !                             "testEnLoc_Ens_bak_"//TRIM(varName)//"_"//memid//"_"//gid, TRIM(varName))
    !   END DO
    ! END DO

  END DO
  PRINT *, "Done ensembles interpolating to multi-grid."

  IF (mpddGlob%isBaseProc()) THEN
    DO i = mgStart, mgEnd
      ! ALLOCATE(EnLoc)
      CALL EnLoc%b_getEnsdata(configFile, XbMG(i, :))
      PRINT *, 'Done EnLoc%b_getEnsdata in grid', i

      CALL EnLoc%b_qrDecomp(i)
      PRINT *, 'Done EnLoc%b_qrDecomp in grid', i

      ! CALL EnLoc%b_jacoEigenT(i)
      ! PRINT *, 'Done EnLoc%b_jacoEigenT in grid', i

      CALL EnLoc%b_destroy()
      PRINT *, 'Done EnLoc%b_destroy in grid', i
      ! DEALLOCATE(EnLoc)
    END DO
  END IF

  DEALLOCATE (XbMG)

  ! Destroy the structures
  CALL mpddGlob%finalize

END PROGRAM Test_QR_Multigrid
