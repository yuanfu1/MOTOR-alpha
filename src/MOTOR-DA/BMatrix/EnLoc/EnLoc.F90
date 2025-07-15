!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA.EnLoc.EnLoc
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Jiongming Pang, Yuanfu Xie
! VERSION           : V 0.0   Calling BKErr_m and restore and resultes.
! HISTORY           :
!   Created by Jiongming Pang (pang.j.m@hotmail.com), 2022/3/29, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!!
MODULE EnLoc_m

  USE kinds_m
  USE BKErr_m
  USE State_m, ONLY: State_t
  USE GetEns_m, ONLY: Ens_t
  USE parameters_m, ONLY: degree2radian, radian2degree
  USE geoTools_m
  USE ncWriteGrid_m
  USE WriteVar_m
  USE YAMLRead_m

  TYPE EnLoc_t
    TYPE(Ens_t) :: Ens
    TYPE(BKErr_t), ALLOCATABLE :: BKErr(:)

  CONTAINS
    PROCEDURE, PUBLIC :: b_getEnsData => sub_b_getEnsData
    PROCEDURE, PUBLIC :: b_qrDecomp => sub_b_qrDecomp
    PROCEDURE, PUBLIC :: b_qrSave => sub_b_qrSave
    PROCEDURE, PUBLIC :: b_qrLoad => sub_b_qrLoad
    PROCEDURE, PUBLIC :: b_qrRinv => sub_b_qrRinv
    PROCEDURE, PUBLIC :: b_FWD => sub_b_FWD
    PROCEDURE, PUBLIC :: b_ADJ => sub_b_ADJ
    PROCEDURE, PUBLIC :: b_destroy => sub_b_destroy
  END TYPE EnLoc_t

  PRIVATE :: sub_b_getEnsData, sub_b_qrDecomp, sub_b_qrSave, sub_b_qrRinv
  PRIVATE :: sub_b_qrLoad, sub_b_FWD, sub_b_ADJ, sub_b_destroy


CONTAINS


  SUBROUTINE sub_b_getEnsData(this, configFile, Xm)
    IMPLICIT NONE
    CLASS(EnLoc_t) :: this
    CHARACTER(LEN=1024), INTENT(IN)     :: configFile
    TYPE(State_t), INTENT(IN), OPTIONAL :: Xm(:)

    ! get ensembles
    IF (PRESENT(Xm)) THEN
      this%Ens = Ens_t(configFile, Xm)
    ELSE
      this%Ens = Ens_t(configFile)
    END IF
    PRINT *, 'DONE getting ensembles.'

  END SUBROUTINE sub_b_getEnsData


  SUBROUTINE sub_b_qrDecomp(this, configFile, Xm, iGrid)
    IMPLICIT NONE
    CLASS(EnLoc_t) :: this
    CHARACTER(LEN=1024), INTENT(IN)       :: configFile
    TYPE(State_t), INTENT(IN)             :: Xm(:)
    INTEGER(i_kind), INTENT(IN), OPTIONAL :: iGrid

    INTEGER(i_kind)           :: i, m, n, status, pass_test, iv, numVar, ifile, stat
    INTEGER(i_kind)           :: numLevel, numHorz, numTime, iens
    REAL(r_kind), ALLOCATABLE :: at(:, :), a(:, :), values(:), vectors(:, :), amin, amax
    REAL(r_kind), ALLOCATABLE :: dataSwap(:, :, :)
    CHARACTER(LEN=1024)       :: outputDir, logFileName, datFileName
    CHARACTER(LEN=20)         :: varName
    CHARACTER(LEN=3)          :: gid
    character(len=1024)       :: msg

    IF (.NOT. this%Ens%dataReady) THEN
      PRINT *, 'ERROR!!! There is no data for QR factorization.'
      STOP
    ELSE
      PRINT *, 'Ensemble data is ready.'
    END IF

    ifile  = yaml_get_var(TRIM(configFile), 'IO', 'output_dir', outputDir)
    numVar = this%Ens%numVar

    IF (Xm(1)%mpddGlob%isBaseProc())  THEN

      ALLOCATE (this%BKErr(numVar))
      DO iv = 1, numVar
        ! Default as passing tests:
        pass_test = 0

        varName = TRIM(this%Ens%EnsData(iv)%varName)
        PRINT *, "dealing with variable:", TRIM(varName)

        IF (PRESENT(iGrid)) THEN
          WRITE (gid, '(A1,I2.2)') "G", iGrid
          logFileName = TRIM(outputDir)//'/BKErr_QRdecomp_'//TRIM(varName)//'_'//gid//'.log'
        ELSE
          logFileName = TRIM(outputDir)//'/BKErr_QRdecomp_'//TRIM(varName)//'.log'
        END IF

        m = this%Ens%EnsData(iv)%numC
        n = this%Ens%EnsData(iv)%ensNum
        numLevel = this%Ens%EnsData(iv)%numLevel
        numHorz  = this%Ens%EnsData(iv)%numHorz
        numTime  = this%Ens%EnsData(iv)%tSlots
        PRINT *, 'm, n:', m, n

        !Open log file for keeping run information:
        OPEN (10, FILE=TRIM(logFileName), status="replace", action="write", iostat=stat, iomsg=msg)
        !if (stat /= 0) then
        !  print *, "state of creating log file: ", stat
        !  print *, "Error creating log file: ", trim(msg)
        !  stop
        !end if
        WRITE (10, *, iostat=stat, iomsg=msg) 'QR factorization for ', TRIM(varName)
        !print *, "state of writing log file: ", stat
        !print *, "Error writing log file: ", trim(msg)

        ! Initialize
        CALL this%BKErr(iv)%init(m, REAL(1.0E-4,kind=r_kind))
        IF (ALLOCATED(a)) DEALLOCATE(a)
        IF (ALLOCATED(at)) DEALLOCATE(at)
        ALLOCATE (a(m, n))
        ALLOCATE (at(m, n))
        
        IF (ALLOCATED(dataSwap)) DEALLOCATE(dataSwap)
        ALLOCATE(dataSwap(numLevel, numHorz, numTime))
        DO iens = 1, n
          CALL Xm(iens)%fields(iv)%sg%aggrGridRealForFieldGrid(Xm(iens)%fields(iv)%data, &
                  dataSwap, [numLevel, numHorz, numTime])
          at(:, iens) = RESHAPE(dataSwap, (/m/))
        END DO
        DEALLOCATE(dataSwap)
        
        DO i = 1, n
          a(:, i) = at(:, i) - SUM(at, 2) / n
        END DO

        a = a / SQRT(m * 1.0D0)

        CALL this%BKErr(iv)%b_qrDecomposition(n, a, 1, status)
        PRINT *, 'done b_qrDecomposition, status:', status

        WRITE (10, 11) this%BKErr(iv)%nz
11      FORMAT('Number of rank: ', i5)

        IF (status .EQ. 0) THEN
          WRITE (10, 12) this%BKErr(iv)%qr_aux
12        FORMAT('qr_aux: ', 1E18.8)

          WRITE (10, 13) this%BKErr(iv)%xe_dr
13        FORMAT('qr_xdr: ', 1E18.8)
        ELSE
          WRITE (10, *) 'No qraux or qr_a is available!'
          pass_test = 1
          PRINT *, 'pass_status:', pass_test
          STOP
        END IF

        CLOSE (10)
        PRINT *, 'pass_status:', pass_test
        DEALLOCATE (a, at)

      END DO

    END IF

  END SUBROUTINE sub_b_qrDecomp


  SUBROUTINE sub_b_qrSave(this, configFile, Xm, iGrid)
    IMPLICIT NONE
    CLASS(EnLoc_t)  :: this
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    TYPE(State_t), INTENT(IN)       :: Xm
    INTEGER(i_kind), OPTIONAL       :: iGrid

    INTEGER(i_kind) :: i, m, n, status, pass_test, iv, ifile
    CHARACTER(LEN=1024)        :: outputDir, datFileName
    CHARACTER(LEN=20)          :: varName
    CHARACTER(LEN=3)           :: gid

    IF (Xm%sg%mpddGlob%isBaseProc()) THEN

      ifile = yaml_get_var(TRIM(configFile), 'IO', 'output_dir', outputDir)
      DO iv = 1, this%Ens%numVar
        varName = TRIM(this%Ens%EnsData(iv)%varName)
        PRINT *, "saving qr decomp of variable:", TRIM(varName)

        IF (PRESENT(iGrid)) THEN
          WRITE (gid, '(A1,I2.2)') "G", iGrid
          datFileName = TRIM(outputDir)//'/BKErr_QRdecomp_'//TRIM(varName)//'_'//gid//'.dat'
        ELSE
          datFileName = TRIM(outputDir)//'/BKErr_QRdecomp_'//TRIM(varName)//'.dat'
        END IF

        CALL this%BKErr(iv)%b_qrDecomp_save(TRIM(datFileName), this%Ens%ensBKDiag, status)
        PRINT *, 'done qr decomp saving'
      END DO

    END IF

  END SUBROUTINE sub_b_qrSave


  SUBROUTINE sub_b_qrLoad(this, configFile, Xm, iGrid)
    IMPLICIT NONE
    CLASS(EnLoc_t)  :: this
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    TYPE(State_t), INTENT(IN)       :: Xm
    INTEGER(i_kind), OPTIONAL       :: iGrid

    INTEGER(i_kind)     :: status, iv, ifle, ifile
    CHARACTER(LEN=1024) :: outputDir, datFileName
    CHARACTER(LEN=20)   :: varName
    CHARACTER(LEN=3)    :: gid

    IF (Xm%sg%mpddGlob%isBaseProc()) THEN

      ifile = yaml_get_var(TRIM(configFile), 'IO', 'output_dir', outputDir)

      IF (.NOT. ALLOCATED(this%BKErr)) ALLOCATE (this%BKErr(this%Ens%numVar))

      DO iv = 1, this%Ens%numVar
        varName = TRIM(this%Ens%EnsData(iv)%varName)
        PRINT *, "dealing with variable:", TRIM(varName)

        IF (PRESENT(iGrid)) THEN
          WRITE (gid, '(A1,I2.2)') "G", iGrid
          datFileName = TRIM(outputDir)//'/BKErr_QRdecomp_'//TRIM(varName)//'_'//gid//'.dat'
        ELSE
          datFileName = TRIM(outputDir)//'/BKErr_QRdecomp_'//TRIM(varName)//'.dat'
        END IF
        PRINT *, 'datFileName:', TRIM(datFileName)

        CALL this%BKErr(iv)%b_qrDecomp_read(TRIM(datFileName), status)
        PRINT *, "done b_qrDecomp_read of var ", TRIM(varName)
      END DO

    END IF

  END SUBROUTINE sub_b_qrLoad


  SUBROUTINE sub_b_qrRinv(this, Xm)
    IMPLICIT NONE
    CLASS(EnLoc_t)  :: this
    TYPE(State_t), INTENT(IN) :: Xm

    INTEGER(i_kind)           :: iv

    IF (Xm%sg%mpddGlob%isBaseProc()) THEN

      DO iv = 1, this%Ens%numVar
        CALL this%BKErr(iv)%b_qrDecomp_Rinv()
      END DO
      PRINT *, "done sub_b_qrRinv"

    END IF

  END SUBROUTINE sub_b_qrRinv


  SUBROUTINE sub_b_FWD(this, configFile, Xm, iGrid)
    IMPLICIT NONE
    CLASS(EnLoc_t)  ::  this
    CHARACTER(LEN=1024), INTENT(IN)       :: configFile
    TYPE(State_t), INTENT(IN)             :: Xm
    INTEGER(i_kind), INTENT(IN), OPTIONAL :: iGrid

    CHARACTER(LEN=1024) :: outputDir, datFileName
    CHARACTER(LEN=20)   :: varName
    CHARACTER(LEN=3)    :: gid
    INTEGER(i_kind)     :: iv, status, ifile

    IF (Xm%sg%mpddGlob%isBaseProc()) THEN

      DO iv = 1, this%Ens%numVar
        CALL this%BKErr(iv)%b_forward()
      END DO
      PRINT *, "done FWD calculation"

      ifile = yaml_get_var(TRIM(configFile), 'IO', 'output_dir', outputDir)

      DO iv = 1, this%Ens%numVar
        status = 0
        varName = TRIM(this%Ens%EnsData(iv)%varName)

        IF (PRESENT(iGrid)) THEN
          WRITE (gid, '(A1,I2.2)') "G", iGrid
          datFileName = TRIM(outputDir)//'/BEC_EnLoc_FWDT_'//TRIM(varName)//'_'//gid
        ELSE
          datFileName = TRIM(outputDir)//'/BEC_EnLoc_FWDT_'//TRIM(varName)
        END IF

        CALL this%BKErr(iv)%b_forwardT_save(TRIM(datFileName), status)
        IF (status .EQ. 0) THEN
          PRINT *, "done FWDT save of var ", TRIM(varName)
        ELSE
          PRINT *, "error in FWDT save of var ", TRIM(varName)
        END IF
      END DO

      IF (status .EQ. 0) THEN
        PRINT *, "done sub_b_FWD"
      ELSE
        PRINT *, "error sub_b_FWD"
      END IF

    END IF

  END SUBROUTINE sub_b_FWD


  SUBROUTINE sub_b_ADJ(this, configFile, Xm, iGrid)
    IMPLICIT NONE
    CLASS(EnLoc_t)  ::  this
    CHARACTER(LEN=1024), INTENT(IN)       :: configFile
    TYPE(State_t), INTENT(IN)             :: Xm
    INTEGER(i_kind), INTENT(IN), OPTIONAL :: iGrid

    CHARACTER(LEN=1024) :: outputDir, datFileName
    CHARACTER(LEN=20)   :: varName
    CHARACTER(LEN=3)    :: gid
    INTEGER(i_kind)     :: iv, status, ifile
    INTEGER(i_kind)     :: count_start, count_end, count_rate
    REAL(r_kind)        :: ts

    IF (Xm%sg%mpddGlob%isBaseProc()) THEN

      CALL system_clock(count_start, count_rate)
      DO iv = 1, this%Ens%numVar
        CALL this%BKErr(iv)%b_adjoint()
      END DO
      CALL system_clock(count_end)
      ts = REAL(count_end - count_start) / REAL(count_rate)
      PRINT *, "Done ADJ calculation, cost ", ts

      ifile = yaml_get_var(TRIM(configFile), 'IO', 'output_dir', outputDir)

      DO iv = 1, this%Ens%numVar
        status = 0
        varName = TRIM(this%Ens%EnsData(iv)%varName)

        IF (PRESENT(iGrid)) THEN
          WRITE (gid, '(A1,I2.2)') "G", iGrid
          datFileName = TRIM(outputDir)//'/BEC_EnLoc_ADJ_'//TRIM(varName)//'_'//gid
        ELSE
          datFileName = TRIM(outputDir)//'/BEC_EnLoc_ADJ_'//TRIM(varName)
        END IF

        CALL system_clock(count_start, count_rate)
        CALL this%BKErr(iv)%b_adjoint_save(TRIM(datFileName), status)
        CALL system_clock(count_end)
        ts = REAL(count_end - count_start) / REAL(count_rate)
        IF (status .EQ. 0) THEN
          PRINT *, "Done ADJ save of var ", TRIM(varName), ', cost ', ts
        ELSE
          PRINT *, "error in ADJ save of var ", TRIM(varName)
        END IF
      END DO

      IF (status .EQ. 0) THEN
        PRINT *, "done sub_b_ADJ"
      ELSE
        PRINT *, "error sub_b_ADJ"
      END IF

    END IF

  END SUBROUTINE sub_b_ADJ


  SUBROUTINE sub_b_destroy(this)
    IMPLICIT NONE
    CLASS(EnLoc_t) :: this

    print *, 'begin EnLoc_t%sub_b_destroy'
    IF (ALLOCATED(this%BKErr)) DEALLOCATE (this%BKErr)
    print *, 'done EnLoc_t%sub_b_destroy'

  END SUBROUTINE sub_b_destroy


END MODULE EnLoc_m
