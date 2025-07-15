!!--------------------------------------------------------------------------------------------------
! PROJECT           : rttov_nl
! AFFILIATION       : GBA-MWF(SIMI)
! AUTOHR(S)         : Yali Wu
! VERSION           : V 0.0
! HISTORY           :
!   Created by Yali Wu, 2021/8/27
! Information       : Add interfaces to call rttov NL..
!                   : Only the clear-sky capability is enabled so far.
!                   : Requirement: rttov V13
! Function          : Calculate J_radiance and grad_J_radiance wrt CV variables

!> @brief
!!
MODULE rttov_diag_out_m
  USE State_m, ONLY: State_t
  USE ObsSet_m, ONLY: ObsSet_t
  USE ObsField_m, ONLY: ObsField_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE kinds_m, ONLY: r_kind, i_kind
  ! the same functionality, but on a single grid point
  USE rttov_nl_sp_m, ONLY: rttov_nl_sp_t
  USE rttov_typedefine_m
  USE YAMLRead_m
  USE mpObs_m, ONLY: mpObs_t
  USE Satellite_utils_m
 
  CONTAINS

  ! Each channel may have different obs numbers, and those not used obs have been removed after ObsThinning
  ! Thus, it may be useful if we retain T/Q/PRES for each channel

  !!! Note to use_kmodel together with write_diag !!!
  SUBROUTINE write_diag_vars(configFile, X, platform_name, inst_name, Y, X2)
    ! USE QVCtl2QV_sv_m, ONLY: QVCtl2QV_sv_t
    ! USE TempCtl2Temp_Exp_m, ONLY: TempCtl2Temp_Exp_t
    ! USE TempCtl2Temp_sv_m, ONLY: TempCtl2Temp_sv_t
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    TYPE(State_t), INTENT(IN) :: X
    CHARACTER(len=*), INTENT(IN) :: platform_name, inst_name
    TYPE(ObsSet_t), INTENT(IN) :: Y
    TYPE(State_t), INTENT(OUT) :: X2

    ! Local variables:
    TYPE(rttov_nl_sp_t) :: rttov_nl_sp
    ! TYPE(QVCtl2QV_sv_t) :: QVCtl2QV
    ! TYPE(TempCtl2Temp_Exp_t) :: TempCtl2Temp_Exp
    ! TYPE(TempCtl2Temp_sv_t) :: TempCtl2Temp_sv
    ! CHARACTER(LEN=20) :: rhov_scale_scheme = 'DomainAvg'
    ! CHARACTER(LEN=20) :: temp_scale_scheme = 'DomainAvg'
    TYPE(State_t) :: X1
    TYPE(rttov_diag_output) :: diag
    TYPE(rttov_setup_input_t)   :: rttov_setup_input
    REAL(r_kind) :: Hx_sp
    INTEGER(i_kind) :: i_inst, n_insts, i_prof, nprofs, nchans
    INTEGER(i_kind) :: l, i, j, k, iz, iobs, ichan, nobs
    INTEGER(i_kind) :: locX(2), ih, it, ivar, ilev, nvars, nvars1, ivartmp
    CHARACTER(LEN=1024) :: satinfo_file
    LOGICAL :: istat
    REAL(r_kind) :: SatAngles(4)
    INTEGER(i_kind) :: ifile
    CHARACTER(len=5) :: ichar
    CHARACTER(len=10) :: diag_name(14)=(/'pres   ','T      ', 'qv     ', 'qc     ', 'qi     ', 'cldfrac', &
                                                   'JacT   ', 'Jacqv  ', 'Jacqc  ', 'Jacqi  ', &  
                                                   'tau    ', 'dtau   ', 'bt     ', 'btclr  '/)
    ! CHARACTER(len=10) :: diag_name(10)=(/'pres   ','T      ','qv     ', 'cldfrac', 'JacT   ', &
    !                                      'Jacqv  ','tau    ','dtau   ', 'bt     ', 'btclr  '/)
    INTEGER(i_kind), ALLOCATABLE :: chan_lists(:), rttov_chan_lists(:)
    INTEGER(i_kind), ALLOCATABLE :: chan_list_used(:)
    INTEGER(i_kind)      :: IdxStart, ifield

    X1 = X
    X2 = X
    DO ivar = LBOUND(X2%fields, 1), UBOUND(X2%fields, 1)
      CALL X2%rmVar(TRIM(X%fields(ivar)%Get_Name()))
    END DO

    ! ifile = yaml_get_var(TRIM(configFile), 'CV_Transform', 'rhov_scale_scheme', rhov_scale_scheme)
    ! ifile = yaml_get_var(TRIM(configFile), 'CV_Transform', 'temp_scale_scheme', temp_scale_scheme)
    ! print *,'inst_name:',inst_name

    ! IF (X1%getVarIdx(TRIM('qvapor_ctl')) .NE. 0) THEN
    !   QVCtl2QV = QVCtl2QV_sv_t(configFile, X)
    ! END IF
    ! IF (X1%getVarIdx(TRIM('temp_ctl')) .NE. 0) THEN
    !   TempCtl2Temp_Exp = TempCtl2Temp_Exp_t(configFile, X)
    !   TempCtl2Temp_sv = TempCtl2Temp_sv_t(configFile, X)
    ! END IF

    ! IF (X1%getVarIdx(TRIM('qvapor_ctl')) .NE. 0) &
    ! CALL QVCtl2QV%fwdNL(X1) ! get qvapor
  
    WHERE (X1%fields(X1%getVarIdx('qvapor'))%data < qv_limit) X1%fields(X1%getVarIdx('qvapor'))%data = qv_limit
    WHERE (X1%fields(X1%getVarIdx('qcloud'))%data < 1.0D-12) X1%fields(X1%getVarIdx('qcloud'))%data = 0.0D0
    WHERE (X1%fields(X1%getVarIdx('qice'))%data < 1.0D-12) X1%fields(X1%getVarIdx('qice'))%data = 0.0D0

    ! IF (X1%getVarIdx(TRIM('temp_ctl')) .NE. 0 ) THEN
    !   IF (TRIM(temp_scale_scheme) .EQ. 'Exp') THEN 
    !     CALL TempCtl2Temp_Exp%fwdNL(X1) ! get temp
    !   ELSE 
    !     CALL TempCtl2Temp_sv%fwdNL(X1) ! get temp
    !   END IF
    ! END IF

    IF (X%sg%isActiveProc()) THEN

      nvars1 = SIZE(diag_name, 1)

      CALL rttov_nl_sp%initialize(configFile, X1, TRIM(inst_name), TRIM(platform_name))

      CALL Get_rttov_chan_info(trim(platform_name)//'-'//trim(inst_name), nchans)
      ALLOCATE (chan_lists(nchans), rttov_chan_lists(nchans))
      CALL Get_rttov_chan_info(trim(platform_name)//'-'//trim(inst_name), nchans, chan_lists, rttov_chan_lists)

      ! background profile, jacobian, tau and dtau are written for each channel
      IdxStart = Y%getObsIdx('tbb')
      print *, 'nchans in diag: ', IdxStart, IdxStart + nchans -1, size(Y%ObsFields,1)
      print *, 'chan_lists: ', chan_lists
      print *, 'rttov_chan_lists: ', rttov_chan_lists
      IF (IdxStart > 0) THEN
        ichan = 0
        DO ifield = IdxStart, IdxStart + nchans-1
          ichan = ichan + 1
          ! PRINT *, 'diag for ', ichan, ' channel'
          nobs = UBOUND(Y%ObsFields(ifield)%values, 1)
          PRINT *, 'ichan = ', ichan, 'nobs = ', nobs

          IF (chan_lists(ichan) <10) THEN
            write(ichar, "(I1)") chan_lists(ichan)
          ELSE
            IF (chan_lists(ichan) <100) write(ichar, "(I2)") chan_lists(ichan)
          END IF

          ivartmp = 0
          DO ivar = nvars1 * (ichan-1) + 1, nvars1 * ichan 
            ivartmp = ivartmp + 1
            CALL X2%addVar(TRIM(diag_name(ivartmp))//'_ch'//TRIM(ichar))
            ! print *,ivar,ivartmp,'X2%fields(ivar)%data:',SIZE(X2%fields(ivar)%data),ichan,TRIM(diag_name(ivartmp))//'_ch'//TRIM(ichar)
            X2%fields(ivar)%data = missing
          END DO
          ! IF (iobs == 1) PRINT *, 'diag for ', ivar, ' var'

          DO iobs = 1, nobs
            ih = Y%ObsFields(ifield)%idx(iobs)%hIdx
            it = Y%ObsFields(ifield)%idx(iobs)%tIdx
            locX = (/ih, it/)
            ! PRINT *, 'locX ', locX,ih,it
            SatAngles(1) = Y%ObsFields(ifield)%ObsAttrSat%zenangle(iobs) 
            SatAngles(2) = Y%ObsFields(ifield)%ObsAttrSat%azangle(iobs)
            SatAngles(3) = Y%ObsFields(ifield)%ObsAttrSat%sunzenangle(iobs) 
            SatAngles(4) = Y%ObsFields(ifield)%ObsAttrSat%sunazangle(iobs) 
            CALL rttov_nl_sp%rttov_nl_sp_simobs(X1, rttov_chan_lists(ichan), locX, SatAngles, Hx_sp, diag)

            ivartmp = 0
            DO ivar = nvars1 * (ichan-1) + 1, nvars1 * ichan 
              ivartmp = ivartmp + 1
              X2%fields(ivar)%data(:, ih, it) = diag%outvar(ivartmp, :)
              ! PRINT *, 'Check diag var values: ', ivar, MAXVAL(X2%fields(ivar)%data(:, ih, it)), MINVAL(X2%fields(ivar)%data(:, ih, it))
            END DO
          END DO
        END DO
      END IF
        
      IF (ALLOCATED(chan_lists)) DEALLOCATE(chan_lists) 
      IF (ALLOCATED(rttov_chan_lists)) DEALLOCATE(rttov_chan_lists)

    END IF

  END SUBROUTINE write_diag_vars

  ! SUBROUTINE write_diag_vars(configFile, X, platform_name, inst_name, Y, X2)
  !   CHARACTER(LEN=1024), INTENT(IN) :: configFile
  !   TYPE(State_t), INTENT(IN) :: X
  !   CHARACTER(len=*), INTENT(IN) :: platform_name, inst_name
  !   TYPE(ObsSet_t), INTENT(IN) :: Y
  !   TYPE(State_t), INTENT(OUT) :: X2

  !   ! Local variables:
  !   TYPE(rttov_nl_sp_t) :: rttov_nl_sp
  !   TYPE(State_t) :: X1
  !   TYPE(rttov_diag_output) :: diag
  !   TYPE(rttov_setup_input_t)   :: rttov_setup_input
  !   REAL(r_kind) :: Hx_sp
  !   INTEGER(i_kind) :: i_inst, n_insts, i_prof, nprofs, nchans
  !   INTEGER(i_kind) :: l, i, j, k, iz, iobs, ichan, nobs
  !   INTEGER(i_kind) :: locX(2), ih, it, ivar, ilev, nvars, ivartmp
  !   CHARACTER(LEN=1024) :: satinfo_file
  !   LOGICAL :: istat
  !   REAL(r_kind) :: SatAngles(4)
  !   INTEGER(i_kind) :: ifile
  !   CHARACTER(len=5) :: ichar
  !   CHARACTER(len=10) :: diag_name(8)=(/'pres ','T    ','qv   ','JacT ','Jacqv','tau  ','dtau ', 'bt   '/)
  !   INTEGER(i_kind), ALLOCATABLE :: chan_lists(:), rttov_chan_lists(:)

  !   X1 = X
  !   X2 = X
  !   DO ivar = LBOUND(X2%fields, 1), UBOUND(X2%fields, 1)
  !     CALL X2%rmVar(TRIM(X%fields(ivar)%Get_Name()))
  !   END DO

  !   IF (X%sg%isActiveProc()) THEN

  !     CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", satinfo_file)

  !     satinfo_file = TRIM(satinfo_file)//"/Satellite/satinfo_"// &
  !                   trim(platform_name)//'-'//trim(inst_name)//".txt"

  !     inquire (file=satinfo_file, exist=istat)
  !     IF (istat) THEN
  !       open (unit=121, file=satinfo_file, status='old')
  !       read (121, *)
  !       read (121, *) nchans
  !       ALLOCATE (chan_lists(nchans))
  !       ALLOCATE (rttov_chan_lists(nchans))

  !       DO ichan = 1, nchans
  !         read (121, *) chan_lists(ichan), rttov_chan_lists(ichan)
  !       END DO
  !       close (121)
  !     ELSE
  !       PRINT *, '----------satinfo file was not found----------'
  !       STOP
  !     END IF

  !     nvars = SIZE(diag_name, 1)
  !     ! profiles and jac profiles are only written once
  !     ! tau and dtau are written for each channel
  !     nvars = nvars + 5 * (nchans -1)
  !     ! PRINT *, 'Out diag nvars = ', nvars
  !     nobs = UBOUND(Y%ObsFields(1)%values, 1)

  !     rttov_nl_sp = rttov_nl_sp_t(configFile, X1, TRIM(inst_name), TRIM(platform_name))

  !     ! bkg profiles are only written once
  !     ! PRINT *, 'nobs = ', nobs
  !     DO ivar = 1, 3
  !       CALL X2%addVar(TRIM(diag_name(ivar)))
  !       X2%fields(ivar)%data = missing
  !     END DO
  !     DO iobs = 1, nobs
  !       ih = Y%ObsFields(1)%idx(iobs)%hIdx
  !       it = Y%ObsFields(1)%idx(iobs)%tIdx
  !       locX = (/ih, it/)
  !       SatAngles(1) = Y%ObsFields(1)%ObsAttrSat%zenangle(iobs)
  !       SatAngles(2) = Y%ObsFields(1)%ObsAttrSat%azangle(iobs)
  !       SatAngles(3) = Y%ObsFields(1)%ObsAttrSat%sunzenangle(iobs)
  !       SatAngles(4) = Y%ObsFields(1)%ObsAttrSat%sunazangle(iobs)
  !       IF (.NOT. ALL(((Y%ObsFields(1)%values) - missing ) < 1.0D0)) &
  !         CALL rttov_nl_sp%rttov_nl_sp_simobs(X1, rttov_chan_lists(1), locX, SatAngles, Hx_sp, diag)

  !       DO ivar = 1, 3
  !         X2%fields(ivar)%data(:, ih, it) = diag%outvar(ivar, :)
  !         IF (iobs == 1) PRINT *, 'diag for ', ivar, ' var'
  !       END DO
  !     END DO

  !     ! jacobian, tau and dtau are written for each channel
  !     DO ichan = 1, nchans
  !       ! PRINT *, 'diag for ', ichan, ' channel'
  !       IF (ALL(((Y%ObsFields(ichan)%values) - missing ) < 1.0D0)) RETURN
  !       write(ichar, "(I1)") ichan
  !       ivartmp = 4
  !       DO ivar = 5 * ichan - 1, 5 * ichan + 3
  !         CALL X2%addVar(TRIM(diag_name(ivartmp))//'_ch'//TRIM(ichar))
  !         X2%fields(ivar)%data = missing
  !         ivartmp = ivartmp + 1
  !       END DO
  !       ! IF (iobs == 1) PRINT *, 'diag for ', ivar, ' var'

  !       DO iobs = 1, nobs
  !         ih = Y%ObsFields(1)%idx(iobs)%hIdx
  !         it = Y%ObsFields(1)%idx(iobs)%tIdx
  !         locX = (/ih, it/)
  !         SatAngles(1) = Y%ObsFields(1)%ObsAttrSat%zenangle(iobs)
  !         SatAngles(2) = Y%ObsFields(1)%ObsAttrSat%azangle(iobs)
  !         SatAngles(3) = Y%ObsFields(1)%ObsAttrSat%sunzenangle(iobs)
  !         SatAngles(4) = Y%ObsFields(1)%ObsAttrSat%sunazangle(iobs)
  !         CALL rttov_nl_sp%rttov_nl_sp_simobs(X1, rttov_chan_lists(ichan), locX, SatAngles, Hx_sp, diag)

  !         ivartmp = 4
  !         DO ivar = 5 * ichan - 1, 5 * ichan + 3
  !           X2%fields(ivar)%data(:, ih, it) = diag%outvar(ivartmp, :)
  !           ! PRINT *, 'Check diag var values: ', ivar, ivartmp, MAXVAL(X2%fields(ivar)%data(:, ih, it)), MINVAL(X2%fields(ivar)%data(:, ih, it))
  !           ivartmp = ivartmp + 1
  !         END DO
  !       END DO
  !     END DO

  !     IF (ALLOCATED(chan_lists)) DEALLOCATE(chan_lists)
  !     IF (ALLOCATED(rttov_chan_lists)) DEALLOCATE(rttov_chan_lists)

  !   END IF

  ! END SUBROUTINE write_diag_vars

END MODULE rttov_diag_out_m
