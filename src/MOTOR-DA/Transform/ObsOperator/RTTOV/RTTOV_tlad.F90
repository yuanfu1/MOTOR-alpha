!!--------------------------------------------------------------------------------------------------
! PROJECT           : rttov_tlad
! AFFILIATION       : GBA-MWF(SIMI)
! AUTOHR(S)         : Yali Wu
! VERSION           : V 0.0
! HISTORY           :
!   Created by Yali Wu, 2021/8/27
! Information       : Add interfaces to call rttov TL/AD..
!                   : Only the clear-sky capability is enabled so far.
!                   : Requirement: rttov V13
! Function          : Calculate J_radiance and grad_J_radiance wrt CV variables
!
!> @brief
!!
MODULE rttov_tlad_m
  USE State_m, ONLY: State_t
  USE State2NC_m
  USE ObsSet_m, ONLY: ObsSet_t
  USE ObsField_m, ONLY: ObsField_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE kinds_m, ONLY: r_kind, i_kind
  ! Self defined rttov utilities and types
  USE rttov_utils_m, ONLY: rttov_utils_t
  USE rttov_tlad_sp_m, ONLY: rttov_tlad_sp_t
  USE rttov_typedefine_m
  USE YAMLRead_m
  USE mpObs_m, ONLY: mpObs_t
  USE parameters_m, ONLY: hPa2Pa
  USE parkind1, ONLY: jpim

  TYPE rttov_tlad_t
    TYPE(State_t), POINTER :: X
    TYPE(ObsSet_t), POINTER :: Y
    TYPE(rttov_utils_t) :: rttov_utils
    TYPE(rttov_tlad_sp_t) :: rttov_tlad_sp
    CHARACTER(LEN=1024) :: configFile
    INTEGER(i_kind), ALLOCATABLE :: rttov_chan_lists(:)
    INTEGER(i_kind) :: nchans, IdxStart
    LOGICAL :: diag_grad = .FALSE.

  CONTAINS
    PROCEDURE, PUBLIC, PASS(this) :: initialize
    FINAL :: destructor_tlad
    PROCEDURE, PUBLIC :: rttov_tl_simobs
    PROCEDURE, PUBLIC :: rttov_ad_simobs
  END TYPE rttov_tlad_t

CONTAINS

SUBROUTINE initialize(this, configFile, X, Y, IdxStart, nchans, rttov_chan_lists, inst_name, platform_name)
    IMPLICIT NONE
    CLASS(rttov_tlad_t) :: this
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    TYPE(State_t), TARGET, INTENT(IN) :: X
    TYPE(ObsSet_t), TARGET, INTENT(IN) :: Y
    INTEGER(i_kind), INTENT(IN) :: IdxStart, nchans
    INTEGER(i_kind), INTENT(IN) :: rttov_chan_lists(:)
    CHARACTER(len=*), INTENT(IN) :: inst_name, platform_name

    this%X => X
    this%Y => Y
    this%configFile = configFile
    this%nchans = nchans
    this%IdxStart = IdxStart
    this%rttov_chan_lists = rttov_chan_lists
    CALL this%rttov_tlad_sp%initialize(configFile, X, TRIM(inst_name), TRIM(platform_name))

  END SUBROUTINE initialize

  SUBROUTINE rttov_tl_simobs(this, X, dX, dhx)
    CLASS(rttov_tlad_t) :: this
    TYPE(State_t), INTENT(IN) :: X
    TYPE(State_t), INTENT(IN) :: dX ! delta x
    TYPE(ObsSet_t), INTENT(INOUT) :: dhx  ! delta_BT

    ! Local variables:
    INTEGER(i_kind) :: i, i_inst, i_prof, nlevels, iobs, ichan, nobs
    INTEGER(i_kind) :: nprofs, ninst, nchans, nsimobs
    INTEGER(kind=jpim)  :: errstatus
    TYPE(rttov_setup_input_t)   :: rttov_setup_input
    REAL(r_kind) :: dhx_sp
    INTEGER(i_kind) :: locX(2), ih, it, ifield
    LOGICAL :: istat
    REAL(r_kind) :: SatAngles(4)

    ! nprofs = X%sg%num_icell ! sg%num_icell or sg%num_cell

    ! DO j = LBOUND(X%fields, 1), UBOUND(X%fields, 1)
    !   PRINT *, 'Check TL X varnames: ', TRIM(X%fields(j)%Get_Name()), MAXVAL(X%fields(j)%data), MINVAL(X%fields(j)%data)
    ! END DO
    ! DO j = LBOUND(dX%fields, 1), UBOUND(dX%fields, 1)
    !   PRINT *, 'Check TL dX varnames: ', TRIM(dX%fields(j)%Get_Name()), MAXVAL(dX%fields(j)%data), MINVAL(dX%fields(j)%data)
    ! END DO

    IF (X%sg%isActiveProc()) THEN

      ! DO ichan = 1, this%nchans
      ichan = 0
      ! DO ifield = 1, SIZE(dhx%ObsFields, 1)
      DO ifield = this%IdxStart, this%IdxStart + this%nchans - 1

        IF (TRIM(dhx%ObsFields(ifield)%Get_Name()) .NE. 'tbb') CYCLE
        ichan = ichan + 1

        IF (ichan .GE. this%nchans + this%IdxStart) PRINT *, 'ERROR: channel number exceeds MAXMA'

        ! nobs = UBOUND(dhx%ObsFields(i_inst)%valueArray(ichan, :), 1)
        nobs = UBOUND(dhx%ObsFields(ifield)%values(:), 1)

        DO iobs = 1, nobs
          ih = dhx%ObsFields(ifield)%idx(iobs)%hIdx
          it = dhx%ObsFields(ifield)%idx(iobs)%tIdx
          locX = (/ih, it/)
          SatAngles(1) = dhx%ObsFields(ifield)%ObsAttrSat%zenangle(iobs)
          SatAngles(2) = dhx%ObsFields(ifield)%ObsAttrSat%azangle(iobs)
          SatAngles(3) = dhx%ObsFields(ifield)%ObsAttrSat%sunzenangle(iobs)
          SatAngles(4) = dhx%ObsFields(ifield)%ObsAttrSat%sunazangle(iobs)
          CALL this%rttov_tlad_sp%rttov_tl_sp_simobs(X, this%rttov_chan_lists(ichan), locX, SatAngles, dX, dhx_sp)
          ! dhx%ObsFields(i_inst)%valueArray(ichan,iobs) = dhx_sp
          dhx%ObsFields(ifield)%values(iobs) = dhx_sp
        END DO

        ! IF (nobs>0) PRINT *, 'TL values = ', MAXVAL(dhx%ObsFields(ichan)%values), MINVAL(dhx%ObsFields(ichan)%values)
      END DO

    END IF
    ! PRINT *, 'END rttov_tl_simobs'
    ! PRINT *, 'TL values =', dhx%ObsFields(2)%values
  END SUBROUTINE rttov_tl_simobs

  SUBROUTINE rttov_ad_simobs(this, X, norm_d, grad_X)
    CLASS(rttov_tlad_t) :: this
    TYPE(State_t), INTENT(IN) :: X
    TYPE(ObsSet_t), INTENT(IN) :: norm_d  !O^(-1).(hx-y); d is innovation
    TYPE(State_t), TARGET, INTENT(INOUT) :: grad_X

    ! Local variables:
    INTEGER(i_kind) :: i, i_inst, i_prof, nlevels, iobs, ichan, nobs
    INTEGER(i_kind) :: nprofs, ninst, nchans, nsimobs
    INTEGER(kind=jpim)  :: errstatus
    TYPE(rttov_setup_input_t)   :: rttov_setup_input
    REAL(r_kind)  :: norm_d_sp, t2, t1
    INTEGER(i_kind) :: locX(2), ih, it, nlevs, ifield
    REAL(r_kind), ALLOCATABLE :: temp_sp(:), qvapor_sp(:), qcloud_sp(:), qice_sp(:)
    REAL(r_kind) :: psfc_sp
    LOGICAL :: istat
    REAL(r_kind) :: SatAngles(4)
    REAL(r_kind), POINTER :: temp(:,:,:), qvapor(:,:,:), qcloud(:,:,:), qice(:,:,:)

    ! DO j = LBOUND(X%fields, 1), UBOUND(X%fields, 1)
    !   PRINT *, 'Check AD X varnames: ', TRIM(X%fields(j)%Get_Name()), MAXVAL(X%fields(j)%data), MINVAL(X%fields(j)%data)
    ! END DO
    ! DO j = LBOUND(grad_X%fields, 1), UBOUND(grad_X%fields, 1)
    !   PRINT *, 'Check AD grad_X varnames: ', TRIM(grad_X%fields(j)%Get_Name()), MAXVAL(grad_X%fields(j)%data), MINVAL(grad_X%fields(j)%data)
    ! END DO

    IF (X%sg%isActiveProc() .AND. norm_d%getObsIdx('tbb') .NE. 0) THEN

      IF (grad_X%getVarIdx('temp') .NE. 0 .AND. grad_X%getVarIdx('qvapor') .NE. 0) THEN

        temp => grad_X%fields(grad_X%getVarIdx('temp'))%data
        qvapor => grad_X%fields(grad_X%getVarIdx('qvapor'))%data
        IF ( this%rttov_tlad_sp%rttov_utils%opts%rt_ir%addclouds ) &
          qcloud => grad_X%fields(grad_X%getVarIdx('qcloud'))%data
        IF ( this%rttov_tlad_sp%rttov_utils%opts%rt_ir%addclouds ) &
          qice => grad_X%fields(grad_X%getVarIdx('qice'))%data
        ! psfc => grad_X%fields(grad_X%getVarIdx('psfc')) )
      
        nlevs = SIZE(temp, 1)
        IF (.NOT. ALLOCATED(temp_sp)) ALLOCATE(temp_sp(nlevs))
        IF (.NOT. ALLOCATED(qvapor_sp)) ALLOCATE(qvapor_sp(nlevs))
        IF (.NOT. ALLOCATED(qcloud_sp)) ALLOCATE(qcloud_sp(nlevs))
        IF (.NOT. ALLOCATED(qice_sp)) ALLOCATE(qice_sp(nlevs))

        temp = ZERO
        qvapor = ZERO
        IF ( this%rttov_tlad_sp%rttov_utils%opts%rt_ir%addclouds ) qcloud = ZERO
        IF ( this%rttov_tlad_sp%rttov_utils%opts%rt_ir%addclouds ) qice = ZERO
        
        ! psfc = ZERO
      
        ! DO ichan = 1, this%nchans
        ichan = 0
        ! DO ifield = 1, SIZE(norm_d%ObsFields, 1)
        DO ifield = this%IdxStart, this%IdxStart + this%nchans - 1

          IF ( TRIM(norm_d%ObsFields(ifield)%Get_Name()) .NE. 'tbb' ) CYCLE
          ichan = ichan + 1

          IF (ichan .GE. this%nchans + this%IdxStart) PRINT *, 'ERROR: channel number exceeds MAXMA'

          ! IF (ALL(((this%Y%ObsFields(ichan)%values) - missing ) < 1.0D0)) RETURN
          ! nobs = UBOUND(norm_d%ObsFields(i_inst)%valueArray(ichan, :), 1)
          nobs = UBOUND(norm_d%ObsFields(ifield)%values(:), 1)
          ! IF (nobs>0) PRINT *, 'check norm_d', MAXVAL(norm_d%ObsFields(ichan)%values), MINVAL(norm_d%ObsFields(ichan)%values)
          DO iobs = 1, nobs
            ih = norm_d%ObsFields(ifield)%idx(iobs)%hIdx
            it = norm_d%ObsFields(ifield)%idx(iobs)%tIdx
            locX = (/ih, it/)
            ! norm_d_sp = norm_d%ObsFields(i_inst)%valueArray(ichan,iobs)
            norm_d_sp = norm_d%ObsFields(ifield)%values(iobs)
            
            temp_sp = ZERO
            qvapor_sp = ZERO
            ! psfc_sp = ZERO
            qcloud_sp = ZERO
            qice_sp = ZERO

            SatAngles(1) = norm_d%ObsFields(ifield)%ObsAttrSat%zenangle(iobs) 
            SatAngles(2) = norm_d%ObsFields(ifield)%ObsAttrSat%azangle(iobs)
            SatAngles(3) = norm_d%ObsFields(ifield)%ObsAttrSat%sunzenangle(iobs)
            SatAngles(4) = norm_d%ObsFields(ifield)%ObsAttrSat%sunazangle(iobs)

            IF (this%diag_grad) norm_d_sp = 1.0D0
            CALL this%rttov_tlad_sp%rttov_ad_sp_simobs(X, this%rttov_chan_lists(ichan), locX, SatAngles, &
            norm_d_sp, temp_sp, qvapor_sp, qcloud_sp, qice_sp)

            temp(:, ih, it) = temp(:, ih, it)+ temp_sp
            qvapor(:, ih, it) = qvapor(:, ih, it)+ qvapor_sp
            IF ( this%rttov_tlad_sp%rttov_utils%opts%rt_ir%addclouds ) qcloud(:, ih, it) = qcloud(:, ih, it)+ qcloud_sp
            IF ( this%rttov_tlad_sp%rttov_utils%opts%rt_ir%addclouds ) qice(:, ih, it) = qice(:, ih, it)+ qice_sp

            ! PRINT *, 'AD values : ', temp%data(:, ih, it), qvapor%data(:, ih, it)
            ! psfc%data(1, ih, it) = psfc%data(1, ih, it)+ psfc_sp *hPa2Pa
            ! ! ! Note: Add_Value is only for scalar, not work for satellite arrays ! ! ! 
            ! CALL psfc%Add_Value(norm_d%ObsFields(ichan)%idx(iobs), psfc_sp)
            ! CALL cpu_time(t2)
            ! PRINT *, 't2-t1 = ', t2-t1
            ! grad_x = grad_x + grad_x_tmp  !!!Time consuming!!!
          END DO 

          ! PRINT *, 'AD Temp values = ', MAXVAL(temp), MINVAL(temp)
          ! PRINT *, 'AD qvapor values = ', MAXVAL(qvapor), MINVAL(qvapor)
        END DO 

        CALL grad_X%sg%ExchangeMatOnHaloForFieldGrid(grad_X%sg%tSlots, grad_X%sg%vLevel, temp)
        CALL grad_X%sg%ExchangeMatOnHaloForFieldGrid(grad_X%sg%tSlots, grad_X%sg%vLevel, qvapor)
        IF ( this%rttov_tlad_sp%rttov_utils%opts%rt_ir%addclouds ) &
        CALL grad_X%sg%ExchangeMatOnHaloForFieldGrid(grad_X%sg%tSlots, grad_X%sg%vLevel, qcloud)
        IF ( this%rttov_tlad_sp%rttov_utils%opts%rt_ir%addclouds ) &
        CALL grad_X%sg%ExchangeMatOnHaloForFieldGrid(grad_X%sg%tSlots, grad_X%sg%vLevel, qice)

        IF (ALLOCATED(temp_sp)) DEALLOCATE(temp_sp)
        IF (ALLOCATED(qvapor_sp)) DEALLOCATE(qvapor_sp)
        IF (ALLOCATED(qcloud_sp)) DEALLOCATE(qcloud_sp)
        IF (ALLOCATED(qice_sp)) DEALLOCATE(qice_sp)
        NULLIFY(qvapor, temp)
        IF ( this%rttov_tlad_sp%rttov_utils%opts%rt_ir%addclouds ) NULLIFY(qcloud, qice)

        IF (this%diag_grad) THEN
          CALL Output_NC_State_AV(grad_X, "/Users/yaliwu/Desktop/MOTOR/MOTOR/output/", &
                                    "RTTOV_grad", .TRUE., .TRUE.)
        END IF                                  

      END IF

    END IF

    ! PRINT *, 'END rttov_ad_simobs'
    ! PRINT *, 'Check qvapor profile of AD: ', MAXVAL(grad_X%fields(grad_X%getVarIdx('qvapor'))%data), MINVAL(grad_X%fields(grad_X%getVarIdx('qvapor'))%data)
    ! PRINT *, 'Check temp profile of AD: ', MAXVAL(grad_X%fields(grad_X%getVarIdx('temp'))%data), MINVAL(grad_X%fields(grad_X%getVarIdx('temp'))%data)
  END SUBROUTINE rttov_ad_simobs

  IMPURE ELEMENTAL SUBROUTINE destructor_tlad(this)
    IMPLICIT NONE
    TYPE(rttov_tlad_t), INTENT(INOUT) :: this

    IF (ALLOCATED(this%rttov_chan_lists)) DEALLOCATE (this%rttov_chan_lists)

  END SUBROUTINE destructor_tlad

END MODULE rttov_tlad_m
