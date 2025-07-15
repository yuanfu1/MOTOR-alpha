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
MODULE rttov_nl_m
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

  TYPE rttov_nl_t
    TYPE(State_t), POINTER :: X
    TYPE(ObsSet_t), POINTER :: Y
    CHARACTER(LEN=1024) :: configFile
    TYPE(rttov_nl_sp_t) :: rttov_nl_sp
    INTEGER(i_kind), ALLOCATABLE :: rttov_chan_lists(:)
    INTEGER(i_kind) :: nchans, IdxStart

  CONTAINS
    PROCEDURE, PUBLIC, PASS(this) :: initialize
    FINAL :: destructor
    PROCEDURE, PUBLIC :: rttov_nl_simobs
  END TYPE rttov_nl_t

CONTAINS

SUBROUTINE initialize(this, configFile, X, Y, IdxStart, nchans, rttov_chan_lists, inst_name, platform_name)
  IMPLICIT NONE
  CLASS(rttov_nl_t) :: this
  CHARACTER(LEN=1024), INTENT(IN) :: configFile
  TYPE(State_t), TARGET, INTENT(IN) :: X
  TYPE(ObsSet_t), TARGET, INTENT(IN) :: Y
  INTEGER(i_kind), INTENT(IN) :: IdxStart, nchans
  INTEGER(i_kind), INTENT(IN) :: rttov_chan_lists(:)
  CHARACTER(len=*), INTENT(IN) :: inst_name, platform_name
  INTEGER(i_kind) :: i

  this%X => X
  this%Y => Y
  this%configFile = configFile
  this%nchans = nchans
  this%IdxStart = IdxStart
  this%rttov_chan_lists = rttov_chan_lists

    CALL this%rttov_nl_sp%initialize(configFile, X, TRIM(inst_name), TRIM(platform_name))

  END SUBROUTINE initialize

  SUBROUTINE rttov_nl_simobs(this, X, Hx)
    CLASS(rttov_nl_t) :: this
    TYPE(State_t), INTENT(IN) :: X
    TYPE(ObsSet_t), INTENT(INOUT) :: Hx

    ! Local variables:
    REAL(r_kind) :: Hx_sp
    INTEGER(i_kind) :: i_inst, n_insts, i_prof, nprofs, nchans
    INTEGER(i_kind) :: l, i, j, k, iz, iobs, ichan, nobs
    INTEGER(i_kind) :: locX(2), ih, it, ifield
    LOGICAL :: istat
    REAL(r_kind) :: SatAngles(4)
    INTEGER(i_kind) :: ifile

    ! PRINT *, 'Total channel numbers: ', SIZE(Hx%ObsFields,1)
    ! PRINT *, 'rttov_chan_lists: ', this%rttov_chan_lists
    ! PRINT *, 'nchans: ', this%nchans

    DO j = LBOUND(X%fields, 1), UBOUND(X%fields, 1)
      ! PRINT *, 'Check NL X varnames: ', TRIM(X%fields(j)%Get_Name()), MAXVAL(X%fields(j)%data), MINVAL(X%fields(j)%data)
      IF ('qvapor' .EQ. TRIM(X%fields(j)%Get_Name()) .AND. (MAXVAL(X%fields(j)%DATA) > 0.055 .OR. MINVAL(X%fields(j)%DATA) < 0.0)) THEN
        PRINT *, '--- Check qvapor values during minimization, STOP --- '
        PRINT *, 'MAX=', MAXVAL(X%fields(j)%DATA), 'MIN=', MINVAL(X%fields(j)%DATA)
        STOP
      END IF
      IF ('temp' .EQ. TRIM(X%fields(j)%Get_Name()) .AND. (MAXVAL(X%fields(j)%DATA) > 350.0 .OR. MINVAL(X%fields(j)%DATA) < 150.0)) THEN
        PRINT *, '--- Check temp values during minimization, STOP --- '
        PRINT *, 'MAX=', MAXVAL(X%fields(j)%DATA), 'MIN=', MINVAL(X%fields(j)%DATA)
        STOP
      END IF
    END DO

    ! nprofs = X%sg%num_icell ! sg%num_icell or sg%num_cell

    ! For Hx%ObsFields(n), n is nchans
    ! PRINT *, 'Check SIZE of Y ', SIZE(Hx%ObsFields, 1)
    IF (X%sg%isActiveProc()) THEN
      ! DO ichan = 1, this%nchans
      ichan = 0
      ! DO ifield = 1, SIZE(Hx%ObsFields, 1)
      DO ifield = this%IdxStart, this%IdxStart + this%nchans - 1

        IF (TRIM(Hx%ObsFields(ifield)%Get_Name()) .NE. 'tbb') CYCLE
        ichan = ichan + 1

        IF (ichan .GE. this%nchans + this%IdxStart) PRINT *, 'ERROR: channel number exceeds MAXMA'

        ! nobs = UBOUND(Hx%ObsFields(i_inst)%valueArray(ichan, :), 1)
        nobs = UBOUND(Hx%ObsFields(ifield)%values(:), 1)
        ! PRINT *,'RTTOV NL nobs = ', nobs

        DO iobs = 1, nobs
          ih = Hx%ObsFields(ifield)%idx(iobs)%hIdx
          it = Hx%ObsFields(ifield)%idx(iobs)%tIdx
          locX = (/ih, it/)
          SatAngles(1) = Hx%ObsFields(ifield)%ObsAttrSat%zenangle(iobs)
          SatAngles(2) = Hx%ObsFields(ifield)%ObsAttrSat%azangle(iobs)
          SatAngles(3) = Hx%ObsFields(ifield)%ObsAttrSat%sunzenangle(iobs)
          SatAngles(4) = Hx%ObsFields(ifield)%ObsAttrSat%sunazangle(iobs)

          CALL this%rttov_nl_sp%rttov_nl_sp_simobs(X, this%rttov_chan_lists(ichan), locX, SatAngles, Hx_sp)
          ! Hx%ObsFields(i_inst)%valueArray(ichan,iobs) = Hx_sp
          Hx%ObsFields(ifield)%values(iobs) = Hx_sp

        END DO

        ! IF (nobs>0) PRINT *, 'HX = ', MAXVAL(Hx%ObsFields(ifield)%values), MINVAL(Hx%ObsFields(ifield)%values)
        ! IF (nobs>0) PRINT *, 'Y = ', MAXVAL(this%Y%ObsFields(ifield)%values), MINVAL(this%Y%ObsFields(ifield)%values)

      END DO

    END IF
    ! PRINT *, 'END rttov_nl_simobs'
  END SUBROUTINE rttov_nl_simobs

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(rttov_nl_t), INTENT(INOUT) :: this

    IF (ALLOCATED(this%rttov_chan_lists)) DEALLOCATE (this%rttov_chan_lists)

  END SUBROUTINE destructor

END MODULE rttov_nl_m
