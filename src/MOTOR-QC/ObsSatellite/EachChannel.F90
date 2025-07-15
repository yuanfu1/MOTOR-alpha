!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-QC.Obs
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Yali Wu
! VERSION           : V 0.0
! HISTORY           :
!   Created by Yali Wu (wuyali@gbamwf.com), 2021/12/17, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module implements a satellite radiance data structure.
MODULE EachChannel_m
  USE kinds_m
  USE parameters_m
  USE SingleGrid_m, ONLY: SingleGrid_t, GridIdx_t
  USE ObsSet_m, ONLY: ObsSet_t
  USE ObsField_m, ONLY: ObsField_t
  USE ObsBase_m, ONLY: ObsBase_t
  USE State_m, ONLY: State_t
  USE geometry_m, ONLY: geometry_t
  USE mpObs_m, ONLY: mpObs_t
  USE YAMLRead_m
  USE FLog_m, ONLY: logger
  USE Satellite_QC_utils_m

  IMPLICIT NONE

  TYPE, EXTENDS(ObsBase_t) :: EachChannel_t
    REAL(r_kind), ALLOCATABLE :: tb_bakAtmodel(:,:)
    REAL(r_kind), ALLOCATABLE :: Thinned_Angles(:,:)
    REAL(r_kind) :: obsData_avg
    INTEGER(i_kind), ALLOCATABLE :: Thinned_Angles_hidx(:), Thinned_Angles_tidx(:)

  CONTAINS
    PROCEDURE, PUBLIC  :: ObsInitial => EachChannel_setup
    PROCEDURE, PUBLIC  :: ObsIngest => EachChannel_read
    PROCEDURE, PUBLIC  :: ObsForward => EachChannel_fwd
    PROCEDURE, PUBLIC  :: ObsTangent => EachChannel_tgt
    PROCEDURE, PUBLIC  :: ObsAdjoint => EachChannel_adj

    PROCEDURE, PUBLIC :: GetForwardValue  ! Used by superObs to calculate the increment of Obs
    PROCEDURE, PUBLIC :: ObsPrepareForSg  ! Used for observations which has to be prepared before thinning, like radar and satellite
    PROCEDURE, NOPASS, PUBLIC :: Field2ObsIsExisted
    PROCEDURE, PUBLIC :: destroy_EachChannel
  END TYPE EachChannel_t

CONTAINS

  SUBROUTINE ObsPrepareForSg(this, X)
    CLASS(EachChannel_t) :: this
    TYPE(State_t) :: X

    ! This is an optional holder for data prepare polymorphic subroutine
    ALLOCATE (this%tb_bakAtmodel(X%sg%num_cell, X%sg%tSlots))
    this%tb_bakAtmodel = ZERO

    CALL Calc_average(this%ObsData(:,1), this%ObsData_avg)
    PRINT *, 'Calc_average: ', this%ObsData_avg
  END SUBROUTINE

  FUNCTION GetForwardValue(this, X, varnametmp, vIdx, hIdx, tIdx, iv) RESULT(bt)
    CLASS(EachChannel_t) :: this
    TYPE(State_t) :: X
    INTEGER(i_kind), INTENT(IN) :: vIdx, hIdx, tIdx
    INTEGER(i_kind), INTENT(IN), OPTIONAL :: iv
    
    REAL(r_kind) :: bt ! For each instrument, each channel
    CHARACTER(LEN=*) :: varnametmp

    IF (iv < 1) RETURN
    ! print *,'GetForwardValue ',vIdx, hIdx, tIdx,iv,X%getVarIdx('tbb'),TRIM(varnametmp),TRIM(this%Source_of_bkg)
    IF (X%getVarIdx('tbb') .NE. 0) THEN
      bt = X%Fields(X%getVarIdx(TRIM(varnametmp)))%data(vIdx, hIdx, tIdx)
    ELSE
      IF  ( TRIM(varnametmp) .EQ. 'tbb' ) THEN 
        bt = this%tb_bakAtmodel(hIdx, tIdx)
      ELSE
        bt = this%obsData_avg
      END IF
      !  print *, 'bt1 = ', bt
    END IF

  END FUNCTION

  FUNCTION Field2ObsIsExisted(X, varnametmp)
    TYPE(State_t) :: X
    LOGICAL :: Field2ObsIsExisted
    CHARACTER(LEN=*) :: varnametmp

    Field2ObsIsExisted = .FALSE.

    ! Add more mandatory variables here
    IF ((X%getVarIdx(TRIM('temp')) .NE. 0) .AND. &
        (X%getVarIdx(TRIM('qvapor')) .NE. 0)) Field2ObsIsExisted = .TRUE.

    IF ((X%getVarIdx(TRIM('tbb')) .NE. 0)) Field2ObsIsExisted = .TRUE.

  END FUNCTION

  SUBROUTINE EachChannel_setup(this, configFile, idxFile)

    IMPLICIT NONE
    CLASS(EachChannel_t) :: this
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    INTEGER, OPTIONAL :: idxFile
    INTEGER(i_kind) :: istatus

    this%configFile = configFile
    this%numVars = 1
    istatus = yaml_get_var(TRIM(this%configFile), 'obs_thinning', 'interpolation', this%interpolation)
    ! print *, 'check interpolation scheme: ', this%interpolation

  END SUBROUTINE EachChannel_setup

  SUBROUTINE EachChannel_read(this, X)
    IMPLICIT NONE
    CLASS(EachChannel_t) :: this
    TYPE(State_t) :: X

  END SUBROUTINE EachChannel_read

  FUNCTION EachChannel_fwd(this, X, O) RESULT(Y)
    IMPLICIT NONE
    CLASS(EachChannel_t) :: this
    TYPE(State_t) :: X
    TYPE(ObsSet_t), INTENT(IN) :: O ! Telling the operator to map X on to this obs structure
    TYPE(ObsSet_t) :: Y

  END FUNCTION EachChannel_fwd

  FUNCTION EachChannel_tgt(this, dX, X) RESULT(Y)
    IMPLICIT NONE
    CLASS(EachChannel_t) :: this
    TYPE(State_t) :: dX, X
    TYPE(ObsSet_t) :: Y

  END FUNCTION EachChannel_tgt

  FUNCTION EachChannel_adj(this, dY, X) RESULT(dX)
    CLASS(EachChannel_t) :: this
    TYPE(State_t) :: dX, X
    TYPE(ObsSet_t) :: dY

  END FUNCTION EachChannel_adj

  SUBROUTINE destroy_EachChannel(this)
    IMPLICIT NONE
    CLASS(EachChannel_t) :: this

    IF (ALLOCATED(this%tb_bakAtmodel)) DEALLOCATE(this%tb_bakAtmodel)
    IF (ALLOCATED(this%Thinned_Angles)) DEALLOCATE (this%Thinned_Angles)
    IF (ALLOCATED(this%Thinned_Angles_hidx)) DEALLOCATE(this%Thinned_Angles_hidx)
    IF (ALLOCATED(this%Thinned_Angles_tidx)) DEALLOCATE(this%Thinned_Angles_tidx)

  END SUBROUTINE destroy_EachChannel

END MODULE EachChannel_m
