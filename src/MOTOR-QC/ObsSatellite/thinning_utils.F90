!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-QC.Thinning_utils
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Yali Wu
! VERSION           : V 0.0
! HISTORY           :
!   Created by Yali Wu (wuyali@gbamwf.com), 2022/06/23, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module implements a satellite radiance data structure.
MODULE Thinning_utils_m
  USE kinds_m
  USE parameters_m
  USE EachChannel_m, ONLY: EachChannel_t
  USE ObsSet_m, ONLY: ObsSet_t
  USE Satellite_utils_m, ONLY: Convert_int2char, findhIdx_all

  IMPLICIT NONE

CONTAINS

  SUBROUTINE SetObsAttrSAT(YAngles, EachChannel)
    IMPLICIT NONE
    TYPE(ObsSet_t), INTENT(IN) :: YAngles
    TYPE(EachChannel_t), INTENT(INOUT) :: EachChannel
    INTEGER(i_kind) :: ifield, numAngles(4), numOUT, IdxOUT(1)

    ! Each ObsField has its ObsAttrSat
    IF (SIZE(YAngles%obsFields) .EQ. 0) RETURN

    DO ifield = 1, 4
      numAngles(ifield) = SIZE(YAngles%obsFields(ifield)%values)
      IF (ifield > 1) THEN
        IF (numAngles(ifield) .NE. numAngles(ifield - 1)) PRINT *, 'NUM of four satellite angles are different: ', numAngles(ifield), numAngles(ifield - 1)
      END IF
    END DO
    ! print *, 'check nums: ', numAngles
    numOUT = MAXVAL(numAngles)
    IdxOUT(1:1) = MAXLOC(numAngles)
    
    IF (.not. ALLOCATED(EachChannel%Thinned_Angles)) ALLOCATE(EachChannel%Thinned_Angles(4, numOUT))
    IF (.not. ALLOCATED(EachChannel%Thinned_Angles_hidx)) ALLOCATE(EachChannel%Thinned_Angles_hidx(numOUT))
    IF (.not. ALLOCATED(EachChannel%Thinned_Angles_tidx)) ALLOCATE(EachChannel%Thinned_Angles_tidx(numOUT))
    
    EachChannel%Thinned_Angles_hidx(:) = YAngles%obsFields(IdxOUT(1))%idx(:)%hIdx
    EachChannel%Thinned_Angles_tidx(:) = YAngles%obsFields(IdxOUT(1))%idx(:)%tIdx ! Four angles should have the same dimension sizes
    EachChannel%Thinned_Angles = missing
    DO ifield = 1, 4  
      EachChannel%Thinned_Angles(ifield, :) = YAngles%obsFields(ifield)%values(:)
    END DO 
 
  END SUBROUTINE SetObsAttrSAT

  SUBROUTINE SetObsAttrSAT_Y(Y, YAngles, X)
    USE State_m, ONLY: State_t
    IMPLICIT NONE
    TYPE(ObsSet_t), INTENT(INOUT) :: Y
    TYPE(ObsSet_t), INTENT(IN) :: YAngles
    TYPE(State_t), INTENT(IN)  :: X
    ! INTEGER(i_kind), INTENT(IN) :: tSlots
    INTEGER(i_kind) :: ichan
    CHARACTER(len=10) :: chan_num
    CHARACTER(len=100) :: new_obstype
    INTEGER(i_kind) :: ih, it, i, iobs, nobs, obsIdx
    INTEGER(i_kind), ALLOCATABLE :: obshIdx(:)
    INTEGER(i_kind), ALLOCATABLE :: IdxArray(:)

    ! Each ObsField has its ObsAttrSat
    ! Each channel has the same ObsAttrSat
    ! For each channel/ivar, all time slots
    
    ALLOCATE(obshIdx(X%sg%tSlots)) ! Max size is tSlots
    obshIdx = 0
    IdxArray = YAngles%obsFields(1)%idx(:)%hIdx

    DO ichan = 1, SIZE(Y%ObsFields)
      ! PRINT *, 'Check if TBB and angles match: ', SIZE(Y%obsFields(ichan)%values, 1), &
      !                                             SIZE(YAngles%obsFields(1)%values, 1)
      ! IF (SIZE(Y%obsFields(ichan)%values, 1) .EQ. 0 .OR. &
      !     SIZE(YAngles%obsFields(1)%values, 1) .EQ. 0) RETURN

      nobs = UBOUND(Y%ObsFields(ichan)%values, 1)
      ALLOCATE(Y%obsFields(ichan)%ObsAttrSat%zenangle(nobs))
      ALLOCATE(Y%obsFields(ichan)%ObsAttrSat%azangle(nobs))
      ALLOCATE(Y%obsFields(ichan)%ObsAttrSat%sunzenangle(nobs))
      ALLOCATE(Y%obsFields(ichan)%ObsAttrSat%sunazangle(nobs))
      ALLOCATE(Y%obsFields(ichan)%ObsAttrSat%latitude(nobs))
      ALLOCATE(Y%obsFields(ichan)%ObsAttrSat%longitude(nobs))
      
      DO iobs = 1, nobs

        ih = Y%ObsFields(ichan)%idx(iobs)%hIdx
        it = Y%ObsFields(ichan)%idx(iobs)%tIdx

        Y%obsFields(ichan)%ObsAttrSat%latitude(iobs) = X%sg%cell_cntr(1,ih)
        Y%obsFields(ichan)%ObsAttrSat%longitude(iobs) = X%sg%cell_cntr(2,ih)

        CALL findhIdx_all(IdxArray, ih, obshIdx)
        DO i = 1, SIZE(obshIdx, 1)
          IF (obshIdx(i) .EQ. 0) CYCLE
          IF (YAngles%obsFields(1)%idx(obshIdx(i))%tIdx .EQ. it) THEN
            obsIdx = obshIdx(i)
          END IF
        END DO

        IF (obsIdx == 0) RETURN
        Y%obsFields(ichan)%ObsAttrSat%zenangle(iobs) = YAngles%obsFields(1)%values(obsIdx)
        Y%obsFields(ichan)%ObsAttrSat%azangle(iobs) = YAngles%obsFields(2)%values(obsIdx)
        Y%obsFields(ichan)%ObsAttrSat%sunzenangle(iobs) = YAngles%obsFields(3)%values(obsIdx)
        Y%obsFields(ichan)%ObsAttrSat%sunazangle(iobs) = YAngles%obsFields(4)%values(obsIdx)

      END DO
      ! Set new_obstype, for output diag convenience
      chan_num = Convert_int2char(ichan)
      new_obstype = TRIM(Y%obsFields(ichan)%Get_ObsType())//'_ch_'//ADJUSTL(chan_num)
      CALL Y%ObsFields(ichan)%Set_ObsType(ADJUSTL(new_obstype))
      ! PRINT *, 'Reset obstype: ', adjustL(new_obstype)
      WHERE (Y%obsFields(ichan)%values > 999.0 .OR. Y%obsFields(ichan)%values < -999.0 ) Y%obsFields(ichan)%values = missing
      PRINT *, 'thinning_utils: Num of OBS of each channel: ', ichan, size(Y%obsFields(ichan)%ObsAttrSat%zenangle)
    END DO
    DEALLOCATE (IdxArray)
    DEALLOCATE (obshIdx)

  END SUBROUTINE SetObsAttrSAT_Y

END MODULE Thinning_utils_m
