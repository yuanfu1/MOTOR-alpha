!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-QC.ObsBase
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Jiongming Pang
! VERSION           : V 0.0
! HISTORY           : 2021-01-20   Created by Jiongming Pang
!
!   Created by Jiongming Pang (pang.j.m@hotmail.com), 2022/01/20, @GBA-MWF, Shenzhen
!
!   Modified by Yuanfu Xie (yuanfu_xie@yahoo.com), 2025/07/22, @GBA-MWF, Shenzhen
!     created a DeepCopyObsField subroutine to copy ObsField_t objects suggested by chatGPT to
!     avoid segmentation fault when concatenating ObsSet_t objects by copying ObsFields using
!     ObsSetAll%ObsFields(indxBeg:indxEnd) = ObsSetList(i)%ObsFields(:)
!     as the idx array is not allocated in ObsSetList(i)%ObsFields(:)
!     and the idx array is not allocated in ObsSetAll%ObsFields(indxBeg:indxEnd)
!
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module provides an abstract data structure of general observation data.
MODULE ObsUtilities_m
  USE kinds_m, ONLY: i_kind, r_kind, r_double
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE State_m, ONLY: State_t
  USE ObsSet_m, ONLY: ObsSet_t
  USE mpObs_m, ONLY: mpObs_t
  USE geoTools_m, ONLY: GeoBox_t
  USE slint, ONLY: slint_init, tgt_grid
  USE ObsField_m, ONLY: ObsField_t
  USE AdvanceTime_m

  IMPLICIT NONE

CONTAINS

  SUBROUTINE ObsConcat_s(ObsSetList, ObsSetAll)
    IMPLICIT NONE
    TYPE(ObsSet_t), INTENT(IN)  :: ObsSetList(:)
    TYPE(ObsSet_t), INTENT(INOUT) :: ObsSetAll

    INTEGER(i_kind) :: numSets, numFields, i, j, indxBeg, indxEnd, num

    numSets = UBOUND(ObsSetList, 1)

    IF (numSets .EQ. 0) THEN
      STOP "ERROR: There is no available ObsSet of ObsSetList in sub ObsConcat_s!"
    END IF

    numFields = 0
    DO i = 1, numSets
      IF (ALLOCATED(ObsSetList(i)%ObsFields)) THEN
        num = UBOUND(ObsSetList(i)%ObsFields, 1)
        numFields = numFields + num
      END IF
    END DO

    ALLOCATE (ObsSetAll%ObsFields(numFields))
    indxBeg = 1
    DO i = 1, numSets
      IF (ALLOCATED(ObsSetList(i)%ObsFields)) THEN
        num = UBOUND(ObsSetList(i)%ObsFields, 1)
        IF (num .EQ. 0) CYCLE
        indxEnd = indxBeg + num - 1
        ! Yuanfu Xie added this line to deep copy ObsField_t objects and turned off the following assignment: 2025/07/22
        ! ObsSetAll%ObsFields(indxBeg:indxEnd) = ObsSetList(i)%ObsFields(:)
        DO j=indxBeg, indxEnd
          CALL DeepCopyObsField(ObsSetAll%ObsFields(j), ObsSetList(i)%ObsFields(j-indxBeg+1))
        END DO
        indxBeg = indxEnd + 1
      END IF
    END DO
  END SUBROUTINE ObsConcat_s

  ! Yuanfu Xie added this subroutine to deep copy ObsField_t objects by modifying the suggested code from chatGPT
  SUBROUTINE DeepCopyObsField(dest, src)
    TYPE(ObsField_t), INTENT(OUT) :: dest
    TYPE(ObsField_t), INTENT(IN)  :: src
    dest%name = src%name
    dest%valType = src%valType
    dest%obsType = src%obsType
    dest%mpObs => src%mpObs
    IF (ALLOCATED(src%idx)) THEN
      ALLOCATE(dest%idx(SIZE(src%idx)))
      ALLOCATE(dest%values(SIZE(src%values)))
      dest%idx = src%idx
      dest%values = src%values
    END IF
    IF (ALLOCATED(src%valueArray)) THEN
      ALLOCATE(dest%valueArray(SIZE(src%valueArray,1), SIZE(src%valueArray,2)))
      dest%valueArray = src%valueArray
    END IF
  END SUBROUTINE

END MODULE ObsUtilities_m
