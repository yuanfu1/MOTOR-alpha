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

    INTEGER(i_kind) :: numSets, numFields, i, indxBeg, indxEnd, num

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
!         WRITE (*, 1) indxBeg, indxEnd, i, ObsSetAll%mpObs%myrank
! 1       FORMAT('Obs concat - indxBeg/indxEnd: ', 2I6, ' at obs Set: ', I2, ' at proc: ', I2)
        ObsSetAll%ObsFields(indxBeg:indxEnd) = ObsSetList(i)%ObsFields(:)
        indxBeg = indxEnd + 1
      END IF
    END DO

  END SUBROUTINE ObsConcat_s

END MODULE ObsUtilities_m
