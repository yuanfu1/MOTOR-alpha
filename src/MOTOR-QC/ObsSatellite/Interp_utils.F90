!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-QC.Interp_utils
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Yali Wu
! VERSION           : V 0.0
! HISTORY           :
!   Created by Yali Wu (wuyali@gbamwf.com), 2022/06/23, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!! RadBase module implements a satellite radiance data structure.
MODULE Interp_utils_m
  USE kinds_m
  USE parameters_m
  USE State_m, ONLY: State_t
  USE EachChannel_m, ONLY: EachChannel_t
  USE Satellite_QC_utils_m
  USE Satellite_utils_m
  USE rttov_nl_sp_m, ONLY: rttov_nl_sp_t

CONTAINS

SUBROUTINE TB_bkgAtmodel(X, OprRTTOV, rttov_chan, EachChannel)
  TYPE(State_t), INTENT(IN) :: X
  TYPE(rttov_nl_sp_t), INTENT(IN) :: OprRTTOV
  INTEGER(i_kind), INTENT(IN) :: rttov_chan
  TYPE(EachChannel_t), INTENT(INOUT) :: EachChannel
  INTEGER(i_kind) :: i, hIdx, tIdx
  REAL(r_kind) :: forward
  INTEGER(i_kind) :: locX(2)
  REAL(r_kind) :: SatAngles(4)
  INTEGER(i_kind), ALLOCATABLE :: obsIdx(:)
  INTEGER(i_kind) :: IdxOUT=0
  REAL(r_kind), ALLOCATABLE :: tb_bakAtmodel(:,:,:)

  ! EachChannel%tb_bakAtmodel = EachChannel%obsData_avg
  EachChannel%tb_bakAtmodel = missing
  ALLOCATE(obsIdx(X%sg%tSlots)) ! Max size is tSlots

  DO tIdx = 1, X%sg%tSlots
    DO hIdx = 1, X%sg%num_icell
    ! DO hIdx = 1, X%sg%num_cell ! either num_cell or num_icell is ok

      locX = (/hIdx, tIdx/)
      obsIdx = 0
      CALL findhIdx_ht(EachChannel%Thinned_Angles_hidx, EachChannel%Thinned_Angles_tidx, hIdx, tIdx, obsIdx(1))
      ! CALL findhIdx_all(EachChannel%Thinned_Angles_hidx, hIdx, obsIdx)
      ! print *, 'check sizes ', SIZE(EachChannel%Thinned_Angles,2), obsIdx

      IdxOUT = obsIdx(1)
      ! IdxOUT = 0
      ! DO i = 1, SIZE(obsIdx, 1)
      !   IF (obsIdx(i) .NE. 0) THEN
      !     IF (EachChannel%Thinned_Angles_tidx(obsIdx(i)) .EQ. tIdx) IdxOUT = obsIdx(i)
      !   END IF
      ! END DO
      ! PRINT *, 'obsIdx = ', obsIdx(1), 'IdxOUT = ', IdxOUT
      ! IF (SIZE(EachChannel%Thinned_Angles,2) .LT. IdxOUT) PRINT *, 'IdxOUT=',IdxOUT, SIZE(EachChannel%Thinned_Angles,2)
      IF (IdxOUT > 0) THEN
        IF (SIZE(EachChannel%Thinned_Angles, 2) .LT. IdxOUT) THEN
          ! forward = EachChannel%obsData_avg
          forward = missing
        ELSE
          SatAngles = EachChannel%Thinned_Angles(1:4, IdxOUT)
          CALL OprRTTOV%rttov_nl_sp_simobs(X, rttov_chan, locX, SatAngles, forward)
        END IF
        EachChannel%tb_bakAtmodel(hIdx, tIdx) = forward
      END IF

    END DO
  END DO
  DEALLOCATE (obsIdx)

  ! Exchange halo to avoid parallel blocks
  ! print *, 'check forward: ', maxval(EachChannel%tb_bakAtmodel), minval(EachChannel%tb_bakAtmodel)
  ALLOCATE(tb_bakAtmodel(X%sg%vLevel, X%sg%num_cell, X%sg%tSlots))
  DO k = 1, X%sg%vLevel
    tb_bakAtmodel(k, :, :) = EachChannel%tb_bakAtmodel(:,:)
  END DO
  CALL X%sg%ExchangeMatOnHaloForFieldGrid(X%sg%tSlots, X%sg%vLevel, tb_bakAtmodel)
  EachChannel%tb_bakAtmodel(:, :) = tb_bakAtmodel(1, :,:)
  DEALLOCATE(tb_bakAtmodel)

  print *, 'finish TB_bkgAtmodel calculation'

END SUBROUTINE TB_bkgAtmodel

END MODULE Interp_utils_m
