!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-QC.ObsUtilities
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Yongjian Huang
! VERSION           : V 0.0
! HISTORY           :
!   Created by Yongjian Huang, 2024/10/31, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

MODULE ObsOutlierDetector_m
  
  USE ObsBase_m, ONLY: ObsBase_t
  USE kinds_m, ONLY: r_kind, i_kind
  USE State_m, ONLY: State_t
  USE ObsSet_m, ONLY: ObsSet_t


  INTEGER(i_kind), PARAMETER :: TYPE_LEN = 20                         !> len of the detector type identifier
  
  TYPE, PUBLIC, ABSTRACT :: ObsOutlierDetector_t
    CHARACTER(len=TYPE_LEN):: det_type

    CONTAINS
      PROCEDURE(detect_outlier), DEFERRED :: detect
      
  END TYPE


  INTERFACE
  SUBROUTINE detect_outlier(this, outlier_type, obs_data, bck_at_obs, upper_bound, lower_bound, outlier_mask)
    IMPORT :: ObsOutlierDetector_t
    IMPORT :: ObsBase_t
    IMPORT  :: State_t
    IMPORT :: r_kind
    IMPORT :: i_kind

    CLASS(ObsOutlierDetector_t)       :: this
    CHARACTER(len=*) , INTENT(IN)         :: outlier_type
    REAL(r_kind),    INTENT(IN)         :: obs_data(:,:)
    REAL(r_kind),    INTENT(IN)         :: bck_at_obs(:,:)
    REAL(r_kind) , INTENT(IN)         :: upper_bound
    REAL(r_kind) , INTENT(IN)         :: lower_bound
    INTEGER(i_kind), ALLOCATABLE, INTENT(INOUT)     :: outlier_mask( :) 

  END SUBROUTINE detect_outlier


END INTERFACE

 




END MODULE ObsOutlierDetector_m