!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA.AuxType
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2021/4/25, @GBA-MWF, Shenzhen
!!-------------------------------------------------------------------------------------------------

!> @brief
!! This type enclosed the essential parameters in constructing types.
MODULE AuxTypeSG_m
  USE geometry_m, ONLY: geometry_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE mpddinfo_sg_m, ONLY: mpddinfo_sg_t
  USE kinds_m, ONLY: i_kind, r_kind

  TYPE :: AuxTypeSG_t
    TYPE(SingleGrid_t), POINTER :: sg
    TYPE(mpddGlob_t), POINTER :: mpddGlob
    TYPE(mpddinfo_sg_t), POINTER :: mpddSub
  CONTAINS
    PROCEDURE, PUBLIC :: aux_initialize
    FINAL :: destructor
    PROCEDURE, PUBLIC, PASS :: EqualTo
    GENERIC :: OPERATOR(.EQ.) => EqualTo
  END TYPE AuxTypeSG_t

CONTAINS

  LOGICAL FUNCTION EqualTo(x1, x2)
    CLASS(AuxTypeSG_t), INTENT(IN) :: x1
    TYPE(AuxTypeSG_t), INTENT(IN) :: x2

    EqualTo = .FALSE.
    IF (x1%sg%id .EQ. x2%sg%id) EqualTo = .TRUE.
  END FUNCTION EqualTo

  SUBROUTINE aux_initialize(this, sg)
    IMPLICIT NONE
    CLASS(AuxTypeSG_t) :: this
    TYPE(SingleGrid_t), TARGET, INTENT(IN) :: sg

    this%sg => sg
    this%mpddGlob => sg%mpddGlob
    this%mpddSub => sg%mpddInfo_sg
  END SUBROUTINE aux_initialize

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(AuxTypeSG_t), INTENT(INOUT) :: this

  END SUBROUTINE destructor
END MODULE AuxTypeSG_m
