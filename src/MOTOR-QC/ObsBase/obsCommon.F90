!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-QC.obsCommon.F90
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           : 2022-11-07   Created by Yuanfu Xie
!
!   Created by Yuanfu Xie (xieyf@gbamwf.com), 2022/11/07, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module provides an abstract data structure of general observation data.
MODULE obsCommon_m
  USE kinds_m, ONLY: i_kind, r_kind, r_double

  TYPE :: obsCommon_t
    CHARACTER(LEN=8), ALLOCATABLE :: varObsSpace(:)
  CONTAINS
    PROCEDURE, PUBLIC :: initialize
    PROCEDURE, PUBLIC :: obs2obsSpace
  END TYPE obsCommon_t

CONTAINS
  SUBROUTINE initialize(this, configFile)
    USE YAMLRead_m, ONLY: yaml_get_var
    IMPLICIT NONE
    CLASS(obsCommon_t) :: this
    CHARACTER(LEN=1024), INTENT(IN) :: configFile

    ! Local variables:
    INTEGER(i_kind) :: i, i_status

    i_status = yaml_get_var(TRIM(configFile), 'obs_thinning', 'varObsSpace', this%varObsSpace)
    IF (i_status .NE. 0) THEN
      WRITE (*, 1) i_status
1     FORMAT('obsCommon: There is no obs space defined, use default')
    ELSE
      WRITE (*, 2) (TRIM(this%varObsSpace(i)), i=1, UBOUND(this%varObsSpace, 1))
2     FORMAT('obsCommon: vars in obs space: ', 50(1X, A))
    END IF

  END SUBROUTINE initialize

  ! Converting observations to var in obs space:
  SUBROUTINE obs2obsSpace(this, num, obs, obsSpace, var)
    IMPLICIT NONE
    CLASS(obsCommon_t) :: this
    CHARACTER(*), INTENT(IN) :: obs
    CHARACTER(*), INTENT(IN) :: obsSpace
    INTEGER(i_kind), INTENT(IN) :: num
    REAL(r_kind), INTENT(INOUT) :: var(num)

    ! Converting:
    SELECT CASE (TRIM(obs))
    CASE ('pres')
      IF (TRIM(obsSpace) .EQ. 'lnp') THEN
        PRINT *, 'The pres is converted to lnp in the obs2obsSpace.'
        var = LOG(var)
      END IF
    END SELECT
  END SUBROUTINE obs2obsSpace
END MODULE obsCommon_m
