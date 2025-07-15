! Description:
!> @file
!!   Subroutines for thread-safe management of logical units
!
!> @brief
!!   Subroutines for thread-safe management of logical units
!!
!
! Copyright:
!    This software was developed within the context of
!    the EUMETSAT Satellite Application Facility on
!    Numerical Weather Prediction (NWP SAF), under the
!    Cooperation Agreement dated 7 December 2016, between
!    EUMETSAT and the Met Office, UK, by one or more partners
!    within the NWP SAF. The partners in the NWP SAF are
!    the Met Office, ECMWF, DWD and MeteoFrance.
!
!    Copyright 2017, EUMETSAT, All Rights Reserved.
!
MODULE rttov_lun
!$ use omp_lib

  USE parkind1,       ONLY : jpim, jplm
  USE rttov_unix_env, ONLY : rttov_exit

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: rttov_get_lun, rttov_put_lun

  LOGICAL(KIND=jplm) :: lunf(10:25)
  LOGICAL(KIND=jplm) :: init = .FALSE.

  INTEGER(KIND=jpim), PARAMETER :: lun_min = LBOUND(lunf,1), &
                                   lun_max = UBOUND(lunf,1)

  CONTAINS

  SUBROUTINE rttov_dump_lun()
    INTEGER(KIND=jpim) :: lun
    !$write(0,'(I3," ")',advance='no') omp_get_thread_num ()
    DO lun = lun_min, lun_max
      WRITE(0, '(L1)', ADVANCE='no') lunf(lun)
    ENDDO
    WRITE(0, *)
  END SUBROUTINE

  SUBROUTINE rttov_get_lun(lun_out)
    INTEGER(KIND=jpim), INTENT(OUT) :: lun_out
    INTEGER(KIND=jpim) :: lun

    lun_out = -1
!$OMP CRITICAL (rttov_lun_lock)
!$OMP FLUSH(lunf)
    IF (.NOT. init) THEN
      lunf = .TRUE.
      init = .TRUE.
    ENDIF
    DO lun = lun_min, lun_max
      IF (lun_out <= 0 .AND. lunf(lun)) THEN
        lunf(lun) = .FALSE.
        lun_out = lun
      ENDIF
    ENDDO
!$OMP FLUSH(lunf)
!$OMP END CRITICAL (rttov_lun_lock)
    ! write(0,*) " get = ", lun_out, omp_get_thread_num ()
  END SUBROUTINE

  SUBROUTINE rttov_put_lun(lun)
    INTEGER(KIND=jpim), intent(in) :: lun
!$OMP CRITICAL (rttov_lun_lock)
!$OMP FLUSH(lunf)
    IF (lunf(lun)) THEN
      WRITE(0, *) lun, "was already free"
      CALL rttov_exit(1_jpim)
    ENDIF
    lunf(lun) = .TRUE.
!$OMP FLUSH(lunf)
!$OMP END CRITICAL (rttov_lun_lock)
    ! write(0,*) " put = ", lun, omp_get_thread_num ()
  END SUBROUTINE

END MODULE
