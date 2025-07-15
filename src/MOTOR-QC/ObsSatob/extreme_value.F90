
MODULE extreme_value_m
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Description:
! Define extreme value
!  extreme value are assigned directly,
!  HISTORY           : Origionally from RCNMP
!   Modified by Yuanfu Xie 2022-11-27 for covering printing info contained
!     in the temp_extreme_values.list file. This file is from GRAPES MESO
!   Modified by Yuanfu Xie 2024-04-16 for changing the file name from .dat
!     to .list as our gitignore sets .dat to exclude this file.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  USE Meteoro_constants

  IMPLICIT NONE

  !TYPE extreme_t
  REAL :: h_min_mid_sus(maxstl), h_min_mid_err(maxstl), h_max_mid_sus(maxstl), h_max_mid_err(maxstl), &
          h_min_hig_sus(maxstl), h_min_hig_err(maxstl), h_max_hig_sus(maxstl), h_max_hig_err(maxstl), &
          v_max_sus(maxstl), T_min_mid_sus(maxstl), T_min_mid_err(maxstl), T_max_mid_sus(maxstl), &
          T_max_mid_err(maxstl), T_min_hig_sus(maxstl), T_min_hig_err(maxstl), T_max_hig_sus(maxstl), &
          T_max_hig_err(maxstl)
  REAL   :: delta(2, 2, 3, maxstl)
  REAL   :: lapse(maxstl)
  REAL   :: dha(maxstl)
  REAL   :: rms_z(5, maxstl), rms_u(5, maxstl), rms_v(5, maxstl), rms_t(5, maxstl), rms_rh(5, maxstl)

  !CONTAINS
  !   PROCEDURE, PUBLIC  :: read_temp_extreme
  !   PROCEDURE, PUBLIC  :: read_temp_inverT
  !   PROCEDURE, PUBLIC  :: read_temp_lapse_rate
  !   PROCEDURE, PUBLIC  :: read_temp_hydro_cirt
  !   PROCEDURE, PUBLIC  :: read_temp_rms
  !END TYPE extreme_t
!-------------------------
!read the extreme value form temp_extreme_value.dat
!----------------------------
CONTAINS

  SUBROUTINE read_temp_extreme
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! PURPOSE:
! History:
! Version   Date         Comment
! 1.0       20121009     coded by tianwh
!
!-++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    IMPLICIT NONE
!CLASS(extreme_t) :: this
    CHARACTER(len=200)     :: filename, satinfo_file

    ! filename='/sources/data/static/'//'temp_extreme_value.list'
    CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", satinfo_file)
    filename = TRIM(satinfo_file)//'/'//'temp_extreme_value.list'
    PRINT *, 'temp_extreme:', filename
    OPEN (12, file=filename, STATUS='OLD')
    !       ACTION ='READ',      &
    !     blocksize=1048576,buffered='yes',buffercount=1 )

    READ (12, *) h_min_mid_sus(:)
    READ (12, *) h_min_mid_err(:)
    READ (12, *) h_max_mid_sus(:)
    READ (12, *) h_max_mid_err(:)
    READ (12, *) h_min_hig_sus(:)
    READ (12, *) h_min_hig_err(:)
    READ (12, *) h_max_hig_sus(:)
    READ (12, *) h_max_hig_err(:)
    READ (12, *) v_max_sus(:)
    READ (12, *) T_min_mid_sus(:)
    READ (12, *) T_min_mid_err(:)
    READ (12, *) T_max_mid_sus(:)
    READ (12, *) T_max_mid_err(:)
    READ (12, *) T_min_hig_sus(:)
    READ (12, *) T_min_hig_err(:)
    READ (12, *) T_max_hig_sus(:)
    READ (12, *) T_max_hig_err(:)
    CLOSE (12)
    ! Yuanfu Xie temporarily covers the following printings as it contains info in the temp_extreme_value.list
    !print *,h_min_mid_sus
    !print *,h_min_mid_err
    RETURN
    ! 99  WRITE (0,'(A)') ' ERROR READING temp_extreme_value.list in extreme_value.F'
    !       stop

  END SUBROUTINE read_temp_extreme

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  SUBROUTINE read_temp_inverT
!------------------------------
! read out the invertion criteration
!-------------------------------
    IMPLICIT NONE
    !CLASS(extreme_t) :: this
    INTEGER :: i, j, k, l
    CHARACTER(len=200)    ::   filename

!--------------------------------------------
! read out the criterion from the data file
!-------------------------------------------
    filename = '/sources/data/static/'//'temp_inverT_value.dat'
    PRINT *, 'temp_inverT:', filename
    ! open(11,file=filename, STATUS = 'OLD',ACTION= 'READ',&
    !          ERR = 9000,blocksize=1048576,buffered='yes',buffercount=1)
    OPEN (11, file=filename, STATUS='OLD')
    DO i = 1, 2
      DO j = 1, 2
        DO k = 1, 3
          READ (11, *, ERR=9000) delta(i, j, k, :)
        END DO
      END DO
    END DO
    CLOSE (11)
    RETURN
9000 WRITE (0, '(A)') ' ERROR READING temp_inverT_value.dat'
    STOP
  END SUBROUTINE read_temp_inverT

!+++++++++++++++++++++++++++++++++++++++++++++++++++++
  SUBROUTINE read_temp_lapse_rate
    IMPLICIT NONE
    !CLASS(extreme_t) :: this
    CHARACTER(len=200)    ::   filename
!------------------------------
! read out the invertion criteration
!-------------------------------
    filename = '/sources/data/static/'//'temp_lapse_rate.dat'
    PRINT *, 'temp_lapse_rate:', filename
!------------------------------------------------

!open(12,file=filename, STATUS = 'OLD',&
!        ACTION= 'READ',ERR =9000,   &
!    blocksize=1048576,buffered='yes',buffercount=1 )
    OPEN (12, file=filename, STATUS='OLD')
    READ (12, *, ERR=9000) lapse(:)
    CLOSE (12)
    RETURN
9000 WRITE (0, '(A)') ' ERROR READING temp_lapse_rate.dat'
    STOP
  END SUBROUTINE read_temp_lapse_rate

!+++++++++++++++++++++++++++++++++++++++++++++++++++++
  SUBROUTINE read_temp_hydro_cirt
    IMPLICIT NONE
    !CLASS(extreme_t) :: this
    CHARACTER(len=200)    ::   filename
!------------------------------
! read out the invertion criteration
!-------------------------------
    filename = '/sources/data/static/'//'temp_hydro_dha.dat'
    PRINT *, 'temp_hydro_dha:', filename
!------------------------------------------------
!open(12,file='temp_hydro_dha.dat', STATUS = 'OLD',&
!            ACTION= 'READ',ERR =9000,  &
!     blocksize=1048576,buffered='yes',buffercount=1)
    OPEN (12, file='temp_hydro_dha.dat', STATUS='OLD')
    READ (12, *, ERR=9000) dha(:)
    CLOSE (12)
    RETURN
9000 WRITE (0, '(A)') ' ERROR  read_temp_hydro_cirt      '

    STOP
  END SUBROUTINE read_temp_hydro_cirt

!----------------------------------------------------------
  SUBROUTINE read_temp_rms
    IMPLICIT NONE
    !CLASS(extreme_t) :: this
    INTEGER :: i, k
    CHARACTER(len=200)    ::   filename
!------------------------------
! read out the error of the o-b
!-------------------------------
    filename = '/sources/data/static/'//'temp_first_guess.dat'
    PRINT *, 'temp_first_guess:', filename
    !------------------------------------------------
!open(12,file=filename, STATUS = 'OLD', &
!          ACTION= 'READ',ERR =9000 ,  &
!     blocksize=1048576,buffered='yes',buffercount=1)
    OPEN (12, file=filename, STATUS='OLD')

    DO i = 1, 5
      DO k = 1, maxstl
        READ (12, *, ERR=9000) rms_z(i, k), rms_u(i, k), rms_v(i, k), rms_t(i, k), rms_rh(i, k)
!   print*,i,k,rms_z(i,k),rms_u(i,k),rms_v(i,k),rms_t(i,k),rms_rh(i,k)
      END DO
    END DO

    CLOSE (12)
    RETURN
9000 WRITE (0, '(A)') ' ERROR READING temp_background_rms '
    STOP
  END SUBROUTINE read_temp_rms

END MODULE extreme_value_m
