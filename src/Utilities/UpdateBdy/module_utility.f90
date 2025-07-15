! Developed by CMA GRAPES_MESO 5.0
! Modified by YaliWu for CMA-GD intaerfaces 20220720

MODULE module_utility

  USE module_variables

CONTAINS
  SUBROUTINE get_ijk(s_we, e_we, s_sn, e_sn, s_vert, e_vert, &
                     ids, ide, jds, jde, kds, kde, &
                     ims, ime, jms, jme, kms, kme, &
                     its, ite, jts, jte, kts, kte)

    IMPLICIT NONE

    INTEGER, INTENT(IN)  :: s_we, e_we, s_sn, e_sn, s_vert, e_vert

    INTEGER, INTENT(OUT) :: ids, ide, jds, jde, kds, kde, &
                            ims, ime, jms, jme, kms, kme, &
                            its, ite, jts, jte, kts, kte

    ids = s_we
    ide = e_we
    jds = s_sn
    jde = e_sn
    kds = s_vert
    kde = e_vert + 1
    ims = s_we
    ime = e_we
    jms = s_sn
    jme = e_sn
    kms = s_vert - 1
    kme = e_vert + 1
    its = s_we
    ite = e_we
    jts = s_sn
    jte = e_sn
    kts = s_vert
    kte = e_vert
    PRINT *, 'ids = ', ids
    PRINT *, 'ide = ', ide
    PRINT *, 'jds = ', jds
    PRINT *, 'jde = ', jde
    PRINT *, 'kts = ', kts
    PRINT *, 'kte = ', kte
    PRINT *, 'kds = ', kds
    PRINT *, 'kde = ', kde
    PRINT *, 'kms = ', kms
    PRINT *, 'kme = ', kme
  END SUBROUTINE get_ijk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE alloc_space_field( &
    ids, ide, jds, jde, kds, kde, &
    ims, ime, jms, jme, kms, kme, &
    its, ite, jts, jte, kts, kte, spec_bdy_width)

    IMPLICIT NONE

    !  Input data.

    INTEGER, INTENT(IN)            :: ids, ide, jds, jde, kds, kde
    INTEGER, INTENT(IN)            :: ims, ime, jms, jme, kms, kme
    INTEGER, INTENT(IN)            :: its, ite, jts, jte, kts, kte
    INTEGER, INTENT(IN)            :: spec_bdy_width

    !  Local data.
    REAL    :: initial_data_value
    INTEGER :: ierr

    initial_data_value = 0.

    ALLOCATE (pi(ims:ime, jms:jme, kms:kme), STAT=ierr)
    IF (ierr .NE. 0) THEN
      STOP 'allocate pi'
    END IF
    pi = initial_data_value
    ALLOCATE (th(ims:ime, jms:jme, kms:kme), STAT=ierr)
    IF (ierr .NE. 0) THEN
      STOP 'allocate th'
    END IF
    th = initial_data_value
    ALLOCATE (u(ims:ime, jms:jme, kms:kme), STAT=ierr)
    IF (ierr .NE. 0) THEN
      STOP 'allocate u'
    END IF
    u = initial_data_value
    ALLOCATE (v(ims:ime, jms:jme, kms:kme), STAT=ierr)
    IF (ierr .NE. 0) THEN
      STOP 'allocate v'
    END IF
    v = initial_data_value
    ALLOCATE (q(ims:ime, jms:jme, kms:kme), STAT=ierr)
    IF (ierr .NE. 0) THEN
      STOP 'allocate q'
    END IF
    q = initial_data_value

!
    ALLOCATE (pi_b(MAX(ide, jde), (kds - 1):kde, spec_bdy_width, 4), STAT=ierr)
    ALLOCATE (pi_b_old(MAX(ide, jde), (kds - 1):kde, spec_bdy_width, 4), STAT=ierr)
    IF (ierr .NE. 0) THEN
      STOP 'allocate pi_b'
    END IF
    pi_b = initial_data_value
    ALLOCATE (pi_b_1(MAX(ide, jde), (kds - 1):kde, spec_bdy_width, 4), STAT=ierr)
    IF (ierr .NE. 0) THEN
      STOP 'allocate pi_b_1'
    END IF
    pi_b_1 = initial_data_value
    ALLOCATE (pi_b_2(MAX(ide, jde), (kds - 1):kde, spec_bdy_width, 4), STAT=ierr)
    IF (ierr .NE. 0) THEN
      STOP 'allocate pi_b_2'
    END IF
    pi_b_2 = initial_data_value

    ALLOCATE (th_b(MAX(ide, jde), (kds - 1):kde, spec_bdy_width, 4), STAT=ierr)
    ALLOCATE (th_b_old(MAX(ide, jde), (kds - 1):kde, spec_bdy_width, 4), STAT=ierr)
    IF (ierr .NE. 0) THEN
      STOP 'allocate th_b'
    END IF
    th_b = initial_data_value
    ALLOCATE (th_b_1(MAX(ide, jde), (kds - 1):kde, spec_bdy_width, 4), STAT=ierr)
    IF (ierr .NE. 0) THEN
      STOP 'allocate th_b_1'
    END IF
    th_b_1 = initial_data_value
    ALLOCATE (th_b_2(MAX(ide, jde), (kds - 1):kde, spec_bdy_width, 4), STAT=ierr)
    IF (ierr .NE. 0) THEN
      STOP 'allocate th_b_2'
    END IF
    th_b_2 = initial_data_value

    ALLOCATE (u_b(MAX(ide, jde), (kds - 1):kde, spec_bdy_width, 4), STAT=ierr)
    ALLOCATE (u_b_old(MAX(ide, jde), (kds - 1):kde, spec_bdy_width, 4), STAT=ierr)
    IF (ierr .NE. 0) THEN
      STOP 'allocate u_b'
    END IF
    u_b = initial_data_value
    ALLOCATE (u_b_1(MAX(ide, jde), (kds - 1):kde, spec_bdy_width, 4), STAT=ierr)
    IF (ierr .NE. 0) THEN
      STOP 'allocate u_b_1'
    END IF
    u_b_1 = initial_data_value
    ALLOCATE (u_b_2(MAX(ide, jde), (kds - 1):kde, spec_bdy_width, 4), STAT=ierr)
    IF (ierr .NE. 0) THEN
      STOP 'allocate u_b_2'
    END IF
    u_b_2 = initial_data_value

    ALLOCATE (v_b(MAX(ide, jde), (kds - 1):kde, spec_bdy_width, 4), STAT=ierr)
    ALLOCATE (v_b_old(MAX(ide, jde), (kds - 1):kde, spec_bdy_width, 4), STAT=ierr)
    IF (ierr .NE. 0) THEN
      STOP 'allocate v_b'
    END IF
    v_b = initial_data_value
    ALLOCATE (v_b_1(MAX(ide, jde), (kds - 1):kde, spec_bdy_width, 4), STAT=ierr)
    IF (ierr .NE. 0) THEN
      STOP 'allocate v_b_1'
    END IF
    v_b_1 = initial_data_value
    ALLOCATE (v_b_2(MAX(ide, jde), (kds - 1):kde, spec_bdy_width, 4), STAT=ierr)
    IF (ierr .NE. 0) THEN
      STOP 'allocate v_b_2'
    END IF
    v_b_2 = initial_data_value

    ALLOCATE (q_b(MAX(ide, jde), (kds - 1):kde, spec_bdy_width, 4), STAT=ierr)
    ALLOCATE (q_b_old(MAX(ide, jde), (kds - 1):kde, spec_bdy_width, 4), STAT=ierr)
    IF (ierr .NE. 0) THEN
      STOP 'allocate q_b'
    END IF
    q_b = initial_data_value
    ALLOCATE (q_b_1(MAX(ide, jde), (kds - 1):kde, spec_bdy_width, 4), STAT=ierr)
    IF (ierr .NE. 0) THEN
      STOP 'allocate q_b_1'
    END IF
    q_b_1 = initial_data_value
    ALLOCATE (q_b_2(MAX(ide, jde), (kds - 1):kde, spec_bdy_width, 4), STAT=ierr)
    IF (ierr .NE. 0) THEN
      STOP 'allocate q_b_2'
    END IF
    q_b_2 = initial_data_value

    ALLOCATE (pi_bt(MAX(ide, jde), (kds - 1):kde, spec_bdy_width, 4), STAT=ierr)
    ALLOCATE (pi_bt_old(MAX(ide, jde), (kds - 1):kde, spec_bdy_width, 4), STAT=ierr)
    IF (ierr .NE. 0) THEN
      STOP 'allocate pi_bt'
    END IF
    pi_bt = initial_data_value
    ALLOCATE (th_bt(MAX(ide, jde), (kds - 1):kde, spec_bdy_width, 4), STAT=ierr)
    ALLOCATE (th_bt_old(MAX(ide, jde), (kds - 1):kde, spec_bdy_width, 4), STAT=ierr)
    IF (ierr .NE. 0) THEN
      STOP 'allocate th_bt'
    END IF
    th_bt = initial_data_value
    ALLOCATE (u_bt(MAX(ide, jde), (kds - 1):kde, spec_bdy_width, 4), STAT=ierr)
    ALLOCATE (u_bt_old(MAX(ide, jde), (kds - 1):kde, spec_bdy_width, 4), STAT=ierr)
    IF (ierr .NE. 0) THEN
      STOP 'allocate u_bt'
    END IF
    u_bt = initial_data_value
    ALLOCATE (v_bt(MAX(ide, jde), (kds - 1):kde, spec_bdy_width, 4), STAT=ierr)
    ALLOCATE (v_bt_old(MAX(ide, jde), (kds - 1):kde, spec_bdy_width, 4), STAT=ierr)
    IF (ierr .NE. 0) THEN
      STOP 'allocate v_bt'
    END IF
    v_bt = initial_data_value
    ALLOCATE (q_bt(MAX(ide, jde), (kds - 1):kde, spec_bdy_width, 4), STAT=ierr)
    ALLOCATE (q_bt_old(MAX(ide, jde), (kds - 1):kde, spec_bdy_width, 4), STAT=ierr)
    IF (ierr .NE. 0) THEN
      STOP 'allocate q_bt'
    END IF
    q_bt = initial_data_value

    RETURN
  END SUBROUTINE alloc_space_field

  SUBROUTINE get_initial_value(pi, th, u, v, q, fid, &
                               ids, ide, jds, jde, kds, kde, &
                               ims, ime, jms, jme, kms, kme, &
                               its, ite, jts, jte, kts, kte)

    IMPLICIT NONE

    INTEGER, INTENT(in)    :: ids, ide, jds, jde, kds, kde, &
                              ims, ime, jms, jme, kms, kme, &
                              its, ite, jts, jte, kts, kte
    INTEGER, INTENT(in)     :: fid

    REAL, DIMENSION(ims:ime, jms:jme, kms:kme), INTENT(OUT)  :: pi, th, u, v, q

    ! local variables
    INTEGER                      :: i, j, k, fid2, kk
    INTEGER, DIMENSION(5)         :: mdate
    REAL, DIMENSION(ids:ide, jds:jde, kms:kme) :: globbuf

    READ (fid) mdate
    PRINT *, 'mdate = ', mdate
    !-------------------

    DO kk = 1, 11
      READ (fid) (((globbuf(i, j, k), i=ids, ide), k=kms, kme), j=jds, jde)
      PRINT *, MAXVAL(globbuf), MINVAL(globbuf)
    END DO

    READ (fid) (((pi(i, j, k), i=ids, ide), k=kms, kme), j=jds, jde)
    PRINT *, 'pi=', MAXVAL(pi), MINVAL(pi)

    READ (fid) (((u(i, j, k), i=ids, ide), k=kms, kme), j=jds, jde)
    READ (fid) (((v(i, j, k), i=ids, ide), k=kms, kme), j=jds, jde)
    PRINT *, 'u=', MAXVAL(u), MINVAL(u)
    PRINT *, 'v=', MAXVAL(v), MINVAL(v)

    DO kk = 1, 2
      READ (fid) (((globbuf(i, j, k), i=ids, ide), k=kms, kme), j=jds, jde)
      PRINT *, MAXVAL(globbuf), MINVAL(globbuf)
    END DO

    READ (fid) (((th(i, j, k), i=ids, ide), k=kms, kme), j=jds, jde)
    READ (fid) (((q(i, j, k), i=ids, ide), k=kms, kme), j=jds, jde)
    PRINT *, 'th=', MAXVAL(th), MINVAL(th)
    PRINT *, 'q=', MAXVAL(q), MINVAL(q)

    RETURN
  END SUBROUTINE get_initial_value

  SUBROUTINE split_date_char(date, century_year, month, day, hour, minute, second, ten_thousandth)

    IMPLICIT NONE

    !  Input data.

    CHARACTER(LEN=24), INTENT(IN) :: date

    !  Output data.

    INTEGER, INTENT(OUT) :: century_year, month, day, hour, minute, second, ten_thousandth

    READ (date, FMT='(    I4)') century_year
    READ (date, FMT='( 5X,I2)') month
    READ (date, FMT='( 8X,I2)') day
    READ (date, FMT='(11X,I2)') hour
    READ (date, FMT='(14X,I2)') minute
    READ (date, FMT='(17X,I2)') second
    READ (date, FMT='(20X,I4)') ten_thousandth

  END SUBROUTINE split_date_char

  SUBROUTINE read_initial_boundary(mdate, u_b, v_b, pi_b, th_b, q_b, &
                                   u_b_old, v_b_old, pi_b_old, th_b_old, q_b_old, &
                                   u_bt, v_bt, pi_bt, th_bt, q_bt, &
                                   interval_seconds, spec_bdy_width, fid, &
                                   ids, ide, jds, jde, kds, kde, &
                                   ims, ime, jms, jme, kms, kme, &
                                   its, ite, jts, jte, kts, kte)

    IMPLICIT NONE
    INTEGER, INTENT(IN) ::      ims, ime, jms, jme, kms, kme, &
                           ids, ide, jds, jde, kds, kde, &
                           its, ite, jts, jte, kts, kte

    ! Arguments
    INTEGER, INTENT(IN)            :: interval_seconds, spec_bdy_width, fid
    INTEGER, DIMENSION(5), INTENT(out)  :: mdate
    REAL, DIMENSION(MAX(ide, jde), (kds - 1):kde, spec_bdy_width, 4), INTENT(OUT) :: u_b, v_b, pi_b, th_b, q_b
    REAL, DIMENSION(MAX(ide, jde), (kds - 1):kde, spec_bdy_width, 4) :: u_b_old, v_b_old, pi_b_old, th_b_old, q_b_old
    !local variables
    REAL, DIMENSION(ims:ime, jms:jme, kms:kme) :: u, v, pi, th, q
    REAL, DIMENSION(MAX(ide, jde), (kds - 1):kde, spec_bdy_width, 4) :: u_bt, v_bt, pi_bt, th_bt, q_bt

    !
    INTEGER     :: i, j, k

    READ (fid) mdate
    PRINT *, 'read_initial_boundary: ', mdate
    READ (fid) (((pi_b(i, k, j, 1), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    READ (fid) (((pi_b(i, k, j, 2), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    READ (fid) (((pi_b(i, k, j, 3), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    READ (fid) (((pi_b(i, k, j, 4), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)

    READ (fid) (((pi_bt(i, k, j, 1), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    READ (fid) (((pi_bt(i, k, j, 2), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    READ (fid) (((pi_bt(i, k, j, 3), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    READ (fid) (((pi_bt(i, k, j, 4), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)

    READ (fid) (((th_b(i, k, j, 1), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    READ (fid) (((th_b(i, k, j, 2), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    READ (fid) (((th_b(i, k, j, 3), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    READ (fid) (((th_b(i, k, j, 4), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    PRINT *, 'check th_b ', MAXVAL(th_b), MINVAL(th_b)
    PRINT *, 'check pi_b ', MAXVAL(pi_b), MINVAL(pi_b)

    READ (fid) (((th_bt(i, k, j, 1), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    READ (fid) (((th_bt(i, k, j, 2), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    READ (fid) (((th_bt(i, k, j, 3), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    READ (fid) (((th_bt(i, k, j, 4), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)

    READ (fid) (((u_b(i, k, j, 1), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    READ (fid) (((u_b(i, k, j, 2), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    READ (fid) (((u_b(i, k, j, 3), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    READ (fid) (((u_b(i, k, j, 4), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)

    READ (fid) (((u_bt(i, k, j, 1), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    READ (fid) (((u_bt(i, k, j, 2), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    READ (fid) (((u_bt(i, k, j, 3), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    READ (fid) (((u_bt(i, k, j, 4), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)

    READ (fid) (((v_b(i, k, j, 1), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    READ (fid) (((v_b(i, k, j, 2), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    READ (fid) (((v_b(i, k, j, 3), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    READ (fid) (((v_b(i, k, j, 4), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)

    READ (fid) (((v_bt(i, k, j, 1), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    READ (fid) (((v_bt(i, k, j, 2), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    READ (fid) (((v_bt(i, k, j, 3), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    READ (fid) (((v_bt(i, k, j, 4), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)

    READ (fid) (((q_b(i, k, j, 1), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    READ (fid) (((q_b(i, k, j, 2), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    READ (fid) (((q_b(i, k, j, 3), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    READ (fid) (((q_b(i, k, j, 4), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)

    READ (fid) (((q_bt(i, k, j, 1), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    READ (fid) (((q_bt(i, k, j, 2), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    READ (fid) (((q_bt(i, k, j, 3), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    READ (fid) (((q_bt(i, k, j, 4), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)

    PRINT *, 'check initial bdy values: th_bt ', MAXVAL(th_bt), MINVAL(th_bt)
    PRINT *, 'check initial bdy values: pi_bt ', MAXVAL(pi_bt), MINVAL(pi_bt)
    PRINT *, 'check initial bdy values: th_bt ', MAXVAL(th_bt), MINVAL(th_bt)
    PRINT *, 'check initial bdy values: u_bt ', MAXVAL(u_bt), MINVAL(u_bt)
    PRINT *, 'check initial bdy values: v_bt ', MAXVAL(v_bt), MINVAL(v_bt)
    PRINT *, 'check initial bdy values: q_bt ', MAXVAL(q_bt), MINVAL(q_bt)

    u_b_old = u_b
    v_b_old = v_b
    q_b_old = q_b
    th_b_old = th_b
    pi_b_old = pi_b

    CALL get_second_time_boundary(pi_b, pi_bt, &
                                  spec_bdy_width, interval_seconds, &
                                  ids, ide, jds, jde, kds, kde, &
                                  ims, ime, jms, jme, kms, kme, &
                                  its, ite, jts, jte, kts - 1, kte + 1)

    CALL get_second_time_boundary(th_b, th_bt, &
                                  spec_bdy_width, interval_seconds, &
                                  ids, ide, jds, jde, kds, kde, &
                                  ims, ime, jms, jme, kms, kme, &
                                  its, ite, jts, jte, kts, kte + 1)

    CALL get_second_time_boundary(u_b, u_bt, &
                                  spec_bdy_width, interval_seconds, &
                                  ids, ide, jds, jde, kds, kde, &
                                  ims, ime, jms, jme, kms, kme, &
                                  its, ite, jts, jte, kts, kte)

    CALL get_second_time_boundary(v_b, v_bt, &
                                  spec_bdy_width, interval_seconds, &
                                  ids, ide, jds, jde, kds, kde, &
                                  ims, ime, jms, jme, kms, kme, &
                                  its, ite, jts, jte, kts, kte)

    CALL get_second_time_boundary(q_b, q_bt, &
                                  spec_bdy_width, interval_seconds, &
                                  ids, ide, jds, jde, kds, kde, &
                                  ims, ime, jms, jme, kms, kme, &
                                  its, ite, jts, jte, kts, kte + 1)

    RETURN

  END SUBROUTINE read_initial_boundary

  SUBROUTINE get_second_time_boundary(field, field_bt, &
                                      spec_bdy_width, interval_seconds, &
                                      ids, ide, jds, jde, kds, kde, &
                                      ims, ime, jms, jme, kms, kme, &
                                      its, ite, jts, jte, kts, kte)

    IMPLICIT NONE
    INTEGER, INTENT(IN) ::      ims, ime, jms, jme, kms, kme, &
                           ids, ide, jds, jde, kds, kde, &
                           its, ite, jts, jte, kts, kte

    INTEGER, INTENT(IN) ::  spec_bdy_width, interval_seconds
    REAL, DIMENSION(MAX(ide, jde), (kds - 1):kde, spec_bdy_width, 4), INTENT(IN)  :: field_bt
    REAL, DIMENSION(MAX(ide, jde), (kds - 1):kde, spec_bdy_width, 4), INTENT(INOUT):: field
    ! local var
    INTEGER       :: i, j, k
    REAL          :: c

    c = float(interval_seconds)
    !
    ! Y-start boundary
    DO j = 1, spec_bdy_width
      DO k = kts, kte
        DO i = ids, ide
          field(i, k, j, P_YSB) = field(i, k, j, P_YSB) + field_bt(i, k, j, P_YSB) * c
        END DO
      END DO
    END DO

    ! Y-end boundary
    DO j = 1, spec_bdy_width
      DO k = kts, kte
        DO i = ids, ide
          field(i, k, j, P_YEB) = field(i, k, j, P_YEB) + field_bt(i, k, j, P_YEB) * c
        END DO
      END DO
    END DO

    ! X-start boundary
    DO i = 1, spec_bdy_width
      DO k = kts, kte
        DO j = jds, jde
          field(j, k, i, P_XSB) = field(j, k, i, P_XSB) + field_bt(j, k, i, P_XSB) * c
        END DO
      END DO
    END DO

    ! X-end boundary
    DO i = 1, spec_bdy_width
      DO k = kts, kte
        DO j = jds, jde
          field(j, k, i, P_XEB) = field(j, k, i, P_XEB) + field_bt(j, k, i, P_XEB) * c
        END DO
      END DO
    END DO

  END SUBROUTINE get_second_time_boundary

  SUBROUTINE get_boundary_value(u, v, pi, th, q, &
                                u_b, v_b, pi_b, th_b, q_b, &
                                spec_bdy_width, &
                                ids, ide, jds, jde, kds, kde, &
                                ims, ime, jms, jme, kms, kme, &
                                its, ite, jts, jte, kts, kte)

    IMPLICIT NONE
    INTEGER, INTENT(IN) ::      ims, ime, jms, jme, kms, kme, &
                           ids, ide, jds, jde, kds, kde, &
                           its, ite, jts, jte, kts, kte

    ! Arguments
    INTEGER, INTENT(IN)            :: spec_bdy_width
    REAL, DIMENSION(ims:ime, jms:jme, kms:kme), INTENT(IN) :: u, v, pi, th, q
    REAL, DIMENSION(MAX(ide, jde), (kds - 1):kde, spec_bdy_width, 4), INTENT(OUT) :: u_b, v_b, pi_b, th_b, q_b

    !

    CALL get_field_boundary_value(pi, pi_b, &
                                  spec_bdy_width, &
                                  ids, ide, jds, jde, kds, kde, &
                                  ims, ime, jms, jme, kms, kme, &
                                  its, ite, jts, jte, kts - 1, kte + 1)

    CALL get_field_boundary_value(th, th_b, &
                                  spec_bdy_width, &
                                  ids, ide, jds, jde, kds, kde, &
                                  ims, ime, jms, jme, kms, kme, &
                                  its, ite, jts, jte, kts, kte + 1)

    CALL get_field_boundary_value(u, u_b, &
                                  spec_bdy_width, &
                                  ids, ide, jds, jde, kds, kde, &
                                  ims, ime, jms, jme, kms, kme, &
                                  its, ite, jts, jte, kts, kte)

    CALL get_field_boundary_value(v, v_b, &
                                  spec_bdy_width, &
                                  ids, ide, jds, jde, kds, kde, &
                                  ims, ime, jms, jme, kms, kme, &
                                  its, ite, jts, jte, kts, kte)

    CALL get_field_boundary_value(q, q_b, &
                                  spec_bdy_width, &
                                  ids, ide, jds, jde, kds, kde, &
                                  ims, ime, jms, jme, kms, kme, &
                                  its, ite, jts, jte, kts, kte + 1)
  END SUBROUTINE get_boundary_value

  SUBROUTINE get_field_boundary_value(field, field_b, &
                                      spec_bdy_width, &
                                      ids, ide, jds, jde, kds, kde, &
                                      ims, ime, jms, jme, kms, kme, &
                                      its, ite, jts, jte, kts, kte)

    IMPLICIT NONE
    INTEGER, INTENT(IN) ::      ims, ime, jms, jme, kms, kme, &
                           ids, ide, jds, jde, kds, kde, &
                           its, ite, jts, jte, kts, kte

    ! Arguments
    INTEGER, INTENT(IN)            :: spec_bdy_width
    REAL, DIMENSION(ims:ime, jms:jme, kms:kme), INTENT(IN) :: field
    REAL, DIMENSION(MAX(ide, jde), (kds - 1):kde, spec_bdy_width, 4), INTENT(OUT) :: field_b
    ! local var
    INTEGER       :: i, ii, j, jj, k

    DO j = 1, spec_bdy_width
      DO k = kts, kte
        DO i = ids, ide
          field_b(i, k, j, P_YSB) = field(i, j, k)
        END DO
      END DO
    END DO

    DO j = jde, jde - spec_bdy_width + 1, -1
      jj = jde - j + 1
      DO k = kts, kte
        DO i = ids, ide
          field_b(i, k, jj, P_YEB) = field(i, j, k)
        END DO
      END DO
    END DO

    DO i = 1, spec_bdy_width
      DO k = kts, kte
        DO j = jds, jde
          field_b(j, k, i, P_XSB) = field(i, j, k)
        END DO
      END DO
    END DO

    DO i = ide, ide - spec_bdy_width + 1, -1
      ii = ide - i + 1
      DO k = kts, kte
        DO j = jds, jde
          field_b(j, k, ii, P_XEB) = field(i, j, k)
        END DO
      END DO
    END DO

  END SUBROUTINE get_field_boundary_value

  SUBROUTINE get_boundary_tendency(spec_bdy_width, &
                                   pi_b_1, th_b_1, u_b_1, v_b_1, q_b_1, &
                                   pi_b_2, th_b_2, u_b_2, v_b_2, q_b_2, &
                                   pi_bt, th_bt, u_bt, v_bt, q_bt, &
                                   interval_seconds, &
                                   ids, ide, jds, jde, kds, kde, &
                                   ims, ime, jms, jme, kms, kme, &
                                   its, ite, jts, jte, kts, kte)

    IMPLICIT NONE
    INTEGER, INTENT(IN) ::      ims, ime, jms, jme, kms, kme, &
                           ids, ide, jds, jde, kds, kde, &
                           its, ite, jts, jte, kts, kte
    INTEGER, INTENT(IN) ::  interval_seconds, spec_bdy_width

    ! Arguments

    REAL, DIMENSION(MAX(ide, jde), (kds - 1):kde, spec_bdy_width, 4), INTENT(IN) :: pi_b_1, th_b_1, u_b_1, v_b_1, q_b_1
    REAL, DIMENSION(MAX(ide, jde), (kds - 1):kde, spec_bdy_width, 4), INTENT(IN) :: pi_b_2, th_b_2, u_b_2, v_b_2, q_b_2
    REAL, DIMENSION(MAX(ide, jde), (kds - 1):kde, spec_bdy_width, 4), INTENT(OUT) :: pi_bt, th_bt, u_bt, v_bt, q_bt

    ! local var
    REAL       :: c

    c = 1./float(interval_seconds)

    CALL get_field_boundary_tendency(pi_b_2, pi_b_1, pi_bt, &
                                     spec_bdy_width, c, &
                                     ids, ide, jds, jde, kds, kde, &
                                     ims, ime, jms, jme, kms, kme, &
                                     its, ite, jts, jte, kts - 1, kte + 1)

    CALL get_field_boundary_tendency(th_b_2, th_b_1, th_bt, &
                                     spec_bdy_width, c, &
                                     ids, ide, jds, jde, kds, kde, &
                                     ims, ime, jms, jme, kms, kme, &
                                     its, ite, jts, jte, kts, kte + 1)
    CALL get_field_boundary_tendency(u_b_2, u_b_1, u_bt, &
                                     spec_bdy_width, c, &
                                     ids, ide, jds, jde, kds, kde, &
                                     ims, ime, jms, jme, kms, kme, &
                                     its, ite, jts, jte, kts, kte)

    CALL get_field_boundary_tendency(v_b_2, v_b_1, v_bt, &
                                     spec_bdy_width, c, &
                                     ids, ide, jds, jde, kds, kde, &
                                     ims, ime, jms, jme, kms, kme, &
                                     its, ite, jts, jte, kts, kte)

    CALL get_field_boundary_tendency(q_b_2, q_b_1, q_bt, &
                                     spec_bdy_width, c, &
                                     ids, ide, jds, jde, kds, kde, &
                                     ims, ime, jms, jme, kms, kme, &
                                     its, ite, jts, jte, kts, kte + 1)

  END SUBROUTINE get_boundary_tendency

  SUBROUTINE get_field_boundary_tendency(field_2, field_1, field_bt, &
                                         spec_bdy_width, c, &
                                         ids, ide, jds, jde, kds, kde, &
                                         ims, ime, jms, jme, kms, kme, &
                                         its, ite, jts, jte, kts, kte)

    IMPLICIT NONE
    INTEGER, INTENT(IN) ::      ims, ime, jms, jme, kms, kme, &
                           ids, ide, jds, jde, kds, kde, &
                           its, ite, jts, jte, kts, kte

    INTEGER, INTENT(IN) ::  spec_bdy_width
    REAL, INTENT(IN) :: c
    REAL, DIMENSION(MAX(ide, jde), (kds - 1):kde, spec_bdy_width, 4), INTENT(IN)  :: field_2, field_1
    REAL, DIMENSION(MAX(ide, jde), (kds - 1):kde, spec_bdy_width, 4), INTENT(OUT)  :: field_bt
    ! local var
    INTEGER       :: i, j, k
    !
    ! Y-start boundary
    DO j = 1, spec_bdy_width
      DO k = kts, kte
        DO i = ids, ide
          field_bt(i, k, j, P_YSB) = (field_2(i, k, j, P_YSB) - field_1(i, k, j, P_YSB)) * c
        END DO
      END DO
    END DO

    ! Y-end boundary
    DO j = 1, spec_bdy_width
      DO k = kts, kte
        DO i = ids, ide
          field_bt(i, k, j, P_YEB) = (field_2(i, k, j, P_YEB) - field_1(i, k, j, P_YEB)) * c
        END DO
      END DO
    END DO

    ! X-start boundary
    DO i = 1, spec_bdy_width
      DO k = kts, kte
        DO j = jds, jde
          field_bt(j, k, i, P_XSB) = (field_2(j, k, i, P_XSB) - field_1(j, k, i, P_XSB)) * c
        END DO
      END DO
    END DO

    ! X-end boundary
    DO i = 1, spec_bdy_width
      DO k = kts, kte
        DO j = jds, jde
          field_bt(j, k, i, P_XEB) = (field_2(j, k, i, P_XEB) - field_1(j, k, i, P_XEB)) * c
        END DO
      END DO
    END DO

  END SUBROUTINE get_field_boundary_tendency

  SUBROUTINE read_boundary_value(pi_b, th_b, u_b, v_b, q_b, &
                                 pi_bt, th_bt, u_bt, v_bt, q_bt, &
                                 mdate, spec_bdy_width, fid, &
                                 ids, ide, jds, jde, kds, kde, &
                                 ims, ime, jms, jme, kms, kme, &
                                 its, ite, jts, jte, kts, kte)

    IMPLICIT NONE
    INTEGER, INTENT(IN) ::      ims, ime, jms, jme, kms, kme, &
                           ids, ide, jds, jde, kds, kde, &
                           its, ite, jts, jte, kts, kte

    ! Arguments
    INTEGER, INTENT(IN)            ::  spec_bdy_width, fid
    INTEGER, DIMENSION(5), INTENT(out)  :: mdate
    REAL, DIMENSION(MAX(ide, jde), (kds - 1):kde, spec_bdy_width, 4), INTENT(OUT) :: u_b, v_b, pi_b, th_b, q_b
    REAL, DIMENSION(MAX(ide, jde), (kds - 1):kde, spec_bdy_width, 4), INTENT(OUT) :: u_bt, v_bt, pi_bt, th_bt, q_bt

    !
    INTEGER    :: i, j, k

    READ (fid) mdate
    PRINT *, 'bdy mdate ', mdate
    READ (fid) (((pi_b(i, k, j, 1), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    READ (fid) (((pi_b(i, k, j, 2), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    READ (fid) (((pi_b(i, k, j, 3), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    READ (fid) (((pi_b(i, k, j, 4), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)

    READ (fid) (((pi_bt(i, k, j, 1), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    READ (fid) (((pi_bt(i, k, j, 2), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    READ (fid) (((pi_bt(i, k, j, 3), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    READ (fid) (((pi_bt(i, k, j, 4), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)

    READ (fid) (((th_b(i, k, j, 1), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    READ (fid) (((th_b(i, k, j, 2), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    READ (fid) (((th_b(i, k, j, 3), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    READ (fid) (((th_b(i, k, j, 4), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)

    READ (fid) (((th_bt(i, k, j, 1), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    READ (fid) (((th_bt(i, k, j, 2), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    READ (fid) (((th_bt(i, k, j, 3), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    READ (fid) (((th_bt(i, k, j, 4), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)

    READ (fid) (((u_b(i, k, j, 1), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    READ (fid) (((u_b(i, k, j, 2), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    READ (fid) (((u_b(i, k, j, 3), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    READ (fid) (((u_b(i, k, j, 4), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)

    READ (fid) (((u_bt(i, k, j, 1), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    READ (fid) (((u_bt(i, k, j, 2), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    READ (fid) (((u_bt(i, k, j, 3), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    READ (fid) (((u_bt(i, k, j, 4), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)

    READ (fid) (((v_b(i, k, j, 1), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    READ (fid) (((v_b(i, k, j, 2), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    READ (fid) (((v_b(i, k, j, 3), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    READ (fid) (((v_b(i, k, j, 4), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)

    READ (fid) (((v_bt(i, k, j, 1), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    READ (fid) (((v_bt(i, k, j, 2), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    READ (fid) (((v_bt(i, k, j, 3), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    READ (fid) (((v_bt(i, k, j, 4), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)

    READ (fid) (((q_b(i, k, j, 1), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    READ (fid) (((q_b(i, k, j, 2), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    READ (fid) (((q_b(i, k, j, 3), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    READ (fid) (((q_b(i, k, j, 4), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)

    READ (fid) (((q_bt(i, k, j, 1), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    READ (fid) (((q_bt(i, k, j, 2), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    READ (fid) (((q_bt(i, k, j, 3), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    READ (fid) (((q_bt(i, k, j, 4), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)

    RETURN

  END SUBROUTINE read_boundary_value

  SUBROUTINE Compute_num_loop(num_loop, start_year, start_month, start_day, start_hour, &
                              end_year, end_month, end_day, end_hour, interval_seconds)

    IMPLICIT NONE

    INTEGER, INTENT(IN)   :: start_year, start_month, start_day, start_hour
    INTEGER, INTENT(IN)   :: end_year, end_month, end_day, end_hour, interval_seconds
    INTEGER, INTENT(OUT)  :: num_loop
    ! variable
    INTEGER  :: year, month, day, hour

    num_loop = 1
    loopa: DO

      CALL Get_Newtime(num_loop, year, month, day, hour, &
                       start_year, start_month, start_day, start_hour, interval_seconds)

      IF ((year == end_year) .AND. (month == end_month) .AND. &
          (day == end_day) .AND. (hour == end_hour)) EXIT

      IF ((year > end_year) .OR. ((year == end_year) .AND. &
                                  (month > end_month)) .OR. ((year == end_year) .AND. &
                                                             (month == end_month) .AND. (day > end_day)) .OR. &
          ((year == end_year) .AND. (month == end_month) .AND. &
           (day == end_day) .AND. (hour > end_hour))) THEN
        PRINT *, "DATE in your namelist is wrong!"
        CALL abort
      END IF
      num_loop = num_loop + 1

    END DO loopa

  END SUBROUTINE Compute_num_loop

  SUBROUTINE Get_Newtime(loop, current_year, current_month, current_day, current_hour, &
                         start_year, start_month, start_day, start_hour, interval_seconds)

    IMPLICIT NONE

    INTEGER, INTENT(IN)    :: loop  ! interval number
    INTEGER, INTENT(IN)    :: start_year, start_month, start_day, &
                              start_hour, interval_seconds
    INTEGER, INTENT(OUT)   :: current_year, current_month, &
                              current_day, current_hour

    !local variable:
    INTEGER :: intv_hour, intv_day

    current_year = start_year
    current_month = start_month
    current_day = start_day

    intv_hour = interval_seconds * (loop - 1) / 3600
    current_hour = start_hour + intv_hour

    intv_day = current_hour / 24
    current_hour = MOD(current_hour, 24)

    IF (intv_day >= 1) THEN
      current_day = current_day + intv_day
      IF ((current_month == 2) .AND. (current_day > 28)) THEN
        IF (MOD(current_year, 4) /= 0) THEN
          current_day = current_day - 28
          current_month = current_month + 1
        ELSEIF (MOD(current_year, 100) == 0 .AND. MOD(current_year, 400) /= 0) THEN
          current_day = current_day - 28
          current_month = current_month + 1
        ELSEIF (current_day > 29) THEN
          current_day = current_day - 29
          current_month = current_month + 1
        END IF
      ELSEIF (current_day > 30) THEN
        SELECT CASE (current_month)
        CASE (4, 6, 9, 11)
          current_day = current_day - 30
          current_month = current_month + 1
        CASE (1, 3, 5, 7, 8, 10)
          IF (current_day > 31) THEN
            current_day = current_day - 31
            current_month = current_month + 1
          END IF
        CASE (12)
          IF (current_day > 31) THEN
            current_day = current_day - 31
            current_month = 1
            current_year = current_year + 1
          END IF
        END SELECT
      END IF
    END IF

  END SUBROUTINE Get_Newtime

  SUBROUTINE write_grapes_var_bdy(pi_b, th_b, u_b, v_b, q_b, &
                                  pi_bt, th_bt, u_bt, v_bt, q_bt, &
                                  spec_bdy_width, mdate, fid, &
                                  ids, ide, jds, jde, kds, kde, &
                                  ims, ime, jms, jme, kms, kme, &
                                  its, ite, jts, jte, kts, kte)

    IMPLICIT NONE
    INTEGER, INTENT(IN) ::      ims, ime, jms, jme, kms, kme, &
                           ids, ide, jds, jde, kds, kde, &
                           its, ite, jts, jte, kts, kte
    INTEGER, INTENT(IN) ::      spec_bdy_width

    REAL, DIMENSION(MAX(ide, jde), (kds - 1):kde, spec_bdy_width, 4), INTENT(IN) :: pi_b, th_b, u_b, v_b, q_b
    REAL, DIMENSION(MAX(ide, jde), (kds - 1):kde, spec_bdy_width, 4), INTENT(IN) :: pi_bt, th_bt, u_bt, v_bt, q_bt

    INTEGER, INTENT(IN)                :: fid
    INTEGER, DIMENSION(5), INTENT(in)   :: mdate

    INTEGER                       ::  i, j, k
    INTEGER, DIMENSION(5)         :: mdate_ivf

    DO i = 1, 4
      mdate_ivf(i) = mdate(i)
    END DO
    mdate_ivf(5) = 0

    WRITE (fid) mdate
    WRITE (fid) (((pi_b(i, k, j, 1), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    WRITE (fid) (((pi_b(i, k, j, 2), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    WRITE (fid) (((pi_b(i, k, j, 3), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    WRITE (fid) (((pi_b(i, k, j, 4), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)

    WRITE (fid) (((pi_bt(i, k, j, 1), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    WRITE (fid) (((pi_bt(i, k, j, 2), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    WRITE (fid) (((pi_bt(i, k, j, 3), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    WRITE (fid) (((pi_bt(i, k, j, 4), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)

    WRITE (fid) (((th_b(i, k, j, 1), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    WRITE (fid) (((th_b(i, k, j, 2), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    WRITE (fid) (((th_b(i, k, j, 3), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    WRITE (fid) (((th_b(i, k, j, 4), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)

    WRITE (fid) (((th_bt(i, k, j, 1), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    WRITE (fid) (((th_bt(i, k, j, 2), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    WRITE (fid) (((th_bt(i, k, j, 3), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    WRITE (fid) (((th_bt(i, k, j, 4), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)

    WRITE (fid) (((u_b(i, k, j, 1), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    WRITE (fid) (((u_b(i, k, j, 2), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    WRITE (fid) (((u_b(i, k, j, 3), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    WRITE (fid) (((u_b(i, k, j, 4), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)

    WRITE (fid) (((u_bt(i, k, j, 1), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    WRITE (fid) (((u_bt(i, k, j, 2), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    WRITE (fid) (((u_bt(i, k, j, 3), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    WRITE (fid) (((u_bt(i, k, j, 4), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)

    WRITE (fid) (((v_b(i, k, j, 1), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    WRITE (fid) (((v_b(i, k, j, 2), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    WRITE (fid) (((v_b(i, k, j, 3), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    WRITE (fid) (((v_b(i, k, j, 4), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)

    WRITE (fid) (((v_bt(i, k, j, 1), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    WRITE (fid) (((v_bt(i, k, j, 2), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    WRITE (fid) (((v_bt(i, k, j, 3), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    WRITE (fid) (((v_bt(i, k, j, 4), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)

    WRITE (fid) (((q_b(i, k, j, 1), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    WRITE (fid) (((q_b(i, k, j, 2), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    WRITE (fid) (((q_b(i, k, j, 3), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    WRITE (fid) (((q_b(i, k, j, 4), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)

    WRITE (fid) (((q_bt(i, k, j, 1), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    WRITE (fid) (((q_bt(i, k, j, 2), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    WRITE (fid) (((q_bt(i, k, j, 3), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)
    WRITE (fid) (((q_bt(i, k, j, 4), i=1, MAX(ide, jde)), k=kms, kme), j=1, spec_bdy_width)

  END SUBROUTINE write_grapes_var_bdy

  SUBROUTINE write_grapes_bdy_grads(fid1, fid2, fid3, fid4, &
                                    pi_b, th_b, u_b, v_b, q_b, &
                                    pi_bt, th_bt, u_bt, v_bt, q_bt, &
                                    spec_bdy_width, mdate, &
                                    ids, ide, jds, jde, kds, kde, &
                                    ims, ime, jms, jme, kms, kme, &
                                    its, ite, jts, jte, kts, kte)

    IMPLICIT NONE
    INTEGER, INTENT(IN) ::      fid1, fid2, fid3, fid4
    INTEGER, INTENT(IN) ::      ims, ime, jms, jme, kms, kme, &
                           ids, ide, jds, jde, kds, kde, &
                           its, ite, jts, jte, kts, kte
    INTEGER, INTENT(IN) ::      spec_bdy_width

    REAL, DIMENSION(MAX(ide, jde), (kds - 1):kde, spec_bdy_width, 4), INTENT(IN) :: pi_b, th_b, u_b, v_b, q_b
    REAL, DIMENSION(MAX(ide, jde), (kds - 1):kde, spec_bdy_width, 4), INTENT(IN) :: pi_bt, th_bt, u_bt, v_bt, q_bt

    !   INTEGER, INTENT(IN)                :: fid
    INTEGER, DIMENSION(4), INTENT(in)   :: mdate

    INTEGER                       ::  i, j, k
    INTEGER, DIMENSION(5)         :: mdate_ivf

    DO k = kms, kme
      WRITE (fid1) ((pi_b(i, k, j, 1), i=1, MAX(ide, jde)), j=1, spec_bdy_width)
    END DO
    DO k = kms, kme
      WRITE (fid1) ((pi_bt(i, k, j, 1), i=1, MAX(ide, jde)), j=1, spec_bdy_width)
    END DO
    DO k = kms, kme
      WRITE (fid1) ((th_b(i, k, j, 1), i=1, MAX(ide, jde)), j=1, spec_bdy_width)
    END DO
    DO k = kms, kme
      WRITE (fid1) ((th_bt(i, k, j, 1), i=1, MAX(ide, jde)), j=1, spec_bdy_width)
    END DO
    DO k = kms, kme
      WRITE (fid1) ((u_b(i, k, j, 1), i=1, MAX(ide, jde)), j=1, spec_bdy_width)
    END DO
    DO k = kms, kme
      WRITE (fid1) ((u_bt(i, k, j, 1), i=1, MAX(ide, jde)), j=1, spec_bdy_width)
    END DO
    DO k = kms, kme
      WRITE (fid1) ((v_b(i, k, j, 1), i=1, MAX(ide, jde)), j=1, spec_bdy_width)
    END DO
    DO k = kms, kme
      WRITE (fid1) ((v_bt(i, k, j, 1), i=1, MAX(ide, jde)), j=1, spec_bdy_width)
    END DO
    DO k = kms, kme
      WRITE (fid1) ((q_b(i, k, j, 1), i=1, MAX(ide, jde)), j=1, spec_bdy_width)
    END DO
    DO k = kms, kme
      WRITE (fid1) ((q_bt(i, k, j, 1), i=1, MAX(ide, jde)), j=1, spec_bdy_width)
    END DO

    DO k = kms, kme
      WRITE (fid2) ((pi_b(i, k, j, 2), i=1, MAX(ide, jde)), j=1, spec_bdy_width)
    END DO
    DO k = kms, kme
      WRITE (fid2) ((pi_bt(i, k, j, 2), i=1, MAX(ide, jde)), j=1, spec_bdy_width)
    END DO
    DO k = kms, kme
      WRITE (fid2) ((th_b(i, k, j, 2), i=1, MAX(ide, jde)), j=1, spec_bdy_width)
    END DO
    DO k = kms, kme
      WRITE (fid2) ((th_bt(i, k, j, 2), i=1, MAX(ide, jde)), j=1, spec_bdy_width)
    END DO
    DO k = kms, kme
      WRITE (fid2) ((u_b(i, k, j, 2), i=1, MAX(ide, jde)), j=1, spec_bdy_width)
    END DO
    DO k = kms, kme
      WRITE (fid2) ((u_bt(i, k, j, 2), i=1, MAX(ide, jde)), j=1, spec_bdy_width)
    END DO
    DO k = kms, kme
      WRITE (fid2) ((v_b(i, k, j, 2), i=1, MAX(ide, jde)), j=1, spec_bdy_width)
    END DO
    DO k = kms, kme
      WRITE (fid2) ((v_bt(i, k, j, 2), i=1, MAX(ide, jde)), j=1, spec_bdy_width)
    END DO
    DO k = kms, kme
      WRITE (fid2) ((q_b(i, k, j, 2), i=1, MAX(ide, jde)), j=1, spec_bdy_width)
    END DO
    DO k = kms, kme
      WRITE (fid2) ((q_bt(i, k, j, 2), i=1, MAX(ide, jde)), j=1, spec_bdy_width)
    END DO

    DO k = kms, kme
      WRITE (fid3) ((pi_b(i, k, j, 3), i=1, MAX(ide, jde)), j=1, spec_bdy_width)
    END DO
    DO k = kms, kme
      WRITE (fid3) ((pi_bt(i, k, j, 3), i=1, MAX(ide, jde)), j=1, spec_bdy_width)
    END DO
    DO k = kms, kme
      WRITE (fid3) ((th_b(i, k, j, 3), i=1, MAX(ide, jde)), j=1, spec_bdy_width)
    END DO
    DO k = kms, kme
      WRITE (fid3) ((th_bt(i, k, j, 3), i=1, MAX(ide, jde)), j=1, spec_bdy_width)
    END DO
    DO k = kms, kme
      WRITE (fid3) ((u_b(i, k, j, 3), i=1, MAX(ide, jde)), j=1, spec_bdy_width)
    END DO
    DO k = kms, kme
      WRITE (fid3) ((u_bt(i, k, j, 3), i=1, MAX(ide, jde)), j=1, spec_bdy_width)
    END DO
    DO k = kms, kme
      WRITE (fid3) ((v_b(i, k, j, 3), i=1, MAX(ide, jde)), j=1, spec_bdy_width)
    END DO
    DO k = kms, kme
      WRITE (fid3) ((v_bt(i, k, j, 3), i=1, MAX(ide, jde)), j=1, spec_bdy_width)
    END DO
    DO k = kms, kme
      WRITE (fid3) ((q_b(i, k, j, 3), i=1, MAX(ide, jde)), j=1, spec_bdy_width)
    END DO
    DO k = kms, kme
      WRITE (fid3) ((q_bt(i, k, j, 3), i=1, MAX(ide, jde)), j=1, spec_bdy_width)
    END DO

    DO k = kms, kme
      WRITE (fid4) ((pi_b(i, k, j, 4), i=1, MAX(ide, jde)), j=1, spec_bdy_width)
    END DO
    DO k = kms, kme
      WRITE (fid4) ((pi_bt(i, k, j, 4), i=1, MAX(ide, jde)), j=1, spec_bdy_width)
    END DO
    DO k = kms, kme
      WRITE (fid4) ((th_b(i, k, j, 4), i=1, MAX(ide, jde)), j=1, spec_bdy_width)
    END DO
    DO k = kms, kme
      WRITE (fid4) ((th_bt(i, k, j, 4), i=1, MAX(ide, jde)), j=1, spec_bdy_width)
    END DO
    DO k = kms, kme
      WRITE (fid4) ((u_b(i, k, j, 4), i=1, MAX(ide, jde)), j=1, spec_bdy_width)
    END DO
    DO k = kms, kme
      WRITE (fid4) ((u_bt(i, k, j, 4), i=1, MAX(ide, jde)), j=1, spec_bdy_width)
    END DO
    DO k = kms, kme
      WRITE (fid4) ((v_b(i, k, j, 4), i=1, MAX(ide, jde)), j=1, spec_bdy_width)
    END DO
    DO k = kms, kme
      WRITE (fid4) ((v_bt(i, k, j, 4), i=1, MAX(ide, jde)), j=1, spec_bdy_width)
    END DO
    DO k = kms, kme
      WRITE (fid4) ((q_b(i, k, j, 4), i=1, MAX(ide, jde)), j=1, spec_bdy_width)
    END DO
    DO k = kms, kme
      WRITE (fid4) ((q_bt(i, k, j, 4), i=1, MAX(ide, jde)), j=1, spec_bdy_width)
    END DO

  END SUBROUTINE write_grapes_bdy_grads

END MODULE

