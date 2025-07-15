! Developed by CMA GRAPES_MESO 5.0
! Modified by YaliWu for CMA-GD intaerfaces 20220720
MODULE UpdateBdy
  USE module_domain, ONLY: domain_t
  USE module_configure, ONLY: grid_config_rec_type
  USE module_variables
  USE module_utility
  USE kinds_m
  USE NMLRead_m
  USE YAMLRead_m
  IMPLICIT NONE

CONTAINS

  SUBROUTINE update_bdy(bdyPathName, anaPathName, nmlPathName)
    !-------------------------------------------------------------------------------
    CHARACTER(LEN=*), INTENT(IN) :: bdyPathName, anaPathName, nmlPathName
    TYPE(domain_t)            :: hgrid
    INTEGER(i_kind) :: s_we, e_we, s_sn, e_sn, s_vert, e_vert
    INTEGER(i_kind) :: spec_bdy_width
    REAL(r_kind)    :: xs_we, ys_sn, xd, yd
    INTEGER(i_kind) :: start_year, start_month, start_day, start_hour
    INTEGER(i_kind) :: interval_seconds
    INTEGER(i_kind) :: end_year, end_month, end_day, end_hour
    CHARACTER(LEN=1024) :: nlFileName, OutputFileName
    INTEGER(i_kind)     :: fid_init_in, fid_bdy_in, fid_out
    INTEGER(i_kind)     :: loop, num_loop
    INTEGER(i_kind)     :: k, i, j, nout, ierr
    LOGICAL             :: readnext
    INTEGER(i_kind), DIMENSION(5)    :: mdate
    CHARACTER(len=64)   :: flm_in, flm_out, flm_in_mpi
    LOGICAL  :: setbdy
    INTEGER                          :: ids, ide, jds, jde, kds, kde, &
                                        ims, ime, jms, jme, kms, kme, &
                                        its, ite, jts, jte, kts, kte
    INTEGER :: ifile
    CHARACTER(LEN=1024) :: INPUT_SI, BDY_SI, BDY_UPDATE
    !----------------------------------------------------------------------------------

    nlFileName = TRIM(nmlPathName)
    OutputFileName = TRIM(anaPathName)
    PRINT *, 'nlFileName: ', TRIM(nlFileName)
    PRINT *, 'OutputFileName: ', TRIM(OutputFileName)

    CALL namelist_read(nlFileName, 's_we', s_we)
    CALL namelist_read(nlFileName, 'e_we', e_we)
    CALL namelist_read(nlFileName, 's_sn', s_sn)
    CALL namelist_read(nlFileName, 'e_sn', e_sn)
    CALL namelist_read(nlFileName, 's_vert', s_vert)
    CALL namelist_read(nlFileName, 'e_vert', e_vert)
    CALL namelist_read(nlFileName, 'spec_bdy_width', spec_bdy_width)
    CALL namelist_read(nlFileName, 'xd', xd)
    CALL namelist_read(nlFileName, 'yd', yd)
    CALL namelist_read(nlFileName, 'start_year', start_year)
    CALL namelist_read(nlFileName, 'start_month', start_month)
    CALL namelist_read(nlFileName, 'start_day', start_day)
    CALL namelist_read(nlFileName, 'start_hour', start_hour)
    CALL namelist_read(nlFileName, 'end_year', end_year)
    CALL namelist_read(nlFileName, 'end_month', end_month)
    CALL namelist_read(nlFileName, 'end_day', end_day)
    CALL namelist_read(nlFileName, 'end_hour', end_hour)
    CALL namelist_read(nlFileName, 'interval_seconds', interval_seconds)
    CALL namelist_read(nlFileName, 'xs_we', xs_we)
    CALL namelist_read(nlFileName, 'ys_sn', ys_sn)

    hgrid = domain_t(TRIM(nlFileName))
    CALL get_ijk(s_we, e_we, s_sn, e_sn, s_vert, e_vert, &
                 ids, ide, jds, jde, kds, kde, &
                 ims, ime, jms, jme, kms, kme, &
                 its, ite, jts, jte, kts, kte)

    ! idn = hgrid%ide - hgrid%ids + 1
    ! jdn = hgrid%jde - hgrid%jds + 1
    ! kdn = hgrid%kde - hgrid%kds + 1
    ! ids = hgrid%ids; ide = hgrid%ide
    ! jds = hgrid%jds; jde = hgrid%jde
    ! kds = hgrid%kds; kde = hgrid%kde

    fid_init_in = 11
    fid_bdy_in = 12
    fid_out = 13

    CALL alloc_space_field(ids, ide, jds, jde, kds, kde, &
                           ims, ime, jms, jme, kms, kme, &
                           its, ite, jts, jte, kts, kte, spec_bdy_width)

    INPUT_SI = TRIM(OutputFileName)
    BDY_SI = TRIM(bdyPathName)
    BDY_UPDATE = TRIM(bdyPathName)//'_update'
    PRINT *, 'grapesinput IN: ', TRIM(INPUT_SI)
    PRINT *, 'grapesbdy   IN: ', TRIM(BDY_SI)
    PRINT *, 'grapesbdy   OUT: ', TRIM(BDY_UPDATE)

    OPEN (fid_init_in, file=TRIM(OutputFileName), form='unformatted', status='unknown', ACCESS='sequential', convert="big_endian")
    OPEN (fid_bdy_in, file=TRIM(BDY_SI), form='unformatted', status='unknown', ACCESS='sequential', convert="big_endian")
    OPEN (fid_out, file=TRIM(BDY_UPDATE), form='unformatted', status='unknown', ACCESS='sequential', convert="big_endian")

    !----output for test------------------------------------------
    !     OPEN(203,file=TRIM(output_dir)//'/bdy1.dat',form='unformatted',ACCESS='sequential', convert="big_endian")
    !     OPEN(204,file=TRIM(output_dir)//'/bdy2.dat',form='unformatted',ACCESS='sequential', convert="big_endian")
    !     OPEN(205,file=TRIM(output_dir)//'/bdy3.dat',form='unformatted',ACCESS='sequential', convert="big_endian")
    !     OPEN(206,file=TRIM(output_dir)//'/bdy4.dat',form='unformatted',ACCESS='sequential', convert="big_endian")
    !
    !     OPEN(207,file=TRIM(output_dir)//'/bdy1_si.dat',form='unformatted',ACCESS='sequential', convert="big_endian")
    !     OPEN(208,file=TRIM(output_dir)//'/bdy2_si.dat',form='unformatted',ACCESS='sequential', convert="big_endian")
    !     OPEN(209,file=TRIM(output_dir)//'/bdy3_si.dat',form='unformatted',ACCESS='sequential', convert="big_endian")
    !     OPEN(300,file=TRIM(output_dir)//'/bdy4_si.dat',form='unformatted',ACCESS='sequential', convert="big_endian")

    !-------------------------------------------------------------------
    CALL compute_num_loop(num_loop, start_year, start_month, &
                          start_day, start_hour, &
                          end_year, end_month, &
                          end_day, end_hour, interval_seconds)

    PRINT *, 'Number of loops: ', num_loop
    !--------------------------------------------------------------
    ! Get grapesinput_new and store 5 full variables( u,v,th,pi,q)
    ! Revised by YaliWu for CMA-GD grapesinput
    CALL get_initial_value(pi, th, u, v, q, fid_init_in, &
                           ids, ide, jds, jde, kds, kde, &
                           ims, ime, jms, jme, kms, kme, &
                           its, ite, jts, jte, kts, kte)

    ! Galculate new bdys (e.g.pi_b of 6 spec_bdy_width) at the initial time (i.e. t1')
    CALL get_boundary_value(u, v, pi, th, q, &
                            u_b, v_b, pi_b, th_b, q_b, &
                            spec_bdy_width, &
                            ids, ide, jds, jde, kds, kde, &
                            ims, ime, jms, jme, kms, kme, &
                            its, ite, jts, jte, kts, kte)

    ! Get grapesbyd_old of pi_b & pi_bt (not stored, overwritten by the second-time pi_b)
    ! Get the second time bdy of pi_b
    ! Revised by YaliWu to get old pi_b & pi_b_t also
    CALL read_initial_boundary(mdate, &
                               u_b_2, v_b_2, pi_b_2, th_b_2, q_b_2, &
                               u_b_old, v_b_old, pi_b_old, th_b_old, q_b_old, &
                               u_bt_old, v_bt_old, pi_bt_old, th_bt_old, q_bt_old, &
                               interval_seconds, spec_bdy_width, fid_bdy_in, &
                               ids, ide, jds, jde, kds, kde, &
                               ims, ime, jms, jme, kms, kme, &
                               its, ite, jts, jte, kts, kte)

    ! Calculate time tendency, e.g., pi_b_t
    CALL get_boundary_tendency(spec_bdy_width, &
                               pi_b, th_b, u_b, v_b, q_b, &
                               pi_b_2, th_b_2, u_b_2, v_b_2, q_b_2, &
                               pi_bt, th_bt, u_bt, v_bt, q_bt, &
                               interval_seconds, &
                               ids, ide, jds, jde, kds, kde, &
                               ims, ime, jms, jme, kms, kme, &
                               its, ite, jts, jte, kts, kte)

    !CALL write_grapes_var_bdy(pi_b_old,th_b,u_b,v_b,q_b,           &
    !                        pi_bt_old,th_bt,u_bt,v_bt,q_bt,       &
    CALL write_grapes_var_bdy(pi_b, th_b, u_b, v_b, q_b, &
                              pi_bt, th_bt, u_bt, v_bt, q_bt, &
                              spec_bdy_width, mdate, fid_out, &
                              ids, ide, jds, jde, kds, kde, &
                              ims, ime, jms, jme, kms, kme, &
                              its, ite, jts, jte, kts, kte)
    ! For plot
    !     CALL write_grapes_bdy_grads( 203, 204, 205, 206, &
    !                                    pi_b,th_b,u_b,v_b,q_b,            &
    !                                    pi_bt,th_bt,u_bt,v_bt,q_bt,       &
    !                                    spec_bdy_width, mdate,        &
    !                                    ids,ide,jds,jde,kds,kde,          &
    !                                    ims,ime,jms,jme,kms,kme,          &
    !                                    its,ite,jts,jte,kts,kte)
    !
    !     CALL write_grapes_bdy_grads( 207, 208, 209, 300, &
    !                                    pi_b_old,th_b_old,u_b_old,v_b_old,q_b_old,            &
    !                                    pi_bt_old,th_bt_old,u_bt_old,v_bt_old,q_bt_old,       &
    !                                    spec_bdy_width, mdate,        &
    !                                    ids,ide,jds,jde,kds,kde,          &
    !                                    ims,ime,jms,jme,kms,kme,          &
    !                                    its,ite,jts,jte,kts,kte)
    !     CLOSE(203)
    !     CLOSE(204)
    !     CLOSE(205)
    !     CLOSE(206)
    !     CLOSE(207)
    !     CLOSE(208)
    !     CLOSE(209)
    !     CLOSE(300)

    PRINT *, 'check updated bdy values: th_bt ', MAXVAL(th_bt), MINVAL(th_bt)
    PRINT *, 'check updated bdy values: pi_bt ', MAXVAL(pi_bt), MINVAL(pi_bt)
    PRINT *, 'check updated bdy values: th_bt ', MAXVAL(th_bt), MINVAL(th_bt)
    PRINT *, 'check updated bdy values: u_bt ', MAXVAL(u_bt), MINVAL(u_bt)
    PRINT *, 'check updated bdy values: v_bt ', MAXVAL(v_bt), MINVAL(v_bt)
    PRINT *, 'check updated bdy values: q_bt ', MAXVAL(q_bt), MINVAL(q_bt)

    DO loop = 2, num_loop - 1
      CALL read_boundary_value(pi_b, th_b, u_b, v_b, q_b, &
                               pi_bt, th_bt, u_bt, v_bt, q_bt, &
                               mdate, spec_bdy_width, fid_bdy_in, &
                               ids, ide, jds, jde, kds, kde, &
                               ims, ime, jms, jme, kms, kme, &
                               its, ite, jts, jte, kts, kte)
      CALL write_grapes_var_bdy(pi_b, th_b, u_b, v_b, q_b, &
                                pi_bt, th_bt, u_bt, v_bt, q_bt, &
                                spec_bdy_width, mdate, fid_out, &
                                ids, ide, jds, jde, kds, kde, &
                                ims, ime, jms, jme, kms, kme, &
                                its, ite, jts, jte, kts, kte)

    END DO

  END SUBROUTINE update_bdy
END MODULE
