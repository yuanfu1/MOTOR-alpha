!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Yuanfu Xie (yuanfu_xie@yahoo.com), 2025/7/3, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!!
!! @author Yuanfu Xie
!! @copyright (C) 2020 GBA-MWF, All rights reserved.
PROGRAM Test_VerticalInterp
  USE kinds_m, ONLY: r_kind, i_kind
  USE verticalInterp_m, ONLY: verticalInterp_t
  USE YAMLRead_m, ONLY: yaml_get_var
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE parameters_m, ONLY: machineEps, pi

  IMPLICIT NONE

  TYPE(mpddGlob_t), TARGET :: mpddGlob
  TYPE(verticalInterp_t) :: vertInterp

  CHARACTER(LEN=1024) :: yamlFile, hybridFile, staticDir
  INTEGER(i_kind) :: istatus, logOpt, nin, nout, ng, num_cell, lvl = 8
  INTEGER(i_kind) :: i, j, k, l, ierr, ndvd(2), im1, ip1
  REAL(r_kind) :: top, z, dz, a, b, am
  REAL(r_kind), PARAMETER :: off = 1.0D0 ! 100.0_r_kind
  REAL(r_kind), ALLOCATABLE :: fin(:), flg(:), vin(:), fout(:), flog(:), vout(:), &
    ana(:), analog(:), err(:), errlog(:)
  REAL(r_kind), ALLOCATABLE :: f2in(:,:), f2lg(:,:), v2in(:,:), f2out(:,:), f2log(:,:), v2out(:,:), &
    ana2(:,:), ana2log(:,:), err2(:), err2log(:)
  REAL(r_kind), ALLOCATABLE :: hybrida(:), hybridb(:), psfc(:)

  ! Get the configFile
  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", staticDir)
  yamlFile = TRIM(staticDir)//"/Application/App_4DVAR.yaml"

  ! Initialize the mpdd
  CALL mpddGlob%initialize()

  ! Initialize the vertical interpolation object
  ! CALL vertInterp%initialize(mpddGlob)

  istatus = 0 ! default pass number, yes
  ! Read in the configuration file
  ! ifile = yaml_get_var(configFile, 'IO', 'output_dir', vertInterp%outputDir)

  top = 30000.0D0

  ! Test case 1:
  nin = 5
  nout = 30
  nin = 4
  nout = 25
  logOpt = 1 ! flag for log interpolation

  ! Analytical vertical profile for comparison:
  ALLOCATE(ana(nout), analog(nout), fout(nout), flog(nout), vout(nout))
  DO i=1,nout
    z = REAL(i-1, r_kind) * top / REAL(nout-1, r_kind)+off
    vout(i) = z
    ana(i) = DSIN(pi * vout(i) / top)  ! Example vertical profile function
    analog(i) = 3.0D0 - 2.0D0 * (z / top)**2  ! Example analytical profile for comparison
  END DO

  ! Multigrid vertical interpolation test
  ALLOCATE(err(lvl), errlog(lvl))
  err = 0.0_r_kind
  errlog = 0.0_r_kind
  DO l=1,lvl

    nin = 2 * (nin-1) + 1

    ! vertical grid points:
    ng = 3 * (nin - 1) + 1

    PRINT*,'Test level', l, 'with vertical grid points:', ng, ' and output points:', nout

    ALLOCATE(fin(ng), flg(ng), vin(ng))

    ! 1 dimensional input and output arrays:
    ! Break the vertical into 2 sections, 1: from bottom to the middle, 2: from middle to the top
    ! The bottom section uses twice grid points as the top section
    dz = top / REAL(4*(nin-1)+1, r_kind)
    DO i=1, ng
      IF (i .LE. 2*(nin-1)+1) THEN
        z = REAL(i-1, r_kind) * dz
      ELSE
        z = 0.5D0 * top + REAL(i - 2*(nin-1) - 1) * 2.0D0 * dz
      END IF
      z = z + off
      vin(i) = z
      fin(i) = DSIN(pi * vin(i) / top)   ! Example vertical profile function

      flg(i) = EXP(3.0D0 - 2.0D0 * (vin(i) / top)**2)  ! Example analytical profile for comparison
    END DO

    ! Perform vertical interpolation test
    PRINT*,'Interpolating from', ng, 'to', nout, 'vertical levels at level', l, fin(1), vin(1)/top
    logOpt = 2
    CALL vertInterp%verticalInterp(flg, vin, flog, vout, istatus, logOpt)
    logOpt = 0
    CALL vertInterp%verticalInterp(fin, vin, fout, vout, istatus, logOpt)

    ! Debug output:
    IF (l .EQ. 8) THEN
      DO i=1, ng
        WRITE(10,2) i, flg(i), vin(i)/top
      END DO
      DO i=1, nout
        WRITE(20,2) i, flog(i) - analog(i), vout(i)/top  !
      END DO
    END IF
2   FORMAT('vlvl', I8, ': f: ', E12.6, ', z: ', E12.6)
    PRINT*,'Interpolation status:', istatus

    ! Calculate error:
    ierr = 0
    err(l) = 0.0_r_kind
    DO i=1, nout
      IF (ABS(fout(i) - ana(i)) .GT. err(l)) THEN
        ierr = i
        err(l) = ABS(fout(i) - ana(i))
      END IF
    END DO
    WRITE(*,5) l, err(l), MAXVAL(ABS(fout - ana)), ierr, nout
5   FORMAT('Error of test 1 at level', I2, ': ', E12.6, ' max error: ', E12.6, ' at ', I3, ' points out of ', I3)
    ! Log error:
    ierr = 0
    errlog(l) = 0.0_r_kind
    DO i=1, nout
      IF (ABS(flog(i) - analog(i)) .GT. errlog(l)) THEN
        ierr = i
        errlog(l) = ABS(flog(i) - analog(i))
      END IF
    END DO
    PRINT*,'Size of Log',size(flog),size(analog),flog(1),analog(1)
    WRITE(*,6) l, errlog(l), MAXVAL(ABS(flog - analog)), ierr, nout
6   FORMAT('Error of test 1 Log(f) at level', I2, ': ', E12.6, ' max error: ', E12.6, ' at ', I3, ' points out of ', I3)

    DEALLOCATE(fin, flg, vin)

  END DO

  DO l=2,lvl
    IF (err(l-1)/err(l) .LT. 3.6) THEN
      istatus = -1  ! Fail if the error does not decrease significantly
      WRITE(*,1) l, err(l-1)/err(l), err(l-1), err(l), mpddGlob%myrank
1     FORMAT('Std test 1 at', I2, ' failed. Err reduction ratio:', E12.4,' with errs:',2E12.4, ' pc: ',I1)
      EXIT
    ELSE
      IF (mpddGlob%isBaseProc()) THEN
        WRITE(*,3) l,  err(l-1)/err(l), err(l-1), err(l), mpddGlob%myrank
3       FORMAT('Std test 1 at', I2, ' passed. Err reduction ratio:', E12.4, ' with errs:',2E12.4, ' pc: ',I1)
      END IF
    END IF
  END DO
  PRINT*,''

  DO l=2,lvl
    ! Log error check:
    IF (errlog(l-1)/errlog(l) .LT. 3.6) THEN
      istatus = -1  ! Fail if the error does not decrease significantly
      WRITE(*,7) l, errlog(l-1)/errlog(l), errlog(l-1), errlog(l), mpddGlob%myrank
7     FORMAT('Log test 1 at', I2, ' failed. Err reduction ratio:', E12.4,' with errs:',2E12.4, ' pc: ',I1)
      EXIT
    ELSE
      IF (mpddGlob%isBaseProc()) THEN
        WRITE(*,8) l,  errlog(l-1)/errlog(l), errlog(l-1), errlog(l), mpddGlob%myrank
8       FORMAT('Log test 1 at', I2, ' passed. Err reduction ratio:', E12.4, ' with errs:',2E12.4, ' pc: ',I1)
      END IF
    END IF
  END DO

  IF (istatus .EQ. 0) THEN
    IF (mpddGlob%isBaseProc()) WRITE(*,*) 'Test passed'
  ELSE
    IF (mpddGlob%isBaseProc()) WRITE(*,*) 'Test failed with failure code:', istatus
  END IF

  DEALLOCATE(ana, analog, fout, flog, vout)

  !*************************************************************************************************
  ! Test case 2: two-dimensional vertical interpolation
  ! read in the configuration file for hybrid pressure levels:
  PRINT*, 'yamlFile:', TRIM(yamlFile)
  k = yaml_get_var(yamlFile, 'DASpace', 'hybridFile', hybridFile)
  IF (k .NE. 0) THEN
    PRINT*, 'Error reading hybrid file from yaml:', hybridFile
    istatus = -1
    RETURN
  END IF
  hybridFile = TRIM(staticDir)//'/'//TRIM(hybridFile)
  PRINT*, 'hybridFile:', TRIM(hybridFile)
  k = yaml_get_var(hybridFile, 'hybrid_levels', 'a', hybrida)
  k = yaml_get_var(hybridFile, 'hybrid_levels', 'b', hybridb)
  DO i=1, SIZE(hybrida)
    PRINT*, 'hybrid level', i, ': a = ', hybrida(i), ', b = ', hybridb(i)
  END DO

  num_cell = 1000
  nin = 5
  nout = 65
  am = 1.0D0    ! Test function sine amplitude
  ! Only test the first 65 output levels of the hybrid levels
  IF (SIZE(hybrida) .LT. nout) THEN
    PRINT*, 'Error: hybrid levels are less than nout:', SIZE(hybrida), nout
    istatus = -3
    GOTO 999
  END IF

  ! Since the finest level is the output, the multilevel test is set as follow:
  lvl = INT(DLOG(REAL((nout-1)/(nin-1), r_kind))/DLOG(2.0D0))-1
  PRINT*,'Test 2 total testing levels: ',lvl

  logOpt = 1 ! flag for log interpolation

  ! Analytical vertical profile for comparison:
  ALLOCATE(ana2(nout,num_cell), ana2log(nout,num_cell), f2out(nout,num_cell), f2log(nout,num_cell), &
    v2out(nout,num_cell), psfc(num_cell))

  ! Example surface pressure in Pa:
  psfc = 101325.0_r_kind  ! Example surface pressure in Pa
  top = hybrida(nout) + hybridb(nout) * psfc(1)  ! Top level based on the last hybrid level

  ! Initialize the vertical profile for each output level
  DO i=1,nout
    ! z = REAL(i-1, r_kind) * top / REAL(nout-1, r_kind)
    DO j=1,num_cell
      z = hybrida(i) + hybridb(i) * psfc(j)     ! Example vertical profile calculation
      IF (j .EQ. 1) WRITE(55,10) i,z, DLOG(z)
10    FORMAT(I3,2E14.6)
      v2out(i,j) = z

      ana2(i,j) = DSIN(am * pi * v2out(i,j) / top)  ! Example vertical profile function
      ana2log(i,j) = 3.0D0 - 2.0D0 * (z / top)**2  ! Example analytical profile for comparison
    END DO
  END DO

  ! Multigrid vertical interpolation test
  ALLOCATE(err2(lvl), err2log(lvl))
  err2 = 0.0_r_kind
  err2log = 0.0_r_kind
  DO l=1,lvl

    nin = 2 * (nin-1) + 1

    ndvd = 0
    IF (nin .GT. nout) THEN
      ndvd(1) = (nin-1)/(nout-1) ! Assume the nin and nou are in a relation of 2 times
    ELSE
      ndvd(2) = (nout-1)/(nin-1) ! Assume the nin and nou are in a relation of 2 times
    END IF

    PRINT*,'Test 2 at level', l, 'with vertical grid points:', nin, ' and output points:', nout,' top:',top

    ALLOCATE(f2in(nin,num_cell), f2lg(nin,num_cell), v2in(nin,num_cell))

    ! 2 dimensional input and output arrays:
    im1 = 1
    DO i=1, nin
      a = 1.0D0; b = 0.0D0
      IF (ndvd(1) .GT. 0) THEN
        im1 = INT((i-1) / ndvd(1)) + 1
        b = REAL(MOD(i-1, ndvd(1)),  r_kind) / REAL(ndvd(1), r_kind)
        a = 1.0D0 - b
      END IF
      IF (ndvd(2) .GT. 0) im1 = ndvd(2)*(i-1) + 1
      ip1 = im1 + 1
      IF (ip1 .GT. nout) ip1 = nout
      DO j=1,num_cell
        z = a * (hybrida(im1) + hybridb(im1) * psfc(j)) + &
            b * (hybrida(ip1) + hybridb(ip1) * psfc(j))  ! Example vertical profile calculation
        v2in(i,j) = off + z
        f2in(i,j) = DSIN(am * pi * v2in(i,j) / top)   ! Example vertical profile function
        f2lg(i,j) = EXP(3.0D0 - 2.0D0 * (v2in(i,j) / top)**2)  ! Example analytical profile for comparison
      END DO
    END DO

    ! Perform vertical interpolation test
    PRINT*,'Test 2 - interpolating from', nin, 'to', nout, 'vertical levels at level', l, f2in(1,1), v2in(1,1)/top
    logOpt = 2
    CALL vertInterp%verticalInterp(f2lg, v2in, f2log, v2out, istatus, logOpt)
    logOpt = 0
    CALL vertInterp%verticalInterp(f2in, v2in, f2out, v2out, istatus, logOpt)
    PRINT*,'Interpolation test 2 status:', istatus

    ! Debug output:
    IF (l .EQ. 2) THEN
      j = 1  ! Example for the first cell
      DO i=1, nin
        ! WRITE(12,22) i, f2lg(i,j), v2in(i,j)/top
        WRITE(12,22) i, f2lg(i,j), v2in(i,j)/top
      END DO
      DO i=1, nout
        !WRITE(22,22) i, f2log(i,j) - ana2log(i,j), v2out(i,j)/top
        WRITE(22,22) i, f2log(i,j)-ana2log(i,j), v2out(i,j)/top  !   
      END DO
    END IF
22  FORMAT('vlvl', I8, ': f: ', E12.6, ', z: ', E12.6)
    IF (l .EQ. 3) THEN
      j = 1  ! Example for the first cell
      DO i=1, nin
        ! WRITE(12,22) i, f2lg(i,j), v2in(i,j)/top
        WRITE(13,22) i, f2in(i,j), v2in(i,j)/top
      END DO
      DO i=1, nout
        !WRITE(22,22) i, f2log(i,j) - ana2log(i,j), v2out(i,j)/top
        WRITE(23,22) i, f2out(i,j)-ana2(i,j), v2out(i,j)/top  !   
      END DO
    END IF
    PRINT*,'Interpolation 2 status:', istatus

    ! Calculate error:
    ierr = 0
    err2(l) = 0.0_r_kind
    DO i=1, nout
      DO j=1, 1 !num_cell
        IF (ABS(f2out(i,j) - ana2(i,j)) .GT. err2(l)) THEN
          ierr = i
          err2(l) = ABS(f2out(i,j) - ana2(i,j))
        END IF
      END DO
    END DO
    WRITE(*,24) l, err2(l), MAXVAL(ABS(f2out - ana2)), ierr, nout
24  FORMAT('Error of test 2 at level', I2, ': ', E12.6, ' max error: ', E12.6, ' at ', I3, ' points out of ', I3)
    
    ! Log error:
    ierr = 0
    err2log(l) = 0.0_r_kind
    DO i=1, nout
      DO j=1, 1 !num_cell
        IF (ABS(f2log(i,j) - ana2log(i,j)) .GT. err2log(l)) THEN
          ierr = i
          err2log(l) = ABS(f2log(i,j) - ana2log(i,j))
        END IF
      END DO
    END DO
    WRITE(*,25) l, err2log(l), MAXVAL(ABS(f2log - ana2log)), ierr, nout
25  FORMAT('Error of test 2 in Log func at level', I2, ': ', E12.6, ' max error: ', E12.6, ' at ', I3, ' points out of ', I3)

    DEALLOCATE(f2in, f2lg, v2in)

  END DO

  PRINT*,'Error reduction in standard test 2'
  DO l=2,lvl
    IF (err2(l-1)/err2(l) .LT. 3.6) THEN
      istatus = -21  ! Fail if the error does not decrease significantly
      WRITE(*,21) l, err2(l-1)/err2(l), err2(l-1), err2(l), mpddGlob%myrank
21     FORMAT('Std test 2 at', I2, ' ** failed **. Err reduction ratio:', E12.4,' with errs:',2E12.4, ' pc: ',I1)
      EXIT
    ELSE
      IF (mpddGlob%isBaseProc()) THEN
        WRITE(*,23) l,  err2(l-1)/err2(l), err2(l-1), err2(l), mpddGlob%myrank
23       FORMAT('Std test 2 at', I2, ' passed. Err reduction ratio:', E12.4, ' with errs:',2E12.4, ' pc: ',I1)
      END IF
    END IF
  END DO
  PRINT*,''
  PRINT*,'Error reduction in Log test 2'

  DO l=2,lvl
    ! Log error check:
    IF (err2log(l-1)/err2log(l) .LT. 3.0 .AND. err2log(l-1) .GT. 1.0D2*machineEps) THEN
      istatus = -22  ! Fail if the error does not decrease significantly
      WRITE(*,27) l, err2log(l-1)/err2log(l), err2log(l-1), err2log(l), machineEps, mpddGlob%myrank
27     FORMAT('Log test 2 at', I2, ' ** failed **. Err reduction ratio:', E12.4,' with errs:',3E12.4, ' pc: ',I1)
      EXIT
    ELSE
      IF (mpddGlob%isBaseProc()) THEN
        WRITE(*,28) l,  err2log(l-1)/err2log(l), err2log(l-1), err2log(l), mpddGlob%myrank
28       FORMAT('Log test 2 at', I2, ' passed. Err reduction ratio:', E12.4, ' with errs:',2E12.4, ' pc: ',I1)
      END IF
    END IF
  END DO

999 CONTINUE

  IF (istatus .EQ. 0) THEN
    IF (mpddGlob%isBaseProc()) WRITE(*,*) 'Test passed'
  ELSE
    IF (mpddGlob%isBaseProc()) WRITE(*,*) 'Test failed with failure code:', istatus
  END IF

  DEALLOCATE(ana2, ana2log, f2out, f2log, v2out, psfc)

  CALL mpddGlob%finalize()  ! Finalize the mpdd

END PROGRAM Test_VerticalInterp