!!--------------------------------------------------------------------------------------------------
! PROJECT           : Satellite_utils
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Yali Wu
! VERSION           : V 0.0
! HISTORY           :
!   Created by Yali Wu (wuyali@gbamwf.com), 2022/06/23, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module implements a satellite radiance data structure.
MODULE Satellite_utils_m
  USE kinds_m
  USE parameters_m
  
CONTAINS

  SUBROUTINE Get_time_from_filename(filename, gmt_time)

    IMPLICIT NONE
    CHARACTER(len=*) :: filename
    INTEGER, INTENT(INOUT) :: gmt_time(:)

    !Local vars
    CHARACTER(len=:), ALLOCATABLE :: myString, lastSegment
    INTEGER :: lastSlashIndex, stringLength
    INTEGER :: i, ierr

    myString = TRIM(filename)
    stringLength = len(myString)

    ! Find the index of the last string
    lastSlashIndex = 0
    DO i = stringLength, 1, -1
      IF (myString(i:i) == '/') THEN
          lastSlashIndex = i
          exit
      END IF
    END DO

    ! abstract the last string
    IF (lastSlashIndex > 0) THEN
      lastSegment = TRIM(myString(lastSlashIndex+1:stringLength))
    ELSE
      lastSegment = TRIM(myString)
    END IF

    ! PRINT *, 'original string:', myString
    ! PRINT *, 'the last string:', lastSegment

    ! START get time from the filename
    READ(lastSegment(20:23), '(i4)', iostat=ierr) gmt_time(1)  ! Year
    READ(lastSegment(24:25), '(i2)', iostat=ierr) gmt_time(2)  ! Month
    READ(lastSegment(26:27), '(i2)', iostat=ierr) gmt_time(3)  ! Day 
    READ(lastSegment(29:30), '(i2)', iostat=ierr) gmt_time(4)  ! Hour
    READ(lastSegment(31:32), '(i2)', iostat=ierr) gmt_time(5)  ! Minute
    gmt_time(6) = 0
    ! PRINT *, 'Get_time_from_filename: ', gmt_time

  END SUBROUTINE Get_time_from_filename

  !> Y=A+BX
  SUBROUTINE Linear_Least_Squares(n, x, y, a, b, square_error)

    USE, INTRINSIC :: iso_fortran_env, only: dp => real64
    INTEGER, INTENT(IN) :: n
    REAL(dp), INTENT(IN) :: x(n), y(n)
    REAL(dp), INTENT(OUT) :: a, b
    REAL(dp), INTENT(OUT), OPTIONAL :: square_error  !! 平方误差

    REAL(dp) :: phi(2, 2), f(2, 1)  !! 最小二乘法系数矩阵和右端项
    INTEGER :: ipiv(2)  !! 工作数组
    INTEGER :: info     !! 返回值

    IF (n == 1) THEN

        a = y(1)
        b = 0.0_dp
        return

    END IF

    phi(1, 1) = n
    phi(2, 1) = sum(x)
    phi(1, 2) = phi(2, 1)
    phi(2, 2) = sum(x**2)

    f(1, 1) = sum(y)
    f(2, 1) = sum(x*y)

    ! A function from BLAS
    CALL dgesv(2, 1, phi, 2, ipiv, f, 2, info)
    a = f(1, 1)
    b = f(2, 1)
    
    IF (PRESENT(square_error)) THEN
        square_error = sum((y - a - b*x)**2)
    END IF

  END SUBROUTINE Linear_Least_Squares

  FUNCTION Convert_int2char(varin) RESULT(varout)
    IMPLICIT NONE
    INTEGER(i_kind)  :: varin 
    CHARACTER(len=10) :: varout
  
    IF (varin .GE. 1000 .AND. varin .LT. 10000) THEN
      write(varout,"(I4)") varin
    ELSEIF (varin .GE. 100) THEN 
      write (varout, "(I3)") varin 
    ELSEIF (varin .GE. 10) THEN 
      write (varout, "(I2)") varin
    ELSEIF (varin .GT. 0) THEN 
      write (varout, "(I1)") varin
    END IF
  
    END FUNCTION Convert_int2char

  SUBROUTINE Get_rttov_chan_info(satinfo, nchans, chan_lists, rttov_chan_lists, ifuse, Err, bias, &
                                  camin, camax, ErrClr, ErrCld, BTlim, only_over_sea, &
                                  lon, reso, col_offset, col_scale, radius_earth_a, radius_earth_b, satellite_h)
  CHARACTER(LEN=1024) :: satinfo_file
  INTEGER(i_kind)     :: ichan
  CHARACTER(len=*), INTENT(IN) :: satinfo
  INTEGER(i_kind), INTENT(OUT) :: nchans
  ! local vars
  INTEGER(i_kind), ALLOCATABLE :: a_chan_lists(:), a_rttov_chan_lists(:), a_ifuse(:)
  REAL(r_kind), ALLOCATABLE :: a_Err(:), a_bias(:), a_camin(:), a_camax(:), a_ErrClr(:), a_ErrCld(:), a_BTlim(:)
  INTEGER(i_kind), ALLOCATABLE  :: a_only_over_sea(:)
  REAL(r_kind) :: a_lon, a_reso
  REAL(r_kind) :: a_col_offset, a_col_scale, a_radius_earth_a, a_radius_earth_b, a_satellite_h

  INTEGER(i_kind), ALLOCATABLE, INTENT(INOUT), OPTIONAL :: chan_lists(:), rttov_chan_lists(:), ifuse(:)
  REAL(r_kind), ALLOCATABLE, INTENT(INOUT), OPTIONAL :: Err(:), bias(:), camin(:), camax(:), ErrClr(:), ErrCld(:), BTlim(:)
  INTEGER(i_kind), ALLOCATABLE, INTENT(INOUT), OPTIONAL :: only_over_sea(:)
  REAL(r_kind), INTENT(INOUT), OPTIONAL :: lon, reso
  REAL(r_kind), INTENT(INOUT), OPTIONAL :: col_offset, col_scale, radius_earth_a, radius_earth_b, satellite_h

  LOGICAL :: istat

  ! Get environments
  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", satinfo_file)
  satinfo_file = TRIM(satinfo_file)//"/Satellite/satinfo_"//TRIM(satinfo)//".txt"
  PRINT *, 'satinfo_file in Satellite_utils = ', TRIM(satinfo_file)

  inquire (file=satinfo_file, exist=istat)
    IF (istat) THEN
      open (unit=121, file=satinfo_file, status='old')
      read (121, *)
      read (121, *) nchans
     
      ALLOCATE (a_chan_lists(nchans), a_rttov_chan_lists(nchans), a_ifuse(nchans), a_only_over_sea(nchans))
      ALLOCATE(a_Err(nchans), a_bias(nchans), a_camin(nchans), a_camax(nchans), a_ErrClr(nchans), a_ErrCld(nchans), a_BTlim(nchans))

      DO ichan = 1, nchans
        read (121, *) a_chan_lists(ichan), a_rttov_chan_lists(ichan), a_ifuse(ichan), &
        a_Err(ichan), a_bias(ichan), a_camin(ichan), a_camax(ichan), &
        a_ErrClr(ichan), a_ErrCld(ichan), a_BTlim(ichan), a_only_over_sea(ichan)
      END DO

      IF (PRESENT(lon)) THEN
        read (121, *)
        read (121, *)
        read (121, *)
        read (121, *) a_lon
        read (121, *)
        read (121, *) a_reso
        read (121, *)
        read (121, *) a_col_offset, a_col_scale, a_radius_earth_a, a_radius_earth_b, a_satellite_h
      ELSE
        a_lon = missing
        a_reso = missing
        a_col_offset = missing
        a_col_scale = missing
        a_radius_earth_a = missing
        a_radius_earth_b = missing
        a_satellite_h = missing
      END IF

      close (121)
    ELSE
      PRINT *, '----------satinfo file was not found in Satellite_utils----------'
      STOP
    END IF

    IF (PRESENT(chan_lists)) chan_lists=a_chan_lists
    IF (PRESENT(rttov_chan_lists)) rttov_chan_lists=a_rttov_chan_lists
    IF (PRESENT(ifuse)) ifuse=a_ifuse
    IF (PRESENT(Err)) Err=a_Err
    IF (PRESENT(bias)) bias=a_bias
    IF (PRESENT(camin)) camin=a_camin
    IF (PRESENT(camax)) camax=a_camax
    IF (PRESENT(ErrClr)) ErrClr=a_ErrClr
    IF (PRESENT(ErrCld)) ErrCld=a_ErrCld
    IF (PRESENT(BTlim)) BTlim=a_BTlim
    IF (PRESENT(only_over_sea)) only_over_sea=a_only_over_sea
    IF (PRESENT(lon)) lon=a_lon
    IF (PRESENT(reso)) reso=a_reso
    IF (PRESENT(col_offset)) col_offset=a_col_offset
    IF (PRESENT(col_scale)) col_scale=a_col_scale
    IF (PRESENT(radius_earth_a)) radius_earth_a=a_radius_earth_a
    IF (PRESENT(radius_earth_b)) radius_earth_b=a_radius_earth_b
    IF (PRESENT(satellite_h)) satellite_h=a_satellite_h

    DEALLOCATE (a_chan_lists, a_rttov_chan_lists, a_ifuse, a_only_over_sea)
    DEALLOCATE(a_Err, a_bias, a_camin, a_camax, a_ErrClr, a_ErrCld, a_BTlim)

  END SUBROUTINE Get_rttov_chan_info
  
  SUBROUTINE Get_rttov_giirs_chan_info(nchans, nchan_used, ifuse, chan_idx, satinfo, chan_lists, rttov_chan_lists)
    CHARACTER(LEN=1024) :: satinfo_file
    INTEGER(i_kind)     :: ichan
    CHARACTER(len=*), INTENT(IN) :: satinfo
    INTEGER(i_kind), INTENT(OUT) :: nchans, nchan_used
    INTEGER(i_kind), ALLOCATABLE, INTENT(OUT) :: ifuse(:), chan_lists(:), rttov_chan_lists(:), chan_idx(:) 
    INTEGER(i_kind), ALLOCATABLE :: chan_lists_used(:)
    INTEGER(i_kind)   :: NCH,II

    LOGICAL :: istat

    ! Get environments
    CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", satinfo_file)
    satinfo_file = TRIM(satinfo_file)//"/Satellite/satinfo_"//TRIM(satinfo)//".txt"

    inquire (file=satinfo_file, exist=istat)
      IF (istat) THEN
        open (unit=121, file=satinfo_file, status='old')
        read (121, *)
        read (121, *) nchans
        ALLOCATE (chan_lists(nchans))
        ALLOCATE (chan_lists_used(nchans))
        ALLOCATE (ifuse(nchans))
        ALLOCATE (chan_idx(nchans))

        NCH = 0
        DO ichan = 1, nchans
          READ (121,*) chan_lists(ichan),II,ifuse(ichan)
          IF(ifuse(ichan).GT.0)THEN                      !! =1: use, / 0: not use
            NCH = NCH + 1
            chan_idx(NCH) = ichan
            chan_lists_used(NCH) = chan_lists(ichan)
            !print *,ichan,NCH,chan_idx(NCH),"chan_lists:",chan_lists(ichan),chan_lists_used(NCH),ifuse(ichan)
          ENDIF 
        END DO
        nchan_used = NCH 
        ALLOCATE (rttov_chan_lists(nchan_used))
        rttov_chan_lists(1:nchan_used) = chan_lists_used(1:nchan_used)
        close (121)
      ELSE
        PRINT *, '----------satinfo file was not found in Satellite_utils----------'
        STOP
      END IF
      DEALLOCATE(chan_lists_used)

  END SUBROUTINE  Get_rttov_giirs_chan_info

  SUBROUTINE findhIdx_all(obsData, value, Idx)
    ! only works for 1D interger obsDatas
    IMPLICIT NONE 
    INTEGER(i_kind), INTENT(IN) :: obsData(:)
    INTEGER(i_kind), INTENT(IN) :: value
    INTEGER(i_kind), INTENT(INOUT) :: Idx(:)

    INTEGER :: nvalues, ndx, i,j 

    nvalues = SIZE(obsData, 1)
    ndx = 0
    DO i = 1, nvalues
      IF (obsData(i) .EQ. value) THEN 
        ndx = ndx + 1
        Idx(ndx) = i
      END IF 
    END DO 

  END SUBROUTINE findhIdx_all

  SUBROUTINE findhIdx_ht(obsData_h, obsData_t, VALUE_h, VALUE_t, Idx)
    ! only works for 1D interger obsDatas
    IMPLICIT NONE
    INTEGER(i_kind), INTENT(IN) :: obsData_h(:), obsData_t(:)
    INTEGER(i_kind), INTENT(IN) :: VALUE_h, VALUE_t
    INTEGER(i_kind), INTENT(INOUT) :: Idx

    INTEGER :: nvalues, ndx, i, j

    nvalues = SIZE(obsData_h, 1)
    DO i = 1, nvalues
      IF (obsData_h(i) .EQ. VALUE_h .AND. obsData_t(i) .EQ. VALUE_t) THEN
        Idx = i
      END IF
    END DO

  END SUBROUTINE findhIdx_ht

END MODULE Satellite_utils_m
