!
! Write Portable PixMap Image File (ppm, pbm, pgm)
! with Fortran
!
!   Charles O'Neill Oct 23, 2009
!    charles.oneill@gmail.com
!    www.caselab.okstate.edu
!
! Copyright (c) 2009 Charles O'Neill
!
! Permission is hereby granted, free of charge, to any person
! obtaining a copy of this software and associated documentation
! files (the "Software"), to deal in the Software without
! restriction, including without limitation the rights to use,
! copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the
! Software is furnished to do so, subject to the following
! conditions:
!
! The above copyright notice and this permission notice shall be
! included in all copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
! EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
! OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
! NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
! HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
! WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
! FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
! OTHER DEALINGS IN THE SOFTWARE.
!

MODULE PPM
  IMPLICIT NONE

CONTAINS

  !--------------------------------------------------------------
  ! Portable PixMap Type 1 (Black and White)
  SUBROUTINE writeppm1Matrix(M, text)
    INTEGER :: M(:, :)
    CHARACTER(len=*) :: text
    INTEGER :: cols, rows
    INTEGER :: i, j

    ! Open File
    OPEN (unit=100, file=TRIM(text)//".pbm", status='unknown')

    ! Write Header and ppm file type
    WRITE (100, '( A )') "P1"
    WRITE (100, '( A )') "# PPM Type 1 File (generated with fortran)"

    ! Write Image Size
    cols = SIZE(M, 2)
    rows = SIZE(M, 1)
    WRITE (100, '( i6, 1x, i6 )') cols, rows

    ! Write Image
    DO i = 1, rows
      DO j = 1, cols
        WRITE (100, '( i1 )', advance='no') M(i, j)
      END DO
      WRITE (100, *) ! Endline
    END DO
  END SUBROUTINE

  !--------------------------------------------------------------
  ! Portable PixMap Type 2 (Grayscale)
  SUBROUTINE writeppm2Matrix(M, text)
    INTEGER :: M(:, :)
    CHARACTER(len=*) :: text
    INTEGER :: cols, rows
    INTEGER :: i, j
    INTEGER :: maxvalue

    ! Open File
    OPEN (unit=100, file=TRIM(text)//".pgm", status='unknown')

    ! Write Header and ppm file type
    WRITE (100, '( A )') "P2"
    WRITE (100, '( A )') "# PPM Type 2 File (generated with fortran)"

    ! Write Image Size
    cols = SIZE(M, 2)
    rows = SIZE(M, 1)
    WRITE (100, '( i6, 1x, i6 )') cols, rows

    ! Write Maximum Value
    maxvalue = MAXVAL(MAXVAL(M, dim=1), dim=1)
    WRITE (100, '( i6 )') maxvalue

    ! Write Image
    DO i = 1, rows
      DO j = 1, cols
        WRITE (100, '( i5,1x )', advance='no') M(i, j)
      END DO
      WRITE (100, *) ! Endline
    END DO
  END SUBROUTINE

  !--------------------------------------------------------------
  ! Portable PixMap Type 3 (RGB color)
  SUBROUTINE writeppm3Matrix(R, G, B, text)
    INTEGER :: R(:, :), G(:, :), B(:, :)
    CHARACTER(len=*) :: text
    INTEGER :: cols, rows
    INTEGER :: i, j
    INTEGER :: maxvalue

    ! Open File
    OPEN (unit=100, file=TRIM(text)//".ppm", status='unknown')

    ! Write Header and ppm file type
    WRITE (100, '( A )') "P3"
    WRITE (100, '( A )') "# PPM Type 2 File (generated with fortran)"

    ! Write Image Size
    cols = SIZE(R, 2)
    rows = SIZE(R, 1)
    WRITE (100, '( i6, 1x, i6 )') cols, rows

    ! Write Maximum Value
    maxvalue = MAX(MAXVAL(MAXVAL(R, dim=1), dim=1) &
                   , MAXVAL(MAXVAL(G, dim=1), dim=1) &
                   , MAXVAL(MAXVAL(B, dim=1), dim=1))
    WRITE (6, *) ' Image maxvalue is ', maxvalue

!   Set brightness scaling in the image at a floor value
!   Image will always have a max >= 85 counts
!   if(maxvalue .lt. 64)then
!       maxvalue = 4 * maxvalue
    IF (maxvalue .LT. 85) THEN
      maxvalue = 3 * maxvalue
    ELSE
      maxvalue = 255
    END IF

    WRITE (6, *) ' Image scaled to   ', maxvalue
    WRITE (100, '( i6 )') maxvalue

    ! Write Image
    DO i = 1, rows
      DO j = 1, cols
        WRITE (100, '( 3(i5,1x) )') R(i, j), G(i, j), B(i, j)
      END DO
    END DO

    CLOSE (100)
  END SUBROUTINE

  !--------------------------------------------------------------
  ! Portable PixMap Type 3 (RGB color)
  SUBROUTINE writeppm3_16bit(R, G, B, text)
    INTEGER :: R(:, :), G(:, :), B(:, :)
    CHARACTER(len=*) :: text
    INTEGER :: cols, rows
    INTEGER :: i, j
    INTEGER :: maxvalue

    ! Open File
    OPEN (unit=100, file=TRIM(text)//".ppm", status='unknown')

    ! Write Header and ppm file type
    WRITE (100, '( A )') "P3"
    WRITE (100, '( A )') "# PPM Type 2 File (generated with fortran)"

    ! Write Image Size
    cols = SIZE(R, 2)
    rows = SIZE(R, 1)
    WRITE (100, '( i6, 1x, i6 )') cols, rows

    ! Write Maximum Value
    maxvalue = MAX(MAXVAL(MAXVAL(R, dim=1), dim=1) &
                   , MAXVAL(MAXVAL(G, dim=1), dim=1) &
                   , MAXVAL(MAXVAL(B, dim=1), dim=1))
    WRITE (6, *) ' Image maxvalue is ', maxvalue
    maxvalue = 65535
    WRITE (6, *) ' Image scaled to   ', maxvalue
    WRITE (100, '( i6 )') maxvalue

    ! Write Image
    DO i = 1, rows
      DO j = 1, cols
        WRITE (100, '( 3(i5,1x) )') R(i, j), G(i, j), B(i, j)
      END DO
    END DO
  END SUBROUTINE

  !--------------------------------------------------------------
  ! Test Module
  SUBROUTINE testPPM
    INTEGER, PARAMETER :: N = 100
    INTEGER :: A(N, N)
    INTEGER :: R(N, N)
    INTEGER :: G(N, N)
    INTEGER :: B(N, N)
    REAL :: AA(N, N)
    INTEGER :: i, j

    ! Show the PixMap Format with a simple case
    OPEN (unit=100, file="test.ppm", status='unknown')
    WRITE (100, '( A )') "P1"
    WRITE (100, '( A )') "# This is an example bitmap"
    WRITE (100, '( A )') "18 10"
    WRITE (100, '( A )') "0 0 0 0 0 0   0 0 0 0 0 0   0 0 0 0 0 0"
    WRITE (100, '( A )') "0 0 1 1 0 0   0 0 1 1 1 0   0 1 0 0 1 0"
    WRITE (100, '( A )') "0 1 0 0 1 0   0 1 0 0 1 0   0 1 0 0 1 0"
    WRITE (100, '( A )') "0 1 0 0 1 0   0 1 0 0 0 0   0 1 0 0 1 0"
    WRITE (100, '( A )') "0 1 0 0 1 0   0 0 1 0 0 0   0 1 0 0 1 0"
    WRITE (100, '( A )') "0 1 0 0 1 0   0 0 0 1 0 0   0 1 0 0 1 0"
    WRITE (100, '( A )') "0 1 0 0 1 0   0 0 0 0 1 0   0 1 0 0 1 0"
    WRITE (100, '( A )') "0 1 0 0 1 0   0 1 0 0 1 0   0 1 0 0 1 0"
    WRITE (100, '( A )') "0 0 1 1 0 0   0 1 1 1 0 0   0 0 1 1 0 0"
    WRITE (100, '( A )') "0 0 0 0 0 0   0 0 0 0 0 0   0 0 0 0 0 0"
    CLOSE (100)

    ! Get a matrix of random numbers
    CALL RANDOM_SEED()
    CALL RANDOM_NUMBER(AA)

    ! Setup and write type 1
    A = AA + 0.5
    CALL writeppm1Matrix(A, "subtest")
    WRITE (*, *) "writeppm1Matrix"

    ! Setup and write type 2
    A = AA * 256
    CALL writeppm2Matrix(A, "subtest2")
    WRITE (*, *) "writeppm2Matrix"

    ! Setup and write type 3
    DO i = 1, N
      DO j = 1, N
        R(i, j) = ABS(MOD(i + j, N / 5))
        G(i, j) = ABS(MOD(i - j, N / 5))
        B(i, j) = ABS(MOD(i * j, N / 5))
      END DO
    END DO
    CALL writeppm3Matrix(R, G, B, "subtest3")
    WRITE (*, *) "writeppm3Matrix"
  END SUBROUTINE

  SUBROUTINE read_ppm(u, img, ncol, iwidth, iheight)

    INTEGER ncol, iwidth, iheight
    INTEGER ncol2, iwidth2, iheight2

    INTEGER img(ncol, iwidth, iheight)
    INTEGER u
    CHARACTER(2) :: sign

    READ (u, '(A2)') sign
    READ (u, *) iwidth2, iheight2
    READ (u, *) ncol2

    WRITE (6, *) sign
    WRITE (6, *) iwidth2, iheight2
    WRITE (6, *) ncol2

    IF (ncol2 .NE. 255 .AND. ncol2 .NE. 65535) THEN
      WRITE (6, *) ' error in read_ppm - ncol2 appears incorrect: ', ncol2
      RETURN
    END IF

    READ (u, *) img

    WRITE (6, *) ' Image has been read in subroutine read_ppm'

  END SUBROUTINE

END MODULE
