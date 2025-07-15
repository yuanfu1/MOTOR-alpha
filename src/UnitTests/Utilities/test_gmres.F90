
!!---------------------------------------------------------------------------------------
! PROJECT           : Utility.test_gmres.F90
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center 
!                     for Monitoring Warning and Forecasting (GBA-MWF)
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : V 0.1
! HISTORY           :
!   Created  by Yuanfu Xie, 2024/03 @GBA-MWF, Shenzhen
!!---------------------------------------------------------------------------------------

!> @brief
!! A unit test of the least norm solver of GMRES.
!! It is mainly used for wider interpolation stencil problem that is underdetermined.
!! This package allows users to solve
!!      Ax = b
!! where A is an asymmetric and rank deficient matrix.
!! @author Yuanfu Xie
!! @copyright (C) 2024 GBA-MWF, All rights reserved.
!! @warning
!! @attention
PROGRAM test_gmres

    USE kinds_m, ONLY : i_kind, r_kind

    IMPLICIT NONE

    INTEGER(i_kind) :: i,j,k,m,n,its,info,icase,ioff(2),nn,mm
    REAL(r_kind) :: res,del,xx,yy
    REAL(r_kind), ALLOCATABLE :: x(:),y(:),z(:),b(:),h(:,:),v(:,:)

    EXTERNAL :: ddot,dotprd,psolve,matvec,matvec_Poisson1d,matvec_Poisson2d,matvec_Poisson2dMixed

    ! Test of call of dot product:
    ! dprd = DDOT(x,b)

    icase = 4   ! 1. Simple linear system test; 2. 1D Poisson; 3. 2D Poisson

    ! Parameter setting:
    SELECT CASE (icase)
    CASE(1)
        n = 3
        m = 3
        its = 3     ! on input: max num its; 
    CASE(2)
        n = 169
        m = 169
    CASE(3)
        n = 6561
        m = 6561
        its = 169     ! on input: max num its; 
    CASE(4)
        n = 1681*2
        m = 1681*2
        its = 30000     ! on input: max num its; 
        ioff(1) = 0
        ioff(2) = n/2

        IF (MOD(n,2) .NE. 0) THEN
            PRINT*,'N must be an even number for a mixed Poisson solver!'
            STOP
        END IF
        IF (MOD(DSQRT(DBLE(n/2)),1.0D0) .NE. 0.0D0) THEN
            PRINT*,'N/2 must be a square of an integer for the solver!'
            STOP
        END IF
        mm = DSQRT(DBLE(n/2))
    CASE DEFAULT
        WRITE(*,1)
1       FORMAT('This test case has not implemented. Check and rerun')
    END SELECT

    ! Allocate memory for x, b, h and v:
    ALLOCATE(x(n),y(n),z(n),b(n),h(m+1,m),v(n,m+1))

    x = 0.1D0
    b = 1.0D0
    SELECT CASE (icase)
    CASE(1)
        b(3) = 2.0D0
    CASE(2)
    CASE(3)
    CASE(4)
        ! Set up boundary conditions:
        ! Assume psi = chi = (x^2 + y^2)/4
        ! u = -psi_y + chi_x = (-y+x)/2
        ! v =  psi_x + chi_y = ( x+y)/2
        ! Top and bottom:
        DO j=1,mm,mm-1
            yy = DBLE(j-1)/DBLE(mm-1)-0.5D0
            DO i=1,mm
                xx = DBLE(i-1)/DBLE(mm-1)-0.5D0
                b(i+mm*(j-1)+ioff(1)) = 0.5D0*(-yy+xx)
                b(i+mm*(j-1)+ioff(2)) = 0.5D0*( yy+xx)
            END DO
        END DO
        ! Left and right:
        DO j=2,mm-1
            yy = DBLE(j-1)/DBLE(mm-1)-0.5D0
            DO i=1,mm,mm-1
                xx = DBLE(i-1)/DBLE(mm-1)-0.5D0
                b(i+mm*(j-1)+ioff(1)) = 0.5D0*(-yy+xx)
                b(i+mm*(j-1)+ioff(2)) = 0.5D0*( yy+xx)
            END DO
        END DO
    END SELECT

    res = 1.0D-5 ! on input: |Ax-b|/|b|<res;
                ! on exit:  residual reached
    del = 1.0D1 ! on input: if(del>0) then the x returned is the hookstep
                ! on exit:  norm of next b predicted by hook
                ! on exit:  number of its taken
    info = 1    ! on input: if(info==1) print* residuals
                ! if(info==2) recalc hookstep with new del
                ! on exit:  0 sucessful, 1 method breakdown, 2 max its

    WRITE(*,2)
2   FORMAT('Calling GMRES...')

    ! Calling GMRES with different matrix-vector multiplier matvec*:
    SELECT CASE (icase)
    CASE(1)
        CALL gmresm(m,n,x,b,matvec,psolve,dotprd,h,v,res,del,its,info)
    CASE(2)
        CALL gmresm(m,n,x,b,matvec_Poisson1d,psolve,dotprd,h,v,res,del,its,info)
        CALL nullProjection1d(n,x,y)
        x = y
        DO i=1,n
            y(i) = 0.5D0*((DBLE(i-1)/DBLE(n-1)-0.0D0)**2)
        END DO
        CALL nullProjection1d(n,y,z)
    CASE(3)
        ! Test:
        ! CALL matvec_Poisson2d(n,x,y)
        CALL gmresm(m,n,x,b,matvec_Poisson2d,psolve,dotprd,h,v,res,del,its,info)
        CALL nullProjection2d(n,x,y)
        x = y
        m = INT(DSQRT(DBLE(n)))
        k = 0
        DO j=1,m
            DO i=1,m
                k = k+1
                y(k) = 0.25D0*((DBLE(i-1)/DBLE(m-1)-0.0D0)**2+(DBLE(j-1)/DBLE(m-1)-0.0D0)**2)
            END DO
        END DO
        CALL nullProjection2d(n,y,z)
    CASE(4)
        ! Test: Streamfunction and velocity potential mixed solver:
        CALL gmresm(m,n,x,b,matvec_Poisson2dMixed,psolve,dotprd,h,v,res,del,its,info)
        CALL nullProjection2dMixed(n,x,y)
        x = y
        mm = INT(DSQRT(DBLE(n/2)))
        k = 0
        DO j=1,mm
            DO i=1,mm
                k = k+1
                y(k+ioff(1)) = 0.25D0*((DBLE(i-1)/DBLE(mm-1)-0.0D0)**2+(DBLE(j-1)/DBLE(mm-1)-0.0D0)**2)
                y(k+ioff(2)) = 0.25D0*((DBLE(i-1)/DBLE(mm-1)-0.0D0)**2+(DBLE(j-1)/DBLE(mm-1)-0.0D0)**2)
            END DO
        END DO
        CALL nullProjection2dMixed(n,y,z)
    END SELECT

    WRITE(*,3) info,its,res,DSQRT(DOT_PRODUCT(x,x)),DSQRT(DOT_PRODUCT(z,z))
3   FORMAT('GMRES: success: ',I2,' Iterations: ',I6,' Residue: ',E12.4,' Solution norms: ',2E12.4)

    ! Save data to plot:
    SELECT CASE (icase)
    CASE(2)
        ! 1D poisson solution:
        DO i=1,n
            WRITE(10,4) i,x(i),z(i),x(i)-z(i)
        END DO
4       FORMAT('Solution ',I4,' : ',5E12.4)
    CASE(3)
        ! 2D Poisson solution:
        m = INT(DSQRT(DBLE(n)))
        k = 0
        DO j=1,m
            yy = DBLE(j-1)/DBLE(m-1)-0.5D0
            DO i=1,m
                xx = DBLE(i-1)/DBLE(m-1)-0.5D0
                k = k+1
                WRITE(11,5) yy,xx,x(k)
            END DO
            WRITE(11,*) ! Space to break the gnuplot lines
        END DO
5       FORMAT(3E14.6)

        k = 0
        DO j=1,m
            yy = DBLE(j-1)/DBLE(m-1)-0.5D0
            DO i=1,m
                xx = DBLE(i-1)/DBLE(m-1)-0.5D0
                k = k+1
                WRITE(12,5) yy,xx,z(k)-x(k)
            END DO
            WRITE(12,*) ! Space to break the gnuplot lines
        END DO
    CASE(4)
        ! 2D Poisson solution:
        mm = INT(DSQRT(DBLE(n/2)))
        k = 0
        DO j=1,mm
            yy = DBLE(j-1)/DBLE(mm-1)-0.5D0
            DO i=1,mm
                xx = DBLE(i-1)/DBLE(mm-1)-0.5D0
                k = k+1
                WRITE(13,6) yy,xx,x(k+ioff(1))
                WRITE(14,6) yy,xx,x(k+ioff(2))
                WRITE(15,6) yy,xx,y(k)
                WRITE(16,6) yy,xx,z(k)
            END DO
            WRITE(13,*) ! Space to break the gnuplot lines
            WRITE(14,*)
            WRITE(15,*)
            WRITE(16,*)
        END DO
6       FORMAT(3E14.6)
    END SELECT

    ! Deallocate memory:
    DEALLOCATE(x,y,z,b,h,v)
END PROGRAM test_gmres

SUBROUTINE matvec(n,x,y)

    USE kinds_m, ONLY : i_kind, r_kind

    INTEGER(i_kind), INTENT(IN) :: n
    REAL(r_kind), INTENT(IN) :: x(n)
    REAL(r_kind), INTENT(OUT) :: y(n)

    ! Test case 1: 2x_1 + x_2 = 1; 3x_1 + 2x_2 = 1; 4x_1 + 3x_2 = 1:
    y(1) = 2.0D0*x(1)+1.0D0*x(2)
    y(2) = 3.0D0*x(1)+2.0D0*x(2)
    y(3) = 4.0D0*x(1)+3.0D0*x(2)
END SUBROUTINE matvec

SUBROUTINE matvec_Poisson1d(n,x,y)

    USE kinds_m, ONLY : i_kind, r_kind

    INTEGER(i_kind), INTENT(IN) :: n
    REAL(r_kind), INTENT(IN) :: x(n)
    REAL(r_kind), INTENT(OUT) :: y(n)

    ! Local variables:
    INTEGER(i_kind) :: i
    REAL(r_kind) :: dx

    dx = 1.0D0/DBLE(n-1)

    DO i=1,n
        IF (isnan(x(i))) THEN
            PRINT*,'NaN encountered: ',i,x(i)
            stop
        END IF
    END DO

    ! Test case 2: Poisson equations:
    y(1) = (2.0D0*x(1)-5.0D0*x(2)  +4.0D0*x(3)  -x(4))/(dx*dx)
    y(n) = (2.0D0*x(n)-5.0D0*x(n-1)+4.0D0*x(n-2)-x(n-3))/(dx*dx)

    DO i=2,n-1
        y(i) = (x(i+1)-2.0D0*x(i)+x(i-1))/(dx*dx)
    END DO
END SUBROUTINE matvec_Poisson1d

SUBROUTINE nullProjection1d(n,x,y)

    USE kinds_m, ONLY : i_kind, r_kind

    INTEGER(i_kind), INTENT(IN) :: n
    REAL(r_kind), INTENT(IN) :: x(n)
    REAL(r_kind), INTENT(OUT) :: y(n)

    ! Local variables:
    INTEGER(i_kind) :: i
    REAL(r_kind) :: z(n),d(n)

    ! Project out the constant:
    d = x-SUM(x)/DBLE(n)

    ! Project out the linear function:
    DO i=1,n
        z(i) = DBLE(i-1)/DBLE(n-1)-0.5D0
    END DO
    y = d - DOT_PRODUCT(d,z)/DOT_PRODUCT(z,z)*z

END SUBROUTINE nullProjection1d

SUBROUTINE matvec_Poisson2d(n,x,y)

    USE kinds_m, ONLY : i_kind, r_kind

    IMPLICIT NONE

    INTEGER(i_kind), INTENT(IN) :: n
    REAL(r_kind), INTENT(IN) :: x(n)
    REAL(r_kind), INTENT(OUT) :: y(n)

    ! Local variables:
    INTEGER(i_kind) :: i,j,m
    REAL(r_kind) :: dxdy

    ! PRINT*,'N = ',n,MOD(DSQRT(DBLE(n)),1.0D0),DSQRT(DBLE(n))
    IF (MOD(DSQRT(DBLE(n)),1.0D0) .NE. 0) THEN
        PRINT*,'The N is not a square of an integer, quit!'
        PRINT*,'Current test is on a uniform grid in x y directions'
        STOP
    END IF

    m = DSQRT(DBLE(n))

    dxdy = 1.0D0/DBLE(m-1)

    ! Debugging:
    DO i=1,n
        IF (ISNAN(x(i))) THEN
            WRITE(*,2) i,n,x(i),dxdy
2           FORMAT('Found input NaN: ',2I5,2E12.4)
            STOP
        END IF
    END DO

    ! Test case 3: 2D Poisson equations:
    ! Bottom left corner:
    y(1) = (2.0D0*x(1)-5.0D0*x(2)+4.0D0*x(3)-x(4))/(dxdy*dxdy)+ &
    (2.0D0*x(1)-5.0D0*x(m+1)+4.0D0*x(2*m+1)-x(3*m+1))/(dxdy*dxdy)
    ! Bottom right corner:
    y(m) = (2.0D0*x(m)-5.0D0*x(m-1)+4.0D0*x(m-2)-x(m-3))/(dxdy*dxdy)+ &
    (2.0D0*x(m)-5.0D0*x(2*m)+4.0D0*x(3*m)-x(4*m))/(dxdy*dxdy)
    ! Top left corner:
    y(n-m+1) = (2.0D0*x(n-m+1)-5.0D0*x(n-m+2)+4.0D0*x(n-m+3)-x(n-m+4))/(dxdy*dxdy)+ &
    (2.0D0*x(n-m+1)-5.0D0*x(n-2*m+1)+4.0D0*x(n-3*m+1)-x(n-4*m+1))/(dxdy*dxdy)
    ! Top right corner:
    y(n) = (2.0D0*x(n)-5.0D0*x(n-1)+4.0D0*x(n-2)-x(n-3))/(dxdy*dxdy)+ &
    (2.0D0*x(n)-5.0D0*x(n-m)+4.0D0*x(n-2*m)-x(n-3*m))/(dxdy*dxdy)

    ! Boundaries:
    ! Bottom boundary:
    DO i=2,m-1
        y(i) = (x(i+1)-2.0D0*x(i)+x(i-1))/(dxdy*dxdy)+ &
        (2.0D0*x(i)-5.0D0*x(i+m)+4.0D0*x(i+2*m)-x(i+3*m))/(dxdy*dxdy)
    END DO
    ! Top boundary:
    DO i=n-m+2,n-1
        y(i) = (x(i+1)-2.0D0*x(i)+x(i-1))/(dxdy*dxdy)+ &
        (2.0D0*x(i)-5.0D0*x(i-m)+4.0D0*x(i-2*m)-x(i-3*m))/(dxdy*dxdy)
    END DO
    ! Left boundary:
    DO i=m+1,n-2*m+1,m
        y(i) = (2.0D0*x(i)-5.0D0*x(i+1)+4.0D0*x(i+2)-x(i+3))/(dxdy*dxdy)+ &
        (x(i+m)-2.0D0*x(i)+x(i-m))/(dxdy*dxdy)
    END DO
    ! Right boundary:
    DO i=2*m,n-m,m
        y(i) = (2.0D0*x(i)-5.0D0*x(i-1)+4.0D0*x(i-2)-x(i-3))/(dxdy*dxdy)+ &
        (x(i+m)-2.0D0*x(i)+x(i-m))/(dxdy*dxdy)
    END DO

    DO j=2,m-1
        DO i=2,m-1
            y(i+m*(j-1)) = (x(i+m*(j-1)+1)-2.0D0*x(i+m*(j-1))+x(i+m*(j-1)-1))/(dxdy*dxdy)+ &
                           (x(i+m*(j-1)+m)-2.0D0*x(i+m*(j-1))+x(i+m*(j-1)-m))/(dxdy*dxdy)
        END DO
    END DO

    ! Debugging:
    DO i=1,n
        IF (ISNAN(y(i))) THEN
            WRITE(*,1) i,y(i),x(i)
1           FORMAT('Found NaN: ',I5,2E12.4)
        END IF
    END DO
END SUBROUTINE matvec_Poisson2d

SUBROUTINE nullProjection2d(n,x,y)

    USE kinds_m, ONLY : i_kind, r_kind

    INTEGER(i_kind), INTENT(IN) :: n
    REAL(r_kind), INTENT(IN) :: x(n)
    REAL(r_kind), INTENT(OUT) :: y(n)

    ! Local variables:
    INTEGER(i_kind) :: i,j,k,m
    REAL(r_kind) :: s(n),z(n),d(n),xx,yy

    ! Project out the constant:
    d = x-SUM(x)/DBLE(n)
    y = d
    PRINT*,'MaxMin values: ',maxval(d),minval(d)
    PRINT*,'MaxMin values: ',maxval(x),minval(x)
    return

    ! Project out the linear function:
    k = 0
    DO j=1,m
        DO i=1,m
            k = k+1
            s(k) = DBLE(i-1)/DBLE(m-1)-0.5D0
            z(k) = DBLE(j-1)/DBLE(m-1)-0.5D0
        END DO
    END DO
    d = y - DOT_PRODUCT(y,s)/DOT_PRODUCT(s,s)*s
    y = d - DOT_PRODUCT(d,z)/DOT_PRODUCT(z,z)*z

END SUBROUTINE nullProjection2d

!!---------------------------------------------------------------------------------------
!  This routine loads the coefficient matrix multiplying a vector x of streamfunction and
!  velocity potential combined solver with mixed boundary conditions:
!   -psi_y + chi_x = u
!    psi_x + chi_y = v
!  Yuanfu Xie implemented as a test 2024-03-25
!!---------------------------------------------------------------------------------------
SUBROUTINE matvec_Poisson2dMixed(nn,x,y)

    USE kinds_m, ONLY : i_kind, r_kind

    IMPLICIT NONE

    INTEGER(i_kind), INTENT(IN) :: nn
    REAL(r_kind), INTENT(IN) :: x(nn)
    REAL(r_kind), INTENT(OUT) :: y(nn)

    ! Local variables:
    INTEGER(i_kind) :: i,j,k,l,n,m,ioff(2)
    REAL(r_kind) :: dxdy,uv(nn)

    ! Check if nn is an even number:
    IF (MOD(nn,2) .NE. 0) THEN
        PRINT*,'This is a test of mixed Poisson solver, nn is required to be an even number: ',nn
        STOP
    ELSE
        n = nn/2
        ioff(1) = 0
        ioff(2) = n
    END IF

    IF (MOD(DSQRT(DBLE(n)),1.0D0) .NE. 0) THEN
        PRINT*,'The N is not a square of an integer, quit!'
        PRINT*,'Current test is on a uniform grid in x y directions'
        STOP
    END IF

    m = DSQRT(DBLE(n))

    dxdy = 1.0D0/DBLE(m-1)

    ! Debugging:
    DO i=1,n
        IF (ISNAN(x(i))) THEN
            WRITE(*,2) i,n,x(i),dxdy
2           FORMAT('Found input NaN: ',2I5,2E12.4)
            STOP
        END IF
    END DO

    ! Test case 4: 2D Mixed Poisson equations:
    ! Bottom left corner: -psi_y + chi_x = u; pis_x + chi_y = v;
    y(1+ioff(1)) = ( 3.0D0*x(1+ioff(1))-4.0D0*x(m+1+ioff(1))+x(2*m+1+ioff(1)))/(2.0D0*dxdy)+ &
                   (-3.0D0*x(1+ioff(2))+4.0D0*x(2  +ioff(2))-x(3    +ioff(2)))/(2.0D0*dxdy)
    y(1+ioff(2)) = (-3.0D0*x(1+ioff(2))+4.0D0*x(m+1+ioff(2))-x(2*m+1+ioff(2)))/(2.0D0*dxdy)+ &
                   (-3.0D0*x(1+ioff(1))+4.0D0*x(2  +ioff(1))-x(3    +ioff(1)))/(2.0D0*dxdy)
    ! Bottom right corner:
    y(m+ioff(1)) = ( 3.0D0*x(m+ioff(1))-4.0D0*x(2*m+ioff(1))+x(3*m+ioff(1)))/(2.0D0*dxdy)+ &
                   ( 3.0D0*x(m+ioff(2))-4.0D0*x(m-1+ioff(2))+x(m-2+ioff(2)))/(2.0D0*dxdy)
    y(m+ioff(2)) = (-3.0D0*x(m+ioff(2))+4.0D0*x(2*m+ioff(2))-x(3*m+ioff(2)))/(2.0D0*dxdy)+ &
                   ( 3.0D0*x(m+ioff(1))-4.0D0*x(m-1+ioff(1))+x(m-2+ioff(1)))/(2.0D0*dxdy)
    ! Top left corner:
    l = n-m+1
    y(l+ioff(1)) = (-3.0D0*x(l+ioff(1))+4.0D0*x(l-m+ioff(1))-x(l-2*m+ioff(1)))/(2.0D0*dxdy)+ &
                   (-3.0D0*x(l+ioff(2))+4.0D0*x(l+1+ioff(2))-x(l+2  +ioff(2)))/(2.0D0*dxdy)
    y(l+ioff(2)) = ( 3.0D0*x(l+ioff(2))-4.0D0*x(l-m+ioff(2))+x(l-2*m+ioff(2)))/(2.0D0*dxdy)+ &
                   (-3.0D0*x(l+ioff(1))+4.0D0*x(l+1+ioff(1))-x(l+2  +ioff(1)))/(2.0D0*dxdy)
    ! Top right corner:
    y(n+ioff(1)) = (-3.0D0*x(n+ioff(1))+4.0D0*x(n-m+ioff(1))-x(n-2*m+ioff(1)))/(2.0D0*dxdy)+ &
                   ( 3.0D0*x(n+ioff(2))-4.0D0*x(n-1+ioff(2))+x(n-2  +ioff(2)))/(2.0D0*dxdy)
    y(n+ioff(2)) = ( 3.0D0*x(n+ioff(2))-4.0D0*x(n-m+ioff(2))+x(n-2*m+ioff(2)))/(2.0D0*dxdy)+ &
                   ( 3.0D0*x(n+ioff(1))-4.0D0*x(n-1+ioff(1))+x(n-2  +ioff(1)))/(2.0D0*dxdy)

    ! Boundaries:
    ! Bottom boundary:
    DO i=2,m-1
        y(i+ioff(1)) = (x(i+1+ioff(2))-x(i-1+ioff(2)))/(2.0*dxdy)+ &
            ( 3.0D0*x(i+ioff(1))-4.0D0*x(i+m+ioff(1))+x(i+2*m+ioff(1)))/(2.0D0*dxdy)
        y(i+ioff(2)) = (x(i+1+ioff(1))-x(i-1+ioff(1)))/(2.0*dxdy)+ &
            (-3.0D0*x(i+ioff(2))+4.0D0*x(i+m+ioff(2))-x(i+2*m+ioff(2)))/(2.0D0*dxdy)
    END DO
    ! Top boundary:
    DO i=n-m+2,n-1
        y(i+ioff(1)) = (x(i+1+ioff(2))-x(i-1+ioff(2)))/(2.0D0*dxdy)+ &
            (-3.0D0*x(i+ioff(1))+4.0D0*x(i-m+ioff(1))-x(i-2*m+ioff(1)))/(2.0D0*dxdy)
        y(i+ioff(2)) = (x(i+1+ioff(1))-x(i-1+ioff(1)))/(2.0D0*dxdy)+ &
            ( 3.0D0*x(i+ioff(2))-4.0D0*x(i-m+ioff(2))+x(i-2*m+ioff(2)))/(2.0D0*dxdy)
    END DO
    ! Left boundary:
    DO i=m+1,n-2*m+1,m
        y(i+ioff(1)) = -(x(i+m+ioff(1))-x(i-m+ioff(1)))/(2.0D0*dxdy)+ &
            (-3.0D0*x(i+ioff(2))+4.0D0*x(i+1+ioff(2))-x(i+2+ioff(2)))/(2.0D0*dxdy)
        y(i+ioff(2)) = (x(i+m+ioff(2))-x(i-m+ioff(2)))/(2.0D0*dxdy)+ &
            (-3.0D0*x(i+ioff(1))+4.0D0*x(i+1+ioff(1))-x(i+2+ioff(1)))/(2.0D0*dxdy)
    END DO
    ! Right boundary:
    DO i=2*m,n-m,m
        y(i+ioff(1)) = -(x(i+m+ioff(1))-x(i-m+ioff(1)))/(2.0D0*dxdy)+ &
            ( 3.0D0*x(i+ioff(2))-4.0D0*x(i-1+ioff(2))+x(i-2+ioff(2)))/(2.0D0*dxdy)
        y(i+ioff(2)) =  (x(i+m+ioff(2))-x(i-m+ioff(2)))/(2.0D0*dxdy)+ &
            ( 3.0D0*x(i+ioff(1))-4.0D0*x(i-1+ioff(1))+x(i-2+ioff(1)))/(2.0D0*dxdy)
    END DO

    ! Interior points:
    DO j=2,m-1
        DO i=2,m-1
            l = i+m*(j-1)
            y(l+ioff(1)) = &
            (x(l+1+ioff(1))-2.0D0*x(l+ioff(1))+x(l-1+ioff(1)))/(dxdy*dxdy)+ &
            (x(l+m+ioff(1))-2.0D0*x(l+ioff(1))+x(l-m+ioff(1)))/(dxdy*dxdy)
            y(l+ioff(2)) = &
            (x(l+1+ioff(2))-2.0D0*x(l+ioff(2))+x(l-1+ioff(2)))/(dxdy*dxdy)+ &
            (x(l+m+ioff(2))-2.0D0*x(l+ioff(2))+x(l-m+ioff(2)))/(dxdy*dxdy)
        END DO
    END DO

    ! Debugging:
    DO i=1,n
        IF (ISNAN(y(i))) THEN
            WRITE(*,1) i,y(i),x(i)
1           FORMAT('Found NaN: ',I5,2E12.4)
        END IF
    END DO
END SUBROUTINE matvec_Poisson2dMixed

SUBROUTINE nullProjection2dMixed(n,x,y)

    USE kinds_m, ONLY : i_kind, r_kind

    INTEGER(i_kind), INTENT(IN) :: n
    REAL(r_kind), INTENT(IN) :: x(n)
    REAL(r_kind), INTENT(OUT) :: y(n)

    ! Local variables:
    INTEGER(i_kind) :: i,j,k,m,nn,ioff(2)
    REAL(r_kind) :: s(n),z(n),d(n),xx,yy

    m = DSQRT(DBLE(n/2))
    ioff(1) = 0
    ioff(2) = m

    ! Project out the constant:
    d(1:n/2) = x(1:n/2)-SUM(x(1:n/2))/DBLE(n/2)
    d(n/2+1:n) = x(n/2+1:n)-SUM(x(n/2+1:n))/DBLE(n/2)
    y = d
    PRINT*,'MaxMin values: ',maxval(d),minval(d)
    PRINT*,'MaxMin values: ',maxval(x),minval(x)

    ! Project out the linear function:
    k = 0
    s = 0.0D0
    z = 0.0D0
    DO j=1,m
        DO i=1,m
            k = k+1
            s(k+ioff(1)) =  DBLE(i-1)/DBLE(m-1)-0.5D0
            z(k+ioff(1)) = -(DBLE(j-1)/DBLE(m-1)-0.5D0)
            s(k+ioff(2)) =  DBLE(j-1)/DBLE(m-1)-0.5D0
            z(k+ioff(2)) =  DBLE(i-1)/DBLE(m-1)-0.5D0
        END DO
    END DO
    d(1:n/2) = y(1:n/2) - DOT_PRODUCT(y(1:n/2),s(1:n/2))/DOT_PRODUCT(s(1:n/2),s(1:n/2))*s(1:n/2)
    y(1:n/2) = d(1:n/2) - DOT_PRODUCT(d(1:n/2),z(1:n/2))/DOT_PRODUCT(z(1:n/2),z(1:n/2))*z(1:n/2)
    d(n/2+1:n) = y(n/2+1:n) - DOT_PRODUCT(y(n/2+1:n),s(n/2+1:n))/DOT_PRODUCT(s(n/2+1:n),s(n/2+1:n))*s(n/2+1:n)
    y(n/2+1:n) = d(n/2+1:n) - DOT_PRODUCT(d(n/2+1:n),z(n/2+1:n))/DOT_PRODUCT(z(n/2+1:n),z(n/2+1:n))*z(n/2+1:n)

END SUBROUTINE nullProjection2dMixed

SUBROUTINE psolve(n,x)

    USE kinds_m, ONLY : i_kind, r_kind

    INTEGER(i_kind),INTENT(IN) :: n
    REAL(r_kind), INTENT(INOUT) :: x(n)

    RETURN
END SUBROUTINE psolve

REAL(r_kind) FUNCTION dotprd(n,a,b)

    USE kinds_m, ONLY : i_kind, r_kind
    
    INTEGER(i_kind), INTENT(IN) :: n
    REAL(r_kind), INTENT(IN) :: a(n),b(n)

    dotprd = 0.0D0
    dotprd = DOT_PRODUCT(a,b)
END FUNCTION dotprd