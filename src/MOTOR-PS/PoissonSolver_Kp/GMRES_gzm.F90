!----------------------------------------------------------------------
! Openpipeflow.org.  If used in your work, please cite
! Willis, A. (2017) SoftwareX 6, 124-127.
! https://doi.org/10.1016/j.softx.2017.05.003 (open access)
!                                      Thanks in advance! Ashley 2019.
!----------------------------------------------------------------------
! solve A x = b for x ;  
! minimise |Ax-b| subject to constraint |x| < delta .
! requires lapack routines dgelsy, dgesvd.
!----------------------------------------------------------------------
! m	  gmres dimension
! n 	  dimension of x
! x	  on input:  guess for x, can be 0
!         on exit:  solution x, subject to constraint if del>0
! b	  input b
! matvec  performs y := A x, call matvec(N,x, y)
! psolve  preconditioner, solve M x_out = x_in, call psolve(N,x)
! dotprd  dot product, d = dotprd(n,a,b)
! h       Hessian matrix,  size (m+1)*m
! v       Krylov subspace, size n*(m+1)
! res	  on input: |Ax-b|/|b|<res; 
!         on exit:  residual reached
! del     on input: if(del>0) then the x returned is the hookstep
!         on exit:  norm of next b predicted by hook
! its	  on input: max num its; 
!         on exit:  number of its taken
! info	  on input: if(info==1) print* residuals
!                   if(info==2) recalc hookstep with new del
! 	  on exit:  0 sucessful, 1 method breakdown, 2 max its
!							A.P.Willis 2008
!----------------------------------------------------------------------

!----------------------------------------------------------------------
! Yuanfu Xie modified this package: 2024-03-22
! See the in line comment!
! In some cases, this GMRES may breakdown. Currently, Yuanfu Xie adopts
! a restart strategy with different initial guess. It seems resolving 
! the breakdown problem but he adds a safeguard in calling this package
! to remind the users the breakdown over 5 different restarts with 5
! different initial guesses.
!
! There are still potential problems causing NaN even though Yuanfu Xie
! fixed one of them here.
!----------------------------------------------------------------------

! Yuanfu Xie modified this package on 2024-12-09 for allowing matvec to
! access gzm model.

 subroutine gmresm_gzm(m,n,x,b,matvec,psolve,dotprd,h,v,res,del,its,info, gzm, scaling)
   USE gzm_m, ONLY: gzm_t
   implicit none
   integer,          intent(in)    :: m
   integer,          intent(in)    :: n
   double precision, intent(inout) :: x(n)
   double precision, intent(in)    :: b(n)
   external                        :: matvec,psolve
   double precision, external      :: dotprd
   double precision, intent(inout) :: h(m+1,m)
   double precision, intent(inout) :: v(n,m+1)
   double precision, intent(inout) :: res
   double precision, intent(inout) :: del
   integer,          intent(inout) :: its
   integer,          intent(inout) :: info
   double precision :: tol,res_,stgn, w(n), z(n)
   double precision :: h_(m+1,m), y(m+1), p(m+1), work(4*m+1)
   integer :: imx, piv(m), rank, i
   double precision, save :: beta
   integer, save :: j
   logical :: done   
   
   ! Yuanfu Xie added an option to pass in gzm model:
   TYPE(gzm_t), OPTIONAL, INTENT(IN) :: gzm
   DOUBLE PRECISION, OPTIONAL, INTENT(IN) :: scaling(2)
   
   ! Yuanfu Xie added a counter of the times of restarts:
   INTEGER :: numRestarts

   numRestarts = -1

   PRINT*,'GMRES_gzm - Debugging 0...',info,m,n
   if(info==2) then
      call hookstep(j,h,m,beta,del, y)
      z = matmul(v(:,1:j),y(1:j))
      call psolve(n, z)
      x = z
      info = 0
      return
   end if	 

   tol = res
   imx = its
   its = 0
   v   = 0d0

   PRINT*,'GMRES_gzm - Debugging 0.1',dot_product(x,x),dot_product(b,b)

 1 continue
   numRestarts = numRestarts+1   ! Count the number of restarts
   res_ = 1d99
   stgn = 1d0 - 1d-14
 
   beta = dsqrt(dotprd(n,x,x)) 
   if(beta==0d0)  w = 0d0

   ! Yuanfu Xie modified to call matvec with gzm option:
   if(beta/=0d0)  THEN
      IF (PRESENT(gzm)) THEN
         call matvec(n,x, w, gzm, scaling)
      ELSE
         call matvec(n,x, w)
      END IF
   END IF
   w = b - w
   beta = dsqrt(dotprd(n,w,w)) 
   v(:,1) = w / beta

   ! Debugging:
   WRITE(*,21) w(2:5),x(2:5)
21  FORMAT('GMRES_gzm - Debugging 1.0 VTV: ',8E14.6)

   IF (numRestarts .EQ. -1) THEN
      WRITE(10) w(1:n/2)
      WRITE(10) w(n/2+1:n)
   END IF
     
   h = 0d0
   do j = 1, m
      its = its + 1
      z = v(:,j)      

      call psolve(n, z)
      ! Yuanfu Xie modified to call matvec with gzm option:
      IF (PRESENT(gzm)) THEN
         call matvec(n,z, w, gzm, scaling)

         ! Debugging:
         WRITE(*,11) dot_product(z,z),dot_product(w,w),j,numRestarts,w(2:5)
      11 FORMAT('Call matvec with gzm: zz: ',E12.4,' ww: ',E12.4,' j numRestarts: ',2I3,' w162: ',4E14.6)

         IF (numRestarts .EQ. 0 .AND. j .EQ. 1) THEN
            WRITE(10) z(1:n/2)
            WRITE(10) z(n/2+1:n)
         END IF
      ELSE
         call matvec(n,z, w)
      END IF

      ! Debugging:
      WRITE(*,22) dot_product(w,w),dot_product(v(:,1),v(:,1)),w(2:5)
22    FORMAT('GMRES_gzm - Debugging 1.2 WTW: ',E14.6,' VTV: ',E14.6,' Ax b: ',4E12.4)

      do i = 1, j
         h(i,j) = dotprd(n,w,v(1,i))
         w = w - h(i,j)*v(:,i)
      end do
      h(j+1,j) = dsqrt(dotprd(n,w,w))

      ! Debugging:
      WRITE(*,23) h(j+1,j),dotprd(n,w,w),w(2:5)
23    FORMAT('GMRES_gzm - Debugging 2 h(j+1,j): ',E14.6,' WTW: ',E14.6,' w2-5: ',4E14.6)

      ! Yuanfu Xie modified this part of v(:,j+1) to avoid dividing zero
      ! 2024-03-22
      IF (h(j+1,j) .NE. 0.0D0) THEN
      v(:,j+1) = w / h(j+1,j)
      ELSE
      v(:,j+1) = w / 0.1D-10
      END IF
          
      p(1) = beta
      p(2:j+1) = 0d0
      h_(1:j+1,1:j) = h(1:j+1,1:j)
      call dgelsy(j+1,j,1,h_,m+1,p,m+1,piv,m,rank,work,4*m+1,i)
      if(i/=0) stop 'gmresm: dgelsy'
      y = p

      p(1:j+1) = - matmul(h(1:j+1,1:j),y(1:j))
      p(1) = p(1) + beta
   
      res = dsqrt(dot_product(p(1:j+1),p(1:j+1)))

      ! Debugging:
      WRITE(*,24) j,p(1),p(j+1),beta
24    FORMAT('GMRES_gzm - Debugging 3 j: ',I4,' P1 P(j+1): '2E14.6,' Beta: ',E14.6)

      if(info==1) print*, 'gmresm: it=', its,' res=', real(res)
      
      done = (res<=tol .or. its==imx .or. res>res_)
      PRINT*,'DONE: ',j,numRestarts,done, res,tol,its,imx, res, res_
      if(done .or. j==m) then
         if(del>0d0)  call hookstep(j,h,m,beta,del, y)
         print*,'After hookstep...'
         z = matmul(v(:,1:j),y(1:j))
         call psolve(n, z)
         x = x + z
         if(its==imx) info = 2
         if(res>res_) info = 1
         if(res<=tol) info = 0
         if(done)     return
         if(del>0d0)  print*, 'gmres: warning! restart affects hookstep'
         goto 1       ! (j==m) restart
      end if
      res_ = res*stgn
   PRINT*,'GMRES_gzm - Debugging 4...',DOT_PRODUCT(x,x)

   end do   
 
 end subroutine gmresm_gzm
 
 
!-----------------------------------------------------------------
! replace y with a vector that generates a hookstep
! c.f. Viswanath (2008) arXiv:0809.1498
!-----------------------------------------------------------------
 subroutine hookstep(j,h,m,beta,del, y)
   implicit none
   integer,          intent(in)    :: j, m
   double precision, intent(in)    :: h(m+1,j), beta
   double precision, intent(inout) :: del
   double precision, intent(out)   :: y(j)
   double precision :: a(j+1,j), s(j), u(j+1,j+1), vt(j,j), work(5*(j+1))
   double precision :: p(j+1), q(j), mu, qn
   integer :: info
   
   a = h(1:j+1,1:j)
   
   call dgesvd('A','A',j+1,j,a,j+1,s,u,j+1,vt,j,work,5*(j+1),info)
   if(info/=0) stop 'hookstep: dgesvd'
   
   p(1:j) = beta * u(1,1:j)   

   mu = max(s(j)*s(j)*1d-6,1d-99)
   qn = 1d99
   do while(qn>del)
      mu = mu * 1.1d0
      q = p(1:j)*s/(mu+s*s)
      qn = dsqrt(dot_product(q,q))
   end do

   y = matmul(q,vt)

   p = - matmul(h(1:j+1,1:j),y(1:j))
   p(1) = p(1) + beta
   del = dsqrt(dot_product(p,p))
 
 end subroutine hookstep