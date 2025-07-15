! Description:
!> @file
!!   Defines various mathematical operations on arrays
!
!> @brief
!!   Defines various mathematical operations on arrays
!!
!! @details
!!   This module contains a number of functions which perform common
!!   mathematical operations on arrays of various sizes. TL/AD versions
!!   are also included where relevant. This allows for optimised code
!!   to be implemented in one place (for example, compiler-specific
!!   optimised library calls).
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
!    Copyright 2016, EUMETSAT, All Rights Reserved.
!
MODULE rttov_math_mod

#ifdef RTTOV_INTEL_MKL
  USE mkl95_blas
  USE mkl95_precision 
#endif

  USE parkind1, ONLY : jprb, jpim, jplm

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: EXPONENTIAL, EXPONENTIAL_TL, EXPONENTIAL_AD, &
            RECIPROCAL, RECIPROCAL_TL, RECIPROCAL_AD, RECIPROCAL_K, &
            DIVIDE, PLANCK, PLANCK_TL, PLANCK_AD, &
            INV_PLANCK, INV_PLANCK_TL, INV_PLANCK_AD, &
            SINTOSEC, SINTOSEC_TL, SINTOSEC_AD, SINTOSEC_K, &
            INVSQRT, INVSQRT_TL, INVSQRT_AD, INVSQRT_K, &
            SQRT_TL, SQRT_AD, SQRT_K

  INTERFACE EXPONENTIAL
    MODULE PROCEDURE EXPONENTIAL_1D
    MODULE PROCEDURE EXPONENTIAL_2D
  END INTERFACE EXPONENTIAL

  INTERFACE EXPONENTIAL_TL
    MODULE PROCEDURE EXPONENTIAL_1D_TL
    MODULE PROCEDURE EXPONENTIAL_2D_TL
  END INTERFACE EXPONENTIAL_TL

  INTERFACE EXPONENTIAL_AD
    MODULE PROCEDURE EXPONENTIAL_1D_AD
    MODULE PROCEDURE EXPONENTIAL_2D_AD
  END INTERFACE EXPONENTIAL_AD

  INTERFACE RECIPROCAL
     MODULE PROCEDURE RECIPROCAL_1D
     MODULE PROCEDURE RECIPROCAL_2D
  END INTERFACE RECIPROCAL

  INTERFACE RECIPROCAL_TL
     MODULE PROCEDURE RECIPROCAL_1D_TL
     MODULE PROCEDURE RECIPROCAL_2D_TL
   END INTERFACE RECIPROCAL_TL

  INTERFACE RECIPROCAL_AD
     MODULE PROCEDURE RECIPROCAL_1D_AD
     MODULE PROCEDURE RECIPROCAL_2D_AD
   END INTERFACE RECIPROCAL_AD

  INTERFACE RECIPROCAL_K
    MODULE PROCEDURE RECIPROCAL_2D_K
  END INTERFACE

  INTERFACE DIVIDE
     MODULE PROCEDURE DIVIDE_1D
     MODULE PROCEDURE DIVIDE_2D
  END INTERFACE

  INTERFACE PLANCK
     MODULE PROCEDURE PLANCK_1D
     MODULE PROCEDURE PLANCK_1D_NU
     MODULE PROCEDURE PLANCK_SCALAR
     MODULE PROCEDURE PLANCK_SCALAR_NU
  END INTERFACE

  INTERFACE PLANCK_TL
    MODULE PROCEDURE PLANCK_1D_TL
    MODULE PROCEDURE PLANCK_SCALAR_TL
    MODULE PROCEDURE PLANCK_SCALAR_NU_TL
  END INTERFACE

  INTERFACE PLANCK_AD
    MODULE PROCEDURE PLANCK_1D_AD
    MODULE PROCEDURE PLANCK_SCALAR_AD
  END INTERFACE

  INTERFACE INV_PLANCK
    MODULE PROCEDURE INV_PLANCK_SCALAR
    MODULE PROCEDURE INV_PLANCK_SCALAR_NU
  END INTERFACE

  INTERFACE INV_PLANCK_TL
    MODULE PROCEDURE INV_PLANCK_SCALAR_TL
  END INTERFACE

  INTERFACE INV_PLANCK_AD
    MODULE PROCEDURE INV_PLANCK_SCALAR_AD
  END INTERFACE

  INTERFACE SINTOSEC
    MODULE PROCEDURE SINTOSEC_2D
  END INTERFACE

  INTERFACE SINTOSEC_TL
    MODULE PROCEDURE SINTOSEC_2D_TL
  END INTERFACE

  INTERFACE SINTOSEC_AD
    MODULE PROCEDURE SINTOSEC_2D_AD
  END INTERFACE

  INTERFACE SINTOSEC_K
    MODULE PROCEDURE SINTOSEC_2D_K
  END INTERFACE

  INTERFACE INVSQRT
    MODULE PROCEDURE INVSQRT_2D
  END INTERFACE

  INTERFACE INVSQRT_TL
    MODULE PROCEDURE INVSQRT_1D_TL
    MODULE PROCEDURE INVSQRT_2D_TL
  END INTERFACE

  INTERFACE INVSQRT_AD
    MODULE PROCEDURE INVSQRT_1D_AD
    MODULE PROCEDURE INVSQRT_2D_AD
  END INTERFACE

  INTERFACE SQRT_TL
    MODULE PROCEDURE SQRT_1D_TL
    MODULE PROCEDURE SQRT_2D_TL
  END INTERFACE

   INTERFACE SQRT_AD
     MODULE PROCEDURE SQRT_1D_AD
     MODULE PROCEDURE SQRT_2D_AD
   END INTERFACE

CONTAINS

#ifdef RTTOV_INTEL_MKL
  INCLUDE 'mkl_vml.f90'
#endif

  SUBROUTINE exponential_1d(in,out)
    REAL(jprb), INTENT(in) :: in(:)
    REAL(jprb), INTENT(out) :: out(:)

#ifdef RTTOV_INTEL_MKL
    CALL vdexp(SIZE(in), in, out)
#elif RTTOV_IBM_MASS
    CALL vexp(out, in, SIZE(in))
#else
    out = EXP(in)
#endif
  END SUBROUTINE exponential_1d

  SUBROUTINE exponential_2d(in,out)
    REAL(jprb), INTENT(in) :: in(:,:)
    REAL(jprb), INTENT(out) :: out(:,:)

#ifdef RTTOV_INTEL_MKL
    CALL vdexp(SIZE(in), in, out)
#elif RTTOV_IBM_MASS
    CALL vexp(out, in, SIZE(in))
#else
    out = EXP(in)
#endif
  END SUBROUTINE exponential_2d

  SUBROUTINE exponential_1d_tl(x, x_tl, y_tl)
    REAL(jprb), INTENT(in) :: x(:), x_tl(:)
    REAL(jprb), INTENT(out) :: y_tl(:)

    y_tl(:) = x_tl(:) * x(:)
  END SUBROUTINE exponential_1d_tl

  SUBROUTINE exponential_2d_tl(x, x_tl, y_tl)
    REAL(jprb), INTENT(in) :: x(:,:), x_tl(:,:)
    REAL(jprb), INTENT(out) :: y_tl(:,:)

    y_tl(:,:) = x_tl(:,:) * x(:,:)
  END SUBROUTINE exponential_2d_tl

  SUBROUTINE exponential_2d_ad(x, x_ad, y_ad, acc)
    REAL(jprb), INTENT(in) :: x(:,:)
    REAL(jprb), INTENT(inout) :: x_ad(:,:), y_ad(:,:)
    LOGICAL(jplm), INTENT(IN) :: acc
    
    IF (acc) THEN
      x_ad = x_ad + y_ad * x
    ELSE
      x_ad = y_ad * x
    ENDIF
  END SUBROUTINE exponential_2d_ad

  SUBROUTINE exponential_1d_ad(x, x_ad, y_ad, acc)
    REAL(jprb), INTENT(in) :: x(:)
    REAL(jprb), INTENT(inout) :: x_ad(:), y_ad(:)
    LOGICAL(jplm), INTENT(IN) :: acc
    
    IF (acc) THEN
      x_ad = x_ad + y_ad * x
    ELSE
      x_ad = y_ad * x
    ENDIF
  END SUBROUTINE exponential_1d_ad

  SUBROUTINE reciprocal_2d(in,out,scale)
    REAL(jprb), INTENT(in) :: in(:,:)
    REAL(jprb), INTENT(inout) :: out(:,:)
    REAL(jprb), INTENT(in), OPTIONAL :: scale

#ifdef RTTOV_INTEL_MKL
    CALL vdinv(SIZE(in), in, out)
#elif RTTOV_IBM_MASS
    CALL vrec(out, in, SIZE(in))
#else
    out = 1._jprb / in
#endif

    IF (PRESENT(scale)) out = out * scale
    
  END SUBROUTINE reciprocal_2d

  SUBROUTINE reciprocal_1d(in,out,scale)
    REAL(jprb), INTENT(in) :: in(:)
    REAL(jprb), INTENT(inout) :: out(:)
    REAL(jprb), INTENT(in), OPTIONAL :: scale

#ifdef RTTOV_INTEL_MKL
    CALL vdinv(SIZE(in), in, out)
#elif RTTOV_IBM_MASS
    CALL vrec(out, in, SIZE(in))
#else
    out = 1._jprb / in
#endif

    IF (PRESENT(scale)) out = out * scale
  END SUBROUTINE reciprocal_1d

  SUBROUTINE reciprocal_1d_tl(x,x_tl,y_tl)
    REAL(jprb), INTENT(in) :: x(:),x_tl(:)
    REAL(jprb), INTENT(out) :: y_tl(:)

    y_tl = - x**2 * x_tl 
  END SUBROUTINE reciprocal_1d_tl

  SUBROUTINE reciprocal_2d_tl(x,x_tl,y_tl)
    REAL(jprb), INTENT(in) :: x(:,:),x_tl(:,:)
    REAL(jprb), INTENT(out) :: y_tl(:,:)

    y_tl = - x**2 * x_tl 
  END SUBROUTINE reciprocal_2d_tl

  SUBROUTINE reciprocal_1d_ad(x,x_ad,y_ad, acc)
    REAL(jprb), INTENT(in) :: x(:), y_ad(:)
    REAL(jprb), INTENT(inout) :: x_ad(:)
    LOGICAL(jplm), INTENT(IN) :: acc
    
    IF (acc) THEN
      x_ad = x_ad - x**2 * y_ad
    ELSE
      x_ad = - x**2 * y_ad
    ENDIF
  END SUBROUTINE reciprocal_1d_ad

  SUBROUTINE reciprocal_2d_ad(x,x_ad,y_ad, acc)
    REAL(jprb), INTENT(in) :: x(:,:), y_ad(:,:)
    REAL(jprb), INTENT(inout) :: x_ad(:,:)
    LOGICAL(jplm), INTENT(IN) :: acc

    IF (acc) THEN
      x_ad = x_ad - x**2 * y_ad
    ELSE
      x_ad = - x**2 * y_ad
    ENDIF
  END SUBROUTINE reciprocal_2d_ad

  SUBROUTINE reciprocal_2d_k(x,x_k,y_k, acc, map, prof_stat, limits)
    REAL(jprb), INTENT(in) :: x(:,:), y_k(:,:)
    REAL(jprb), INTENT(inout) :: x_k(:,:)
    LOGICAL(jplm), INTENT(IN) :: acc
    INTEGER(jpim), INTENT(IN) :: map(:,:)
    INTEGER(jpim), INTENT(IN) :: prof_stat
    INTEGER(jpim), INTENT(IN), OPTIONAL :: limits(:,:)

    INTEGER(jpim) :: loc, lb, ub, prof, i, ul, ll
    REAL(jprb) :: z(SIZE(x(:,1)))
    
    loc = 1

    IF (prof_stat >= 0) THEN
      DO prof = 1, MAXVAL(map(:,1))
        lb = loc
        ub = lb
        DO i = lb + 1, SIZE(map(:,1))
          IF (map(i,1) .NE. prof) THEN
            ub = i - 1
            loc = i
            EXIT
          ENDIF
          ub = SIZE(map(:,1))
        ENDDO

        z(:) = -x(:,prof) * x(:,prof)
        
        IF (PRESENT(limits)) THEN
          IF (acc) THEN
            DO i=lb,ub
              ul = limits(2,map(i,2))
              ll = limits(1,map(i,2))
              x_k(ll:ul,i) = x_k(ll:ul,i) + z(ll:ul) * y_k(ll:ul,i)
            ENDDO
          ELSE
            DO i=lb,ub
              ul = limits(2,map(i,2))
              ll = limits(1,map(i,2))
              x_k(ll:ul,i) = z(ll:ul) * y_k(ll:ul,i)
            ENDDO
          ENDIF
        ELSE
          IF (acc) THEN
            DO i=lb,ub
              x_k(:,i) = x_k(:,i) + z(:) * y_k(:,i)
            ENDDO
          ELSE
            DO i=lb,ub
              x_k(:,i) = z(:) * y_k(:,i)
            ENDDO
          ENDIF
        ENDIF
      ENDDO
      
    ELSE IF (prof_stat .EQ. -1) THEN
      DO i = 1, SIZE(y_k(1,:))
        prof = map(i,1)
        IF (acc) THEN
          x_k(:,i) = x_k(:,i) - x(:, prof)**2 * y_k(:,i)
        ELSE
          x_k(:,i) = -x(:, prof)**2 * y_k(:,i)
        ENDIF
      ENDDO
    ENDIF
  END SUBROUTINE reciprocal_2d_k

  SUBROUTINE DIVIDE_1D(a, b, out, acc)
    REAL(jprb), INTENT(in) :: a(:), b(:)
    REAL(jprb), INTENT(inout) :: out(:)
    LOGICAL(jplm), INTENT(in) :: acc

#ifdef RTTOV_INTEL_MKL
    REAL(jprb) :: z(SIZE(a))
    IF (acc) THEN
      CALL vddiv(SIZE(a), a, b, z)
      CALL vdadd(SIZE(a), out, z, out)
    ELSE
      CALL vddiv(SIZE(a), a, b, out)
    ENDIF
#elif RTTOV_IBM_MASS
    REAL(jprb) :: z(SIZE(a))
    IF (acc) THEN
      CALL vdiv(z, a, b, SIZE(a))
      out = out + z
    ELSE
      CALL vdiv(out, a, b, SIZE(a))
    ENDIF
#else
    IF (acc) THEN
      out = out + (a / b)
    ELSE
      out = a / b
    ENDIF
#endif
    
  END SUBROUTINE DIVIDE_1D

  SUBROUTINE DIVIDE_2D(a, b, out, acc)
    REAL(jprb), INTENT(in) :: a(:,:), b(:,:)
    REAL(jprb), INTENT(inout) :: out(:,:)
    LOGICAL(jplm), INTENT(in) :: acc

#ifdef RTTOV_INTEL_MKL
    REAL(jprb) :: z(SIZE(a(:,1)), SIZE(a(1,:)))
    IF (acc) THEN
      CALL vddiv(SIZE(a), a, b, z)
      CALL vdadd(SIZE(a), out, z, out)
    ELSE
      CALL vddiv(SIZE(a), a, b, out)
    ENDIF
#elif RTTOV_IBM_MASS
    REAL(jprb) :: z(SIZE(a(:,1)), SIZE(a(1,:)))
    IF (acc) THEN
      CALL vdiv(z, a, b, SIZE(a))
      out = out + z
    ELSE
      CALL vdiv(out, a, b, SIZE(a))
    ENDIF
#else
    IF (acc) THEN
      out = out + (a / b)
    ELSE
      out = a / b
    ENDIF
#endif
    
  END SUBROUTINE DIVIDE_2D

  SUBROUTINE INVSQRT_2d(in, out)
    REAL(jprb), INTENT(in) :: in(:,:)
    REAL(jprb), INTENT(out) :: out(:,:)

#ifdef RTTOV_INTEL_MKL
    CALL vdinvsqrt(SIZE(in), in, out)
#elif RTTOV_IBM_MASS
    CALL vrsqrt(out,in,SIZE(in)) ! check these
#else
    out(:,:) = 1._jprb / SQRT(in(:,:))
#endif

  END SUBROUTINE INVSQRT_2d

  SUBROUTINE INVSQRT_2D_TL(x, x_tl, y_tl)
    REAL(jprb), INTENT(in) :: x(:,:), x_tl(:,:)
    REAL(jprb), INTENT(out) :: y_tl(:,:)

    y_tl = -0.5_jprb * x**3_jpim * x_tl
  END SUBROUTINE INVSQRT_2D_TL

  SUBROUTINE INVSQRT_1D_TL(x, x_tl, y_tl)
    REAL(jprb), INTENT(in) :: x(:), x_tl(:)
    REAL(jprb), INTENT(out) :: y_tl(:)

    y_tl = -0.5_jprb * x**3_jpim * x_tl
  END SUBROUTINE INVSQRT_1D_TL

  SUBROUTINE INVSQRT_1D_AD(x, x_ad, y_ad, acc)
    REAL(jprb), INTENT(in) :: x(:), y_ad(:)
    REAL(jprb), INTENT(inout) :: x_ad(:)
    LOGICAL(jplm), INTENT(in) :: acc

    IF (acc) THEN
      x_ad = x_ad - 0.5_jprb * x**3_jpim * y_ad
    ELSE
      x_ad = -0.5_jprb * x**3_jpim * y_ad
    ENDIF
  END SUBROUTINE INVSQRT_1D_AD

  SUBROUTINE INVSQRT_2D_AD(x, x_ad, y_ad, acc)
    REAL(jprb), INTENT(in) :: x(:,:), y_ad(:,:)
    REAL(jprb), INTENT(inout) :: x_ad(:,:)
    LOGICAL(jplm), INTENT(in) :: acc

    IF (acc) THEN
      x_ad = x_ad - 0.5_jprb * x**3_jpim * y_ad
    ELSE
      x_ad = -0.5_jprb * x**3_jpim * y_ad
    ENDIF
  END SUBROUTINE INVSQRT_2D_AD

  SUBROUTINE invsqrt_K(x,x_k,y_k, acc, map, prof_stat, limits)
    REAL(jprb), INTENT(in) :: x(:,:), y_k(:,:)
    REAL(jprb), INTENT(inout) :: x_k(:,:)
    LOGICAL(jplm), INTENT(IN) :: acc
    INTEGER(jpim), INTENT(IN) :: map(:,:)
    INTEGER(jpim), INTENT(IN) :: prof_stat
    INTEGER(jpim), INTENT(IN), OPTIONAL :: limits(:,:)

    INTEGER(jpim) :: loc, lb, ub, prof, i, ll, ul
    REAL(jprb) :: z(SIZE(x(:,1)))
    
    loc = 1

    IF (prof_stat >= 0) THEN
      DO prof = 1, MAXVAL(map(:,1))
        lb = loc
        ub = lb
        DO i = lb + 1, SIZE(map(:,1))
          IF (map(i,1) .NE. prof) THEN
            ub = i - 1
            loc = i
            EXIT
          ENDIF
          ub = SIZE(map(:,1))
        ENDDO

        z(:) = - 0.5_jprb * x(:,prof)**3_jpim
       
        IF (PRESENT(limits)) THEN
          IF (acc) THEN
            DO i=lb,ub
              ul = limits(2,map(i,2))
              ll = limits(1,map(i,2))
              x_k(ll:ul,i) = x_k(ll:ul,i) + z(ll:ul) * y_k(ll:ul,i)
            ENDDO
          ELSE
            DO i=lb,ub
              ul = limits(2,map(i,2))
              ll = limits(1,map(i,2))
              x_k(ll:ul,i) = z(ll:ul) * y_k(ll:ul,i)
            ENDDO
          ENDIF
        ELSE
          IF (acc) THEN
            DO i=lb,ub
              x_k(:,i) = x_k(:,i) + z(:) * y_k(:,i)
            ENDDO
          ELSE
            DO i=lb,ub
              x_k(:,i) = z(:) * y_k(:,i)
            ENDDO
          ENDIF
        ENDIF
      ENDDO
      
    ELSE IF (prof_stat .EQ. -1) THEN
      DO i = 1, SIZE(y_k(1,:))
        prof = map(i,1)
        IF (acc) THEN
          x_k(:,i) = x_k(:,i) - 0.5_jprb * x(:,prof)**3_jpim * y_k(:,i)
        ELSE
          x_k(:,i) = - 0.5_jprb * x(:,prof)**3_jpim * y_k(:,i)
        ENDIF
      ENDDO
    ENDIF
  END SUBROUTINE invsqrt_K

  SUBROUTINE SQRT_1D_TL(x_sqrt, x_tl, y_tl, x_rsqrt)
    REAL(jprb), INTENT(in) :: x_sqrt(:), x_tl(:)
    REAL(jprb), INTENT(out) :: y_tl(:)
    REAL(jprb), INTENT(in), OPTIONAL :: x_rsqrt(:)

    IF(PRESENT(x_rsqrt)) THEN
      y_tl = 0.5_jprb * x_tl * x_rsqrt
    ELSE
      CALL divide_1d(0.5_jprb * x_tl, x_sqrt, y_tl, acc = .FALSE._jplm)
    ENDIF
  END SUBROUTINE SQRT_1D_TL
  
  SUBROUTINE SQRT_2D_TL(x_sqrt, x_tl, y_tl, x_rsqrt)
    REAL(jprb), INTENT(in) :: x_sqrt(:,:), x_tl(:,:)
    REAL(jprb), INTENT(out) :: y_tl(:,:)
    REAL(jprb), INTENT(in), OPTIONAL :: x_rsqrt(:,:)

    IF(PRESENT(x_rsqrt)) THEN
      y_tl = 0.5_jprb * x_tl * x_rsqrt
    ELSE
      CALL divide_2d(0.5_jprb * x_tl, x_sqrt, y_tl, acc = .FALSE._jplm)
    ENDIF
  END SUBROUTINE SQRT_2D_TL

  SUBROUTINE SQRT_1D_AD(x_sqrt, x_ad, y_ad, acc, x_rsqrt)
    REAL(jprb), INTENT(in) :: x_sqrt(:), y_ad(:)
    REAL(jprb), INTENT(inout) :: x_ad(:)
    LOGICAL(jplm), INTENT(in) :: acc
    REAL(jprb), INTENT(in), OPTIONAL :: x_rsqrt(:)

    IF(PRESENT(x_rsqrt)) THEN
      IF(acc) THEN
        x_ad = x_ad + 0.5_jprb * y_ad * x_rsqrt
      ELSE
        x_ad = 0.5_jprb * y_ad * x_rsqrt
      ENDIF
    ELSE
      CALL divide_1d(0.5_jprb * y_ad, x_sqrt, x_ad, acc = acc)
    ENDIF
  END SUBROUTINE SQRT_1D_AD

  SUBROUTINE SQRT_2D_AD(x_sqrt, x_ad, y_ad, acc, x_rsqrt)
    REAL(jprb), INTENT(in) :: x_sqrt(:,:), y_ad(:,:)
    REAL(jprb), INTENT(inout) :: x_ad(:,:)
    LOGICAL(jplm), INTENT(in) :: acc
    REAL(jprb), INTENT(in), OPTIONAL :: x_rsqrt(:,:)

    IF(PRESENT(x_rsqrt)) THEN
      IF(acc) THEN
        x_ad = x_ad + 0.5_jprb * y_ad * x_rsqrt
      ELSE
        x_ad = 0.5_jprb * y_ad * x_rsqrt
      ENDIF
    ELSE
      CALL divide_2d(0.5_jprb * y_ad, x_sqrt, x_ad, acc = acc)
    ENDIF
  END SUBROUTINE SQRT_2D_AD

  SUBROUTINE sqrt_K(x,x_k,y_k, acc, map, prof_stat, limits)
    REAL(jprb), INTENT(in) :: x(:,:), y_k(:,:)
    REAL(jprb), INTENT(inout) :: x_k(:,:)
    LOGICAL(jplm), INTENT(IN) :: acc
    INTEGER(jpim), INTENT(IN) :: map(:,:)
    INTEGER(jpim), INTENT(IN) :: prof_stat
    INTEGER(jpim), INTENT(IN), OPTIONAL :: limits(:,:)

    INTEGER(jpim) :: loc, lb, ub, prof, i, ll, ul
    REAL(jprb) :: z(SIZE(x(:,1)))
    
    loc = 1

    IF (prof_stat >= 0) THEN
      DO prof = 1, MAXVAL(map(:,1))
        lb = loc
        ub = lb
        DO i = lb + 1, SIZE(map(:,1))
          IF (map(i,1) .NE. prof) THEN
            ub = i - 1
            loc = i
            EXIT
          ENDIF
          ub = SIZE(map(:,1))
        ENDDO
        
        CALL reciprocal(x(:,prof), z(:))
      
        IF (PRESENT(limits)) THEN
          IF (acc) THEN
            DO i=lb,ub
              ul = limits(2,map(i,2))
              ll = limits(1,map(i,2))
              x_k(ll:ul,i) = x_k(ll:ul,i) + 0.5_jprb * z(ll:ul) * y_k(ll:ul,i)
            ENDDO
          ELSE
            DO i=lb,ub
              ul = limits(2,map(i,2))
              ll = limits(1,map(i,2))
              x_k(ll:ul,i) = 0.5_jprb * z(ll:ul) * y_k(ll:ul,i)
            ENDDO
          ENDIF
        ELSE
          IF (acc) THEN
            DO i=lb,ub
              x_k(:,i) = x_k(:,i) + 0.5_jprb * z(:) * y_k(:,i)
            ENDDO
          ELSE
            DO i=lb,ub
              x_k(:,i) = 0.5_jprb * z(:) * y_k(:,i)
            ENDDO
          ENDIF
        ENDIF
      ENDDO
      
    ELSE IF (prof_stat .EQ. -1) THEN
      DO i = 1, SIZE(y_k(1,:))
        prof = map(i,1)

        CALL reciprocal(x(:,prof), z(:))

        IF (acc) THEN
          x_k(:,i) = x_k(:,i) + 0.5_jprb * z(:) * y_k(:,i)
        ELSE
          x_k(:,i) = 0.5_jprb * z(:) * y_k(:,i)
        ENDIF
      ENDDO
    ENDIF
  END SUBROUTINE sqrt_K

  SUBROUTINE sintosec_2D(x, y)
    REAL(jprb), INTENT(in) :: x(:,:)
    REAL(jprb), INTENT(out) :: y(:,:)

#ifdef RTTOV_INTEL_MKL
    CALL vdsqr(SIZE(x),x,y)
    CALL vdinvsqrt(SIZE(x),1._jprb - y, y)
#elif RTTOV_IBM_MASS
    y = x**2
    CALL vrsqrt(y,1._jprb - y,SIZE(x)) ! check these
#else
    y(:,:) = 1._jprb / SQRT((1._jprb - x(:,:)**2))
#endif
    
  END SUBROUTINE SINTOSEC_2D

  SUBROUTINE sintosec_2D_tl(x, x_tl, y, y_tl)
    REAL(jprb), INTENT(in) :: x(:,:), x_tl(:,:), y(:,:)
    REAL(jprb), INTENT(out) :: y_tl(:,:)

    y_tl(:,:) = y**3 * x * x_tl
    
  END SUBROUTINE SINTOSEC_2D_TL

  SUBROUTINE sintosec_2D_ad(x, x_ad, y, y_ad, acc)
    REAL(jprb), INTENT(in) :: x(:,:), y(:,:), y_ad(:,:)
    REAL(jprb), INTENT(inout) :: x_ad(:,:)
    LOGICAL(jplm), INTENT(IN) :: acc
    
    IF (acc) THEN
      x_ad = x_ad + y ** 3 * x * y_ad
    ELSE
      x_ad = y ** 3 * x * y_ad
    ENDIF
    
  END SUBROUTINE SINTOSEC_2D_ad

  SUBROUTINE sintosec_2D_k(x, x_k, y, y_k, acc, map, prof_stat)
    REAL(jprb), INTENT(in) :: x(:,:), y(:,:), y_k(:,:)
    REAL(jprb), INTENT(inout) :: x_k(:,:)
    LOGICAL(jplm), INTENT(IN) :: acc
    INTEGER(jpim), INTENT(IN) :: map(:,:)
    INTEGER(jpim), INTENT(IN) :: prof_stat

    INTEGER(jpim) :: loc, lb, ub, prof, i
    REAL(jprb) :: z(SIZE(y(:,1)))

    loc = 1

    IF (prof_stat >= 0) THEN
      DO prof = 1, MAXVAL(map(:,1))
        lb = loc
        ub = lb
        DO i = lb + 1, SIZE(map(:,1))
          IF (map(i,1) .NE. prof) THEN
            ub = i - 1
            loc = i
            EXIT
          ENDIF
          ub = SIZE(map(:,1))
        ENDDO

        z(:) = y(:,prof)**3 * x(:,prof)
        
        IF (acc) THEN
          DO i=lb,ub
            x_k(:,i) = x_k(:,i) + y_k(:,i) * z(:)
          ENDDO
!          x_k(:,lb:ub) = x_k(:,lb:ub) + y_k(:,lb:ub) * &
!            SPREAD(y(:,prof)**3 * x(:,prof),2 ,ub - lb + 1) 
        ELSE
          DO i=lb,ub
            x_k(:,i) = y_k(:,i) * z(:)
          ENDDO
!          x_k(:,lb:ub) = y_k(:,lb:ub) * &
!            SPREAD(y(:,prof)**3 * x(:,prof),2 ,ub - lb + 1) 
        ENDIF
      ENDDO
      
    ELSEIF (prof_stat == -1) THEN
      DO i = 1, SIZE(y_k(1,:))
        prof = map(i,1)
        IF (acc) THEN
          x_k(:,i) = x_k(:,i) + y(:, prof)**3 * x(:, prof) * y_k(:,i)
        ELSE
          x_k(:,i) = y(:, prof)**3 * x(:, prof) * y_k(:,i)
        ENDIF
      ENDDO
    ENDIF
  END SUBROUTINE Sintosec_2D_k

  SUBROUTINE PLANCK_SCALAR(c1, c2, t, out)

    REAL(jprb), INTENT(in) :: c1, c2
    REAL(jprb), INTENT(in) :: t
    REAL(jprb), INTENT(out) :: out

    out = c1 / (EXP(c2 / t) - 1._jprb)

  END SUBROUTINE PLANCK_SCALAR

  SUBROUTINE PLANCK_SCALAR_NU(nu, c1, c2, t, out)

    REAL(jprb), INTENT(in) :: nu
    REAL(jprb), INTENT(in) :: c1, c2
    REAL(jprb), INTENT(in) :: t
    REAL(jprb), INTENT(out) :: out

    CALL PLANCK_SCALAR(c1 * nu**3, c2 * nu, t, out)

  END SUBROUTINE PLANCK_SCALAR_NU

  SUBROUTINE PLANCK_1D(c1, c2, t, out)

    REAL(jprb), INTENT(in) :: c1, c2
    REAL(jprb), INTENT(in) :: t(:)
    REAL(jprb), INTENT(out) :: out(:)

#ifdef RTTOV_INTEL_MKL
    REAL(jprb) :: z(SIZE(t))
    CALL vdinv(SIZE(t), t, z)
    CALL scal(z, c2)
    CALL vdexp(SIZE(t), z, z)
    CALL vdinv(SIZE(t), z - 1._jprb, out)
    CALL scal(out, c1)
#elif RTTOV_IBM_MASS
    REAL(jprb) :: z(SIZE(t))
    CALL vrec(z, t, SIZE(t))
    CALL dscal(SIZE(t), c2, z, 1)
    CALL vexp(z, z, SIZE(t))
    CALL vrec(out, z - 1._jprb, SIZE(t))
    CALL dscal(SIZE(t), c1, out, 1)
#else
    out = c1 / (EXP(c2 / t(:)) - 1._jprb)
#endif

  END SUBROUTINE PLANCK_1D

  SUBROUTINE PLANCK_1D_NU(nu, c1, c2, t, out)

    REAL(jprb), INTENT(in) :: nu
    REAL(jprb), INTENT(in) :: c1, c2
    REAL(jprb), INTENT(in) :: t(:)
    REAL(jprb), INTENT(out) :: out(:)

    CALL PLANCK_1D(c1 * nu**3, c2 * nu, t, out)

  END SUBROUTINE PLANCK_1D_NU

  SUBROUTINE PLANCK_SCALAR_TL(c1, c2, t, t_tl, B, B_tl)
    
    REAL(jprb), INTENT(in) :: c1, c2
    REAL(jprb), INTENT(in) :: B
    REAL(jprb), INTENT(in) :: t, t_tl
    REAL(jprb), INTENT(out) :: B_tl

    B_tl = ((c2 * B * (c1 + B)) / (c1 * t**2)) * t_tl

  END SUBROUTINE PLANCK_SCALAR_TL

  SUBROUTINE PLANCK_SCALAR_NU_TL(nu, c1, c2, t, t_tl, B, B_tl)

    REAL(jprb), INTENT(in) :: nu
    REAL(jprb), INTENT(in) :: c1, c2
    REAL(jprb), INTENT(in) :: B
    REAL(jprb), INTENT(in) :: t, t_tl
    REAL(jprb), INTENT(out) :: B_tl

    CALL PLANCK_SCALAR_TL(c1 * nu**3, c2 * nu, t, t_tl, B, B_tl)

  END SUBROUTINE PLANCK_SCALAR_NU_TL

  SUBROUTINE PLANCK_1D_TL(c1, c2, t, t_tl, B, B_tl)
    REAL(jprb), INTENT(in) :: c1, c2
    REAL(jprb), INTENT(in) :: t(:), t_tl(:)
    REAL(jprb), INTENT(in) :: B(:)
    REAL(jprb), INTENT(out) :: B_tl(:)

#ifdef RTTOV_INTEL_MKL
    REAL(jprb) :: z(SIZE(t)), y(SIZE(t))
    CALL vdsqr(SIZE(t), t, y)
    CALL vddiv(SIZE(t), B, y, z)
    CALL vdmul(SIZE(t), z, (c2 + (c2/c1)*B), y)
    CALL vdmul(SIZE(t), y, t_tl, B_tl)
#elif RTTOV_IBM_MASS
    REAL(jprb) :: z(SIZE(t)), y(SIZE(t))
    y = t**2
    CALL vdiv(z, B, y, SIZE(t))    
    B_tl = (c2 + (c2/c1)*B) * z * t_tl
#else
    B_tl = ((c2 * B * (c1 + B)) / (c1 * t**2)) * t_tl
#endif

  END SUBROUTINE PLANCK_1D_TL

  SUBROUTINE PLANCK_SCALAR_AD(c1, c2, t, t_ad, B, B_ad, acc)
    
    REAL(jprb), INTENT(in) :: c1, c2
    REAL(jprb), INTENT(in) :: B
    REAL(jprb), INTENT(in) :: t
    REAL(jprb), INTENT(inout) :: t_ad
    REAL(jprb), INTENT(in) :: B_ad
    LOGICAL(jplm), INTENT(IN) :: acc

    IF (acc) THEN
      t_ad = t_ad + ((c2 * B * (c1 + B)) / (c1 * t**2)) * B_ad
    ELSE
      t_ad = ((c2 * B * (c1 + B)) / (c1 * t**2)) * B_ad
    ENDIF

  END SUBROUTINE PLANCK_SCALAR_AD

 SUBROUTINE PLANCK_1D_AD(c1, c2, t, t_ad, B, B_ad, acc)
    REAL(jprb), INTENT(in) :: c1, c2
    REAL(jprb), INTENT(in) :: t(:)
    REAL(jprb), INTENT(inout) :: t_ad(:)
    REAL(jprb), INTENT(in) :: B(:)
    REAL(jprb), INTENT(in) :: B_ad(:)
    LOGICAL(jplm), INTENT(IN) :: acc

#ifdef RTTOV_INTEL_MKL
    REAL(jprb) :: z(SIZE(t)), y(SIZE(t))
    CALL vdsqr(SIZE(t), t, y)
    CALL vddiv(SIZE(t), B, y, z)
    CALL vdmul(SIZE(t), z, (c2 + (c2/c1)*B), y)
    IF (acc) THEN
      CALL vdmul(SIZE(t), y, B_ad, t_ad)
    ELSE
      CALL vdmul(SIZE(t), y, B_ad, z)
      CALL axpy(z, t_ad) ! t_ad = t_ad + z
    ENDIF

#elif RTTOV_IBM_MASS
    REAL(jprb) :: z(SIZE(t)), y(SIZE(t))
    y = t**2
    CALL vdiv(z, B, y, SIZE(t))
    IF (acc) THEN
      t_ad = t_ad + (c2 + (c2/c1)*B) * z * B_ad
    ELSE
      t_ad = (c2 + (c2/c1)*B) * z * B_ad
    ENDIF
#else
    IF (acc) THEN
      t_ad = t_ad + ((c2 * B * (c1 + B)) / (c1 * t**2)) * B_ad
    ELSE
      t_ad = ((c2 * B * (c1 + B)) / (c1 * t**2)) * B_ad
    ENDIF
#endif
  END SUBROUTINE PLANCK_1D_AD

  SUBROUTINE INV_PLANCK_SCALAR_NU(nu, c1, c2, r, out)

    REAL(jprb), INTENT(in) :: nu
    REAL(jprb), INTENT(in) :: c1, c2
    REAL(jprb), INTENT(in) :: r
    REAL(jprb), INTENT(out) :: out

    CALL INV_PLANCK_SCALAR(c1 * nu**3, c2 * nu, r, out)

  END SUBROUTINE INV_PLANCK_SCALAR_NU

  SUBROUTINE INV_PLANCK_SCALAR(c1, c2, r, out)

    REAL(jprb), INTENT(in) :: c1, c2
    REAL(jprb), INTENT(in) :: r
    REAL(jprb), INTENT(out) :: out

    out = c2 / LOG(c1 / r + 1._jprb)

  END SUBROUTINE INV_PLANCK_SCALAR

  SUBROUTINE INV_PLANCK_SCALAR_TL(c1, c2, r, r_tl, t, t_tl)

    REAL(jprb), INTENT(in) :: c1, c2
    REAL(jprb), INTENT(in) :: r, r_tl
    REAL(jprb), INTENT(in) :: t
    REAL(jprb), INTENT(out) :: t_tl

    t_tl = c1 * t**2 / (c2 * r * (r + c1)) * r_tl

  END SUBROUTINE INV_PLANCK_SCALAR_TL

  SUBROUTINE INV_PLANCK_SCALAR_AD(c1, c2, r, r_ad, t, t_ad, acc)

    REAL(jprb), INTENT(in) :: c1, c2
    REAL(jprb), INTENT(in) :: r
    REAL(jprb), INTENT(inout) :: r_ad
    REAL(jprb), INTENT(in) :: t
    REAL(jprb), INTENT(in) :: t_ad
    LOGICAL(jplm), INTENT(in) :: acc

    IF (acc) THEN
      r_ad = r_ad + c1 * t**2 / (c2 * r * (r + c1)) * t_ad
    ELSE
      r_ad = c1 * t**2 / (c2 * r * (r + c1)) * t_ad
    ENDIF
  END SUBROUTINE INV_PLANCK_SCALAR_AD

END MODULE rttov_math_mod
