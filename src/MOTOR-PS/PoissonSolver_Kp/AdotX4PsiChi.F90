!!--------------------------------------------------------------------------------------------------
! PROJECT         : MOTOR-DA.AdotX4PsiChi.F90
! AFFILIATION     : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring 
!                   Warning and Forecasting (GBA-MWF) Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Yuanfu Xie, on 2024-12-05
!   Updated by Yuanfu Xie, on 2025-03-05  modified uvVelocityOnInterior routine to output wind velocity at
!     cell centers
!!--------------------------------------------------------------------------------------------------
!> @brief
!! This routine is required by GMRES solver for solving streamfunction and velocity potential from 
!!  vorticity and divergence with u and v as boundary conditions.
!! @copyright (C) 2024 GBA-MWF, All rights reserved.
!! @note
!! @warning
!! @attention
SUBROUTINE AdotX4PsiChi_s(n,x,b, gzm, scaling)  ! b = Ax where A is the coefficient matrix for PsiChi solver
  USE kinds_m, ONLY: i_kind, r_kind
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE gzm_m, ONLY: gzm_t

  IMPLICIT NONE

  INTEGER(i_kind), INTENT(IN) :: n    ! n = total grid cells of Psi and Chi
  REAL(r_kind), INTENT(IN) :: x(n)
  REAL(r_kind), INTENT(OUT) :: b(n)
  TYPE(gzm_t), INTENT(IN) :: gzm
  REAL(r_kind), INTENT(IN) :: scaling(2)

  ! Local variables:
  INTEGER(i_kind) :: ic,j,iq,is,half
  REAL(r_kind), ALLOCATABLE :: u(:),v(:),norm_value(:),tang_value(:)
  REAL(r_kind), ALLOCATABLE :: psi(:,:),chi(:,:),vor(:,:),div(:,:),uwind(:,:),vwind(:,:)

  b = 0.0D0

  ! Check the dimensions:
  IF (n .NE. 2*gzm%sg%vLevel*gzm%sg%num_cell) THEN
    WRITE(*,1) n,gzm%sg%vLevel,gzm%sg%num_cell
1   FORMAT('AdotX4PsiChi_s: GMRES requires n: ',I8,' to be the total cell numbers: ',I3,I8)
    STOP
  END IF
  half = gzm%sg%vLevel * gzm%sg%num_cell

  ! Allocate memory for psi and chi:
  ALLOCATE(psi(gzm%sg%vLevel,gzm%sg%num_cell), chi(gzm%sg%vLevel,gzm%sg%num_cell), &
    uwind(gzm%sg%vLevel,gzm%sg%num_cell), &
    vwind(gzm%sg%vLevel,gzm%sg%num_cell), &
    vor(gzm%sg%vLevel,gzm%sg%num_cell), div(gzm%sg%vLevel,gzm%sg%num_cell))
  ALLOCATE(u(gzm%sg%vLevel),v(gzm%sg%vLevel),norm_value(gzm%sg%vLevel),tang_value(gzm%sg%vLevel))

  ! Initialize:
  psi = RESHAPE(x(1:n/2), (/gzm%sg%vLevel,gzm%sg%num_cell/))
  chi = RESHAPE(x(n/2+1:n), (/gzm%sg%vLevel,gzm%sg%num_cell/))

  ! Scaling psi and chi to physical space:
  psi = psi/scaling(1)
  chi = chi/scaling(1)

  ! x is psi and chi at all vertical levels and horizontal grid:
  CALL gzm%Laplacia(psi,vor)
  b(1                              :  gzm%sg%vLevel*gzm%sg%num_cell) = &
    RESHAPE(vor,(/gzm%sg%vLevel*gzm%sg%num_cell/))
  CALL gzm%Laplacia(chi,div)
  b(gzm%sg%vLevel*gzm%sg%num_cell+1:2*gzm%sg%vLevel*gzm%sg%num_cell) = &
    RESHAPE(div,(/gzm%sg%vLevel*gzm%sg%num_cell/))
  
  ! Scaling the right hand side of the interior points:
  b = b/scaling(1)
  WRITE(*,111) dot_product(b,b),scaling(1),gzm%sg%vLevel
111 FORMAT('B psi at 1:2: ',E14.6,' scaling: ',E14.6,' v-g lvls: ',I3)

  ! Get the boundary values of u and v:
  CALL uvVelocityOnInterior(psi,chi,uwind,vwind,gzm%sg)

  ! Pass the boundary conditions:
  ! DO ic=1,gzm%sg%num_cell
  !   ! For boundary points ONLY, no corner points:
  !   IF (gzm%sg%cell_type(ic) .EQ. 1) THEN
  !     DO j=1,UBOUND(gzm%sg%cell_adnb,1)
  !       ! Skip those cells outside the domain:
  !       IF (gzm%sg%cell_adnb(j,ic) .NE. 0) THEN
  !         IF (gzm%sg%cell_type(gzm%sg%cell_adnb(j,ic)) .EQ. 0) THEN
  !           ! Found an edge divides the boundary:
  !           b(gzm%sg%vLevel*(ic-1)+1     :gzm%sg%vLevel*ic     ) = uwind(:,ic)*scaling(2)
  !           b(gzm%sg%vLevel*(ic-1)+1+half:gzm%sg%vLevel*ic+half) = vwind(:,ic)*scaling(2)
  !         END IF
  !       END IF
  !     END DO
  !   END IF
  ! END DO
  print*,'BB norm: ',DOT_PRODUCT(b,b),DOT_PRODUCT(x,x)
END SUBROUTINE AdotX4PsiChi_s

SUBROUTINE uvVelocityOnInterior(psi,chi,uwind,vwind,sg)
  USE kinds_m, ONLY: i_kind, r_kind
  USE SingleGrid_m, ONLY: SingleGrid_t
  IMPLICIT NONE
  TYPE(singleGrid_t), INTENT(IN) :: sg
  REAL(r_kind), INTENT(IN) :: psi(sg%vlevel,sg%num_cell),chi(sg%vlevel,sg%num_cell)
  REAL(r_kind), INTENT(OUT) :: uwind(sg%vlevel,sg%num_cell),vwind(sg%vlevel,sg%num_cell)

  ! Local variables:
  INTEGER(i_kind) :: ic,j,iq,is,nedge
  REAL(r_kind), ALLOCATABLE :: u(:),v(:),norm_value(:),tang_value(:)

  ALLOCATE(u(sg%vLevel),v(sg%vLevel),norm_value(sg%vLevel),tang_value(sg%vLevel))

  ! Boundary conditions:
  uwind = 0.0D0; vwind = 0.0D0
  DO ic=1,sg%num_cell
    ! For all points except halos:
    IF (sg%cell_type(ic) .LE. 2) THEN
      nedge = 0
      DO j=1,UBOUND(sg%cell_adnb,1)
        ! Calculate the velocity at each edge:
        nedge = nedge + 1
        u = 0.0D0
        v = 0.0D0
        DO iq = 1, sg%numQuadPerEdge
          norm_value = 0.0D0
          tang_value = 0.0D0
          DO is = 1, UBOUND(sg%edge_stcl, 1)
            norm_value = norm_value + &
              chi(:, sg%edge_stcl(is, j, ic)) * sg%coef_norm(is, iq, j, ic)
            tang_value = tang_value + &
              psi(:, sg%edge_stcl(is, j, ic)) * sg%coef_norm(is, iq, j, ic)
          END DO
          u = u + sg%coef_gl(iq) * norm_value         ! u_n =  chi_n
          v = v + sg%coef_gl(iq) * tang_value         ! v_t =  psi_n

          norm_value = 0.0D0
          tang_value = 0.0D0
          DO is = 1, UBOUND(sg%edge_stcl, 1)
            norm_value = norm_value + &
              psi(:, sg%edge_stcl(is, j, ic)) * sg%coef_tang(is, iq, j, ic)
            tang_value = tang_value + &
              chi(:, sg%edge_stcl(is, j, ic)) * sg%coef_tang(is, iq, j, ic)
          END DO
          u = u - sg%coef_gl(iq) * norm_value         ! u_n += -psi_t
          v = v + sg%coef_gl(iq) * tang_value         ! v_t +=  chi_t

        END DO
        u = u/DBLE(sg%numQuadPerEdge) ! Averaged normal velocity
        v = v/DBLE(sg%numQuadPerEdge) ! Averaged tangential velocity

        ! Replace the b value with BC values:
        uwind(:,ic) = uwind(:,ic) + u*sg%edgeNorm2(1,j,ic) + v*sg%edgeTang2(1,j,ic)
        vwind(:,ic) = vwind(:,ic) + u*sg%edgeNorm2(2,j,ic) + v*sg%edgeTang2(2,j,ic)
      END DO
      ! Average the velocity at the cell center:
      uwind(:,ic) = uwind(:,ic)/DBLE(nedge)
      vwind(:,ic) = vwind(:,ic)/DBLE(nedge)
    END IF
  END DO
END SUBROUTINE uvVelocityOnInterior

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