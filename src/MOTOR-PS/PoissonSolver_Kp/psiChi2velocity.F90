!!--------------------------------------------------------------------------------------------------
! PROJECT         : MOTOR-PS.PoissonSolver_Kp.psiChi2velocity.F90
! AFFILIATION     : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring 
!                   Warning and Forecasting (GBA-MWF) Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Yuanfu Xie, on 2024-12-05
!!--------------------------------------------------------------------------------------------------
!> @brief
!! This routine computes the normal and tangential velocities at cell edges by an option of users.
!! @copyright (C) 2024 GBA-MWF, All rights reserved.
!! @note
!! @warning
!! @attention
SUBROUTINE psiChi2velocity(psi,chi,vnorm,vtang,sg)
  USE kinds_m, ONLY: i_kind, r_kind
  USE parameters_m, ONLY: EarthRadius
  USE SingleGrid_m, ONLY: SingleGrid_t

  IMPLICIT NONE

  TYPE(singleGrid_t), INTENT(IN) :: sg
  REAL(r_kind), INTENT(IN) :: &
    psi(sg%vlevel,sg%num_cell),chi(sg%vlevel,sg%num_cell)
  REAL(r_kind), INTENT(OUT) :: &
    vnorm(sg%vlevel,sg%numQuadPerEdge,sg%num_edge,sg%num_cell), &
    vtang(sg%vlevel,sg%numQuadPerEdge,sg%num_edge,sg%num_cell)

  ! Local variables:
  INTEGER(i_kind) :: ic,ie,iq,is

  ! vn = psi_t + chi_n; vt = -psi_n + chi_t:
  ! Note that the actually physical boundary cells are those interior points sharing
  ! an edge with its neighbor that has cell_type = 1, or 2: Yuanfu Xie 2025-02-14
  ! Here the normal and tangential velocity are computed for all cells, including the
  ! fititious points:
  vnorm = 0.0D0; vtang = 0.0D0
  DO ic=1,sg%num_icell
    DO ie=1,sg%num_edge
      DO iq=1,sg%numQuadPerEdge
        DO is = 1, UBOUND(sg%edge_stcl, 1)
          vnorm(:,iq,ie,ic) = vnorm(:,iq,ie,ic) + &
            psi(:, sg%edge_stcl(is, ie, ic)) * sg%coef_tang(is, iq, ie, ic) + &
            chi(:, sg%edge_stcl(is, ie, ic)) * sg%coef_norm(is, iq, ie, ic)
          vtang(:,iq,ie,ic) = vtang(:,iq,ie,ic) - &
            psi(:, sg%edge_stcl(is, ie, ic)) * sg%coef_norm(is, iq, ie, ic) + &
            chi(:, sg%edge_stcl(is, ie, ic)) * sg%coef_tang(is, iq, ie, ic)
        END DO
      END DO
    END DO
  END DO
END SUBROUTINE psiChi2velocity

SUBROUTINE checkPsiChi2UV_s(psi, chi, u, v, sg, errMax)
  USE kinds_m, ONLY: i_kind, r_kind
  USE SingleGrid_m, ONLY: SingleGrid_t

  IMPLICIT NONE
  TYPE(SingleGrid_t), INTENT(IN) :: sg
  REAL(r_kind), INTENT(IN) :: psi(sg%vLevel,sg%num_cell),chi(sg%vLevel,sg%num_cell)
  REAL(r_kind), INTENT(IN) :: u(sg%vLevel,sg%num_cell),v(sg%vLevel,sg%num_cell)
  REAL(r_kind), INTENT(OUT) :: errMax

  ! Local variables:
  INTEGER(i_kind) :: imx,i
  REAL(r_kind), ALLOCATABLE :: uwind(:,:), vwind(:,:)

  ! Allocate memory for computed wind:
  ALLOCATE(uwind(sg%vLevel,sg%num_cell),vwind(sg%vLevel,sg%num_cell))
  ! Convert to U and V:
  CALL uvVelocityOnInterior(psi,chi,uwind,vwind,sg)

  ! Find the max error in V component only (extension is needed for other test cases):
  imx = 0
  errMax = 0.0D0
  DO i=1,sg%num_cell
    IF (sg%cell_type(i) .LE. 2) THEN
      IF (ABS(v(1,i)-vwind(1,i)) .GT. errMax) THEN
        imx = i
        errMax = ABS(v(1,i)-vwind(1,i))
      END IF
    END IF
  END DO
  WRITE(*,1) errMax,imx,v(1,imx),vwind(1,imx)
1     FORMAT('V computed - anaV: ',E14.6,' at:',I6,' true/computed: ',2E14.6)
  PRINT*,'Calculated U component: ',MAXVAL(ABS(uwind(1,:)))
END SUBROUTINE checkPsiChi2UV_s