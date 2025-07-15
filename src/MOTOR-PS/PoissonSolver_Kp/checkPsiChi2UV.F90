
!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA.Analytic.F90
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring 
!                     Warning and Forecasting (GBA-MWF) Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : V 0.0
! Development docs  : "A Mixed Solver for Stream Function and Velocity Potential with Velocity Components
!                     as Boundary Conditions and Its Adjoint Solver.docs"
! HISTORY           :
!   Created by Yuanfu Xie on 2024-12-09
!   Updated by Yuanfu Xie on 2025-03-05
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This routine constructs analytic functions of Psi, Chi and voricity and divergence with U and V at
!  cell centers for testing a latitude and longitude grid.
!! @copyright (C) 2020 GBA-MWF, All rights reserved.
! @note
! @warning
! Check if the operators from Psi and Chi to U V are correct:
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