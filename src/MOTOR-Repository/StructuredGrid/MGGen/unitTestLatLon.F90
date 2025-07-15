!!--------------------------------------------------------------------------------------------------
! PROJECT           : Multi-Grid generation
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : Beta 0.0
! HISTORY           :
!   Created by Yuanfu Xie (yuanfu_xie@yahoo.com), 2020/05/29, @GBA-MWF, Shenzhen
!   Modified by ..... (???@gmail.com), YYYY/MM/DD, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!!===================================================================
!> @brief
!! # Multi-Grid Latlon Generation unit test
!!
!!  *This module sets a unit test for a regional latlon grid*
!!
!! @author Yuanfu Xie
!! @copyright (C) 2020 GBA-MWF, All rights reserved.
!! @test This is a beta testing
!!
!!===================================================================
SUBROUTINE unitTestLatLon(configFile, PASS)
  USE kinds_m, ONLY: i_kind, r_kind
  USE mgGen_m, ONLY: mgGen_t, mgGenLatlon_t
  USE parameters_m, ONLY: degree2radian, EarthRadius

  IMPLICIT NONE

  LOGICAL, INTENT(OUT) :: PASS
  CHARACTER(LEN=1024), INTENT(IN) :: configFile

  ! Local variables:
  INTEGER(i_kind) :: nv, nc, numgrd(2), i, j, k, ig, imx(20), jmx(20), var(20), itmx(20), immx(20)
  REAL(r_kind) :: domain(2, 2), emx(20), etmx(20), emmx(20), &
                  utrue, vtrue, uc, vc, vorticity, derivativ
  TYPE(mgGenLatlon_t) :: grid
  REAL(r_kind), ALLOCATABLE :: stream(:), velpot(:), ucomp(:, :), vcomp(:), &
                               vort(:), divg(:), utemp(:), vtemp(:)

  ! Get the latlon grid information:
  CALL getDomain_latlon(configFile, nv, nc, numgrd, domain)

  CALL grid%mlgrid%read_mlgrid(.FALSE.)
  grid%mglvls = grid%mlgrid%num_lvls
  grid%mgstts = grid%mlgrid%lvl_stts
  grid%mgends = grid%mlgrid%lvl_ends

  ! Allocate memory for testing functions:
  ALLOCATE (stream(nc), velpot(nc), ucomp(nc, grid%mgstts:grid%mgends), vcomp(nc), &
            vort(nc), divg(nc), utemp(nc), vtemp(nc))

  IF (grid%mglvls .LT. 4) THEN
    PRINT *, 'Unit test requires at least 5 levels of a multigrid...', grid%mglvls
    STOP
  END IF

  ! Calculate the maximum errors at each level:
  imx = 0
  jmx = 0
  emx = 0.0D0
  var = 0
  itmx = 0
  etmx = 0.0D0
  immx = 0
  emmx = 0.0D0
  DO ig = grid%mgstts, grid%mgends
    ! Calculate the stream function, velocity and vorticity of a simple polynomial:
    DO i = 1, grid%mlgrid%Params(ig)%num_cell
      stream(i) = grid%mlgrid%Geoqty(ig)%cell_cntr(1, i)**3 * EarthRadius + &
                  grid%mlgrid%Geoqty(ig)%cell_cntr(2, i)**3 * EarthRadius
      velpot(i) = grid%mlgrid%Geoqty(ig)%cell_cntr(1, i)**3 * EarthRadius + &
                  grid%mlgrid%Geoqty(ig)%cell_cntr(2, i)**3 * EarthRadius
      ucomp(i, ig) = -3.0D0 * grid%mlgrid%Geoqty(ig)%cell_cntr(1, i)**2 + &
                     3.0D0 * grid%mlgrid%Geoqty(ig)%cell_cntr(2, i)**2 &
                     / DCOS(grid%mlgrid%Geoqty(ig)%cell_cntr(1, i))
      vcomp(i) = 3.0D0 * grid%mlgrid%Geoqty(ig)%cell_cntr(1, i)**2 + &
                 3.0D0 * grid%mlgrid%Geoqty(ig)%cell_cntr(2, i)**2 &
                 / DCOS(grid%mlgrid%Geoqty(ig)%cell_cntr(1, i))
      vort(i) = 3.0D0 * (2.0D0 * grid%mlgrid%Geoqty(ig)%cell_cntr(1, i)**1 - &
                         grid%mlgrid%Geoqty(ig)%cell_cntr(1, i)**2 * &
                         DTAN(grid%mlgrid%Geoqty(ig)%cell_cntr(1, i)) + &
                         2.0D0 * grid%mlgrid%Geoqty(ig)%cell_cntr(2, i)**1 &
                         / DCOS(grid%mlgrid%Geoqty(ig)%cell_cntr(1, i)) &
                         / DCOS(grid%mlgrid%Geoqty(ig)%cell_cntr(1, i))) &
                / EarthRadius
    END DO

    DO i = 1, grid%mlgrid%Params(ig)%num_cell
      ! Check if it encounts an boundary cell:
      IF (MINVAL(grid%mlgrid%Geoqty(ig)%edge_stcl(:, :, i)) .LE. 0) THEN
        DO j = 1, grid%mlgrid%Params(ig)%numVrtxPerCell
          ! Found a boundary edge:
          IF (grid%mlgrid%Geoqty(ig)%edge_stcl(1, j, i) .LT. 0) THEN
            utemp(i) = 0.0D0
            vtemp(i) = 0.0D0
            utrue = 0.0D0
            vtrue = 0.0D0
            uc = 0.0D0
            vc = 0.0D0
            IF (MINVAL(ABS(grid%mlgrid%Geoqty(ig)%edge_stcl(:, j, i))) .GT. 0) THEN
            DO k = 1, grid%mlgrid%Params(ig)%numStclPerEdge
              ! Assuming mid_edge quadrature is used: second index
              uc = uc - stream(ABS(grid%mlgrid%Geoqty(ig)%edge_stcl(k, j, i))) * &
                   grid%mlgrid%Geoqty(ig)%coef_tgnt(k, 1, j, i) + &
                   velpot(ABS(grid%mlgrid%Geoqty(ig)%edge_stcl(k, j, i))) * &
                   grid%mlgrid%Geoqty(ig)%coef_norm(k, 1, j, i)
              vc = vc + stream(ABS(grid%mlgrid%Geoqty(ig)%edge_stcl(k, j, i))) * &
                   grid%mlgrid%Geoqty(ig)%coef_norm(k, 1, j, i) + &
                   velpot(ABS(grid%mlgrid%Geoqty(ig)%edge_stcl(k, j, i))) * &
                   grid%mlgrid%Geoqty(ig)%coef_tgnt(k, 1, j, i)

              utrue = utrue + ucomp(ABS(grid%mlgrid%Geoqty(ig)%edge_stcl(k, j, i)), ig) * &
                      grid%mlgrid%Geoqty(ig)%coef_func(k, 1, j, i)
              vtrue = vtrue + vcomp(ABS(grid%mlgrid%Geoqty(ig)%edge_stcl(k, j, i))) * &
                      grid%mlgrid%Geoqty(ig)%coef_func(k, 1, j, i)
            END DO
            END IF
            utemp(i) = grid%mlgrid%Geoqty(ig)%norm_vct2(1, j, i) * uc + &
                       grid%mlgrid%Geoqty(ig)%tgnt_vct2(1, j, i) * vc
            vtemp(i) = grid%mlgrid%Geoqty(ig)%norm_vct2(2, j, i) * uc + &
                       grid%mlgrid%Geoqty(ig)%tgnt_vct2(2, j, i) * vc
            IF (ABS(utemp(i) - utrue) .GT. emx(ig)) THEN
              imx(ig) = i
              jmx(ig) = j
              var(ig) = 1
              emx(ig) = ABS(utemp(imx(ig)) - utrue)
              !write(*,11)ig,i,j,utemp(i),utrue
11            FORMAT('U computed ', I2, I8, I2, 2F14.6)
              !stop
            END IF
            IF (ABS(vtemp(i) - vtrue) .GT. emx(ig)) THEN
              imx(ig) = i
              jmx(ig) = j
              emx(ig) = ABS(vtemp(imx(ig)) - vtrue)
              var(ig) = 2
              !write(*,12)ig,i,j,vtemp(i),vtrue
12            FORMAT('V computed ', I2, I8, I2, 2F14.6)
              !stop
            END IF
          END IF
        END DO
      ELSE
        vorticity = 0.0D0
        DO j = 1, grid%mlgrid%Params(ig)%numVrtxPerCell
          derivativ = 0.0D0
          DO k = 1, grid%mlgrid%Params(ig)%numStclPerEdge
            derivativ = derivativ + &
                        stream(grid%mlgrid%Geoqty(ig)%edge_stcl(k, j, i)) * &
                        grid%mlgrid%Geoqty(ig)%coef_norm(k, 1, j, i)
          END DO
          vorticity = vorticity + derivativ * grid%mlgrid%Geoqty(ig)%edge_lnth(j, i)
        END DO
        vorticity = vorticity / grid%mlgrid%Geoqty(ig)%cell_area(i)
        IF (ABS(vorticity - vort(i)) .GT. etmx(ig)) THEN
          itmx(ig) = i
          etmx(ig) = ABS(vorticity - vort(i))
        END IF
      END IF
    END DO
    !print*,'Max error at boundaries: ',emx(ig),imx(ig),jmx(ig),var(ig),ig
    !print*,'Max error at interiors : ',etmx(ig),itmx(ig),ig

    ! Check the multigrid prolongations:
    IF (ig .GT. grid%mgstts) THEN
      DO i = 1, grid%mlgrid%Params(ig)%num_cell
        uc = 0.0D0
        DO k = 1, grid%mlgrid%Params(ig)%num_chld
          uc = uc + &
               ucomp(grid%mlgrid%GeoQty(ig)%cTof_stcl(k, i), ig - 1) * &
               grid%mlgrid%GeoQty(ig)%cTof_coef(k, i)
        END DO
        IF (ABS(uc - ucomp(i, ig)) .GT. emmx(ig)) THEN
          immx(ig) = i
          emmx(ig) = ABS(uc - ucomp(i, ig))
          !write(*,13)ig,i,uc,ucomp(i,ig)
          !write(*,14) grid%mlgrid%GeoQty(ig)%cTof_stcl(1:4,i)
13        FORMAT('U interpolated ', I2, I8, 2F14.6)
14        FORMAT('Child cells:', 4I8)
        END IF
      END DO
      !print*,'Max error in multigrid intp : ',emmx(ig),immx(ig),ig
    END IF
  END DO

  ! Check error reductions:
  DO ig = grid%mgstts + 1, grid%mgends
    IF (emx(ig - 1) / emx(ig) .LT. 3.5) PASS = .FALSE.
    PRINT *, 'Boundary error reduction: ', emx(ig - 1) / emx(ig), ' at level:', ig
  END DO
  DO ig = grid%mgstts + 1, grid%mgends
    IF (etmx(ig - 1) / etmx(ig) .LT. 3.0) PASS = .FALSE.
    PRINT *, 'Interior error reduction: ', etmx(ig - 1) / etmx(ig), ' at level:', ig
  END DO

  ! Check the error reductions between multigrid prolongation:
  DO ig = grid%mgstts + 2, grid%mgends
    IF (emmx(ig - 1) / emmx(ig) .LT. 3.1) PASS = .FALSE.
    PRINT *, 'Multigrid error reduction: ', emmx(ig - 1) / emmx(ig), ' at level:', ig
  END DO

  DEALLOCATE (stream, velpot, ucomp, vcomp, vort, divg, utemp, vtemp)

END SUBROUTINE unitTestLatLon
