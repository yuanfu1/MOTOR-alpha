!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Jia Wang
! VERSION           : V 0.0
! HISTORY           :
!   Created by Jia Wang (), 2022/3/20, @GBA-MWF, Shenzhen
!   Modified by Yuanfu Xie (yuanfu_xie@yahoo.com), 2022/07/27, @GBA-MWF, Shenzhen
!     for changing the balance from incremental scheme to a full nonlinear constraints.
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module contains the UV2W Object which is tranform the variables from control variable space to
!! model space.
MODULE GeosBal_m
  USE State_m, ONLY: State_t
  USE Field_m, ONLY: Field_t
  USE SingleGrid_m, ONLY: SingleGrid_t, GridIdx_t
  USE kinds_m, ONLY: i_kind, r_kind, r_double
  USE ObsSet_m, ONLY: ObsSet_t
  USE M2MBase_m, ONLY: M2MBase_t
  USE parameters_m, ONLY: Coriolis, Omega, dry_air_gas_const

  USE YAMLRead_m
  USE State2NC_m

  IMPLICIT NONE

  TYPE, EXTENDS(M2MBase_t) :: GeosBal_t
    TYPE(State_t), POINTER  :: X

    INTEGER(i_kind) :: geobalBeg = 0, geobalEnd = 0, iq = 0, ip = 0, ir = 0, it = 0, iu = 0, iv = 0, ipr = 0
    INTEGER(i_kind) :: ng(3) = 0
    REAL(r_kind) :: weight = 0.0D0  ! Weak Geostrophic balance weight

  CONTAINS
    PROCEDURE, PUBLIC, PASS(this) :: initial

    PROCEDURE, PUBLIC, PASS(this) :: transFwdNonLinear_opr
    PROCEDURE, PUBLIC, PASS(this) :: transFwdTanLinear_opr
    PROCEDURE, PUBLIC, PASS(this) :: transAdjMultiply_opr

    PROCEDURE, PUBLIC, PASS(this) :: transFwdNonLinear
    PROCEDURE, PUBLIC, PASS(this) :: transFwdTanLinear
    PROCEDURE, PUBLIC, PASS(this) :: transAdjMultiply

    PROCEDURE :: fwdNL_opr => transFwdNonLinear_opr
    PROCEDURE :: fwdTL_opr => transFwdTanLinear_opr
    PROCEDURE :: adjMul_opr => transAdjMultiply_opr

    PROCEDURE :: fwdNL => transFwdNonLinear
    PROCEDURE :: fwdTL => transFwdTanLinear
    PROCEDURE :: adjMul => transAdjMultiply

    ! Weak constraint procedures:
    PROCEDURE, PUBLIC :: dLogP
    PROCEDURE, PUBLIC :: JcGeos
    PROCEDURE, PUBLIC :: dJcDUV

    FINAL :: destructor
  END TYPE GeosBal_t

CONTAINS

  ! Use of initial to replace constructor as Intel compiler does not support constructor well (Qin Zilong)
  SUBROUTINE initial(this, configFile, X)
    IMPLICIT NONE
    CLASS(GeosBal_t) :: this
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    TYPE(State_t), TARGET, INTENT(IN) :: X

    ! Local variables:
    INTEGER(i_kind) :: istatus
    REAL(r_kind), ALLOCATABLE :: weights(:)

    this%X => X

    ! Select the requested control:
    this%ip = X%getVarIdx(TRIM('pres_ctl'))
    this%ir = X%getVarIdx(TRIM('rho_ctl'))
    this%iq = X%getVarIdx(TRIM('qvapor_ctl'))
    this%ipr = 0
    IF (this%ip .GT. 0) THEN
      this%ipr = this%ip
    ELSE
      this%ipr = this%ir
    END IF
    this%it = X%getVarIdx(TRIM('temp'))
    this%iu = X%getVarIdx(TRIM('uwnd'))
    this%iv = X%getVarIdx(TRIM('vwnd'))

    ! Currently, we implement controls of (rho or pres), temp, uwnd and vwnd:
    IF (this%ipr .EQ. 0 .OR. this%it .EQ. 0 .OR. &
        this%iu .EQ. 0 .OR. this%iv .EQ. 0 .OR. this%iq .EQ. 0) THEN
      IF (X%sg%isBaseProc()) THEN
        PRINT *, ''
        WRITE (*, 1)
1       FORMAT('#========================================================================================#')
        WRITE (*, 2) this%ip, this%ir, this%it, this%iu, this%iv
2       FORMAT('GeosBal - requires (pres or rho) temp, u and v must be present for Geostrophic balance', 5I2, /, &
               '  the current analysis is no Geostrophic balance applied')
        WRITE (*, 2)
        PRINT *, ''
      END IF
    END IF

    ! Set grid dimensions:
    this%ng(1) = UBOUND(X%fields(1)%DATA, 1)
    this%ng(2) = UBOUND(X%fields(1)%DATA, 2)
    this%ng(3) = UBOUND(X%fields(1)%DATA, 3)

    ! Get the levels of Geostrophic balance applied:
    IF (yaml_get_var(configFile, 'modelState', 'geobalBeg', this%geobalBeg) /= 0) STOP
    IF (yaml_get_var(configFile, 'modelState', 'geobalEnd', this%geobalEnd) /= 0) STOP
    IF (istatus .NE. 0) THEN
      this%geobalBeg = 0
      this%geobalEnd = 0
    END IF
    IF (X%sg%isBaseProc()) THEN
      IF (this%geobalBeg .LE. X%sg%gLevel .AND. this%geobalEnd .GE. X%sg%gLevel) THEN
        WRITE (*, 4) this%geobalBeg, this%geobalEnd
4       FORMAT('Geostrophic balance is applied between: ', 2I3)
      ELSE
        WRITE (*, 5) this%geobalBeg, this%geobalEnd
5       FORMAT('Geostrophic balance is not applied!', 2I3)
      END IF
    END IF

    ! Get the weak Geostrophic balance weights:
    this%weight = 0.0D0
    IF (this%geobalBeg .GT. 0 .AND. this%geobalEnd .GE. this%geobalBeg) THEN
      IF (yaml_get_var(configFile, 'modelState', 'geobalWeights', weights) /= 0) STOP
      IF (istatus .EQ. 0 .AND. X%sg%gLevel .LE. this%geobalEnd .AND. &
          X%sg%gLevel .GE. this%geobalBeg) THEN
        this%weight = weights(X%sg%gLevel - this%geobalBeg + 1)
        WRITE (*, 6) this%weight, X%sg%gLevel
6       FORMAT('GeosBal initial weak constraint weight: ', D12.4, ' at G:', I2)
      END IF
      DEALLOCATE (weights)
    END IF
  END SUBROUTINE initial

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(GeosBal_t), INTENT(INOUT) :: this

  END SUBROUTINE destructor

  ! Nonlinear forward operator:
  SUBROUTINE transFwdNonLinear(this, X)
    IMPLICIT NONE
    CLASS(GeosBal_t) :: this
    TYPE(State_t), INTENT(INOUT) :: X

    ! Local variables:
    INTEGER(i_kind) :: i, j, k, iedge, istcl !,ir,it,ip,iu,iv,ipr
    REAL(r_kind) :: rf
    REAL(r_kind), DIMENSION(this%X%sg%tSlots)  :: tx, ty, rx, ry, px, py
    REAL(r_kind), DIMENSION(this%X%sg%vLevel, this%X%sg%tSlots) :: tz, rz, pz
    REAL(r_kind) :: rnorm(this%X%sg%tSlots, UBOUND(X%sg%edge_stcl, 2)), &  ! For each edge of a cell
                    rtang(this%X%sg%tSlots, UBOUND(X%sg%edge_stcl, 2)), &  ! For each edge of a cell
                    tnorm(this%X%sg%tSlots, UBOUND(X%sg%edge_stcl, 2)), &  ! For each edge of a cell
                    ttang(this%X%sg%tSlots, UBOUND(X%sg%edge_stcl, 2)), &  ! For each edge of a cell
                    pnorm(this%X%sg%tSlots, UBOUND(X%sg%edge_stcl, 2)), &  ! For each edge of a cell
                    ptang(this%X%sg%tSlots, UBOUND(X%sg%edge_stcl, 2))
    TYPE(State_t) :: OX ! Save the output variable temporarily

    ! Notice that this routine does not use X at all

    ! Initialize output:
    OX = OX%zeroCopy() ! initialize OX as zeros

    ASSOCIATE (temp => X%Fields(this%it)%DATA, &
               uwnd => OX%Fields(this%iu)%DATA, &
               vwnd => OX%Fields(this%iv)%DATA, &
               prro => X%Fields(this%ipr)%DATA, &
               sg => X%sg, coef => X%sg%coef_sigma)

      BLOCK
        ! Applied to the interior vertical levels as terrain following involves derivatives in Z
        DO i = 2, sg%vLevel - 1
          DO j = 1, sg%num_icell
            ! Coriolis:
            rf = dry_air_gas_const / (2.0D0 * Omega * DSIN(X%sg%cell_cntr(1, j)))

            ! Intermediate variable initial:
            tx = 0.0D0; rx = 0.0D0; ty = 0.0D0; ry = 0.0D0; tz = 0.0D0; rz = 0.0D0

            ! Calculate the norm and tang derivatives:
            rnorm = 0.0D0
            rtang = 0.0D0
            tnorm = 0.0D0
            ttang = 0.0D0
            pnorm = 0.0D0
            ptang = 0.0D0
            DO iedge = 1, UBOUND(sg%edge_stcl, 2)
              DO istcl = 1, UBOUND(sg%edge_stcl, 1)
                IF (sg%edge_stcl(istcl, iedge, j) .GT. 0) THEN
                  rnorm(:, iedge) = rnorm(:, iedge) + &
                                    sg%coef_norm(istcl, 1, iedge, j) * prro(i, sg%edge_stcl(istcl, iedge, j), :)
                  rtang(:, iedge) = rtang(:, iedge) + &
                                    sg%coef_tang(istcl, 1, iedge, j) * prro(i, sg%edge_stcl(istcl, iedge, j), :)
                  tnorm(:, iedge) = tnorm(:, iedge) + &
                                    sg%coef_norm(istcl, 1, iedge, j) * temp(i, sg%edge_stcl(istcl, iedge, j), :)
                  ttang(:, iedge) = ttang(:, iedge) + &
                                    sg%coef_tang(istcl, 1, iedge, j) * temp(i, sg%edge_stcl(istcl, iedge, j), :)
                END IF
              END DO
              rx = rx + sg%edgeNorm2(1, iedge, j) * rnorm(:, iedge) + &
                   sg%edgeTang2(1, iedge, j) * rtang(:, iedge)
              ry = ry + sg%edgeNorm2(2, iedge, j) * rnorm(:, iedge) + &
                   sg%edgeTang2(2, iedge, j) * rtang(:, iedge)
              tx = tx + sg%edgeNorm2(1, iedge, j) * tnorm(:, iedge) + &
                   sg%edgeTang2(1, iedge, j) * ttang(:, iedge)
              ty = ty + sg%edgeNorm2(2, iedge, j) * tnorm(:, iedge) + &
                   sg%edgeTang2(2, iedge, j) * ttang(:, iedge)

            END DO
            rx = rx / DBLE(UBOUND(sg%edge_stcl, 2))
            ry = ry / DBLE(UBOUND(sg%edge_stcl, 2))
            tx = tx / DBLE(UBOUND(sg%edge_stcl, 2))
            ty = ty / DBLE(UBOUND(sg%edge_stcl, 2))

            rz(i, :) = prro(i - 1, j, :) * coef(1, i) + prro(i, j, :) * coef(2, i) + prro(i + 1, j, :) * coef(3, i)
            tz(i, :) = temp(i - 1, j, :) * coef(1, i) + temp(i, j, :) * coef(2, i) + temp(i + 1, j, :) * coef(3, i)

            ! Always choose pressure_ctl if it is available
            IF (this%ip .GT. 0) THEN
              uwnd(i, j, :) = -rf * temp(i, j, :) * (ry + rz(i, :) * sg%sigmay(i, j))
              vwnd(i, j, :) = rf * temp(i, j, :) * (rx + rz(i, :) * sg%sigmax(i, j))
            ELSE ! ir > 0
              uwnd(i, j, :) = -rf * (temp(i, j, :) / prro(i, j, :) * (ry + rz(i, :) * sg%sigmay(i, j)) - &
                                     (ty + tz(i, :) * sg%sigmay(i, j)))
              vwnd(i, j, :) = rf * (temp(i, j, :) / prro(i, j, :) * (rx + rz(i, :) * sg%sigmax(i, j)) + &
                                    (tx + tz(i, :) * sg%sigmax(i, j)))
            END IF
          END DO
        END DO
      END BLOCK

    END ASSOCIATE

    X = OX

    CALL X%exHalo()
  END SUBROUTINE transFwdNonLinear

  FUNCTION transFwdNonLinear_opr(this, X) RESULT(X1)
    IMPLICIT NONE
    CLASS(GeosBal_t) :: this
    TYPE(State_t), INTENT(IN) :: X
    TYPE(State_t) :: X1
    INTEGER(i_kind) :: i, j, k

    X1 = X%zeroCopy()
    IF (this%iu .NE. 0 .AND. this%iv .NE. 0) THEN
      X1%fields(this%ipr) = X%fields(this%ipr)
      X1%fields(this%it) = X%fields(this%it)
      CALL this%transFwdNonLinear(X1)
      X1%fields(this%ipr)%DATA = 0.0D0
      X1%fields(this%it)%DATA = 0.0D0
    END IF
  END FUNCTION transFwdNonLinear_opr

  !> @brief by Yuanfu Xie 2022-07-27
  !! This routine contains the GeosBal Object which Jacobian multiples
  !! increment dX from temp/rho to uwnd and vwnd + function value at
  !! at a given point X if it is present; otherwise, simply J*dX
  !!    H*dx: d(temp, rho) -> d(uwnd,vwnd)
  !! Since the geostrophic balance is nonlinear, it cannot directly
  !! call its fully nonliear operator to replace this routine.
  SUBROUTINE transFwdTanLinear(this, dX, X)

!#define TRACE_DERIVATIVES

    IMPLICIT NONE
    CLASS(GeosBal_t) :: this
    TYPE(State_t), INTENT(INOUT) :: dX
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X

    ! Local variables:
    INTEGER(i_kind) :: i, j, k, iedge, istcl !,ip,ir,it,iu,iv,ipr
    REAL(r_kind) :: rf
    ! All intermediate variables should be the same dimensions as the controls:
    REAL(r_kind), ALLOCATABLE :: rx(:, :, :), ry(:, :, :), rz(:, :, :), &
                                 drx(:, :, :), dry(:, :, :), drz(:, :, :), dtx(:, :, :), dty(:, :, :), dtz(:, :, :)
    REAL(r_kind), ALLOCATABLE :: &
      xrnorm(:, :, :), xrtang(:, :, :), xtnorm(:, :, :), xttang(:, :, :), & ! For each edge of a cell
      drnorm(:, :, :), drtang(:, :, :), dtnorm(:, :, :), dttang(:, :, :)    ! For each edge of a cell
    TYPE(State_t) :: OX ! Save the output variable temporarily

    ! Debugging: 2022-09-11:
    INTEGER(i_kind) :: ilvl, icel, igrd
    REAL(r_kind) :: totalCoef, totaln, totalt

    ilvl = 2; icel = 7; igrd = 12 ! Check vertical grid and horizontal grid

    ! Notice that this routine requires X:
    IF (.NOT. PRESENT(X)) THEN
      WRITE (*, 1)
1     FORMAT('GeoBal - transFwdTanLinear: GeoBal is a nonlinear constraint and', /, &
             'it requires X is present. Check and rerun!')
      STOP
    END IF

    ! Initialize output:
    OX = dX%zeroCopy() ! initialize OX as zeros

    ASSOCIATE (temp => X%Fields(this%it)%DATA, &
               dtem => dX%Fields(this%it)%DATA, &
               uwnd => OX%Fields(this%iu)%DATA, &
               vwnd => OX%Fields(this%iv)%DATA, &
               prro => X%Fields(this%ipr)%DATA, &
               dpro => dX%Fields(this%ipr)%DATA, &
               sg => X%sg, coef => X%sg%coef_sigma)

      uwnd = 0.0D0
      vwnd = 0.0D0

      ! Allocate memory for the intermediate variables:
      ALLOCATE (rx(sg%vLevel, sg%num_cell, sg%tSlots), ry(sg%vLevel, sg%num_cell, sg%tSlots), &
                rz(sg%vLevel, sg%num_cell, sg%tSlots), &
                drx(sg%vLevel, sg%num_cell, sg%tSlots), dry(sg%vLevel, sg%num_cell, sg%tSlots), &
                drz(sg%vLevel, sg%num_cell, sg%tSlots), dtx(sg%vLevel, sg%num_cell, sg%tSlots), &
                dty(sg%vLevel, sg%num_cell, sg%tSlots), dtz(sg%vLevel, sg%num_cell, sg%tSlots))
      ALLOCATE (xrnorm(sg%vLevel, sg%num_cell, sg%tSlots), xrtang(sg%vLevel, sg%num_cell, sg%tSlots), &
                xtnorm(sg%vLevel, sg%num_cell, sg%tSlots), xttang(sg%vLevel, sg%num_cell, sg%tSlots), &
                drnorm(sg%vLevel, sg%num_cell, sg%tSlots), drtang(sg%vLevel, sg%num_cell, sg%tSlots), &
                dtnorm(sg%vLevel, sg%num_cell, sg%tSlots), dttang(sg%vLevel, sg%num_cell, sg%tSlots))

      BLOCK

        ! Debugging: 2022-09-12:
        totalCoef = 0.0D0 ! total coef multiplied to a given grid function

        DO j = 1, sg%num_icell

          ! Apply this constraint on interior points only:
          IF (sg%cell_type(j) .NE. 0) CYCLE

          ! Coriolis:
          rf = dry_air_gas_const / (2.0D0 * Omega * DSIN(dX%sg%cell_cntr(1, j)))

          ! Derivatives for each cell:
          rx = 0.0D0; ry = 0.0D0; rz = 0.0D0
          drx = 0.0D0; dry = 0.0D0; drz = 0.0D0
          dtx = 0.0D0; dty = 0.0D0; dtz = 0.0D0

          DO iedge = 1, UBOUND(sg%edge_stcl, 2)
            ! Calculate the norm and tang derivatives for each edge:
            xrnorm = 0.0D0; xrtang = 0.0D0
            xtnorm = 0.0D0; xttang = 0.0D0
            drnorm = 0.0D0; drtang = 0.0D0
            dtnorm = 0.0D0; dttang = 0.0D0

            ! Totaln totalt initial:
            totaln = 0.0D0; totalt = 0.0D0

            DO istcl = 1, UBOUND(sg%edge_stcl, 1)

              IF (sg%edge_stcl(istcl, iedge, j) .GT. 0) THEN
                xrnorm(:, j, :) = xrnorm(:, j, :) + &
                                  sg%coef_norm(istcl, 1, iedge, j) * prro(:, sg%edge_stcl(istcl, iedge, j), :)
                xrtang(:, j, :) = xrtang(:, j, :) + &
                                  sg%coef_tang(istcl, 1, iedge, j) * prro(:, sg%edge_stcl(istcl, iedge, j), :)
                xtnorm(:, j, :) = xtnorm(:, j, :) + &
                                  sg%coef_norm(istcl, 1, iedge, j) * temp(:, sg%edge_stcl(istcl, iedge, j), :)
                xttang(:, j, :) = xttang(:, j, :) + &
                                  sg%coef_tang(istcl, 1, iedge, j) * temp(:, sg%edge_stcl(istcl, iedge, j), :)
                drnorm(:, j, :) = drnorm(:, j, :) + &
                                  sg%coef_norm(istcl, 1, iedge, j) * dpro(:, sg%edge_stcl(istcl, iedge, j), :)
                drtang(:, j, :) = drtang(:, j, :) + &
                                  sg%coef_tang(istcl, 1, iedge, j) * dpro(:, sg%edge_stcl(istcl, iedge, j), :)
                dtnorm(:, j, :) = dtnorm(:, j, :) + &
                                  sg%coef_norm(istcl, 1, iedge, j) * dtem(:, sg%edge_stcl(istcl, iedge, j), :)
                dttang(:, j, :) = dttang(:, j, :) + &
                                  sg%coef_tang(istcl, 1, iedge, j) * dtem(:, sg%edge_stcl(istcl, iedge, j), :)

#ifdef TRACE_DERIVATIVES
                ! Debugging 2022-09-08
                IF (j .EQ. icel .AND. sg%edge_stcl(istcl, iedge, j) .EQ. igrd) THEN
                  totaln = totaln + sg%coef_norm(istcl, 1, iedge, j)
                  totalt = totalt + sg%coef_tang(istcl, 1, iedge, j)
                  WRITE (*, 32) j, istcl, iedge, totaln, totalt, &
                    sg%coef_norm(istcl, 1, iedge, j), sg%coef_tang(istcl, 1, iedge, j), &
                    dtem(ilvl, sg%edge_stcl(istcl, iedge, j), 1), &
                    dX%mpddGlob%myrank
32                FORMAT('CCC: ', I2, ' se', 2I2, ' totalN/T', 2D13.5, ' coe', 2D13.5, ' t', D13.5, ' pc', I2)
                  !call sg%mpddInfo_sg%barrier
                END IF
                ! Debugging 2022-09-08
                IF (j .EQ. icel .OR. j .EQ. icel + 5) THEN
                  WRITE (*, 31) j, iedge, UBOUND(sg%edge_stcl, 1), &
                    sg%edge_stcl(istcl, iedge, j), & ! dtem(i,sg%edge_stcl(istcl,iedge,j),1), &
                    xrnorm(ilvl, j, 1), xrtang(ilvl, j, 1), prro(ilvl, sg%edge_stcl(istcl, iedge, j), 1), sg%coef_tang(istcl, 1, iedge, j), &
                    sg%edgeNorm2(:, iedge, j), &
                    sg%edgeTang2(:, iedge, j), dX%mpddGlob%myrank
31                FORMAT('NNTinfo: ', I2, ' edge', I2, ' Ns/stcl', I2, I3, ' dnT', 4D13.5, ' NV', 2F5.1, ' TV', 2F5.1, ' pc', I2)
                  WRITE (*, 36) j, iedge, sg%edge_stcl(istcl, iedge, j), &
                    sg%coef_norm(istcl, 1, iedge, j), sg%coef_tang(istcl, 1, iedge, j), &
                    dtem(ilvl, sg%edge_stcl(istcl, iedge, j), 1), &
                    sg%mpddInfo_sg%myrank
36                FORMAT('NNTinfo 1:', I2, ' edge', I2, ' stcl', I3, ' coef', 2D12.4, ' temp', D12.4, ' pc', I2)
                END IF
#endif

              END IF
            END DO

            rx = rx + sg%edgeNorm2(1, iedge, j) * xrnorm + &
                 sg%edgeTang2(1, iedge, j) * xrtang
            ry = ry + sg%edgeNorm2(2, iedge, j) * xrnorm + &
                 sg%edgeTang2(2, iedge, j) * xrtang
            drx = drx + sg%edgeNorm2(1, iedge, j) * drnorm + &
                  sg%edgeTang2(1, iedge, j) * drtang
            dry = dry + sg%edgeNorm2(2, iedge, j) * drnorm + &
                  sg%edgeTang2(2, iedge, j) * drtang
            dtx = dtx + sg%edgeNorm2(1, iedge, j) * dtnorm + &
                  sg%edgeTang2(1, iedge, j) * dttang
            dty = dty + sg%edgeNorm2(2, iedge, j) * dtnorm + &
                  sg%edgeTang2(2, iedge, j) * dttang

#ifdef TRACE_DERIVATIVES
            ! Debugging: 2022-09-12:
            IF (j .EQ. icel) THEN
              PRINT *, 'TotalCoef before: ', totalCoef, sg%mpddInfo_sg%myrank
              totalCoef = totalCoef + totaln * sg%edgeNorm2(2, iedge, j) + totalt * sg%edgeTang2(2, iedge, j)
              WRITE (*, 44) i, j, iedge, totalCoef, totaln, totalt, sg%edgeNorm2(2, iedge, j), sg%edgeTang2(2, iedge, j), &
                sg%mpddInfo_sg%myrank
44            FORMAT('Total project', 3I3, ' total', 3D12.4, ' edgeNT', 2F5.1, ' pc', I2)
            END IF
#endif

          END DO  ! End of edge loop

          ! Average:
          rx = rx / DBLE(UBOUND(sg%edge_stcl, 2))
          ry = ry / DBLE(UBOUND(sg%edge_stcl, 2))
          drx = drx / DBLE(UBOUND(sg%edge_stcl, 2))
          dry = dry / DBLE(UBOUND(sg%edge_stcl, 2))
          dtx = dtx / DBLE(UBOUND(sg%edge_stcl, 2))
          dty = dty / DBLE(UBOUND(sg%edge_stcl, 2))

          ! Applied to the interior vertical levels as terrain following involves derivatives in Z
          DO i = 2, sg%vLevel - 1

            ! Arbitrary vertical grid:
            rz(i, j, :) = prro(i - 1, j, :) * coef(1, i) + prro(i, j, :) * coef(2, i) + prro(i + 1, j, :) * coef(3, i)
            drz(i, j, :) = dpro(i - 1, j, :) * coef(1, i) + dpro(i, j, :) * coef(2, i) + dpro(i + 1, j, :) * coef(3, i)
            dtz(i, j, :) = dtem(i - 1, j, :) * coef(1, i) + dtem(i, j, :) * coef(2, i) + dtem(i + 1, j, :) * coef(3, i)

            ! (u,v) = (bar{u},bar{v}) + J * (du, dv):
            IF (this%ip .GT. 0) THEN ! Always choose pressure_ctl if it is available:
              uwnd(i, j, :) = -rf * ((ry(i, j, :) + rz(i, j, :) * sg%sigmay(i, j)) * dtem(i, j, :) + &
                                     temp(i, j, :) * (dry(i, j, :) + drz(i, j, :) * sg%sigmay(i, j)))
              vwnd(i, j, :) = rf * ((rx(i, j, :) + rz(i, j, :) * sg%sigmax(i, j)) * dtem(i, j, :) + &
                                    temp(i, j, :) * (drx(i, j, :) + drz(i, j, :) * sg%sigmax(i, j)))
            ELSE
              uwnd(i, j, :) = &
                -rf * ((dtem(i, j, :) / prro(i, j, :) - temp(i, j, :) / prro(i, j, :)**2 * dpro(i, j, :)) * (ry(i, j, :) + rz(i, j, :) * sg%sigmay(i, j)) &
                       + (dty(i, j, :) + dtz(i, j, :) * sg%sigmay(i, j)) + temp(i, j, :) / prro(i, j, :) * (dry(i, j, :) + drz(i, j, :) * sg%sigmay(i, j)))
              vwnd(i, j, :) = &
                rf * ((dtem(i, j, :) / prro(i, j, :) - temp(i, j, :) / prro(i, j, :)**2 * dpro(i, j, :)) * (rx(i, j, :) + rz(i, j, :) * sg%sigmax(i, j)) &
                      + (dtx(i, j, :) + dtz(i, j, :) * sg%sigmax(i, j)) + temp(i, j, :) / prro(i, j, :) * (drx(i, j, :) + drz(i, j, :) * sg%sigmax(i, j)))
            END IF

#ifdef TRACE_DERIVATIVES
            ! Debugging 2022-09-08:
            IF ((j .EQ. icel .OR. j .EQ. icel + 5) .AND. i .EQ. ilvl) THEN
              !totalCoef = -rf*totalCoef/DBLE(UBOUND(sg%edge_stcl,2))
              WRITE (*, 24) uwnd(i, j, 2), (ry(i, j, 1) + rz(i, j, 1) * sg%sigmay(j)), (dty(i, j, 1) + dtz(i, j, 1) * sg%sigmay(j)), rf * dty(i, j, 1), &
                (dry(i, j, 1) + drz(i, j, 1) * sg%sigmay(j)), -rf * totalCoef / DBLE(UBOUND(sg%edge_stcl, 2)), temp(i, j, 1), sg%mpddInfo_sg%myrank
24            FORMAT('Uwnd', 5D12.4, ' u totalCoef:', D12.4, ' t:', D12.4, ' pc', I2)
              WRITE (*, 23) vwnd(i, j, 2), rx(i, j, 1), sg%sigmax(j), dtx(i, j, 1), rf, &
                this%ip, this%ir, dtem(i, j, 1), sg%mpddInfo_sg%myrank
              ! write(*,23) vwnd(i,j,2),(rx(i,j,1) + rz(i,j,1)*sg%sigmax(j)),(dtx(i,j,1)+dtz(i,j,1)*sg%sigmax(j)),rf*dtx(i,j,1), &
              !   (drx(i,j,1)+drz(i,j,1)*sg%sigmax(j)),this%ip,this%ir,dtem(i,j,1),sg%mpddInfo_sg%myrank
23            FORMAT('Vwnd', 5D12.4, ' ctl idx ip:', I2, ' ir:', I2, ' dt:', D12.4, ' pc', I2)
              WRITE (*, 22) sg%mpddInfo_sg%myrank, sg%cell_type
22            FORMAT('cellTypes: pc', I2, ' types: ', 50I2)
              WRITE (*, 21) i, j, uwnd(i, j, 1), OX%fields(dX%getVarIdx(TRIM('uwnd')))%DATA(ilvl, icel, 1), sg%cell_type(j), rf, dty(i, j, 1), dX%mpddGlob%myrank
21            FORMAT('DTY/DTZ at v/h: ', I2, I4, ' uv', 2D12.5, ' type ', I1, ' val', 2D13.6, ' pc', I2)
              !call dX%mpddSub%barrier
              !call flush
              !stop
            END IF
#endif

          END DO
        END DO
      END BLOCK

    END ASSOCIATE

    DEALLOCATE (rx, ry, rz, drx, dry, drz, dtx, dty, dtz)
    DEALLOCATE (xrnorm, xrtang, xtnorm, xttang, drnorm, drtang, dtnorm, dttang)

    !dX = dX%zeroCopy()
    dX = OX

    CALL dX%exHalo()

#ifdef TRACE_DERIVATIVES
    ! Debugging 2022-09-08:
    WRITE (*, 111) this%iu, OX%getVarIdx(TRIM('uwnd')), dX%fields(dX%getVarIdx(TRIM('uwnd')))%DATA(ilvl, icel, 1), &
      OX%fields(OX%getVarIdx(TRIM('uwnd')))%DATA(ilvl, icel + 1, 1), dX%mpddSub%myrank
111 FORMAT('dX OX FWD after halo exchange: ', 2I2, 2D12.4, ' pc', I2)
#endif
  END SUBROUTINE transFwdTanLinear

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION transFwdTanLinear_opr(this, dX, X) RESULT(X1)
    IMPLICIT NONE
    CLASS(GeosBal_t) :: this
    TYPE(State_t), INTENT(IN) :: dX
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X
    TYPE(State_t) :: X1
    INTEGER(i_kind) :: i, j, k

    X1 = X%zeroCopy()
    IF (this%it .NE. 0 .AND. this%ipr .NE. 0) THEN
      X1%fields(this%it) = dX%fields(this%it)
      X1%fields(this%ipr) = dX%fields(this%ipr)
      CALL this%transFwdTanLinear(X1, X)
      X1%fields(this%it)%DATA = 0.0D0
      X1%fields(this%ipr)%DATA = 0.0D0
    END IF
  END FUNCTION transFwdTanLinear_opr

  !>
  !! @brief:
  !! This routine calculates the adjoint operator transposely multiply a state
  !! in dX. For a nonlinear operator, the state of X is required:
  !!
  !! The geostrophic balance: F(temp,rho) = (u,v) and P = rho R T
  !! u = -P_y/f/rho
  !! v =  P_x/f/rho
  !
  SUBROUTINE transAdjMultiply(this, dX, X)

!#define TRACE_DERIVATIVES

    IMPLICIT NONE
    CLASS(GeosBal_t) :: this
    TYPE(State_t), INTENT(INOUT) :: dX
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X

    ! Local variables:
    INTEGER(i_kind) :: i, j, k, iedge, istcl !,ip,ir,it,iu,iv,ipr
    REAL(r_kind), DIMENSION(dX%sg%vLevel, this%X%sg%tSlots)  :: rx, ry, rz
    REAL(r_kind) :: rf
    REAL(r_kind) :: rnorm(dX%sg%vLevel, dX%sg%tSlots), &  ! For each edge of a cell
                    rtang(dX%sg%vLevel, dX%sg%tSlots), &  ! For each edge of a cell
                    tnorm(dX%sg%vLevel, dX%sg%tSlots), &  ! For each edge of a cell
                    ttang(dX%sg%vLevel, dX%sg%tSlots)     ! For each edge of a cell
    REAL(r_kind), ALLOCATABLE :: &
      drx(:), dry(:), dtx(:), dty(:), &
      drnorm(:), drtang(:), & ! For each edge of a cell, save the interpolation coeffs D(*)
      dtnorm(:), dttang(:)    ! For each edge of a cell
    TYPE(State_t) :: OX

    ! Debugging: 2022-09-11:
    INTEGER(i_kind) :: ilvl, icel, igrd

    ilvl = 1; icel = 7; igrd = 12 ! Check vertical grid and horizontal grid

    IF (.NOT. PRESENT(X)) THEN
      WRITE (*, 1)
1     FORMAT('GeoBal - transAdjMultiply: GeoBal is a nonlinear constraint and', /, &
             'it requires X is present. Check and rerun!')
      STOP
    END IF

    ! Initialize output:
    OX = dX%zeroCopy()

    ASSOCIATE (temp => X%Fields(this%it)%DATA, &
               dtem => OX%Fields(this%it)%DATA, &
               uwnd => dX%Fields(this%iu)%DATA, &
               vwnd => dX%Fields(this%iv)%DATA, &
               prro => X%Fields(this%ipr)%DATA, &
               dpro => OX%Fields(this%ipr)%DATA, &
               sg => X%sg, coef => X%sg%coef_sigma)

      ! Initialize output:
      dtem = 0.0D0
      dpro = 0.0D0

#ifdef TRACE_DERIVATIVES
      ! Debugging 2022-09-11:
      WRITE (*, 3) uwnd(ilvl, icel, 1), vwnd(ilvl, icel, 1), sg%mpddInfo_sg%myrank
3     FORMAT('delta-u: ', 2D12.4, ' pc', I2)
#endif

      ! The adjoint will map du dv values in 1:num_icell to dt dr in 1:num_cell:
      ALLOCATE (drnorm(sg%num_cell), drtang(sg%num_cell), &
                dtnorm(sg%num_cell), dttang(sg%num_cell), &
                drx(sg%num_cell), dry(sg%num_cell), &
                dtx(sg%num_cell), dty(sg%num_cell))

      BLOCK

        ! Applied to the interior vertical levels as terrain following involves derivatives in Z
        DO j = 1, sg%num_icell

          ! Apply this constraint on interior points only:
          IF (sg%cell_type(j) .NE. 0) CYCLE

          ! Coriolis:
          rf = dry_air_gas_const / (2.0D0 * Omega * DSIN(sg%cell_cntr(1, j)))

          rx = 0.0D0; ry = 0.0D0; rz = 0.0D0; 
          drx = 0.0D0; dry = 0.0D0
          dtx = 0.0D0; dty = 0.0D0

          DO iedge = 1, UBOUND(sg%edge_stcl, 2)

            ! Calculate the norm and tang derivatives for each edge:
            rnorm = 0.0D0; rtang = 0.0D0
            tnorm = 0.0D0; ttang = 0.0D0

            ! Initialize the adjoint for each edge:
            drnorm = 0.0D0; drtang = 0.0D0
            dtnorm = 0.0D0; dttang = 0.0D0

            DO istcl = 1, UBOUND(sg%edge_stcl, 1)

              IF (sg%edge_stcl(istcl, iedge, j) .GT. 0) THEN
                rnorm = rnorm + &
                        sg%coef_norm(istcl, 1, iedge, j) * prro(:, sg%edge_stcl(istcl, iedge, j), :)
                rtang = rtang + &
                        sg%coef_tang(istcl, 1, iedge, j) * prro(:, sg%edge_stcl(istcl, iedge, j), :)
                tnorm = tnorm + &
                        sg%coef_norm(istcl, 1, iedge, j) * temp(:, sg%edge_stcl(istcl, iedge, j), :)
                ttang = ttang + &
                        sg%coef_tang(istcl, 1, iedge, j) * temp(:, sg%edge_stcl(istcl, iedge, j), :)

                ! D(*) coefficients in forward operator: accumulate them for the adjoint
                drnorm(sg%edge_stcl(istcl, iedge, j)) = drnorm(sg%edge_stcl(istcl, iedge, j)) + &
                                                        sg%coef_norm(istcl, 1, iedge, j)
                drtang(sg%edge_stcl(istcl, iedge, j)) = drtang(sg%edge_stcl(istcl, iedge, j)) + &
                                                        sg%coef_tang(istcl, 1, iedge, j)
                dtnorm(sg%edge_stcl(istcl, iedge, j)) = dtnorm(sg%edge_stcl(istcl, iedge, j)) + &
                                                        sg%coef_norm(istcl, 1, iedge, j)
                dttang(sg%edge_stcl(istcl, iedge, j)) = dttang(sg%edge_stcl(istcl, iedge, j)) + &
                                                        sg%coef_tang(istcl, 1, iedge, j)

#ifdef TRACE_DERIVATIVES
                ! Debugging 2022-09-08
                IF (j .EQ. icel) THEN
                  WRITE (*, 31) j, iedge, UBOUND(sg%edge_stcl, 1), &
                    sg%edge_stcl(istcl, iedge, j), & !temp(i,sg%edge_stcl(istcl,iedge,j),1),&
                    dtnorm(sg%edge_stcl(istcl, iedge, j)), &
                    dttang(sg%edge_stcl(istcl, iedge, j)), &
                    sg%edgeNorm2(:, iedge, j), &
                    sg%edgeTang2(:, iedge, j), dX%mpddGlob%myrank
31                FORMAT('AD NTinfo:', I2, ' edge', I2, ' Ns/stcl', I2, I3, ' dnT', 2D13.5, ' NV', 2F5.1, ' TV', 2F5.1, ' pc', I2)
                  WRITE (*, 36) j, iedge, sg%edge_stcl(istcl, iedge, j), &
                    sg%coef_norm(istcl, 1, iedge, j), sg%coef_tang(istcl, 1, iedge, j), &
                    temp(ilvl, sg%edge_stcl(istcl, iedge, j), 1), &
                    sg%mpddInfo_sg%myrank
36                FORMAT('AD NTinfo 1:', I2, ' edge', I2, ' stcl', I3, ' coef', 2D12.4, ' temp', D12.4, ' pc', I2)
                END IF
#endif
              END IF
            END DO

            drnorm = drnorm / DBLE(UBOUND(sg%edge_stcl, 2))
            drtang = drtang / DBLE(UBOUND(sg%edge_stcl, 2))
            dtnorm = dtnorm / DBLE(UBOUND(sg%edge_stcl, 2))
            dttang = dttang / DBLE(UBOUND(sg%edge_stcl, 2))

            rx = rx + sg%edgeNorm2(1, iedge, j) * rnorm + &
                 sg%edgeTang2(1, iedge, j) * rtang
            ry = ry + sg%edgeNorm2(2, iedge, j) * rnorm + &
                 sg%edgeTang2(2, iedge, j) * rtang

            drx = drx + sg%edgeNorm2(1, iedge, j) * drnorm + &
                  sg%edgeTang2(1, iedge, j) * drtang
            dry = dry + sg%edgeNorm2(2, iedge, j) * drnorm + &
                  sg%edgeTang2(2, iedge, j) * drtang
            dtx = dtx + sg%edgeNorm2(1, iedge, j) * dtnorm + &
                  sg%edgeTang2(1, iedge, j) * dttang
            dty = dty + sg%edgeNorm2(2, iedge, j) * dtnorm + &
                  sg%edgeTang2(2, iedge, j) * dttang

#ifdef TRACE_DERIVATIVES
            ! Debugging 2022-9-13:
            IF (j .EQ. icel) &
              WRITE (*, 555) sg%edgeNorm2(2, iedge, j), sg%edgeTang2(2, iedge, j), &
              dtnorm(2), dttang(2), sg%mpddInfo_sg%myrank
555         FORMAT('DTY-', 2D12.4, ' NT2:', 2D12.4, ' pc', I2)
#endif

          END DO  ! End of edge loop

          rx = rx / DBLE(UBOUND(sg%edge_stcl, 2))
          ry = ry / DBLE(UBOUND(sg%edge_stcl, 2))

          ! Horizontal derivative terms: D_x and D_y
          DO i = 2, sg%vLevel - 1
            DO k = 1, sg%tSlots
              IF (this%ip .LE. 0) THEN
                dtem(i, :, k) = dtem(i, :, k) - rf * dty * uwnd(i, j, k) + rf * dtx * vwnd(i, j, k)
                dpro(i, :, k) = dpro(i, :, k) - rf * dry * temp(i, j, k) / prro(i, j, k) * uwnd(i, j, k) + &
                                rf * drx * temp(i, j, k) / prro(i, j, k) * vwnd(i, j, k)
              ELSE ! Note: no horizontal terms of delta u/v in dtem
                dpro(i, :, k) = dpro(i, :, k) - rf * dry * temp(i, j, k) * uwnd(i, j, k) + &
                                rf * drx * temp(i, j, k) * vwnd(i, j, k)
              END IF
            END DO
          END DO

#ifdef TRACE_DERIVATIVES
          ! Debugging:  2022-09-11:
          IF (j .EQ. icel) THEN
            WRITE (*, 45) j, dtem(ilvl, igrd, 1), dpro(ilvl, igrd, 1), rf, dty(igrd), uwnd(ilvl, j, 1), sg%mpddInfo_sg%myrank
45          FORMAT('Adjoint T:', I3, ' dt:', 2D12.4, ' r/f', D12.4, ' dty/uwnd', 2D12.4, ' pc', I2)
          END IF
#endif

          ! Calculate the vertical adjoint:
          DO i = 2, sg%vLevel - 1

            rz(i, :) = prro(i - 1, j, :) * coef(1, i) + prro(i, j, :) * coef(2, i) + prro(i + 1, j, :) * coef(3, i)

            ! For temp prro controls:
            IF (this%ip .LE. 0) THEN
              ! No perturbation cross terms: diagonal
              dtem(i, j, :) = dtem(i, j, :) &
                              - rf / prro(i, j, :) * (ry(i, :) + rz(i, :) * sg%sigmay(i, j)) * uwnd(i, j, :) &
                              + rf / prro(i, j, :) * (rx(i, :) + rz(i, :) * sg%sigmax(i, j)) * vwnd(i, j, :)
              dpro(i, j, :) = dpro(i, j, :) &
                              + rf * temp(i, j, :) / prro(i, j, :)**2 * (ry(i, :) + rz(i, :) * sg%sigmay(i, j)) * uwnd(i, j, :) &
                              - rf * temp(i, j, :) / prro(i, j, :)**2 * (rx(i, :) + rz(i, :) * sg%sigmax(i, j)) * vwnd(i, j, :)

              ! Verical derivative terms: D_sigma
              dtem(i + 1, j, :) = dtem(i + 1, j, :) - rf * coef(3, i) * sg%sigmay(i, j) * uwnd(i, j, :)
              dtem(i, j, :) = dtem(i, j, :) - rf * coef(2, i) * sg%sigmay(i, j) * uwnd(i, j, :)
              dtem(i - 1, j, :) = dtem(i - 1, j, :) - rf * coef(1, i) * sg%sigmay(i, j) * uwnd(i, j, :)
              dtem(i + 1, j, :) = dtem(i + 1, j, :) + rf * coef(3, i) * sg%sigmax(i, j) * vwnd(i, j, :)
              dtem(i, j, :) = dtem(i, j, :) + rf * coef(2, i) * sg%sigmax(i, j) * vwnd(i, j, :)
              dtem(i - 1, j, :) = dtem(i - 1, j, :) + rf * coef(1, i) * sg%sigmax(i, j) * vwnd(i, j, :)

              dpro(i + 1, j, :) = dpro(i + 1, j, :) - rf * temp(i, j, :) / prro(i, j, :) * coef(3, i) * sg%sigmay(i, j) * uwnd(i, j, :)
              dpro(i, j, :) = dpro(i, j, :) - rf * temp(i, j, :) / prro(i, j, :) * coef(2, i) * sg%sigmay(i, j) * uwnd(i, j, :)
              dpro(i - 1, j, :) = dpro(i - 1, j, :) - rf * temp(i, j, :) / prro(i, j, :) * coef(1, i) * sg%sigmay(i, j) * uwnd(i, j, :)
              dpro(i + 1, j, :) = dpro(i + 1, j, :) + rf * temp(i, j, :) / prro(i, j, :) * coef(3, i) * sg%sigmax(i, j) * vwnd(i, j, :)
              dpro(i, j, :) = dpro(i, j, :) + rf * temp(i, j, :) / prro(i, j, :) * coef(2, i) * sg%sigmax(i, j) * vwnd(i, j, :)
              dpro(i - 1, j, :) = dpro(i - 1, j, :) + rf * temp(i, j, :) / prro(i, j, :) * coef(1, i) * sg%sigmax(i, j) * vwnd(i, j, :)
            ELSE ! The following blanks used to compare with the rho control equations:
              dtem(i, j, :) = dtem(i, j, :) &
                              - rf * (ry(i, :) + rz(i, :) * sg%sigmay(i, j)) * uwnd(i, j, :) &
                              + rf * (rx(i, :) + rz(i, :) * sg%sigmax(i, j)) * vwnd(i, j, :)
              ! Verical derivative terms: D_sigma
              dpro(i + 1, j, :) = dpro(i + 1, j, :) - rf * temp(i, j, :) * coef(3, i) * sg%sigmay(i, j) * uwnd(i, j, :)
              dpro(i, j, :) = dpro(i, j, :) - rf * temp(i, j, :) * coef(2, i) * sg%sigmay(i, j) * uwnd(i, j, :)
              dpro(i - 1, j, :) = dpro(i - 1, j, :) - rf * temp(i, j, :) * coef(1, i) * sg%sigmay(i, j) * uwnd(i, j, :)
              dpro(i + 1, j, :) = dpro(i + 1, j, :) + rf * temp(i, j, :) * coef(3, i) * sg%sigmax(i, j) * vwnd(i, j, :)
              dpro(i, j, :) = dpro(i, j, :) + rf * temp(i, j, :) * coef(2, i) * sg%sigmax(i, j) * vwnd(i, j, :)
              dpro(i - 1, j, :) = dpro(i - 1, j, :) + rf * temp(i, j, :) * coef(1, i) * sg%sigmax(i, j) * vwnd(i, j, :)
            END IF

          END DO
        END DO

#ifdef TRACE_DERIVATIVES
        ! Debugging: 2022-09-11:
        WRITE (*, 4) MAXVAL(ABS(dtem)), MAXVAL(ABS(dpro)), MAXVAL(ABS(sg%sigmay)), sg%mpddInfo_sg%myrank
4       FORMAT('Adjoint tem and prro: ', 2D12.4, ' max sigmay', D12.4, ' pc', I2)
#endif

      END BLOCK

      ! Deallocate:
      DEALLOCATE (drnorm, drtang, dtnorm, dttang, drx, dry, dtx, dty)
    END ASSOCIATE

#ifdef TRACE_DERIVATIVES
    ! Debugging: 2022-09-11:
    WRITE (*, 5) (MAXVAL(ABS(ox%fields(i)%DATA)), i=1, 4), dx%sg%mpddInfo_sg%myrank
5   FORMAT('Adjoint OX: ', 4D12.4, ' pc', I2)
#endif

    dX = OX

    CALL dX%exHaloRevSum()

  END SUBROUTINE transAdjMultiply

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Adjoint operator
  FUNCTION transAdjMultiply_opr(this, dX, X) RESULT(dX1)
    IMPLICIT NONE
    CLASS(GeosBal_t) :: this
    TYPE(State_t), INTENT(IN) :: dX
    TYPE(State_t) :: dX1
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X

    IF (.NOT. PRESENT(X)) THEN
      WRITE (*, 1)
1     FORMAT('GeoBal - transAdjMultiply_opr: GeoBal is a nonlinear constraint and', /, &
             'it requires X is present. Check and rerun!')
      STOP
    END IF

    dX1 = dX%zeroCopy()
    IF ((dX%getVarIdx(TRIM('temp')) .NE. 0) .AND. &
        (dX%getVarIdx(TRIM('rho')) .NE. 0) .AND. &
        (dX%getVarIdx(TRIM('uwnd')) .NE. 0) .AND. &
        (dX%getVarIdx(TRIM('vwnd')) .NE. 0)) THEN
      dX1%fields(dX1%getVarIdx('uwnd')) = dX%fields(dX%getVarIdx('uwnd'))
      dX1%fields(dX1%getVarIdx('vwnd')) = dX%fields(dX%getVarIdx('vwnd'))
      CALL this%transAdjMultiply(dX1, X)
    END IF
  END FUNCTION transAdjMultiply_opr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE dLogP(this, X, logPx, logPy)
    USE conversions_m, ONLY: hydroPressure, hydroPressTop
    USE State2NC_m
    IMPLICIT NONE

    CLASS(GeosBal_t) :: this
    TYPE(State_t), INTENT(IN) :: X
    REAL(r_kind), INTENT(OUT) :: &
      logPx(this%ng(1), this%ng(2), this%ng(3)), &
      logPy(this%ng(1), this%ng(2), this%ng(3))

    ! Local variables:
    INTEGER(i_kind) :: iedge, istcl, i, j, k
    REAL(r_kind) :: rf
    REAL(r_kind), ALLOCATABLE :: hydro(:, :, :), hydrT(:, :, :), pnorm(:, :, :), ptang(:, :, :), logPz(:, :, :)

    IF (this%iu .LE. 0 .OR. this%iv .LE. 0 .OR. this%iq .LE. 0 &
        .OR. this%it .LE. 0 .OR. this%ip .LE. 0) THEN
      WRITE (*, 1) this%iu, this%iv, this%iq, this%it, this%ip
1     FORMAT('dLogP - GeosBal weak constraint is not applicable: indices of u/v/q/t/p: ', 5I2)
      RETURN
    END IF

    ALLOCATE (hydro(this%ng(1), this%ng(2), this%ng(3)), &
              hydrT(this%ng(1), this%ng(2), this%ng(3)), &
              pnorm(this%ng(1), this%ng(2), this%ng(3)), &
              ptang(this%ng(1), this%ng(2), this%ng(3)), &
              logPz(this%ng(1), this%ng(2), this%ng(3)))

    ASSOCIATE (u => X%fields(this%iu)%DATA, v => X%fields(this%iv)%DATA, &
               q => X%fields(this%iq)%DATA, p => X%fields(this%ip)%DATA, &
               t => X%fields(this%it)%DATA, &
               sg => X%sg, coef => X%sg%coef_sigma, h => X%sg%zHght)
      logPx = 0.0D0; logPy = 0.0D0; hydro = 0.0D0
      pnorm = 0.0D0; ptang = 0.0D0; logPz = 0.0D0

      ! Hydro pressure:
      hydro(1, :, :) = p(1, :, :) ! Surface pressure
      PRINT *, 'HydroPressure range: ', MINVAL(p(1, :, :)), MAXVAL(p(1, :, :))
      DO k = 1, this%ng(3)
        DO j = 1, this%ng(2)
          CALL hydroPressure(this%ng(1), h(:, j), t(:, j, k), q(:, j, k), hydro(:, j, k))
        END DO
      END DO
      hydrT(this%ng(1), :, :) = p(this%ng(1), :, :) ! Top pressure
      DO k = 1, this%ng(3)
        DO j = 1, this%ng(2)
          CALL hydroPressTop(this%ng(1), h(:, j), t(:, j, k), q(:, j, k), hydrT(:, j, k))
        END DO
      END DO
      PRINT *, 'HydroPressure at top range: ', MINVAL(p(this%ng(1), :, :)), MAXVAL(p(this%ng(1), :, :))

      ! pressure gradient terms: interior points only
      DO j = 1, sg%num_icell
        rf = dry_air_gas_const / (2.0D0 * Omega * DSIN(sg%cell_cntr(1, j)))

        DO iedge = 1, UBOUND(sg%edge_stcl, 2)
          DO istcl = 1, UBOUND(sg%edge_stcl, 1)
            ! Valid interpolation points:
            IF (sg%edge_stcl(istcl, iedge, j) .GT. 0) THEN
              pnorm(:, j, :) = pnorm(:, j, :) + &
                               sg%coef_norm(istcl, 1, iedge, j) * LOG(hydro(:, sg%edge_stcl(istcl, iedge, j), :))
              ptang(:, j, :) = ptang(:, j, :) + &
                               sg%coef_tang(istcl, 1, iedge, j) * LOG(hydro(:, sg%edge_stcl(istcl, iedge, j), :))
            END IF
          END DO
          logPx(:, j, :) = logPx(:, j, :) + sg%edgeNorm2(1, iedge, j) * pnorm(:, j, :) + &
                           sg%edgeTang2(1, iedge, j) * ptang(:, j, :)
          logPy(:, j, :) = logPy(:, j, :) + sg%edgeNorm2(2, iedge, j) * pnorm(:, j, :) + &
                           sg%edgeTang2(2, iedge, j) * ptang(:, j, :)
        END DO
        ! Average:
        logPx(:, j, :) = logPx(:, j, :) / DBLE(UBOUND(sg%edge_stcl, 2))
        logPy(:, j, :) = logPy(:, j, :) / DBLE(UBOUND(sg%edge_stcl, 2))

        ! The vertical derivatives:
        DO i = 2, sg%vLevel - 1
          ! Arbitrary vertical grid:
          logPz(i, j, :) = LOG(hydro(i - 1, j, :)) * coef(1, i) + LOG(hydro(i, j, :)) * coef(2, i) + LOG(hydro(i + 1, j, :)) * coef(3, i)
        END DO

        DO k = 1, this%ng(3)
          logPx(:, j, k) = logPx(:, j, k) - logPz(:, j, k) * sg%sigmax(:, j)
          logPy(:, j, k) = logPy(:, j, k) - logPz(:, j, k) * sg%sigmay(:, j)
        END DO
      END DO

      PRINT *, 'ln(P)_xy: ', MINVAL(logPx), MAXVAL(logPx), MINVAL(logPy), MAXVAL(logPy)
      PRINT *, 'ln(P)_z : ', MINVAL(logPz), MAXVAL(logPz), MINVAL(sg%sigmax), MAXVAL(sg%sigmax), MINVAL(sg%sigmay), MAXVAL(sg%sigmay)
      PRINT *, 'coef: norm/tang: ', MINVAL(sg%coef_norm(:, :, :, :)), MAXVAL(sg%coef_norm(:, :, :, :)), &
        MINVAL(sg%coef_tang(:, :, :, :)), MAXVAL(sg%coef_tang(:, :, :, :))

      BLOCK
        TYPE(state_t) :: DD

        DD = X%zeroCopy()
        CALL DD%fields(1)%set_name('hydroP')
        DD%fields(1)%DATA = hydro !LOG(hydro)
        CALL DD%fields(2)%set_name('logPx')
        DO j = 1, sg%num_cell
          DD%fields(2)%DATA(:, j, :) = dry_air_gas_const / &
                                       (2.0D0 * Omega * DSIN(sg%cell_cntr(1, j))) * (1.0D0 + 0.608 * q(:, j, :)) * t(:, j, :) * logPx(:, j, :) !logPx
        END DO

        CALL DD%fields(3)%set_name('logPy')

        DO j = 1, sg%num_cell
          DD%fields(3)%DATA(:, j, :) = dry_air_gas_const / &
                                       (2.0D0 * Omega * DSIN(sg%cell_cntr(1, j))) * (1.0D0 + 0.608 * q(:, j, :)) * t(:, j, :) * logPy(:, j, :) !logPy
        END DO
        CALL DD%fields(4)%set_name('LogT')
        DD%fields(4)%DATA = p !hydrT !LOG(hydrT)
        CALL DD%fields(5)%set_name('Topo')
        DO j = 1, sg%num_icell
          DO k = 1, this%ng(3)
            DD%fields(5)%DATA(1, j, k) = sg%topo(j) !logPz(:,j,k)*sg%sigmax(:,j)
          END DO
        END DO
        CALL Output_NC_State_AV(DD, "/Users/yjzx-xieyuanfu/developments/da/motor_training/test/MOTOR/output", &
                                "XIE"//"_LogP", .TRUE., .TRUE.)
      END BLOCK
    END ASSOCIATE

    DEALLOCATE (hydro, hydrT, pnorm, ptang)
  END SUBROUTINE dLogP

  REAL(r_kind) FUNCTION JcGeos(this, X) RESULT(geos)
    USE conversions_m, ONLY: hydroPressure
    IMPLICIT NONE

    CLASS(GeosBal_t) :: this
    TYPE(State_t), INTENT(IN) :: X

    ! Local variables:
    INTEGER(i_kind) :: j, k, iu, iv, ip, ir, it, ng(3), iedge, istcl
    REAL(r_kind) :: rf
    REAL(r_kind), ALLOCATABLE :: logPx(:, :, :), logPy(:, :, :), hydro(:, :, :)

    geos = 0.0D0
    IF (this%iu .GT. 0 .AND. this%iv .GT. 0 .AND. &
        this%ip .GT. 0 .AND. this%iq .GT. 0) THEN

    END IF

  END FUNCTION JcGeos

! Weak constraint procedures:
  SUBROUTINE dJcDUV(this, X, Y, iflag)
    USE conversions_m, ONLY: hydroPressure, dHydroPs
    IMPLICIT NONE

    CLASS(GeosBal_t) :: this
    TYPE(State_t), INTENT(IN) :: X
    TYPE(State_t), INTENT(OUT) :: Y
    INTEGER(i_kind), INTENT(IN) :: iflag ! 0 no dPs calculated; 1 calculate dPs

    ! Local variables:
    INTEGER(i_kind) :: j, k, iu, iv, iq, ip, it, ng(3), iedge, istcl
    REAL(r_kind) :: rf
    REAL(r_kind), ALLOCATABLE :: logPx(:, :, :), logPy(:, :, :), hydro(:, :, :), dPdPs(:, :, :), pnorm(:, :, :), ptang(:, :, :)

    Y = X%zeroCopy()

    iu = X%getVarIdx('uwnd')
    iv = X%getVarIdx('vwnd')
    iq = X%getVarIdx('qvapor')
    it = X%getVarIdx('temp')
    ip = X%getVarIdx('pres')

    IF (iu .LE. 0 .OR. iv .LE. 0 .OR. iq .LE. 0 .OR. it .LE. 0 .OR. ip .LE. 0) THEN
      WRITE (*, 1) iu, iv, iq, it, ip
1     FORMAT('GeosBal weak constraint is not applicable: indices of u/v/q/t/p: ', 5I2)
      RETURN
    END IF

    ALLOCATE (logPx(this%ng(1), this%ng(2), this%ng(3)), &
              logPy(this%ng(1), this%ng(2), this%ng(3)), &
              hydro(this%ng(1), this%ng(2), this%ng(3)), &
              dPdPs(this%ng(1), this%ng(2), this%ng(3)), &
              pnorm(this%ng(1), this%ng(2), this%ng(3)), &
              ptang(this%ng(1), this%ng(2), this%ng(3)))

    ASSOCIATE (u => X%fields(iu)%DATA, v => X%fields(iv)%DATA, &
               r => X%fields(iq)%DATA, p => X%fields(ip)%DATA, &
               t => X%fields(it)%DATA, &
               sg => X%sg, coef => X%sg%coef_sigma, h => X%sg%zHght)
      logPx = 0.0D0; logPy = 0.0D0; hydro = 0.0D0
      dPdPs = 0.0D0; pnorm = 0.0D0; ptang = 0.0D0

      ! Calculate the log pressure derivatives:
      PRINT *, 'calling dlogP...'
      CALL this%dLogP(X, logPx, logPy)

      ! Hydro pressure and dP/dPs:
      IF (iflag .EQ. 1) THEN
        DO k = 1, ng(3)
          DO j = 1, ng(2)
            CALL dHydroPs(ng(1), h(:, j), t(:, j, k), r(:, j, k), dPdPs(:, j, k))
          END DO
        END DO
      END IF

      ! pressure gradient terms: interior points only

      PRINT *, 'pressureGrad: ', &
        MINVAL(dry_air_gas_const / (2.0D0 * Omega * DSIN(sg%cell_cntr(1, :))) * (1.0D0 + 0.608 * r(2, :, 1)) * t(2, :, 1) * logPy(2, :, 1)), &
        MAXVAL(dry_air_gas_const / (2.0D0 * Omega * DSIN(sg%cell_cntr(1, :))) * (1.0D0 + 0.608 * r(2, :, 1)) * t(2, :, 1) * logPy(2, :, 1)), &
        MINVAL(dry_air_gas_const / (2.0D0 * Omega * DSIN(sg%cell_cntr(1, :))) * (1.0D0 + 0.608 * r(2, :, 1)) * t(2, :, 1) * logPx(2, :, 1)), &
        MAXVAL(dry_air_gas_const / (2.0D0 * Omega * DSIN(sg%cell_cntr(1, :))) * (1.0D0 + 0.608 * r(2, :, 1)) * t(2, :, 1) * logPx(2, :, 1))
      PRINT *, 'temperaturelogP: ', &
        MINVAL(t(1, :, 1) * logPy(1, :, 1)), MAXVAL(t(1, :, 1) * logPy(1, :, 1)), &
        MINVAL(dry_air_gas_const / (2.0D0 * Omega * DSIN(sg%cell_cntr(1, :)))), &
        MAXVAL(dry_air_gas_const / (2.0D0 * Omega * DSIN(sg%cell_cntr(1, :)))), dry_air_gas_const, Omega, &
        MINVAL(DSIN(sg%cell_cntr(1, :))), MAXVAL(DSIN(sg%cell_cntr(1, :))), &
        MINVAL(sg%cell_dist), MAXVAL(sg%cell_dist)
      STOP
      DO j = 1, sg%num_icell
        rf = dry_air_gas_const / (2.0D0 * Omega * DSIN(sg%cell_cntr(1, j)))
        Y%fields(iu)%DATA(:, j, :) = X%fields(iu)%DATA(:, j, :) + &
                                     rf * (1.0D0 + 0.608 * r(:, j, :)) * t(:, j, :) * logPy(:, j, :)
        Y%fields(iv)%DATA(:, j, :) = X%fields(iv)%DATA(:, j, :) - &
                                     rf * (1.0D0 + 0.608 * r(:, j, :)) * t(:, j, :) * logPx(:, j, :)

        IF (j .EQ. 938) PRINT *, 'Vwnd NaN: 938', Y%fields(iv)%DATA(1, j, 1)
        IF (iflag .EQ. 1) THEN
          pnorm(:, j, :) = 0.0D0
          ptang(:, j, :) = 0.0D0
          DO iedge = 1, UBOUND(sg%edge_stcl, 2)
            DO istcl = 1, UBOUND(sg%edge_stcl, 1)
              IF (sg%edge_stcl(istcl, iedge, j) .GT. 0) THEN
              END IF
            END DO
          END DO
          pnorm(:, j, :) = &
            Y%fields(iu)%DATA(:, j, :) * rf * t(:, j, :)
          ptang(:, j, :) = &
            Y%fields(iv)%DATA(:, j, :) * rf * t(:, j, :)

          DO iedge = 1, UBOUND(sg%edge_stcl, 2)
            DO istcl = 1, UBOUND(sg%edge_stcl, 1)
              IF (sg%edge_stcl(istcl, iedge, j) .GT. 0) THEN
              END IF
            END DO
          END DO
        END IF
      END DO
    END ASSOCIATE

    DEALLOCATE (logPx, logPy, dPdPs, pnorm, ptang)

  END SUBROUTINE dJcDUV

END MODULE GeosBal_m
