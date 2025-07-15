!----------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Centers
!                     for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Jiongming Pang
! VERSION           : V 0.0
! HISTORY           :
!   Created by Jiongming Pang (pang.j.m@hotmail.com), 2021-12-1, @GBA-MWF, Shenzhen
!   Modified by Zilong Qin (zilong.qin@gmail.com), 2022/6/25, @GBA-MWF, Shenzhen
!----------------------------------------------------------------------------------------
MODULE conversions_m
  USE kinds_m
  IMPLICIT NONE

  ! a interface for converting wind speed and direction to u, v wind for 2-d data
  ! added by Ting Shu @2023-02-20
  PRIVATE :: wind_to_uv0d, wind_to_uv1d, wind_to_uv2d
  INTERFACE wind_to_uv
    MODULE PROCEDURE :: wind_to_uv0d, wind_to_uv1d, wind_to_uv2d
  END INTERFACE wind_to_uv

CONTAINS

  PURE ELEMENTAL FUNCTION RH_to_qvapor(RH, temp, pres) RESULT(qvapor)
    USE parameters_m, ONLY: KelvinsT0
    REAL(r_kind), INTENT(IN)  ::  RH
    REAL(r_kind), INTENT(IN)  ::  temp
    REAL(r_kind), INTENT(IN)  ::  pres
    REAL(r_kind)  ::  qvapor

    REAL(r_kind) :: Es, E

    Es = 6.1094D0 * EXP(17.625D0 * (temp - KelvinsT0) / (temp - KelvinsT0 + 243.04D0))
    E = RH / 100.D0 * Es * 100.0D0
    qvapor = 0.622D0 * E / (pres - E)

  END FUNCTION

  PURE ELEMENTAL FUNCTION Saturated_Vapor(temp, pres) RESULT(qvapor)
    USE parameters_m, ONLY: KelvinsT0
    REAL(r_kind), INTENT(IN)  ::  temp
    REAL(r_kind), INTENT(IN)  ::  pres
    REAL(r_kind)  ::  qvapor

    REAL(r_kind) :: Es, E

    Es = 6.1094D0 * EXP(17.625D0 * (temp - KelvinsT0) / (temp - KelvinsT0 + 243.04D0))
    E = Es * 100.0D0
    qvapor = 0.622D0 * E / (pres - E)

  END FUNCTION

  PURE ELEMENTAL FUNCTION Td_to_qvapor(Td, pres) RESULT(qvapor)
    USE parameters_m, ONLY: KelvinsT0
    REAL(r_kind), INTENT(IN)  ::  Td    ! In Kelvin
    REAL(r_kind), INTENT(IN)  ::  pres  ! In Pascal
    REAL(r_kind)  ::  qvapor

    REAL(r_kind) :: E

    E = 6.1094D0 * EXP(17.625D0 * (Td - KelvinsT0) / (Td - KelvinsT0 + 243.04D0))
    qvapor = 0.622D0 * E * 100.0D0 / (pres - E * 100.0D0)
  END FUNCTION

  SUBROUTINE wind_to_uv0d(windSp, windDir, u, v)
    REAL(r_kind), INTENT(IN)  ::  windSp
    REAL(r_kind), INTENT(IN)  ::  windDir
    REAL(r_kind), INTENT(INOUT)  ::  u
    REAL(r_kind), INTENT(INOUT)  ::  v

    REAL(r_kind), PARAMETER :: PI = 3.1415926535897932
    u = -SIN(windDir * PI / 180.0D0) * windSp
    v = -COS(windDir * PI / 180.0D0) * windSp

  END SUBROUTINE wind_to_uv0d

  ! convert wind speed and direction to u, v wind for 2-d data
  ! added by Ting Shu @2023-02-20
  SUBROUTINE wind_to_uv1d(windSp, windDir, u, v)
    REAL(r_kind), INTENT(IN)  ::  windSp(:)
    REAL(r_kind), INTENT(IN)  ::  windDir(:)
    REAL(r_kind), INTENT(INOUT)  ::  u(:)
    REAL(r_kind), INTENT(INOUT)  ::  v(:)

    REAL(r_kind), PARAMETER :: PI = 3.1415926535897932

    u = -SIN(windDir * PI / 180.0D0) * windSp
    v = -COS(windDir * PI / 180.0D0) * windSp

  END SUBROUTINE wind_to_uv1d

  ! convert wind speed and direction to u, v wind for 2-d data
  ! added by Ting Shu @2023-02-20
  SUBROUTINE wind_to_uv2d(windSp, windDir, u, v)
    REAL(r_kind), INTENT(IN)  ::  windSp(:, :)
    REAL(r_kind), INTENT(IN)  ::  windDir(:, :)
    REAL(r_kind), INTENT(INOUT)  ::  u(:, :)
    REAL(r_kind), INTENT(INOUT)  ::  v(:, :)

    REAL(r_kind), PARAMETER :: PI = 3.1415926535897932

    u = -SIN(windDir * PI / 180.0D0) * windSp
    v = -COS(windDir * PI / 180.0D0) * windSp

  END SUBROUTINE wind_to_uv2d

  PURE ELEMENTAL SUBROUTINE uv_to_wind(u, v, windSp, windDir)
    REAL(r_kind), INTENT(IN) :: u
    REAL(r_kind), INTENT(IN) :: v
    REAL(r_kind), INTENT(INOUT) :: windSp
    REAL(r_kind), INTENT(INOUT) :: windDir

    REAL(r_kind), PARAMETER :: PI = 3.1415926535897932
    REAL(r_kind) :: tmp
    ! INTEGER :: i

    ! DO i=1,SIZE(windSp)
    windSp = SQRT(u**2 + v**2)

    tmp = 270.0D0 - ATAN2(v, u) * 180.0D0 / PI
    windDir = MOD(tmp, 360.0D0)
    ! END DO

  END SUBROUTINE uv_to_wind

  ! Yuanfu Xie adopts Jiongming's geopotential height code to calculate hydroHeight from ground to top:
  !   INPUT:
  !     n is the number of vertical grid points;
  !     ps is the pressure to start the integration from bottom
  !     p is the pressure field;
  !     t is the temperature field;
  !     r is the water vapor field;
  !
  !   OUTPUT:
  !     h is the hydrostatic height;
  !
  SUBROUTINE hydroHeightUp(n, p, t, r, h)
    USE parameters_m, ONLY: g, dry_air_gas_const
    IMPLICIT NONE
    INTEGER(i_kind), INTENT(IN) :: n
    REAL(r_kind), INTENT(IN) :: p(n), t(n), r(n) ! pls is the pressure to start integrating from bottom
    REAL(r_kind), INTENT(INOUT) :: h(n)

    ! Local variables:
    INTEGER(i_kind) :: k
    REAL(r_kind) :: tv(n), deltaLnP

    ! Integrating from the bottom up:
    DO k = 1, n - 1
      tv(k) = 0.5D0 * (t(k + 1) * (r(k + 1) * 0.608 + 1.0D0) + t(k) * (r(k) * 0.608 + 1.0D0))

      deltaLnP = LOG(p(k)) - LOG(p(k + 1)) ! Assume pressure decreases with height
      h(k + 1) = h(k) + deltaLnP / g * dry_air_gas_const * tv(k)
    END DO

  END SUBROUTINE hydroHeightUp

  ! Yuanfu Xie adopts Jiongming's geopotential height code to calculate hydroHeight from ground to top:
  !   INPUT:
  !     n is the number of vertical grid points;
  !     ps is the pressure to start the integration from bottom
  !     p is the pressure field;
  !     t is the temperature field;
  !     r is the water vapor field;
  !
  !   OUTPUT:
  !     h is the hydrostatic height;
  !
  SUBROUTINE hydroHeightDown(n, p, t, r, h)
    USE parameters_m, ONLY: g, dry_air_gas_const
    IMPLICIT NONE
    INTEGER(i_kind), INTENT(IN) :: n
    REAL(r_kind), INTENT(IN) :: p(n), t(n), r(n) ! pls is the pressure to start integrating from bottom
    REAL(r_kind), INTENT(INOUT) :: h(n)

    ! Local variables:
    INTEGER(i_kind) :: k
    REAL(r_kind) :: tv(n), deltaLnP

    ! Integrating from the bottom up:
    DO k = n, 2, -1
      tv(k) = 0.5D0 * (t(k - 1) * (r(k - 1) * 0.608 + 1.0D0) + t(k) * (r(k) * 0.608 + 1.0D0))

      deltaLnP = LOG(p(k - 1)) - LOG(p(k)) ! Assume pressure decreases with height
      h(k - 1) = h(k) - deltaLnP / g * dry_air_gas_const * tv(k)
    END DO

  END SUBROUTINE hydroHeightDown

  ! Yuanfu Xie adopts Jiongming's geopotential height code to develop hydroHeight:
  !   INPUT:
  !     n is the number of vertical grid points;
  !     pl is the pressure to start the integration from toward the bottom
  !     p is the pressure field;
  !     t is the temperature field;
  !     r is the water vapor field;
  !
  !   OUTPUT:
  !     h is the hydrostatic height;
  !
  SUBROUTINE hydroHeight(n, pl, p, t, r, h)
    USE parameters_m, ONLY: g, dry_air_gas_const
    IMPLICIT NONE
    INTEGER(i_kind), INTENT(IN) :: n
    REAL(r_kind), INTENT(IN) :: pl, p(n), t(n), r(n) ! pl is the pressure to start integrating from toward bottom
    REAL(r_kind), INTENT(INOUT) :: h(n)

    ! Local variables:
    INTEGER(i_kind) :: lvl, k
    REAL(r_kind) :: tv(n), dp

    ! Find lvl to start integrating: p is stored from bottom to top

    !!! Untested and still need to work on this!!!

    ! DO k=1,n
    !   gph(:,:,index_begPress:nz) = gpht(:,:,index_begPress:nz)
    ! do k=1,index_begPress-1
    !   !kd=nz-k
    !   kd=index_begPress-k
    !   read(lev(kd+1),*) ph
    !   read(lev(kd),*) pl
    !   ttmp = (tv(:,:,kd+1) + tv(:,:,kd)) / 2.0d0

    !   do i=1,nx
    !     do j=1,ny
    !       deltaP = log(pl) - log(ph)
    !       dgph = deltaP / 9.80665 * 287.0 * ttmp(i,j)
    !       gph(i,j,kd) = gph(i,j,kd+1) - dgph
    !     end do
    !   end do
    ! end do
    ! END DO

  END SUBROUTINE hydroHeight

  ! Yuanfu Xie adds a hydrostatic pressure calculation along with its adjoints:
  !   INPUT:
  !     n is the number of vertical grid points;
  !     o is the option of gradients in terms of surface pressure or vertical pressure
  !     h is the height field;
  !     t is the temperature field;
  !     r is the water vapor field;
  !     p(1) contains surface pressure as input;
  !
  !   OUTPUT:
  !     p(2:n) hydrostatic pressure;
  !
  SUBROUTINE hydroPressure(n, h, t, r, p)
    USE parameters_m, ONLY: g, dry_air_gas_const
    IMPLICIT NONE
    INTEGER(i_kind), INTENT(IN) :: n
    REAL(r_kind), INTENT(IN) :: h(n), t(n), r(n)
    REAL(r_kind), INTENT(INOUT) :: p(n)

    ! Local variables:
    INTEGER(i_kinD) :: i

    PRINT *, 'hydroPressure: ', n, p(1)
    DO i = 2, n
      p(i) = p(i - 1) * EXP(-(h(i) - h(i - 1)) * g / (dry_air_gas_const * &
                                                      0.5D0 * (t(i) * (r(i) * 0.608 + 1.0D0) + &
                                                               t(i - 1) * (r(i - 1) * 0.608 + 1.0D0))))
    END DO
  END SUBROUTINE hydroPressure

  ! Yuanfu Xie adds a hydrostatic pressure calculation along with its adjoints:
  !   INPUT:
  !     n is the number of vertical grid points;
  !     o is the option of gradients in terms of surface pressure or vertical pressure
  !     h is the height field;
  !     t is the temperature field;
  !     r is the water vapor field;
  !     p(1) contains surface pressure as input;
  !
  !   OUTPUT:
  !     p(2:n) hydrostatic pressure;
  !
  SUBROUTINE hydroPressTop(n, h, t, r, p)
    USE parameters_m, ONLY: g, dry_air_gas_const
    IMPLICIT NONE
    INTEGER(i_kind), INTENT(IN) :: n
    REAL(r_kind), INTENT(IN) :: h(n), t(n), r(n)
    REAL(r_kind), INTENT(INOUT) :: p(n)

    ! Local variables:
    INTEGER(i_kinD) :: i

    DO i = n - 1, 1, -1
      p(i) = p(i + 1) / EXP(-(h(i + 1) - h(i)) * g / (dry_air_gas_const * &
                                                      0.5D0 * (t(i + 1) * (r(i + 1) * 0.608 + 1.0D0) + &
                                                               t(i) * (r(i) * 0.608 + 1.0D0))))
    END DO
  END SUBROUTINE hydroPressTop

  ! Calculate partial P/partial Ps:
  SUBROUTINE dHydroPs(n, h, t, r, dPs)
    USE parameters_m, ONLY: g, dry_air_gas_const
    IMPLICIT NONE
    INTEGER(i_kind), INTENT(IN) :: n
    REAL(r_kind), INTENT(IN) :: h(n), t(n), r(n)
    REAL(r_kind), INTENT(OUT) :: dPs(n)

    ! Local variables:
    INTEGER(i_kind) :: i

    dPs = 0.0D0

    DO i = 2, n
      dPs(i) = dPs(i - 1) - (h(i) - h(i - 1)) * g / (dry_air_gas_const * &
                                                     0.5D0 * (t(i) * (r(i) * 0.608 + 1.0D0) + &
                                                              t(i - 1) * (r(i - 1) * 0.608 + 1.0D0)))
    END DO
    dPs = EXP(dPs)
  END SUBROUTINE dHydroPs

  ! Convert the surface pressure to the sea level pressure
  PURE ELEMENTAL FUNCTION psfc_to_psl(psfc, h, tsfc, tmsl) RESULT(psl)
    USE parameters_m, ONLY: g, dry_air_gas_const
    REAL(r_kind), INTENT(IN) :: psfc
    REAL(r_kind), INTENT(IN) :: h
    REAL(r_kind), INTENT(IN) :: tsfc
    REAL(r_kind), INTENT(IN) :: tmsl
    REAL(r_kind) :: psl

    psl = psfc * EXP(2.0D0 * g * h / (dry_air_gas_const * (tsfc + tmsl)))

  END FUNCTION

  ! Convert the sea level pressure to the surface pressure
  PURE ELEMENTAL FUNCTION psl_to_psfc(psl, h, tsfc, tmsl) RESULT(psfc)
    USE parameters_m, ONLY: g, dry_air_gas_const
    REAL(r_kind), INTENT(IN) :: psl
    REAL(r_kind), INTENT(IN) :: h
    REAL(r_kind), INTENT(IN) :: tsfc
    REAL(r_kind), INTENT(IN) :: tmsl
    REAL(r_kind) :: psfc

    psfc = psl * EXP(-2.0D0 * g * h / (dry_air_gas_const * (tsfc + tmsl)))

  END FUNCTION

END MODULE conversions_m
