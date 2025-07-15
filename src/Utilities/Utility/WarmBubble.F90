MODULE WarmBubble_m
  USE SphericDisCal_m
  USE parameters_m, ONLY: pi, r_d, cp, g, k_d
  USE kinds_m, ONLY: i_kind, r_kind

CONTAINS
  SUBROUTINE GenWarmBubble(lat, lon, z, theta, rho, pres)
    IMPLICIT NONE

    ! 输入参数
    REAL(r_kind), DIMENSION(:, :), INTENT(in) :: z  !z为高度（米）
    REAL(r_kind), DIMENSION(:), INTENT(in) :: lat, lon  ! 纬度和经度（度），z为高度（米）

    ! 输出
    REAL(r_kind), DIMENSION(:, :), INTENT(out) :: theta, rho, pres

    ! 局部变量
    REAL(r_kind), PARAMETER :: theta_bar = 300.0D0  ! 势温（K）
    REAL(r_kind), PARAMETER :: theta_c = 0.5D0      ! K
    REAL(r_kind), PARAMETER :: r_c = 250.0D0        ! 米
    REAL(r_kind), ALLOCATABLE :: r(:, :)
    REAL(r_kind) :: lat_c, lon_c, z_c, ik, P0
    REAL(r_kind), ALLOCATABLE :: theta_prime(:, :), dpres(:, :), presb(:, :)
    INTEGER(i_kind) :: i, k, sz, sx, lev

    sz = SIZE(z, dim=1)
    sx = SIZE(z, dim=2)
    P0 = 1.0D5

    ! 中心坐标（需要根据具体情况定义）
    lat_c = (MAXVAL(lat) + MINVAL(lat)) / 2.0D0  ! 中心纬度
    lon_c = (MAXVAL(lon) + MINVAL(lon)) / 2.0D0  ! 中心经度
    z_c = 350.0  ! 中心高度（米）

    ! 使用spatial_distance函数计算r
    ALLOCATE (r(sz, sx), theta_prime(sz, sx), &
              presb(sz, sx), dpres(sz, sx))
    r = 0.0D0
    CALL SpatialDistance(lat, lon, z, lat_c, lon_c, z_c, r)

    PRINT *, 'max r is: ', MAXVAL(r)

    ! 基于r计算theta_prime
    DO i = 1, sx
      DO k = 1, sz
        IF (r(k, i) .GT. r_c) THEN
          theta_prime(k, i) = 0.0D0
        ELSE
          theta_prime(k, i) = (theta_c / 2.0) * (1.0 + COS(pi * r(k, i) / r_c))
        END IF
      END DO
    END DO
    theta = theta_prime + theta_bar

    presb = P0 * (1.0D0 - g / cp / theta_bar * z)**(cp / r_d)
    dpres = P0 * (1.0D0 - g / cp / theta_bar * z)**(cp / r_d - 1.0D0) * g * z / theta_bar**2.0D0 / r_d * theta_prime
    DO k = 1, sz
      IF (dabs(z(k, 1) - z_c) .LT. 0.1D0) lev = k
    END DO
    dpres(1:lev - 1, :) = -dpres(1:lev - 1, :)
    pres = presb! + dpres
    ! pres = P0**k_d
    ! PRINT *, 'max of pres in gen bubble is: ', maxval(pres)
    ! ik = 1.0D0
    ! DO k = 1, sz
    !    PRINT *, 'k is: ', k
    !    pres(k, :) = pres(1, :)*(1.0D0 - g/cp*SUM(1.0D0/theta(1:k, :), dim=1)/ik*z(k, 1))
    !    PRINT *, 'z(k, 1) is: ', z(k, 1)
    !    ik = ik + 1
    ! END DO
    ! pres = pres**(1/k_d)
    PRINT *, 'max of pres after warm bubble adjust in gen bubble is: ', MAXVAL(pres)

    rho = pres**(1 - k_d) * P0**k_d / r_d / theta

    DEALLOCATE (r, theta_prime, &
                presb, dpres)

  END SUBROUTINE GenWarmBubble

END MODULE
