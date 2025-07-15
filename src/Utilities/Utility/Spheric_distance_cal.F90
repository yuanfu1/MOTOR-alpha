MODULE SphericDisCal_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE parameters_m, ONLY: EarthRadius, pi

  INTERFACE SpatialDistance
    MODULE PROCEDURE spatial_distance_0d
    MODULE PROCEDURE spatial_distance_1d
    MODULE PROCEDURE spatial_distance_2d

  END INTERFACE

CONTAINS

  SUBROUTINE spatial_distance_0d(lat1, lon1, h1, lat2, lon2, h2, distance)
    IMPLICIT NONE

    REAL(r_kind), INTENT(in) :: lat1, lon1, h1  ! 第一点的纬度、经度（度）和高度（米）
    REAL(r_kind), INTENT(in) :: lat2, lon2, h2  ! 第二点的纬度、经度（度）和高度（米）
    REAL(r_kind), INTENT(out) :: distance  ! 距离（米）

    ! 局部变量
    REAL(r_kind), PARAMETER :: deg_to_rad = pi / 180.0D0
    REAL(r_kind) :: x1, y1, z1, x2, y2, z2
    REAL(r_kind) :: rlat1, rlat2, rlon1, rlon2

    ! 将度转换为弧度
    rlat1 = lat1 * deg_to_rad
    rlat2 = lat2 * deg_to_rad
    rlon1 = lon1 * deg_to_rad
    rlon2 = lon2 * deg_to_rad

    ! 计算第一点的笛卡尔坐标
    x1 = (EarthRadius + h1) * dcos(rlat1) * dcos(rlon1)
    y1 = (EarthRadius + h1) * dcos(rlat1) * dsin(rlon1)
    z1 = (EarthRadius + h1) * dsin(rlat1)

    ! 计算第二点的笛卡尔坐标
    x2 = (EarthRadius + h2) * dcos(rlat2) * dcos(rlon2)
    y2 = (EarthRadius + h2) * dcos(rlat2) * dsin(rlon2)
    z2 = (EarthRadius + h2) * dsin(rlat2)

    ! 计算两点之间的直线距离
    distance = SQRT((x2 - x1)**2.0D0 + (y2 - y1)**2.0D0 + (z2 - z1)**2.0D0)
  END SUBROUTINE spatial_distance_0d

  SUBROUTINE spatial_distance_1d(lat1, lon1, h1, lat2, lon2, h2, distance)
    IMPLICIT NONE

    REAL(r_kind), DIMENSION(:), INTENT(in) :: lat1, lon1, h1  ! 第一点的纬度、经度（度）和高度（米）
    REAL(r_kind), INTENT(in) :: lat2, lon2, h2  ! 第二点的纬度、经度（度）和高度（米）
    REAL(r_kind), DIMENSION(:), INTENT(out) :: distance(:)  ! 距离（米）

    ! 局部变量
    REAL(r_kind), PARAMETER :: deg_to_rad = pi / 180.0D0
    REAL(r_kind), DIMENSION(:), ALLOCATABLE :: x1, y1, z1
    REAL(r_kind) :: x2, y2, z2
    REAL(r_kind), DIMENSION(:), ALLOCATABLE :: rlat1, rlon1
    REAL(r_kind) :: rlat2, rlon2

    ! 将度转换为弧度
    rlat1 = lat1 * deg_to_rad
    rlat2 = lat2 * deg_to_rad
    rlon1 = lon1 * deg_to_rad
    rlon2 = lon2 * deg_to_rad

    ! 计算第一点的笛卡尔坐标
    x1 = (EarthRadius + h1) * dcos(rlat1) * dcos(rlon1)
    y1 = (EarthRadius + h1) * dcos(rlat1) * dsin(rlon1)
    z1 = (EarthRadius + h1) * dsin(rlat1)

    ! 计算第二点的笛卡尔坐标
    x2 = (EarthRadius + h2) * dcos(rlat2) * dcos(rlon2)
    y2 = (EarthRadius + h2) * dcos(rlat2) * dsin(rlon2)
    z2 = (EarthRadius + h2) * dsin(rlat2)

    ! 计算两点之间的直线距离
    distance = SQRT((x2 - x1)**2.0D0 + (y2 - y1)**2.0D0 + (z2 - z1)**2.0D0)
    !   deallocate (rlat1, lon1, rlat2, lon2, x1, y1, z1, x2, y2, z2)
  END SUBROUTINE spatial_distance_1d

  SUBROUTINE spatial_distance_2d(lat1, lon1, h1, lat2, lon2, h2, distance)
    IMPLICIT NONE

    REAL(r_kind), DIMENSION(:), INTENT(in) :: lat1, lon1  ! 第一点的纬度、经度（度）
    REAL(r_kind), DIMENSION(:, :), INTENT(in) :: h1  ! 第一点的高度（米）
    REAL(r_kind), INTENT(in) :: lat2, lon2, h2  ! 第二点的纬度、经度（度）和高度（米）
    REAL(r_kind), DIMENSION(:, :), INTENT(out) :: distance  ! 距离（米）

    ! 局部变量
    REAL(r_kind), PARAMETER :: deg_to_rad = 1.0D0 !pi/180.0D0
    REAL(r_kind), DIMENSION(:, :), ALLOCATABLE :: x1, y1, z1
    REAL(r_kind) :: x2, y2, z2
    REAL(r_kind), DIMENSION(:), ALLOCATABLE :: rlat1, rlon1
    REAL(r_kind) :: rlat2, rlon2
    INTEGER(i_kind) :: k

    ! 将度转换为弧度
    rlat1 = lat1 * deg_to_rad
    rlat2 = lat2 * deg_to_rad
    rlon1 = lon1 * deg_to_rad
    rlon2 = lon2 * deg_to_rad

    ! 计算第一点的笛卡尔坐标
    PRINT *, 'size of distance is: ', SIZE(distance, dim=1), SIZE(distance, dim=2)
    ALLOCATE (x1(SIZE(distance, dim=1), SIZE(distance, dim=2)), &
              y1(SIZE(distance, dim=1), SIZE(distance, dim=2)), &
              z1(SIZE(distance, dim=1), SIZE(distance, dim=2)))

    DO k = 1, SIZE(x1, dim=1)
      x1(k, :) = (EarthRadius + h1(k, :)) * dcos(rlat1) * dcos(rlon1)
      y1(k, :) = (EarthRadius + h1(k, :)) * dcos(rlat1) * dsin(rlon1)
      z1(k, :) = (EarthRadius + h1(k, :)) * dsin(rlat1)
    END DO

    ! 计算第二点的笛卡尔坐标
    x2 = (EarthRadius + h2) * dcos(rlat2) * dcos(rlon2)
    y2 = (EarthRadius + h2) * dcos(rlat2) * dsin(rlon2)
    z2 = (EarthRadius + h2) * dsin(rlat2)

    ! 计算两点之间的直线距离

    distance = SQRT((x2 - x1)**2.0D0 + (y2 - y1)**2.0D0 + (z2 - z1)**2.0D0)

    !   deallocate (rlat1, lon1, rlat2, lon2, x1, y1, z1, x2, y2, z2)
  END SUBROUTINE spatial_distance_2d

END MODULE
