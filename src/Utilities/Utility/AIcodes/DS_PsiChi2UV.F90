program limited_area_wind_decomposition
  implicit none
  integer, parameter :: nlat = 100, nlon = 100  ! 纬度和经度网格大小
  real(8), parameter :: pi = 3.14159265358979323846
  real(8), parameter :: R = 6371000.0           ! 地球半径 (单位：米)
  real(8) :: u(nlat, nlon), v(nlat, nlon)       ! 风场分量 (u, v)
  real(8) :: psi(nlat, nlon), chi(nlat, nlon)   ! 流函数和势函数
  real(8) :: lat(nlat), lon(nlon)               ! 纬度和经度
  real(8) :: lat_start, lat_end, lon_start, lon_end  ! 有限区域的经纬度范围
  real(8) :: dlat, dlon                         ! 纬度和经度间隔
  integer :: i, j

  ! 设置有限区域的经纬度范围
  lat_start = 20.0 * pi / 180.0  ! 起始纬度 (20°N)
  lat_end = 50.0 * pi / 180.0    ! 结束纬度 (50°N)
  lon_start = 100.0 * pi / 180.0 ! 起始经度 (100°E)
  lon_end = 150.0 * pi / 180.0   ! 结束经度 (150°E)

  ! 初始化网格
  dlat = (lat_end - lat_start) / (nlat - 1)  ! 纬度间隔 (弧度)
  dlon = (lon_end - lon_start) / (nlon - 1)  ! 经度间隔 (弧度)
  do i = 1, nlat
      lat(i) = lat_start + (i-1)*dlat  ! 纬度从 lat_start 到 lat_end
  end do
  do j = 1, nlon
      lon(j) = lon_start + (j-1)*dlon  ! 经度从 lon_start 到 lon_end
  end do

  ! 读取或生成风场数据 (u, v)
  call generate_wind_field(u, v, lat, lon, nlat, nlon)

  ! 计算流函数和势函数
  call compute_stream_function(u, v, lat, lon, nlat, nlon, psi)
  call compute_potential_function(u, v, lat, lon, nlat, nlon, chi)

  ! 输出结果
  print *, "流函数和势函数计算完成！"
  print *, "流函数 psi(1,1) = ", psi(1,1)
  print *, "势函数 chi(1,1) = ", chi(1,1)

contains

  ! 生成示例风场数据
  subroutine generate_wind_field(u, v, lat, lon, nlat, nlon)
      integer, intent(in) :: nlat, nlon
      real(8), intent(out) :: u(nlat, nlon), v(nlat, nlon)
      real(8), intent(in) :: lat(nlat), lon(nlon)
      integer :: i, j

      do i = 1, nlat
          do j = 1, nlon
              u(i,j) = sin(lat(i)) * cos(lon(j))  ! 示例 u 分量
              v(i,j) = cos(lat(i)) * sin(lon(j))  ! 示例 v 分量
          end do
      end do
  end subroutine generate_wind_field

  ! 计算流函数
  subroutine compute_stream_function(u, v, lat, lon, nlat, nlon, psi)
      integer, intent(in) :: nlat, nlon
      real(8), intent(in) :: u(nlat, nlon), v(nlat, nlon)
      real(8), intent(in) :: lat(nlat), lon(nlon)
      real(8), intent(out) :: psi(nlat, nlon)
      real(8) :: dlat, dlon, integral
      integer :: i, j

      dlat = lat(2) - lat(1)
      dlon = lon(2) - lon(1)

      ! 通过积分计算流函数
      do i = 1, nlat
          do j = 1, nlon
              integral = sum(v(1:i, j)) * dlat
              psi(i,j) = -R * integral
          end do
      end do
  end subroutine compute_stream_function

  ! 计算势函数
  subroutine compute_potential_function(u, v, lat, lon, nlat, nlon, chi)
      integer, intent(in) :: nlat, nlon
      real(8), intent(in) :: u(nlat, nlon), v(nlat, nlon)
      real(8), intent(in) :: lat(nlat), lon(nlon)
      real(8), intent(out) :: chi(nlat, nlon)
      real(8) :: dlat, dlon, integral
      integer :: i, j

      dlat = lat(2) - lat(1)
      dlon = lon(2) - lon(1)

      ! 通过积分计算势函数
      do i = 1, nlat
          do j = 1, nlon
              integral = sum(u(i, 1:j)) * dlon
              chi(i,j) = R * integral
          end do
      end do
  end subroutine compute_potential_function

end program limited_area_wind_decomposition
