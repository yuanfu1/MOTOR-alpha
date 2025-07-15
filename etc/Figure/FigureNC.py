# Created by Zilong Qin(zilong.qin@gmail.com), 2021/7/8, @ GBA-MWF, Shenzhen
import netCDF4 as nc
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

filename = '/Users/qzl/source/MOTOR/input/220602_0730/output/testMO_ana_temp_G10.nc'   # .nc文件名
f = nc.Dataset(filename)

ua = np.array(f['temp'][:])  # 转化为np.array数组
lat = np.array(f['lat'][:])  # 转化为np.array数组
lon = np.array(f['lon'][:])  # 转化为np.array数组
level = 4

# fig = plt.figure()
# # ax = fig.add_subplot(111, projection='3d')

# x, y = np.meshgrid(lon, lat)

# # # Plot the surface
# # ax.plot_surface(x, y, ua[level, :, :], cmap=cm.coolwarm, antialiased=True)
# # plt.show()
# print(ua.shape)

# ax = plt.axes(projection='3d')
# # ax.contour3D(x, y, ua[0, :, :, 0], 50, cmap=cm.coolwarm)
# ax.plot_surface(x, y, ua[0, :, :, 0], cmap=cm.coolwarm, antialiased=True)
# ax.set_xlabel('x')
# ax.set_ylabel('y')
# ax.set_zlabel('z')
# plt.show()

fig, (ax0) = plt.subplots(1, 1)

c = ax0.pcolor(lon, lat, ua[0, :, :, 0], cmap='jet', vmin=0, vmax=5)
ax0.set_title('Without /h^2, scale factor = 1.0')

ax0.set_xlabel('longitude')
ax0.set_ylabel('latitude')

fig.tight_layout()
plt.show()

