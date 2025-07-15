import xarray as xr
import numpy as np
import sys
import matplotlib.pyplot as plt

root = sys.argv[1]
glvl = sys.argv[2].upper()

bcg_ds = xr.open_dataset(f'{root}/MOTOR-3DVar_bak_{glvl}.nc')
# bcg_ds = xr.open_dataset(f'{root}/MOTOR-3DVar_ana_{glvl}.nc.bks4')
old_ds = xr.open_dataset(f'{root}/ctl-T10.9-new/MOTOR-3DVar_ana_{glvl}.nc')
# old_ds = xr.open_dataset(f'{root}/ctl-T10.9-NewRadar/MOTOR-3DVar_ana_{glvl}.nc')
# old_ds = xr.open_dataset(f'{root}/MOTOR-3DVar_ana_G06.nc.o700.1p6.qc5.rnLDV100')
# old_ds = xr.open_dataset(f'{root}/MOTOR-3DVar_ana_G06.nc.oldRadar')
# old_ds = xr.open_dataset(f'{root}/ctl-T10.9-NewRadarNoCons/MOTOR-3DVar_ana_{glvl}.nc')
# 
new_ds = xr.open_dataset(f'{root}/MOTOR-3DVar_ana_{glvl}.nc')
# amo_ds = xr.open_dataset(f'{root}/MOTOR-3DVar_AMO_{glvl}.nc')
era5_ds = xr.open_dataset('/Users/qzl/sources/MOTOR/input/2024110300_T10p1/verify/2024110300.nc')
lats = old_ds.lat
lons = old_ds.lon
times = old_ds.time
press = old_ds.pres
alt = old_ds.height

ea5_gd_u_bg_intp = era5_ds.u.interp(latitude=lats, longitude=lons, valid_time=times[-1], pressure_level=press[-1,:,:,:]/100.0)
ea5_gd_v_bg_intp = era5_ds.v.interp(latitude=lats, longitude=lons, valid_time=times[-1], pressure_level=press[-1,:,:,:]/100.0)

print(press)
print(ea5_gd_u_bg_intp.shape)

era5_ws = (ea5_gd_u_bg_intp**2 + ea5_gd_v_bg_intp**2)**0.5
old_ws = (old_ds.uwnd[-1,:,:,:]**2 + old_ds.vwnd[-1,:,:,:]**2)**0.5
new_ws = (new_ds.uwnd[-1,:,:,:]**2 + new_ds.vwnd[-1,:,:,:]**2)**0.5
bcg_wa = (bcg_ds.uwnd[-1,:,:,:]**2 + bcg_ds.vwnd[-1,:,:,:]**2)**0.5

old_mae = np.nanmean(abs(old_ws - era5_ws), axis=(1, 2))
new_mae = np.nanmean(abs(new_ws - era5_ws), axis=(1, 2))
bcg_mae = np.nanmean(abs(bcg_wa - era5_ws), axis=(1, 2))

# Root Mean Square Error
old_rms = np.nanmean((old_ws - era5_ws)**2, axis=(1, 2))**0.5
new_rms = np.nanmean((new_ws - era5_ws)**2, axis=(1, 2))**0.5
bcg_rms = np.nanmean((bcg_wa - era5_ws)**2, axis=(1, 2))**0.5

# me
old_me = np.nanmean(old_ws - era5_ws, axis=(1, 2))
new_me = np.nanmean(new_ws - era5_ws, axis=(1, 2))
bcg_me = np.nanmean(bcg_wa - era5_ws, axis=(1, 2))

old_ds.close()
new_ds.close()
era5_ds.close()
bcg_ds.close()

x = np.nanmean(press.values[-1,:,:,:]/100.0, axis=(1, 2))

# 创建画布和子图
plt.figure(figsize=(10, 5))  # 设置画布大小

# 第一个子图
plt.subplot(1, 3, 2)  # 1行2列，第1个子图
plt.plot(old_mae, x, color='r', label='old')
plt.plot(new_mae, x, color='b', label='new')
plt.plot(bcg_mae, x, color='g', label='bcg')
plt.title('MAE')
plt.xlabel('wind speed mae(m/s)')
plt.ylabel('pres(hPa)')
plt.legend()
plt.gca().invert_yaxis()
# grid on
plt.grid()

# 第二个子图
plt.subplot(1, 3, 3)  # 1行2列，第2个子图
plt.plot(old_rms, x, color='r', label='old')
plt.plot(new_rms, x, color='b', label='new')
plt.plot(bcg_rms, x, color='g', label='bcg')
plt.title('RMSE')
plt.xlabel('wind speed rmse(m/s)')
# plt.ylabel('pres')
plt.legend()
plt.gca().invert_yaxis()
plt.grid()

# 第三个子图
plt.subplot(1, 3, 1)  # 1行2列，第2个子图
plt.plot(old_me, x, color='r', label='old')
plt.plot(new_me, x, color='b', label='new')
plt.plot(bcg_me, x, color='g', label='bcg')
plt.title('ME')
plt.xlabel('wind speed me(m/s)')
# plt.ylabel('pres')
plt.legend()
plt.gca().invert_yaxis()
plt.grid()


# 显示图形
plt.tight_layout()  # 自动调整子图间距
plt.show()


