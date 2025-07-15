import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

# 读取数据
ds = xr.open_dataset("MOTOR-3DVar_innovation_G06.nc")
temperature = ds["fy4_1_agri_ch_2_tbb"][-1,0,:,:].values.flatten()
temperature = temperature[~np.isnan(temperature)]
print(np.mean(temperature))
print(np.size(temperature))

# 计算 KDE
kde = gaussian_kde(temperature, bw_method=0.4)
x = np.linspace(temperature.min(), temperature.max(), 100)
pdf = kde(x)

# 绘图
plt.figure(figsize=(10, 6))

# 绘制直方图（可选）
plt.hist(
    temperature,
    bins=40,
    density=True,
    alpha=0.3,
    color="lightblue",
    edgecolor="gray"
)

# 绘制 KDE 曲线
plt.plot(x, pdf, color="navy", linewidth=2, label="PDF")

# 添加标注
plt.axvline(
    np.mean(temperature),
    color="red",
    linestyle="--",
    label=f"Mean = {np.mean(temperature):.2f}K"
)

plt.title("PDF", fontsize=14)
plt.xlabel("TBB (K)", fontsize=12)
plt.ylabel("Density", fontsize=12)
plt.legend()
plt.grid(alpha=0.2)
plt.show()