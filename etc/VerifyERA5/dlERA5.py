"""
ERA5数据下载脚本
本脚本用于从ECMWF的CDS API下载ERA5再分析数据, 用于验证MOTOR-DA的同化效果。

主要功能:
1. 下载指定日期时间的ERA5再分析数据
2. 支持下载的变量包括:
   - 位势高度(geopotential)
   - 比湿(specific_humidity) 
   - 温度(temperature)
   - U风场(u_component_of_wind)
   - V风场(v_component_of_wind)
   - 垂直速度(vertical_velocity)
3. 数据区域为中国及周边地区(50°N-0°N, 90°E-135°E)
4. 支持37个标准气压层

使用说明:
1. 需要先配置CDS API密钥
2. 调用download_era5_data函数下载数据
3. 参数date_time为datetime对象,指定下载数据的日期时间
4. 参数filePath为数据保存路径

作者: Zilong Qin
日期: 2025-04-09
版本: 1.0
"""

import cdsapi
from datetime import datetime
import pandas as pd

def download_era5_data(date_time, filePath):
    """
    下载ERA5再分析数据
    
    参数:
        date_time: datetime类型，指定要下载数据的日期和时间
    """
    # 从datetime对象提取年月日和时间
    year = date_time.strftime("%Y")
    month = date_time.strftime("%m")
    day = date_time.strftime("%d")
    time = date_time.strftime("%H:%M")
    
    dataset = "reanalysis-era5-pressure-levels"
    request = {
        "product_type": ["reanalysis"],
        "variable": [
            "geopotential",
            "specific_humidity",
            "temperature",
            "u_component_of_wind",
            "v_component_of_wind",
            "vertical_velocity"
        ],
        "year": [year],
        "month": [month],
        "day": [day],
        "time": [time],
        "pressure_level": [
            "1", "2", "3",
            "5", "7", "10",
            "20", "30", "50",
            "70", "100", "125",
            "150", "175", "200",
            "225", "250", "300",
            "350", "400", "450",
            "500", "550", "600",
            "650", "700", "750",
            "775", "800", "825",
            "850", "875", "900",
            "925", "950", "975",
            "1000"
        ],
        "data_format": "netcdf",
        "download_format": "unarchived",
        "area": [50, 90, 0, 135]
    }

    output_file = f"{filePath}/era5_{year}{month}{day}_{time.replace(':', '')}.nc"
    
    client = cdsapi.Client()
    client.retrieve(dataset, request).download(output_file)
    
    return output_file


# 示例使用方法
if __name__ == "__main__":
    # 下载2024年6月1日12:00的数据
    start_time = datetime(2025, 4, 27, 12, 0)
    end_time = datetime(2025, 4, 27, 12, 0)
    
    filePath = '/home/mgq/proj/MOTOR/era5'
    
    for target_date in pd.date_range(start_time, end_time, freq="12H"):
        file_path = download_era5_data(target_date, filePath)
        print(f"数据已下载到: {file_path}")