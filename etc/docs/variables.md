# Variable names & units

All the variables that will be used in MOTOR-DA must be specified in `template.nl: modelState: varList`. If not, you'll see an error message about the failure. 

### varList specified in template.nl
Note: obtain your variables by matching `type`, `name`, and `data` of the `Field` type.

Control variables in the 3D analysis
``` 
(1)  pres, "air_pressure", unit: Pa
(2)  temp, "air_temperature", unit: K
(3)  uwnd, "x component of 3D wind", unit: m/s
(4)  vwnd, "y component of 3D wind", unit: m/s
(5)  wwnd, "z component of 3D wind", unit: m/s
(6)  qvapor, "specific_humidity of water vapor", unit: kg/kg 
(7)  qcloud, "specific_humidity of cloud water", unit: kg/kg 
(8)  qrain, "specific_humidity of rain water", unit: kg/kg 
(9)  qice, "specific_humidity of cloud ice", unit: kg/kg 
(10)  qsnow, "specific_humidity of snow", unit: kg/kg 
(11)  qgraup, "specific_humidity of graupel", unit: kg/kg 
(12)  rho, "density of dry air", unit: kg/m3
(13)  rhov, "density of vapor", unit: kg/m3
(14)  rhor, "density of rain water", unit: kg/m3
(15)  rhoc, "density of cloud water", unit: kg/m3
(16)  rhoi, "density of ice water", unit: kg/m3
(17)  rhoh, "density of hail", unit: kg/m3
(18)  rhog, "density of graupel", unit: kg/m3
```

Variables used in the 3D analysis, but they are not control variables
```
(1)  qrain, "specific_humidity of rain water", unit: kg/kg 
```

Control variables in the 2D surface analysis
```
(1)  psfc, "surface_pressure", unit: Pa
(2)  tsfc, "surface_temperature", unit: K
(3)  t2m, "temperature_at_two_meters_above_surface", unit: K
(4)  q2m, "specific_humidity_at_two_meters_above_surface", unit: kg/kg
(5)  u10m, "uwind_at_10m", unit: m/s
(6)  v10m, "vwind_at_10m", unit: m/s
```

Parameters which are not belonging to the variables 
```
(1) ztopo, "terrain_height", unit: m
(2) tskin, "skin_temperature", unit: K
(3) seaice, "seaice_fraction", unit: 1
(4) landmask, "landmask", unit: 1
(5) sfctype, "surface_type", unit: 1
(6) vegtype, "veg_type", unit: 1
(7) landtype, "land_type", unit: 1
(8) leaf, "leaf_area_index", unit: 1
(9) tsoil, "soil_temperature", unit: K
(10) qsoil, "soil_moisture", unit: kg/kg
(11) soiltype, "soil_type", unit: 1
```

Observation variables
```
(1) ref, 'Reflectivity in the radar.', type: RADAR, unit: dbZ or mm6/mm3
(2) rwnd, 'Radial wind in the radar.', type: RADAR, unit: m/s
```