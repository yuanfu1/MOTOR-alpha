close all;
clear all;

%% Read the analysis field of both lnp and pres_ctl results
lnp_ana  = ncread('/Users/qzl/sources/MOTOR/input/220527_0000_3km/output/MOTOR-3DVar_ana_G07.nc','lnp');
pres_ctl_ana  = ncread('/Users/qzl/sources/MOTOR/input/220527_0000_3km/output/Alpha-0.1.3/MOTOR-3DVar_ana_G07.nc','pres_ctl');
lat = ncread('/Users/qzl/sources/MOTOR/input/220527_0000_3km/output/Alpha-0.1.3/MOTOR-3DVar_ana_G07.nc','lat');
lon =  ncread('/Users/qzl/sources/MOTOR/input/220527_0000_3km/output/Alpha-0.1.3/MOTOR-3DVar_ana_G07.nc','lon');

%% Read the background field
pres_bak  = ncread('/Users/qzl/sources/MOTOR/input/220527_0000_3km/output/Alpha-0.1.3/MOTOR-3DVar_bak_G07.nc','pres');

pres_ctl_ana = pres_ctl_ana(:,:,1)*100;
pres_lnp_ana = exp(lnp_ana(:,:,1));
pres_bak = pres_bak(:,:,1);

[lat, lon] = meshgrid(lat, lon);
% 
% figure(1)
% pcolor(lon, lat, pres_ctl_ana-pres_bak);
% shading flat;
% 
% figure(2)
% pcolor(lon, lat, pres_lnp_ana-pres_bak);
% shading flat;

%% Figure
figure(3)
pcolor(lon, lat, pres_ctl_ana-pres_lnp_ana);
shading flat;
% caxis([-200 200]);
colormap(jet);
colorbar;

title("pres\_ctl - lnp\_1.3");