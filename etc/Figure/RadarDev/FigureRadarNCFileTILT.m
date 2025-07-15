clear all;
close all;

radialAzim  = ncread('../../input/Obs/Z_RADR_I_Z9010_20200105080000_O_DOR_SA_CAP_FMT.nc','radialAzim');

for i = 1:1:11
    plot(unwrap(radialAzim(:,i)),'LineWidth',2); hold on;
end