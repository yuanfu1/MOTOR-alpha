% close all;
clear all;


latRef  = ncread('../../../output/testRadar.nc','latRef');
lonRef  = ncread('../../../output/testRadar.nc','lonRef');
altRef  = ncread('../../../output/testRadar.nc','altRef');
timRef  = ncread('../../../output/testRadar.nc','timRef');
valRef  = ncread('../../../output/testRadar.nc','valRef');

figure()
scatter3(lonRef,latRef,altRef,5,valRef);

radarMap = [
    255  255  255;
    0  255  255;
    0  157  255;
    0    0  255;
    9  130  175;
    0  255    0;
    8  175   20;
    255  214    0;
    255  152    0;
    255    0    0;
    221    0   27;
    188    0   54;
    121    0  109;
    121   51  160;
    195  163  212;
    255  255  255];
radarMap = uint8(radarMap);

colormap(radarMap);
caxis([0 80]);