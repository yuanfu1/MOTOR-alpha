% Created by Zilong Qin (zilong.qin@gmail.com), 2022/2/18, @GBA-MWF, Shenzhen

close all;
clear all;


latRef  = ncread('../../output/testRadar.nc','latRef');
timRef  = ncread('../../output/testRadar.nc','timRef');
lonRef  = ncread('../../output/testRadar.nc','lonRef');
valRef  = ncread('../../output/testRadar.nc','valRef');
altRef  = ncread('../../output/testRadar.nc','altRef');


scatter3(latRef, lonRef, altRef, 3, valRef);