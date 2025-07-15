% close all;
clear all;

close all;
clear all;

%%
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

filenname = '../../../input/Obs/Z_RADR_I_Z9200_20220217060001_O_DOR_SAD_CAP_FMT.nc';
%%
radialAzim  = ncread(filenname,'azimuthV');
radialElev  = ncread(filenname,'elevationV');
distanceZ  = ncread(filenname,'distanceV')/1000;
Z  = ncread(filenname,'RadialVelocity');

%%
% Z(Z == 0) = 66;
Z(Z< 0) = Z(Z<0)+256;
% Z = (Z-66)/2;
Z(Z==0) = NaN;
Z = Z*0.5-64.5;

%%
figure()
for elev = 1:9
    figure(elev);
    for i = 1:1:365
        scatter(distanceZ*sind(radialAzim(i, elev)), distanceZ*cosd(radialAzim(i, elev)), 3, Z(:, i, elev), 'filled'); hold on;
    end

    subtitle(['elevation ',num2str(elev)]);
    axis equal;


    box on;
    grid on;
    hold off
    ccc = 1;
% % %     colormap(radarMap);
%     caxis([0 80]);
end



