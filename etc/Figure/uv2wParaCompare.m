%%
close all;
close all;

wwnd_s  = ncread('/Users/qzl/source/MOTOR/output/Test_UV2WParaTest_serial.nc','wwnd');
wwnd_p  = ncread('/Users/qzl/source/MOTOR/output/Test_UV2WParaTest_parallel.nc','wwnd');
wwnd_diff = wwnd_s-wwnd_p;

uwnd_s  = ncread('/Users/qzl/source/MOTOR/output/Test_UV2WParaTest_serial.nc','uwnd');
uwnd_p  = ncread('/Users/qzl/source/MOTOR/output/Test_UV2WParaTest_parallel.nc','uwnd');
uwnd_diff = uwnd_s-uwnd_p;

vwnd_s  = ncread('/Users/qzl/source/MOTOR/output/Test_UV2WParaTest_serial.nc','vwnd');
vwnd_p  = ncread('/Users/qzl/source/MOTOR/output/Test_UV2WParaTest_parallel.nc','vwnd');
vwnd_diff = vwnd_s-vwnd_p;

height_s  = ncread('/Users/qzl/source/MOTOR/output/Test_UV2WParaTest_serial.nc','height');
height_p  = ncread('/Users/qzl/source/MOTOR/output/Test_UV2WParaTest_parallel.nc','height');
height_diff = height_s-height_p;


sz = size(wwnd_diff);
level = 2;

figure('position', [100,100,1000,800])

gc = subplot(2,2,1)
pcolor(reshape(uwnd_diff(level,:,:),sz(2),sz(3)))
max(max(max(uwnd_diff)))
min(min(min(uwnd_diff)))
colormap(jet);
colorbar;
shading flat;
axis equal;
title('uwnd')

subplot(2,2,2)
pcolor(reshape(vwnd_diff(level,:,:),sz(2),sz(3)))
max(max(max(vwnd_diff)))
min(min(min(vwnd_diff)))
colormap(jet);
colorbar;
shading flat;
axis equal;
title('vwnd')

subplot(2,2,3)
pcolor(reshape(wwnd_diff(level,:,:),sz(2),sz(3)))
max(max(max(wwnd_diff)))
min(min(min(wwnd_diff)))

colormap(jet);
colorbar;
shading flat;
axis equal;
title('wwnd')

subplot(2,2,4)
pcolor(reshape(height_diff(level,:,:),sz(2),sz(3)))
max(max(max(height_diff)))
min(min(min(height_diff)))
colormap(jet);
colorbar;
shading flat;
axis equal;
title('height')

