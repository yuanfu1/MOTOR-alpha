close all;
clear all;

for ichan=2:7
 varname=['fy4_1_agri_ch_',num2str(ichan,'%01d'), '_tbb']
 vardata = ncread('D:\YJZX\MOTOR\MOTOR\MOTOR\output\MOTOR-NoJc_inv_BC_G08.nc',varname);
 data = vardata(vardata~=0);
 % hist(data,100);
 [a,b] = hist(data,100);
 stdout=std(data);
 meanout=mean(data);
 bar(b, log10(a));
 title (['BC\_fy4\_1\_agri\_ch\_',num2str(ichan,'%01d'), '\_tbb',', std = ', num2str(stdout, '%4.2f'), ', mean = ', num2str(meanout, '%4.2f')]);
 saveas(gcf, ['D:\YJZX\tmp\BC_',varname], 'png');
 close;
end

for ichan=2:7
 varname=['fy4_1_agri_ch_',num2str(ichan,'%01d'), '_tbb']
 vardata = ncread('D:\YJZX\MOTOR\MOTOR\MOTOR\output\MOTOR-obs_noBC_G08.nc',varname);
 data1 = vardata(vardata~=0);
 vardata = ncread('D:\YJZX\MOTOR\MOTOR\MOTOR\output\MOTOR-NoJc_bak_G08.nc',varname);
 data2 = vardata(vardata~=0);
 % hist(data,100);
 data=data1-data2;
 [a,b] = hist(data,100);
 stdout=std(data);
 meanout=mean(data);
 bar(b, log10(a));
 title (['noBC\_fy4\_1\_agri\_ch\_',num2str(ichan,'%01d'), '\_tbb',', std = ', num2str(stdout, '%4.2f'), ', mean = ', num2str(meanout, '%4.2f')]);
 saveas(gcf, ['D:\YJZX\tmp\noBC_',varname], 'png');
 close;
 [a,b] = hist(data-meanout,100);
 bar(b, log10(a));
 stdout=std(data-meanout);
 meanout=mean(data-meanout);
 title (['noBC1\_fy4\_1\_agri\_ch\_',num2str(ichan,'%01d'), '\_tbb',', std = ', num2str(stdout, '%4.2f'), ', mean = ', num2str(meanout, '%4.2f')]);
 saveas(gcf, ['D:\YJZX\tmp\noBC1_',varname], 'png');
 close;
 [a,b] = hist(data1,100);
 bar(b, log10(a));
 title (['OBS\_fy4\_1\_agri\_ch\_',num2str(ichan,'%01d'), '\_tbb']);
 saveas(gcf, ['D:\YJZX\tmp\OBS_',varname], 'png');
 close
 [a,b] = hist(data2,100);
 bar(b, log10(a));
 title (['BAK\_fy4\_1\_agri\_ch\_',num2str(ichan,'%01d'), '\_tbb',]);
 saveas(gcf, ['D:\YJZX\tmp\BAK_',varname], 'png');
 close;
end

varname=['tbb01']
vardata = ncread('D:\YJZX\MOTOR\MOTOR\MOTOR\output\Raw_FY4_AGRI.nc',varname);
data = vardata(vardata~=0);
% hist(data,100);
bar(b, log10(a));
title (['Raw\_tbb\_ch8']);
saveas(gcf, ['D:\YJZX\tmp\Raw_tbb_ch8_',varname], 'png');
close;