close all;
clear all;


filename = '/MOTOR-SfcAna_obsThinned_G06.nc';
var = 'temp';
pres1  = ncread(['/Users/qzl/source/MOTOR/input/220602_0730/output/1/', filename] ,var);
pres4  = ncread(['/Users/qzl/source/MOTOR/input/220602_0730/output/4/', filename] ,var);

presdiff = pres1-pres4;

pcolor(presdiff(:,:,1,3))

max(max(max(max(presdiff(:,:,:,:)))))
min(min(min(min(presdiff(:,:,:,:)))))

size(presdiff)