close all;
clear all;

filename = '/Users/qzl/sources/MOTOR/input/220527_0000_3km/output/MOTOR-3DVar_bak_G07.nc';

vardata = ncread(filename,'qvapor');
lat = ncread(filename,'lat');
lon = ncread(filename,'lon');

[lon,lat] = meshgrid(lon,lat); 
data = vardata(:,:,1,1);

data = double(data-mean(data));
pcolor(lon,lat,data.');
shading flat;
colormap("jet");

% 在data中随机抽取200个点,将这200个点的值设置为0
% 为了方便观察,将这200个点的位置打印出来

dataObs = zeros(size(data));
dataObsIdx = zeros(size(data));

for i = 1:200
    x = randi([1, size(data,1)]);
    y = randi([1, size(data,2)]);
    dataObs(x,y) = data(x,y);
    dataObsIdx(x,y) = 1;
    fprintf("x=%d, y=%d\n", x, y);
end

% 对data和dataObs进行滤波
[b, a] = butter(3, 0.03, 'low');

dataFilt = filtfilt(b, a, data);
dataFilt = filtfilt(b, a, dataFilt.');
dataFilt = dataFilt.';

dataObs = filtfilt(b, a, dataObs*20);
dataObs = filtfilt(b, a, dataObs.');
dataObs = dataObs.';

figure(1)
pcolor(lon,lat,dataFilt.');
shading flat;
colormap("jet");
caxis([-4e-3,4e-3]);

figure(2)
pcolor(lon,lat,dataObs.');
shading flat;
colormap("jet");
caxis([-4e-3,4e-3]);






