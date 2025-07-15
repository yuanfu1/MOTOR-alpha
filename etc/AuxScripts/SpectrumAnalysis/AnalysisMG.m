close all;
clear all;

filepath = '/Users/qzl/sources/MOTOR/input/220527_0000_3km/output';
% 文件名列表
filenames = {
    [filepath '/MOTOR-3DVar_diff_G03.nc'],
    [filepath '/MOTOR-3DVar_diff_G04.nc'],
    [filepath '/MOTOR-3DVar_diff_G05.nc'],
    [filepath '/MOTOR-3DVar_diff_G06.nc'],
    [filepath '/MOTOR-3DVar_diff_G07.nc'],
    };
tIdx = 'end';
vIdx = 40;

% 创建一个figure
figure;
set(gcf, 'Position', [100, 100, 900, 900]);

% 初始化颜色限制
clim = [Inf, -Inf];

% 计算所有文件的颜色限制
for i = 1:length(filenames)
    vardata = ncread(filenames{i}, 'temp');
    clim(1) = min(min(clim(1), min(vardata(:, :, vIdx, end))));
    clim(2) = max(max(clim(2), max(vardata(:, :, vIdx, end))));
end
clim = [-2, 3];

% 绘制5张子图
for i = 1:5
    % 读取第一个文件的数据以获取lat和lon
    lat = ncread(filenames{i}, 'lat');
    lon = ncread(filenames{i}, 'lon');
    [lon, lat] = meshgrid(lon, lat);

    vardata = ncread(filenames{i}, 'temp');
    data = vardata(:, :, vIdx, end);
    subplot(3, 3, i);
    pcolor(lon, lat, data.');
    shading flat;
    colormap("jet");
    caxis(clim); % 设置相同的颜色限制
    title(['G0' num2str(i+2)]);
    axis equal;
    colorbar;

end

%% 对比分析这5个文件的频谱
%% 对数据每一行做fft，然后做一个平均，归一化频率

dx = 2000;
% 初始化频谱
subplot(3, 3, [7 8 9]); % 将频谱图放在最后一行
for i = 1:5
    % 读取第一个文件的数据以获取lat和lon
    lat = ncread(filenames{i}, 'lat');
    lon = ncread(filenames{i}, 'lon');
    [lon, lat] = meshgrid(lon, lat);

    vardata = ncread(filenames{i}, 'temp');
    data = vardata(1:end-1, 1:end-1, vIdx, end);

    dx = dx/2;
    Wd = 100000;
    Fs = 1/dx;
    % 对每一行做fft, 然后平均
    L = size(data, 2);
    dataFft = zeros(size(data, 1), L);
    for j = 1:size(data, 1)
        dataFft(j, :) = abs(fft(data(j, :))/L);
    end
    dataFft = mean(dataFft, 1);

    f = Fs/L*(0:(L/2));
    dataFft = dataFft(1:L/2+1);

    plot(f, dataFft, 'LineWidth',2);
    hold on;
    set(gca, 'YScale', 'log')
end
title('Frequency Spectrum');
xlabel('Frequency (Hz)');
ylabel('Amplitude');
hold off;
