% 读取日志文件
logFile = '/Users/qzl/sources/MOTOR/build/1.log';
fid = fopen(logFile, 'r');
if fid == -1
    error('无法打开日志文件');
end

% 初始化变量
levels = {'G3', 'G4', 'G5', 'G6', 'G7'};
jo_values = containers.Map();
jb_values = containers.Map();
for i = 1:length(levels)
    jo_values(levels{i}) = [];
    jb_values(levels{i}) = [];
end

% 读取日志文件内容
while ~feof(fid)
    line = fgetl(fid);
    for i = 1:length(levels)
        if contains(line, levels{i})
            % processLevel(fid, levels{i}, jo_values(levels{i}), jb_values(levels{i}));

                line = fgetl(fid);
    if contains(line, 'Jb') && contains(line, 'Jo')
        disp(['Processing line: ', line]);
        tokens = regexp(line, 'J-G: Jb\s+(\d+\.\d+D[+-]\d+)\s+Jo\s+(\d+\.\d+D[+-]\d+)', 'tokens');
        if ~isempty(tokens)
            jb_value = str2double(strrep(tokens{1}{1}, 'D', 'E'));
            jo_value = str2double(strrep(tokens{1}{2}, 'D', 'E'));
            jb_values(levels{i}) = [jb_values(levels{i}); jb_value];
            jo_values(levels{i}) = [jo_values(levels{i}); jo_value];
        end
    end
        end
    end
end
fclose(fid);

% 绘制 Jo 和 Jb 的值变化
figure;
for i = 1:length(levels)
    plotLevel(i, levels{i}, jo_values(levels{i}), jb_values(levels{i}));
end
set(gcf, 'Position', [100, 100, 900, 900]);

function plotLevel(position, level, jo_values, jb_values)
    subplot(3, 2, position);
    plot(jo_values, '-o');
    hold on;
    plot(jb_values, '-x');
    xlabel('迭代次数');
    ylabel('值');
    legend('Jo', 'Jb');
    title([level, ' 上 Jo 和 Jb 的值变化']);
    grid on;
end