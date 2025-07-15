close all;
clear all

% 绘制第一个子图
subplot(3,1,1)
lamda = 0:0.0001:0.2;
k = 2*pi*lamda;

H = k.^2;

plot(lamda, H); hold on;
plot(lamda, 4*H);

ylim([0 1]);
% xlim([0 1e5]);

%%  
% 绘制第二个子图
subplot(3,2,3)
N = 3001;
t = 1:1:N;
func = zeros(1, N);
func((N-1)/2+1) = 10;
func((N-1)/4+1) = 10;
func((N-1)*3/4+1) = 10;


plot(t, func); hold on;
xlim([1500 2500]);
xlim([(N-1)/2+1-300 (N-1)/2+1+300]);
title('Initial signal');

% 绘制第三个子图
subplot(3,2,4)

% 进行50次迭代
func = ApplyLaplaceFilter(func, N, 50);
plot(t, (func), 'LineWidth', 2); hold on;

% 进行100次迭代
func = ApplyLaplaceFilter(func, N, 100);
plot(t, (func), 'LineWidth', 2); hold on;

% 进行100次迭代
func = ApplyLaplaceFilter(func, N, 150);
plot(t, (func), 'LineWidth', 2); hold on;

% 进行100次迭代
func = ApplyLaplaceFilter(func, N, 250);
plot(t, (func), 'LineWidth', 2); hold on;

legend('50 iterations', '150 iterations', '250 iterations', '500 iterations');

xlim([(N-1)/2+1-100 (N-1)/2+1+100]);
title('After Laplace filter');


%%  
% 绘制第二个子图
subplot(3,2,5)
N = 3001;
t = 1:1:N;
func = zeros(1, N);

for i = 2:19
    func((N-1)*i/20+1) = 10 - i;
    func((N-1)*i/20+1) = 10 - i;
end


plot(t, func); hold on;
% xlim([(N-1)/2+1-300 (N-1)/2+1+300]);
title('Initial signal');

% 绘制第三个子图
subplot(3,2,6)

% 进行50次迭代
func = ApplyLaplaceFilter(func, N, 50);
plot(t, (func), 'LineWidth', 2); hold on;

% 进行100次迭代
func = ApplyLaplaceFilter(func, N, 100);
plot(t, (func), 'LineWidth', 2); hold on;

% 进行100次迭代
func = ApplyLaplaceFilter(func, N, 150);
plot(t, (func), 'LineWidth', 2); hold on;

% 进行100次迭代
func = ApplyLaplaceFilter(func, N, 250);
plot(t, (func), 'LineWidth', 2); hold on;

legend('50 iterations', '150 iterations', '250 iterations', '500 iterations');

% xlim([(N-1)/2+1-100 (N-1)/2+1+100]);
title('After Laplace filter');

% 定义一个函数来执行FFT和IFFT操作
function func = ApplyLaplaceFilter(func, N, iterations)
    for i = 1:iterations
        fft_func = fft(func);
% for i = 2:20
%     func((N-1)*i/20+1) = 10 - i;
%     func((N-1)*i/20+1) = 10 - i;
% end

        f = 1/N:1/N:1-1/N;
        f = [0, f];
        omega = 2*pi*f;
        fft_func(2:501) = omega(2:501).*omega(2:501).*omega(2:501).*omega(2:501).*fft_func(2:501);
        for i = (N-1)/2+1:N
            fft_func(i) = conj(fft_func(N+1-i+1));
        end
        laps_func = ifft(fft_func);

        % laps_func = func*-2;
        % laps_func = laps_func + circshift(func, 1)+circshift(func, -1);

        func = func - laps_func;
    end
end