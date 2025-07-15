% Created by Zilong Qin (zilong.qin@gmail.com), 2022/2/14, @GBA-MWF, Shenzhen
a1 = 44300;
a2 = 5.256;

H = 0:100:30e3;
P = 100000*(1-H/a1).^(a2);

S1 = 10.^(-12./(log10(P))+2.4);

% plot(P, H);
plot(S1, H, 'LineWidth',3);
grid on;
set(gca,'FontSize',20)
set(gca,'FontName','times')