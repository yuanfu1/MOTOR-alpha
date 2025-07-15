% Created by Zilong Qin (zilong.qin@gmail.com), 2020/12/22, @GBA-MWF, Shenzhen

close all;
clear all;

d(1) = 2^(0.5)/2;
d(2) = 3/2*2^(0.5);
d(3) = (2.5)^(0.5);
d(4) = (2.5)^(0.5);

sum_coef = sum(1./(d.^2));

c = 1./(d.^2)./sum_coef