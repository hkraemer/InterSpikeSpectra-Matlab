clc, clear

%%

clear, clc

x = 0.7;
a=3.55;
for i = 2:200;x(i)=a*x(i-1)*(1-x(i-1));end
[spectrum, rho] = inter_spike_spectrum(x(40:end), 'threshold', 0.99);
area(spectrum)

figure
subplot(211)
plot(1:length(x),x)
grid on
subplot(212)
plot(1:length(x(40:end)),spectrum)
grid on

%%
clear, clc

N=200;
period1 = 22;
period2 = 6;
period3 = 9;
s = zeros(1,N);
s(3:period1:end) = 1;
s(1:period2:end) = 1;
s(4:period3:end) = 1;

s = s + (0.05.*randn(1,N));

[spectrum, rho] = inter_spike_spectrum(s, 'threshold', 0.9);
[spectrum2, rho2] = inter_spike_spectrum(s);
[spectrum3, rho3] = inter_spike_spectrum(s, 'threshold', 0.99);


figure
subplot(231)
plot(1:length(s),s)
grid on
subplot(232)
plot(1:length(s),s)
grid on
subplot(233)
plot(1:length(s),s)
grid on
subplot(234)
plot(1:length(s),spectrum)
title("threshold 0.9")
grid on
subplot(235)
plot(1:length(s),spectrum2)
title("threshold 0.95 (Default)")
grid on
subplot(236)
plot(1:length(s),spectrum3)
title("threshold 0.99")
grid on