clc, clear


N = 300;
M = ceil(N/2);
period1 = 3;
period2 = 7;
% period3 = 13;
% period4 = 30;
s = zeros(1,N);
s(2:period1:end) = 1;
s(1:period2:end) = 1;
%s(9:period3:end) = 1;
%s(8:period4:end) = 1;

rng(1)
s = s + (0.00005.*randn(1,N));

threshold = 0.9;

tic
[spectrum1, rho1] = inter_spike_spectrum(s, 'threshold', threshold);
toc

tic
[spectrum2, rho2] = inter_spike_spectrum(s, 'threshold', threshold, 'method', "STLS");
toc


figure
subplot(311)
plot(1:length(s),s)
grid on
subplot(312)
plot(1:M,spectrum1)
grid on
subplot(313)
plot(1:M,spectrum2)
grid on

%%


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