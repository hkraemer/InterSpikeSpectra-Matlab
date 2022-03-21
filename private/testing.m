clc, clear

N = 300;
M = ceil(N/2);
period1 = 3;
period2 = 7;
period3 = 15;
% period4 = 30;
s = zeros(1,N);
s(1:period1:end) = 1;
s(3:period2:end) = 1;
%s(9:period3:end) = 1;
%s(8:period4:end) = 1;
% s1 = zeros(1,N);
% s2 = zeros(1,N);
% s1(period1:period1:end) = 1;
% s2(period2:period2:end) = 1;
% s = s1 + s2;

rng(1)
s = s + (0.001.*randn(1,N));

threshold = 0.9;
alpha = 0.8;

tic
[spectrum1, rho1] = inter_spike_spectrum(s, 'threshold', threshold, 'logit', false, 'Alpha', alpha);
toc

tic
[spectrum2, rho2] = inter_spike_spectrum(s, 'threshold', threshold, 'logit', true, 'Alpha', alpha);
toc

tic
[spectrum3, rho3] = inter_spike_spectrum(s, 'threshold', threshold, 'method', "STLS");
toc


figure
subplot(411)
plot(1:length(s),s)
title("Spike train input")
grid on
subplot(412)
plot(1:M,spectrum1)
title("LASSO normal")
grid on
subplot(413)
plot(1:M,spectrum2)
title("LASSO logit")
grid on
subplot(414)
plot(1:M,spectrum3)
title("STLS")
grid on

%%


%%
clear, clc
rng(1234)
N = 50;
period1 = 13;
period2 = 8;
test_tauRR2 = create_single_basis_function(N, period1);
test_tauRR3 = create_single_basis_function(N, period2);
test_tauRR = 0.7 .* test_tauRR2(3,:);
test_tauRR_ = 0.8 .* test_tauRR3(7,:);
test_tauRR = test_tauRR + test_tauRR_;
test_tauRR = test_tauRR + 0.05.*randn(1,N);

%%
threshold = 0.99;
[spectrum, rho] = inter_spike_spectrum(test_tauRR, 'logit', false);
[maxis, max_idx] = findpeaks(spectrum);
t_idx = maxis > 0.1;
peak_idxs = max_idx(t_idx);

assert(length(max_idx) == 3)
assert(abs(rho - threshold) < 1e-3)
assert(length(peak_idxs) == 2)
assert(peak_idxs(1) == period2)
assert(peak_idxs(2) == period1)

[spectrum, rho] = inter_spike_spectrum(test_tauRR, 'logit', false, 'Alpha', 0.5);
[maxis, max_idx] = findpeaks(spectrum);
assert(length(max_idx) == 4)
assert(max_idx(1)==period2)
assert(max_idx(2)==period1)
assert(max_idx(3)==16)
assert(max_idx(4)==19)
