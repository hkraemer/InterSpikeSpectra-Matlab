% Test InterSpikeSpectra functionality

clear, clc

cd ../private/
% preconditions
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

cd ../

%% Test 1: 2 Spikes + measurement noise
threshold = 0.85;
tol = 1e-3;
[spectrum, rho] = inter_spike_spectrum(test_tauRR, 'threshold', threshold, 'tol', tol, 'logit', false);
[maxis, max_idx] = findpeaks(spectrum);
assert(abs(rho - threshold) < tol)
assert(sum(max_idx == [period2, period1]) == 2)
assert(0.513 < maxis(1) && maxis(1) < 0.514)
assert(0.486 < maxis(2) && maxis(2) < 0.487)

%% Lasso & Elastic Net
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

[spectrum, ~] = inter_spike_spectrum(test_tauRR, 'logit', false, 'Alpha', 0.5);
[~, max_idx] = findpeaks(spectrum);
assert(length(max_idx) == 4)
assert(max_idx(1)==period2)
assert(max_idx(2)==period1)
assert(max_idx(3)==16)
assert(max_idx(4)==19)


%% Test 2: randomized peak heights
rng(1234)
threshold = 0.95;
tol = 1e-3;
maxcycles = 20;
test_tauRR = abs(randn(1,N)) .* test_tauRR2(1,:);
[spectrum, ~] = inter_spike_spectrum(test_tauRR, 'threshold', threshold, 'tol', tol, 'max_iter', maxcycles, 'logit', false);

[maxis, max_idx] = findpeaks(spectrum);
t_idx = maxis > 0.000001;
peak_idxs = max_idx(t_idx);
assert(sum(rem(peak_idxs, period1)) == 0)

%% Test 3: Perfect spike train
s = zeros(1,100);
period1 = 3;
s(2:period1:end) = 1;

[spectrum, rho] = inter_spike_spectrum(s, 'logit', false);

assert(max(spectrum)>0.999)
assert(rho+ 1e-3 >= 1)
[~, max_idx] = findpeaks(spectrum);
assert(length(max_idx) == 1)
assert(max_idx(1) == period1)

period2 = 11;
s(5:period2:end) = .8;

tol = 1e-4;
threshold = .99;
[spectrum, rho] = inter_spike_spectrum(s, 'threshold', threshold, 'tol', tol, 'logit', false);
assert(abs(rho - threshold) < tol)
[maxis, max_idx] = findpeaks(spectrum);
assert(length(max_idx) == 2)
assert(max_idx(1) == period1)
assert(max_idx(2) == period2 || max_idx(2) / max_idx(1) == period2)
assert(maxis(2) > maxis(1))

%% Test 4: Random time series
rng(1234)
tol = 1e-4;
maxcycles = 20;
s = randn(1,50);

threshold = 0.995;
[spectrum1, rho] = inter_spike_spectrum(s, 'threshold', threshold, 'tol', tol, 'max_iter', maxcycles, 'logit', false);
[maxis, ~] = findpeaks(spectrum1);
numpeaks1 = length(maxis);
assert(abs(rho - threshold) <= tol)

threshold = 0.85;
[spectrum2, rho] = inter_spike_spectrum(s, 'threshold', threshold, 'logit', false);
[maxis, ~] = findpeaks(spectrum2);
numpeaks2 = length(maxis);
assert(abs(rho - threshold) <= 1e-3)
assert(numpeaks2>numpeaks1)
assert(numpeaks2 == 8)
assert(numpeaks1 == 4)
assert(0.1013 < maxis(1) && maxis(1) < 0.1014)
assert(0.098 < maxis(2) && maxis(2) < 0.099)

%% Test 5: Compare normal to logistic regression

N = 300;
M = ceil(N/2);
period1 = 3;
period2 = 7;
s = zeros(1,N);
s(period1:period1:end) = 1;
s(period2:period2:end) = 1;

threshold = 0.99;
[spectrum1, rho1] = inter_spike_spectrum(s, 'threshold', threshold);
[spectrum2, rho2] = inter_spike_spectrum(s, 'threshold', threshold, 'logit', false);

[peaks1, peaks1_idx] = findpeaks(spectrum1);
[peaks2, peaks2_idx] = findpeaks(spectrum2);

assert(abs(rho1 - threshold) <= 1e-3)
assert(abs(rho2 - threshold) <= 1e-3)
assert(0.5633 < peaks1(1) && peaks1(1) < 0.5634)
assert(0.4366 < peaks1(2) && peaks1(2) < 0.4367)
assert(0.332 < peaks2(1) && peaks2(1) < 0.334)
assert(0.665 < peaks2(2) && peaks2(2) < 0.667)
assert(sum(peaks1_idx == [3,7])==2)
assert(sum(peaks2_idx == [3,21])==2)

threshold = 0.9;
[spectrum1b, rho1b] = inter_spike_spectrum(s, 'threshold', threshold);
[spectrum2b, rho2b] = inter_spike_spectrum(s, 'threshold', threshold, 'logit', false);

[peaks1, peaks1_idx] = findpeaks(spectrum1b);
[peaks2, peaks2_idx] = findpeaks(spectrum2b);

assert(abs(rho1b - threshold) <= 1e-3)
assert(abs(rho2b - threshold) <= 1e-3)
assert(0.5633 < peaks1(1) && peaks1(1) < 0.5634)
assert(0.4366 < peaks1(2) && peaks1(2) < 0.4367)
assert(0.5633 < peaks2(1) && peaks2(1) < 0.5634)
assert(0.4366 < peaks2(2) && peaks2(2) < 0.4367)
assert(sum(peaks1_idx == [3,7])==2)
assert(sum(peaks2_idx == [3,7])==2)
