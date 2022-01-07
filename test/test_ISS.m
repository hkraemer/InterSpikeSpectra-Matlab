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
tol = 1e-2;
[~, rho] = inter_spike_spectrum(test_tauRR, threshold, tol);
assert(abs(rho - threshold) < tol)

threshold = 0.99;
[spectrum, rho] = inter_spike_spectrum(test_tauRR, threshold, tol);

[maxis, max_idx] = findpeaks(spectrum);
t_idx = maxis > 0.1;
peak_idxs = max_idx(t_idx);

assert(abs(rho - threshold) < tol)
assert(length(peak_idxs) == 2)
assert(peak_idxs(1) == period2)
assert(peak_idxs(2) == period1)

%% Test 2: randomized peak heights
rng(1234)
threshold = 0.95;
tol = 1e-1;
maxcycles = 100;
test_tauRR = abs(randn(1,N)) .* test_tauRR2(1,:);
[spectrum, ~] = inter_spike_spectrum(test_tauRR, threshold, tol, maxcycles);

[maxis, max_idx] = findpeaks(spectrum);
t_idx = maxis > 0.001;
peak_idxs = max_idx(t_idx);

assert(sum(rem(peak_idxs, period1)) == 0)

%% Test 3: Random time series
rng(1234)
tol = 1e-4;
maxcycles = 100;
s = randn(1,50);

threshold = 0.995;
[spectrum1, rho] = inter_spike_spectrum(s, threshold, tol, maxcycles);
[maxis, ~] = findpeaks(spectrum1);
numpeaks1 = length(maxis);
assert(abs(rho - threshold) <= tol)

threshold = 0.85;
[spectrum2, rho] = inter_spike_spectrum(s, threshold);
[maxis, ~] = findpeaks(spectrum2);
numpeaks2 = length(maxis);
assert(abs(rho - threshold) <= 1e-3)
assert(numpeaks2<numpeaks1)
