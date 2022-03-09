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
[spectrum, rho] = inter_spike_spectrum(test_tauRR, 'threshold', threshold, 'tol', tol);
[maxis, max_idx] = findpeaks(spectrum);
assert(abs(rho - threshold) < tol)
assert(sum(max_idx == [period2, period1]) == 2)
assert(0.513 < maxis(1) && maxis(1) < 0.514)
assert(0.486 < maxis(2) && maxis(2) < 0.487)

%%
threshold = 0.99;
[spectrum, rho] = inter_spike_spectrum(test_tauRR);

[maxis, max_idx] = findpeaks(spectrum);
t_idx = maxis > 0.1;
peak_idxs = max_idx(t_idx);

assert(abs(rho - threshold) < 1e-3)
assert(length(peak_idxs) == 2)
assert(peak_idxs(1) == period2)
assert(peak_idxs(2) == period1)

%% Test 2: randomized peak heights
rng(1234)
threshold = 0.95;
tol = 1e-3;
maxcycles = 20;
test_tauRR = abs(randn(1,N)) .* test_tauRR2(1,:);
[spectrum, ~] = inter_spike_spectrum(test_tauRR, 'threshold', threshold, 'tol', tol, 'max_iter', maxcycles);

[maxis, max_idx] = findpeaks(spectrum);
t_idx = maxis > 0.000001;
peak_idxs = max_idx(t_idx);
assert(sum(rem(peak_idxs, period1)) == 0)

%% Test 3: Perfect spike train
s = zeros(1,100);
period1 = 3;
s(2:period1:end) = 1;

[spectrum, rho] = inter_spike_spectrum(s);

assert(max(spectrum)>0.999)
assert(rho+ 1e-3 >= 1)
[~, max_idx] = findpeaks(spectrum);
assert(length(max_idx) == 1)
assert(max_idx(1) == period1)

period2 = 11;
s(5:period2:end) = .8;

tol = 1e-4;
threshold = .99;
[spectrum, rho] = inter_spike_spectrum(s, 'threshold', threshold, 'tol', tol);
assert(abs(rho - threshold) < tol)
[maxis, max_idx] = findpeaks(spectrum);
assert(length(max_idx) == 2)
assert(max_idx(1) == period1)
assert(max_idx(2) == period2 || max_idx(2) / max_idx(1) == period2)
assert(maxis(1) > maxis(2))

%% Test 4: Random time series
rng(1234)
tol = 1e-4;
maxcycles = 20;
s = randn(1,50);

threshold = 0.995;
[spectrum1, rho] = inter_spike_spectrum(s, 'threshold', threshold, 'tol', tol, 'max_iter', maxcycles);
[maxis, ~] = findpeaks(spectrum1);
numpeaks1 = length(maxis);
assert(abs(rho - threshold) <= tol)

threshold = 0.85;
[spectrum2, rho] = inter_spike_spectrum(s, 'threshold', threshold);
[maxis, ~] = findpeaks(spectrum2);
numpeaks2 = length(maxis);
assert(abs(rho - threshold) <= 1e-3)
assert(numpeaks2>numpeaks1)
assert(numpeaks2 == 7)
assert(numpeaks1 == 5)
assert(0.133 < maxis(1) && maxis(1) < 0.134)
assert(0.129 < maxis(2) && maxis(2) < 0.13)

