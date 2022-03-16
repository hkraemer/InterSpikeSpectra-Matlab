function spectrum = pool_frequencies(y,N)
y = abs(y);
M = ceil(N/2);
% pool the same frequencies
spectrum = zeros(1,M);
cnt = 1;
for i = 1:M
    spectrum(i) = sum(y(cnt:cnt+i-1));
    cnt = cnt + i;
end
end