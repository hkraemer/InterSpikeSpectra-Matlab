function spectrum = pool_frequencies(y,N)
y = abs(y);
M = ceil(N/2);
% pool the same frequencies
spectrum = zeros(1,M);
cnt = 1;
for i = 1:M
    occs = sum(y(cnt:cnt+i-1)>0);
    if occs>0
        spectrum(i) = sum(y(cnt:cnt+i-1))/occs;
    else
        spectrum(i) = 0;
    end
    cnt = cnt + i;
end
end