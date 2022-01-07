function spectrum = pool_frequencies(y,N)
% pool the same frequencies
spectrum = zeros(1,N);
cnt = 1;
for i = 1:N
    occs = sum(y(cnt:cnt+i-1)>0);
    if occs>0
        spectrum(i) = sum(y(cnt:cnt+i-1))/occs;
    else
        spectrum(i) = 0;
    end
    cnt = cnt + i;
end

end