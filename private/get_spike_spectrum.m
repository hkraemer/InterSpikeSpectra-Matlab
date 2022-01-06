function spectrum = get_spike_spectrum(s, lambda)
    assert(sum(s<0)==0, "Please supply a spike train input, which consists of only positive values.")
    assert(sum(s>1)==0, "Please supply a spike train input, which consists of only positive values.")

    N = length(s);
    % get set of basis functions
    Theta = generate_basis_functions(N)';
    % make LASSO regression
    y = lasso(Theta,s,'Lambda',lambda,'Intercept',false);
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