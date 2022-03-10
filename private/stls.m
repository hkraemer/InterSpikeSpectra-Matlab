function coeffs = stls(s, Theta, lambda)
iterations = 10;
% Sequential Thresholded Least Squares sparse regression method
coeffs = pinv(Theta)*s;     % initial guess least squares
%coeffs = Theta\s;  % initial guess least squares
for k = 1:iterations
    smallinds = (abs(coeffs)<lambda);  % find small coefficients 
    coeffs(smallinds) = 0;  % threshold these coeffs
    biginds = ~smallinds;
    % regress onto remaining terms
    coeffs(biginds) = pinv(Theta(:,biginds))*s;
    %coeffs(biginds) = Theta(:,biginds)\s;

end
end