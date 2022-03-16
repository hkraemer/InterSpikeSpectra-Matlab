function [lambda, lambda_min, y_min, rho_min] = find_lambda_max_logit(s, Theta, lambda_step, rho_thres)
% estimate a maximum λ-value for which the correlation coefficient form the re-
% generated signal and the input signal falls below `ρ_thres` or is NaN.
lambda = 0;
lambda_min = 0;
rho_min = 1;
y_min = zeros(size(Theta,2),1);
for i = 1:10000
    % make the regression with specific lambda
    y = lassoglm(Theta, s, 'binomial', 'Lambda', lambda);
    % check whether the regenerated signal matches with the given threshold
    rr = corrcoef(regenerate_signal_logit(Theta, y), s);

    if i == 1
        y_min(:) = y;
        rho_min = rr(2);
    end

    if rr(2) > rho_thres
        lambda_min = lambda;
        y_min(:) = y;
        rho_min = rr(2);
    elseif isnan(rr(2)) || rr(2) <= rho_thres
        break
    end
    lambda = lambda + lambda_step;
end
end