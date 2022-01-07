function [spectrum, rho] = compute_spectrum_according_to_threshold(s, Theta, threshold, tol, maxlambdas)

% initial Lambda (start with least squares solution)
lambda_i = 0;
% initial Lambda-step
lambda_step = 0.1;

pos = true;                         % indicates whether lambda_step is positive or not
upper = false;                      % indicates whether an upper limit for lambda has been reached
cond1 = false;
flag = true;
i = 0;
while flag
    i = i + 1;
    if i == maxlambdas + 1
        warning("Algorithm stopped due to maximum number of lambdas were tried without convergence. Please increase `tol`-input OR increae `maxlambdas` and if this does not help `correlation_threshold` must be higher.")
        break
    end
   
    lambda_f = lambda_i + lambda_step;  % try new lambda
    lambda_i = lambda_f;                % update lambda for the next run
    
    % make the regression with specific lambda
    y = lasso(Theta, s, 'Lambda', lambda_f);
    % check whether the regenerated signal matches with the given threshold
    rr = corrcoef(regenerate_signal(Theta, y), s);

    % check whether upper limit is reached or convergence is fullfilled
    if isnan(rr(2)) || (rr(2) - threshold) < 0
        upper = true;
    end
    if (rr(2) - threshold) == 0 || abs(rr(2) - threshold) < tol
        rho = rr(2);
        spectrum = pool_frequencies(y, length(s));
        flag = false;
        continue
    end
    
    % alter lamba-steps
    if ~isnan(rr(2)) && ~upper
        continue
    elseif (rr(2) - threshold) < 0 || isnan(rr(2))
        upper = true;
        if cond1
            pos = true;
        else
            pos = false;
        end
        cond1 = true; 
    elseif (rr(2) - threshold) > 0
        if cond1
            pos = false;
        else
            pos = true;
        end
        cond1 = false; 
    elseif isnan(rr(2)) && ~upper  
        upper = true;
        pos = false;
        cond1 = true;
    end
    if pos
        lambda_step = lambda_step/2; % split the interval
    else
        lambda_step = (lambda_step*(-1))/2; % split the interval
    end
end
end