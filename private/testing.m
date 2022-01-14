clc, clear

% rng(1234)
% N = 50;
% period1 = 13;
% period2 = 8;
% test_tauRR2 = create_single_basis_function(N, period1);
% test_tauRR3 = create_single_basis_function(N, period2);
% test_tauRR = 0.7 .* test_tauRR2(3,:);
% test_tauRR_ = 0.8 .* test_tauRR3(7,:);
% test_tauRR = test_tauRR + test_tauRR_;
% test_tauRR = test_tauRR + 0.05.*randn(1,N);

N = 100;
s = zeros(1,N);
period1 = 20;
s(2:period1:end) = 1;

%s = test_tauRR;

%%
% normalize time series
s = (s - mean(s)) ./ std(s);
s = s - min(s);
s = s ./ max(s);

N = length(s);
% get set of basis functions
Theta = generate_basis_functions(N)';

%%

plot(s)

%% 
clc
lambda_step = 0.00001;
lambda = 0;
lambda_min = 0;
rho_min = 1;
y_min = zeros(1,size(Theta,2));


figure
for i = 1:30
    % make the regression with specific lambda
    y = lasso(Theta, s, 'Lambda', lambda);
    % check whether the regenerated signal matches with the given threshold
    reg = regenerate_signal(Theta, y);
    spectrum = pool_frequencies(y, length(s));
    rr = corrcoef(spectrum, spec_i);
        rrr = rr(2);   
    subplot(5,6,i)
    %plot(1:N,spectrum)
    plot(1:N,spectrum)

    grid on
    
    title(strcat("\lambda=",num2str(lambda),"{}","{}","\rho=",num2str(rrr)))
    
%     if i == 1
%         y_min(:) = y;
%         rho_min = rr(2);
%     end
% 
%     if rr(2) > rho_thres
%         lambda_min = lambda;
%         y_min(:) = y;
%         rho_min = rr(2);
%     elseif isnan(rr(2)) || rr(2) <= rho_thres
%         break
%     end
    lambda = lambda + lambda_step;
end


%%

%%
clc
[spectrum, rho] = inter_spike_spectrum(s, 'max_iter',20);
rho
figure
subplot(211)
plot(1:N,s)
grid on
subplot(212)
plot(1:N,spectrum)
grid on

%%
spec_i = spectrum;

%%


