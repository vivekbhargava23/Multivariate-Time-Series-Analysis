function [B_hat,residual_covariance_hat, t_ratio,U_hat] = VAR_est(y,p,intercept)
%VAR_EST Summary of this function goes here
%{
INPUTS - 
y - (T+p) x k - matrix of y values
p - number of lags in the VAR p process
intercept - dummy value - 1 to include intercept and 0 for no intercept
%}


T = size(y,1) - p;

for i = 1:T
    Zt = [];
    for j = 1:p
        
        Zt = [Zt, y(p+i-j,:)];
        
    end
    
    Z(:,i) = Zt;
end


Y = y(p+1:p+T,:)';

% Accomodating for the intercept
if intercept == 1
    row_one = ones(1,T);
    Z = [row_one;Z];
end

% Estimating the parameter matrix B

B_hat = Y * Z' * inv (Z*Z');


% Estimating residual covariance matrix
U_hat = Y - B_hat * Z;

DOF = T - size(Z,1);
residual_covariance_hat = (1/DOF) * U_hat * U_hat';

% Calculating t-ratios
Diag1 = diag(kron (inv(Z*Z'),residual_covariance_hat));
std1 = sqrt(Diag1);
vecB = B_hat(:);
t_ratio = vecB ./ std1;

end

