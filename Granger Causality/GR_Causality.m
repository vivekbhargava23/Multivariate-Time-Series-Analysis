function [Wald_statistic,p_wald, F_statistic,p_F_statistic] = GR_Causality(y,p,intercept, vector)
%{   
y -  k x (T+p) - matrix of y values
p - number of lags in the VAR p process
vector - contains1 for z's and 0 for x's
%}

[K2,~] = size (y);
T2 = size(y,2) - p;


if intercept == 1
    interc = 1;
else
    interc = 0;
end

total_columns = K2*p + interc;

zt_index = find(vector==1);
xt_index = find(vector==0);

zt = y(zt_index,:);
xt = y(xt_index,:);

N_z = size(zt_index,2);
N_x = size(xt_index,2);

N = N_z * N_x * p;


%{
Idea to make C matrix 

Following the order - 
1. Going equation by equation for each Zt
    2. In each zt equation - then going by A1,2,3.....= VAR(p=1,2,....) coefficient matrix
        3. for Zt, and A1, finding the column of xt
%}

count = 0;
for i = 1:N_z
    
    for j = 1:p
        
        % creating a matrix with dimensions same as B, 
        % the idea is to replicate B matix and put coefficients as 1
        % corrseponding to restrictions that we want to test
        
        % and then we use vec(B_reference) would form one row of C matrix
        for k = 1:N_x
            
            count = count +1; % Count will later help to make the C matrix, which stores one row for each iteration
            
            B_reference = zeros(K2,K2*p + interc); 
            index_restriction_coefficient = interc + (j-1)*K2 + N_z + k;
            B_reference(i,index_restriction_coefficient) = 1;
            % the above B_ref matrix just stores 1 for coefficient for
            % which we need to test
            
            vecB = B_reference(:);
            C(count,:) = vecB;
            
        end
        
    end
    
end



[B_hat,residual_cov_hat,~,ZZ_prime] = VAR_est(y',p,intercept);

vec_B_hat = B_hat(:);
First_term = (C * vec_B_hat)';
Second_term = inv(C * kron(inv(ZZ_prime),residual_cov_hat) * C');
Third_term = C * vec_B_hat;


Wald_statistic = First_term * Second_term * Third_term;
F_statistic = Wald_statistic / N;


p_wald = 1 - chi2cdf(Wald_statistic, N);
p_F_statistic = 1- fcdf(F_statistic,N, (K2 * T2 - p*(K2^2) - K2));










            
            
            
            
            











end

