function [lags] = FEVD1(B_hat_intercept_removed,lags,B_Cholesky,SigmaU_hat,K)
%FEVD1 Summary of this function goes here
%   Detailed explanation goes here

B1 = B_Cholesky;

phi_0 = eye(K);
phi_B_columnB1 = phi_0 * B1;
phi_B_columnB1_shock1 = phi_B_columnB1(:,1);
phi_B_columnB1_shock2 = phi_B_columnB1(:,2);
phi_B_columnB1_shock3 = phi_B_columnB1(:,3);
phi_B_columnB1_shock4 = phi_B_columnB1(:,4);


for i=1+1:lags % because first column is the h=0 lag so we increment the loop by 1
    
    h=i-1;
    [phi_i] = generate_phi_i(B_hat_intercept_removed,h);
    phi_B_columnB1_shock1(:,i) = phi_i * B1(:,1);
    phi_B_columnB1_shock2(:,i) = phi_i * B1(:,2);
    phi_B_columnB1_shock3(:,i) = phi_i * B1(:,3);
    phi_B_columnB1_shock4(:,i) = phi_i * B1(:,4);
    
end



% Calculating FEVD


for h=1:lags-1
    MSE_forecast = Forecast_MSE(B_hat_intercept_removed,SigmaU_hat,h);
    % [~,MSE_forecast] = Forecasting(Y_t',p,h,intercept)
    
        
        
        for k = 1:K
            
            
            FEVD_matrix(k,h,1) = (sum(phi_B_columnB1_shock1(k,1:h) .* phi_B_columnB1_shock1(k,1:h)))/MSE_forecast(k,k);
            FEVD_matrix(k,h,2) = (sum(phi_B_columnB1_shock2(k,1:h) .* phi_B_columnB1_shock2(k,1:h)))/MSE_forecast(k,k);
            FEVD_matrix(k,h,3) = (sum(phi_B_columnB1_shock3(k,1:h) .* phi_B_columnB1_shock3(k,1:h)))/MSE_forecast(k,k);
            FEVD_matrix(k,h,4) = (sum(phi_B_columnB1_shock4(k,1:h) .* phi_B_columnB1_shock4(k,1:h)))/MSE_forecast(k,k);
            

        end
        
end


sublot1_data = [FEVD_matrix(1,:,1);FEVD_matrix(1,:,2);FEVD_matrix(1,:,3);FEVD_matrix(1,:,4)];
sublot2_data = [FEVD_matrix(2,:,1);FEVD_matrix(2,:,2);FEVD_matrix(2,:,3);FEVD_matrix(2,:,4)];
sublot3_data = [FEVD_matrix(3,:,1);FEVD_matrix(3,:,2);FEVD_matrix(3,:,3);FEVD_matrix(3,:,4)];
sublot4_data = [FEVD_matrix(4,:,1);FEVD_matrix(4,:,2);FEVD_matrix(4,:,3);FEVD_matrix(4,:,4)];

Fontsize1 = 6;

figure
subplot(4,1,1)
bar(1:lags-1,sublot1_data,'stacked')
ylim([0 1])
lgd1 = legend('FFR','TP1','CPI','IP');
lgd1.FontSize = Fontsize1;
title('FFR')

subplot(4,1,2)
bar(1:lags-1,sublot2_data,'stacked')
ylim([0 1])
lgd2 = legend('FFR','TP1','CPI','IP');
lgd2.FontSize = Fontsize1;
title('TP1')

subplot(4,1,3)
bar(1:lags-1,sublot3_data,'stacked')
ylim([0 1])
lgd3 = legend('FFR','TP1','CPI','IP');
lgd3.FontSize = Fontsize1;
title('CPI')

subplot(4,1,4)
bar(1:lags-1,sublot4_data,'stacked')
ylim([0 1])
lgd4 = legend('FFR','TP1','CPI','IP');
lgd4.FontSize = Fontsize1;
title('IP')
sgtitle('FEVD - B Model (Cholesky)')







%{
sublot1_data = [FEVD_matrix(1,:,1);FEVD_matrix(1,:,2);FEVD_matrix(1,:,3);FEVD_matrix(1,:,4)]
sublot2_data = [FEVD_matrix(2,:,1);FEVD_matrix(2,:,2);FEVD_matrix(2,:,3);FEVD_matrix(2,:,4)]
sublot3_data = [FEVD_matrix(3,:,1);FEVD_matrix(3,:,2);FEVD_matrix(3,:,3);FEVD_matrix(3,:,4)]
sublot4_data = [FEVD_matrix(4,:,1);FEVD_matrix(4,:,2);FEVD_matrix(4,:,3);FEVD_matrix(4,:,4)]

figure
subplot(4,1,1)
bar(1:lags-1,sublot1_data,'stacked')

subplot(4,1,2)
bar(1:lags-1,sublot2_data,'stacked')

subplot(4,1,3)
bar(1:lags-1,sublot3_data,'stacked')

subplot(4,1,4)
bar(1:lags-1,sublot4_data,'stacked')
%}


end

