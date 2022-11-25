function [hat_h_matrix,fq_fun] = Iteration_fun(A,P_phi_matrix,bar_PHI_matrix,Pr_matrix,vH_matrix,G_matrix,hH_matrix,f_matrix,Pb_matrix )
hat_h_matrix = zeros(1,4,4,6);

%%
% 更新等效信道
P_Phi_SUM_matrix = zeros(10,10,5);
for i = 1:5 % 第i个IRS
    for j = 1:10 %第j个相移矩阵以及概率
        P_Phi_SUM_matrix(:,:,i) = P_Phi_SUM_matrix(:,:,i) + P_phi_matrix(j,i) * bar_PHI_matrix(:,:,j,i);
    end
end 

for k = 1:6
    for b = 1:4      
        P_V_P_G = zeros(1,4);
        for i = 1:5
            P_V_P_G = P_V_P_G + Pr_matrix(1,i) * vH_matrix(:,:,i,k) * P_Phi_SUM_matrix(:,:,i) * G_matrix(:,:,b,i);
        end
        hat_h_matrix(:,:,b,k) = hH_matrix(:,:,b,k) + P_V_P_G;   
    end
end


%%
%更新信干噪比
%更新gamma

SINR_Fenzi = zeros(1,6);
SINR = zeros(1,6);
T_iter = zeros(1,6);
gamma = zeros(1,6);

for k = 1:6
    SINR_fenzi = 0;
    for b = 1:4
        SINR_fenzi = SINR_fenzi + Pb_matrix(1,b) * hat_h_matrix(:,:,b,k) * f_matrix(:,:,b,k);
    end
    SINR_Fenzi(k) = ( abs( SINR_fenzi ) )^2;
    
    SINR_fenmu = zeros(1,6);
    aa = zeros(1,6);
    for j = 1:6
        if j ~= k
            a = 0;
            for b = 1:4
                a = a +  Pb_matrix(1,b) * hat_h_matrix(:,:,b,k) * f_matrix(:,:,b,k);
            end
            aa(j) = ( abs( a ) )^2;
        else
            aa(j) = 0;
        end
        SINR_fenmu(k) = SINR_fenmu(k) + aa(j);
    end
    SINR_fenmu(k) = SINR_fenmu(k) + A.sigma(7,k)^2;
    
    SINR(k) = SINR_Fenzi(k)./SINR_fenmu(k); 
    T_iter(k) = SINR(k);
    gamma(k) = T_iter(k); 
end

%%
% 更新变换矩阵

tilde_h_matrix = zeros(1,4,4,6);
tilde_f_matrix = zeros(4,1,4,6);
for b = 1:4   
    for k = 1:6
        tilde_h_matrix(:,:,b,k) = hat_h_matrix(:,:,b,k);
        tilde_f_matrix(:,:,b,k) = f_matrix(:,:,b,k) * Pb_matrix(1,b);
    end
end


%%
%更新Z以及fq
[ fq_fun ] = Fq_fun( tilde_f_matrix, tilde_h_matrix, A, gamma );
