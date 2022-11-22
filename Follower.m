function [ fun_follower ] = Follower( A, Pr_matrix, fq_fun )


%初始化
% SINR_Fenzi = zeros(1,6);
% tilde_h_matrix = zeros(1,4,4,6);
% tilde_f_matrix = zeros(4,1,4,6);
% P_phi_matrix = [ones(1,5);zeros(9,5)];
% 
% Pr_matrix = ones(1,5);
% SINR = zeros(1,6);
% %%
% %SINR
% for k = 1:6
%     SINR_fenzi = 0;
%     for b = 1:4
%         SINR_fenzi = SINR_fenzi + abs( hat_h_matrix(:,:,b,k) * f_matrix(:,:,b,k) );
%     end
%     SINR_Fenzi(k) = ( abs( SINR_fenzi ) )^2;
%     
%     SINR_fenmu = zeros(1,6);
%     
%     aa = zeros(1,6);
%     for j = 1:6
%         if j ~= k
%             a = 0;
%             for b = 1:4
%                 a = a +  abs( hat_h_matrix(:,:,b,k) * f_matrix(:,:,b,k) );
%             end
%             aa(j) = ( abs( a ) )^2;
%         else
%             aa(j) = 0;
%         end
%         SINR_fenmu(k) = SINR_fenmu(k) + aa(j);
%     end
%     SINR_fenmu(k) = SINR_fenmu(k) + A.sigma(7,k)^2;
%     
%     SINR(k) = SINR_Fenzi(k)./SINR_fenmu(k); 
%     
% end
% 
% 
% 
% %%
% %复合矩阵h,f
% for b = 1:4   
%     for k = 1:6
%         tilde_h_matrix(:,:,b,k) = hat_h_matrix(:,:,b,k);
%         tilde_f_matrix(:,:,b,k) = f_matrix(:,:,b,k) * Pb_matrix(1,b);
%     end
% end
% 
% [ fq_fun ] = Fq_fun( tilde_f_matrix, tilde_h_matrix, A, F );
% [P_phi_matrix] = F1_P_phi( P_phi_matrix, fq_fun );
% [ Pr_matrix ] = F2_Pr( Pr_matrix, A, fq_fun );

a1 = 0;
a3 = 0;

for r = 1:5%R
    a1 = a1 + Pr_matrix(r) * A.W_r;
end
sum_PWr = a1;
for b = 1:4
    a2 = 0;
    for r = 1:5
        a2 = a2 + Pr_matrix(r) * A.W_bh;
    end
    a3 = a3 + a2;
end
sum_PPWbh = a3;

fun_follower = ( -1 ) * fq_fun + A.alpha * sum_PWr + A.beta * sum_PPWbh




