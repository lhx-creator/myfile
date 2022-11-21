%function [ ] = Follower( )
clc,clear;
SINR_Fenzi = zeros(1,6);
tilde_h_matrix = zeros(1,4,4,6);
tilde_f_matrix = zeros(4,1,4,6);

%%
%SINR
for k = 1:6
    SINR_fenzi = 0;
    for b = 1:4
        SINR_fenzi = SINR_fenzi + abs( hat_h_matrix(:,:,b,k) * f_matrix(:,:,b,k) );
    end
    SINR_Fenzi(k) = ( abs( SINR_fenzi ) )^2;
    
    SINR_fenmu = zeros(1,6);
    
    aa = zeros(1,6);
    for j = 1:6
        if j ~= k
            a = 0;
            for b = 1:4
                a = a +  abs( hat_h_matrix(:,:,b,k) * f_matrix(:,:,b,k) );
            end
            aa(j) = ( abs( a ) )^2;
        else
            aa(j) = 0;
        end
        SINR_fenmu(k) = SINR_fenmu(k) + aa(j);
    end
    SINR_fenmu(k) = SINR_fenmu(k) + A.sigma(7,k)^2;
    
    SINR(k) = SINR_Fenzi(k)./SINR_fenmu(k); 
end


%%
for b = 1:4   
    for k = 1:6
        tilde_h_matrix(:,:,b,k) = hat_h_matrix(:,:,b,k);
        tilde_f_matrix(:,:,b,k) = f_matrix(:,:,b,k);
    end
end
%for i = 1:L.B * A.K
 for i = 1:24
    hH_matrix(:,:,i) = eval(['hH', num2str(i)]); %hH
 end





