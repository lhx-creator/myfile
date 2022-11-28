function [Pb_matrix,f_matrix,miu,fq_fun,T_initial,gamma,Z ] = Initialization( A,hat_h_matrix,f_matrix )


%%
%Leader
%
miu = 0.5;

Pb_matrix = rand(1,A.B);


F_matrix = zeros(A.Nt * A.B * A.K,1);
H_matrix = zeros(A.B+1,A.Nt * A.B * A.K);


for i = 1:A.B * A.K
    F_matrix((i-1)*A.Nt+1:A.Nt*i,1) = f_matrix(:,:,i);
end
for i = 1:A.K%K 行数
        for j = 1:A.K %K
            for k = 1:A.B
                H_matrix(i,(i-1)*A.Nt+1:A.Nt*i) = hat_h_matrix(:,:,k,j) * Pb_matrix(k);
            end
        end
end


%%
%   Follower

%初始化 T_initial(k) Z(k)
SINR_Fenzi = zeros(1,A.K);
tilde_h_matrix = zeros(1,A.Nt,A.B,A.K);
tilde_f_matrix = zeros(A.Nt,1,A.B,A.K);


SINR = zeros(1,A.K);
T_initial = zeros(1,A.K);
gamma = zeros(1,A.K);
Z = zeros(1,A.K);


%SINR
for k = 1:A.K
    SINR_fenzi = 0;
    for b = 1:A.B
        SINR_fenzi = SINR_fenzi + abs( hat_h_matrix(:,:,b,k) * f_matrix(:,:,b,k) );
    end
    SINR_Fenzi(k) = ( abs( SINR_fenzi ) )^2;
    
    SINR_fenmu = zeros(1,A.K);
    
    aa = zeros(1,A.K);
    for j = 1:A.K
        if j ~= k
            a = 0;
            for b = 1:A.B
                a = a +  abs( hat_h_matrix(:,:,b,k) * f_matrix(:,:,b,k) );
            end
            aa(j) = ( abs( a ) )^2;
        else
            aa(j) = 0;
        end
        SINR_fenmu(k) = SINR_fenmu(k) + aa(j);
    end
    SINR_fenmu(k) = SINR_fenmu(k) + A.sigma(A.K+1,k)^2;
    
    SINR(k) = SINR_Fenzi(k)./SINR_fenmu(k); 
    T_initial(k) = SINR(k);
    gamma(k) = T_initial(k);
    
end

%计算初始化的Z
for b = 1:A.B   
    for k = 1:A.K
        tilde_h_matrix(:,:,b,k) = hat_h_matrix(:,:,b,k);
        tilde_f_matrix(:,:,b,k) = f_matrix(:,:,b,k) * Pb_matrix(1,b);
    end
end

tilde_fhH = zeros(1,A.K);
tilde_fh = zeros(1,A.K);
sum_fhhf = zeros(1,A.K);

for k = 1:A.K  
    for b = 1:A.B%B
        tilde_fhH(k) = tilde_fhH(k) + ctranspose( tilde_f_matrix(:,:,b,k) ) * ctranspose( tilde_h_matrix(:,:,b,k) );
    end
end

for k = 1:A.K  
    for b = 1:A.B%B
        tilde_fh(k) = tilde_fh(k) + tilde_h_matrix(:,:,b,k) * tilde_f_matrix(:,:,b,k);
    end
end
for g = 1:A.K
    for k = 1:A.K
    sum_fhhf(g) = sum_fhhf(g) + tilde_fh(k) * tilde_fhH(k);
    end
    sum_fhhf(g) = sum_fhhf(g) + A.sigma(A.K+1,g)^2;
end

[ fq_fun ] = Fq_fun( tilde_f_matrix, tilde_h_matrix, A, gamma);