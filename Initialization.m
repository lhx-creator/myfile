function [Pb_matrix,f_matrix,Pr_matrix,miu,HF_sigma,fq_fun,T_initial,gamma,Z ] = Initialization(A, F, hat_h_matrix )


%%
%Leader
%
miu = 0.5;
%Pb_matrix = ones(1,L.B);
Pb_matrix = rand(1,4);
%f_matrix = zeros(10,10,L.B,A.K);%L.B A.K
f_matrix = zeros(4,1,4,6);

%for m = 1:L.B * A.K%L.B    
for m = 1:24    
        eval(['f',num2str(m), ' =rand(4,1)']);
end
%for i = 1:L.B * A.K
 for i = 1:24
    f_matrix(:,:,i)= eval(['f', num2str(i)]);
 end

% [A,L,F] = Parameter;
% [hH_matrix,vH_matrix,G_matrix,bar_PHI_matrix,P_matrix,hat_h_matrix] = Channel_model(A,L,F);

HF = zeros(1,6);
for i = 1:6 %K
    for j = 1:4
        HkFk = 0;
        HkFk = HkFk + Pb_matrix(j) * abs( hat_h_matrix(:,:,j,i) * f_matrix(:,:,j,i) ); %复信道取模运算
    end
    HF(i) = HkFk;
end
HF = [HF,0]';
HF_total = [HF,HF,HF,HF,HF,HF];
HF_sigma = HF_total + A.sigma;




%%
%   Follower

%初始化 T_initial(k) Z(k)
SINR_Fenzi = zeros(1,6);
tilde_h_matrix = zeros(1,4,4,6);
tilde_f_matrix = zeros(4,1,4,6);
%P_phi_matrix = [ones(1,5);zeros(9,5)];
Pr_matrix = ones(1,5);
SINR = zeros(1,6);
T_initial = zeros(1,6);
gamma = zeros(1,6);
Z = zeros(1,6);
% F.gamma = [1,1,1,1,1,1];
% F.Z = [1,1,1,1,1,1];

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
    T_initial(k) = SINR(k);
    gamma(k) = T_initial(k);
    
end

%计算初始化的Z
for b = 1:4   
    for k = 1:6
        tilde_h_matrix(:,:,b,k) = hat_h_matrix(:,:,b,k);
        tilde_f_matrix(:,:,b,k) = f_matrix(:,:,b,k) * Pb_matrix(1,b);
    end
end

tilde_fhH = zeros(1,6);
tilde_fh = zeros(1,6);
sum_fhhf = zeros(1,6);

for k = 1:6  
    for b = 1:4%B
        tilde_fhH(k) = tilde_fhH(k) + ctranspose( tilde_f_matrix(:,:,b,k) ) * ctranspose( tilde_h_matrix(:,:,b,k) );
    end
end

for k = 1:6  
    for b = 1:4%B
        tilde_fh(k) = tilde_fh(k) + tilde_h_matrix(:,:,b,k) * tilde_f_matrix(:,:,b,k);
    end
end
for g = 1:6
    for k = 1:6
    sum_fhhf(g) = sum_fhhf(g) + tilde_fh(k) * tilde_fhH(k);
    end
    sum_fhhf(g) = sum_fhhf(g) + A.sigma(g)^2;
end

for k = 1:6
    Z(k) = sqrt(1 + gamma(k) ) * tilde_fh(k) / sum_fhhf(g);
end


[ fq_fun ] = Fq_fun( tilde_f_matrix, tilde_h_matrix, A, gamma, Z );