function [hH_matrix,vH_matrix,G_matrix,bar_PHI_matrix,P_phi_matrix,hat_h_matrix,Pr_matrix] = Channel_model(A)

hH_matrix = zeros(1,A.Nt,A.B,A.K);
vH_matrix = zeros(1,A.N,A.R,A.K);
G_matrix = zeros(A.N,A.Nt,A.B,A.R);
bar_PHI_matrix = zeros(A.N,A.N,A.Nphi,A.R);
P_phi_matrix = zeros(A.N,A.N,A.Nphi,A.R);%相移的概率三维矩阵
hat_h_matrix = zeros(1,A.Nt,A.B,A.K);

%%

%BS-UE
%for m = 1:L.B * A.K%L.B    
for m = 1:A.B * A.K    
        eval(['hH',num2str(m), ' =rand(1,A.Nt)']);
end
%for i = 1:L.B * A.K
 for i = 1:A.B * A.K
    hH_matrix(:,:,i) = eval(['hH', num2str(i)]); %hH
 end

%% 
Pr_matrix = rand(1,A.R);%IRS开关概率

%%

%IRS-UE
for m = 1:A.K * A.R    
        eval(['vH',num2str(m), ' =rand(1,A.N)']);
end

 for i = 1:A.K * A.R
    vH_matrix(:,:,i) = eval(['vH', num2str(i)]); %vH
 end

 %%
 
 %BS-IRS
for m = 1:A.B * A.R    
        eval(['G',num2str(m), ' =rand(A.N,A.Nt)']);
end

 for i = 1: A.B * A.R
    G_matrix(:,:,i) = eval(['G', num2str(i)]); %G
 end
 
 %%
 
 %相移矩阵 10种
 for m = 1:A.Nphi * A.R    
        eval(['phi',num2str(m), ' =diag(rand(1,A.N))']);  %做50个IRS相移矩阵
end

 for i = 1:A.Nphi * A.R
    bar_PHI_matrix(:,:,i) = eval(['phi', num2str(i)]); %PHI 首先是每个IRS有10个，然后有5个IRS
 end
 
 %概率模型 10个
for m = 1:A.Nphi  
        a = rand(1,A.Nphi);
        a = a./sum(a);
        eval(['P_phi',num2str(m), ' =diag(a)']);
end

for j = 1:A.R
     for i = 1:A.Nphi
        P_phi_matrix(:,:,i,j) = eval(['P_phi', num2str(i)]); %  (10,5)的概率模型
     end
end

%%
%选择相移矩阵部分

%第一种做法
P_Phi_SUM_matrix = zeros(A.N,A.N,A.R);
for i = 1:A.R % 第i个IRS
    aa = eval(['P_phi',num2str(i)]);
    aaa = zeros(A.N,A.N);
    for j = 1:A.N %第j个相移矩阵以及概率
        aaa = aaa + aa(j,j) * bar_PHI_matrix(:,:,j,i);
    end
    P_Phi_SUM_matrix(:,:,i) = aaa;
end  

%第二种做法
% P_Phi_SUM_matrix = zeros(10,10,5);
% for i = 1:5 % 第i个IRS
%     aa = eval(['P_phi',num2str(i)]);
%     aa = diag(aa);
%     [p,j] = max(aa);
%     aaa = p * bar_PHI_matrix(:,:,j,i);
%     P_Phi_SUM_matrix(:,:,i) = aaa;
% end  



%%


%等效信道集合
%BS-UE
%for m = 1:L.B * A.K%L.B    
for k = 1:A.K
    for b = 1:A.B      

        P_V_P_G = zeros(1,A.B);
        for i = 1:A.R
            P_V_P_G = P_V_P_G + Pr_matrix(1,i) * vH_matrix(:,:,i,k) * P_Phi_SUM_matrix(:,:,i) * G_matrix(:,:,b,i);
        end
        hat_h = eval(['hH',num2str((k-1)*A.B+b)]) + P_V_P_G;
        eval(['hat_h',num2str((k-1)*A.B+b), ' =hat_h']);
        hat_h_matrix(:,:,(k-1)*A.B+b) = eval(['hat_h',num2str((k-1)*A.B+b)]);        
    end
end

        
