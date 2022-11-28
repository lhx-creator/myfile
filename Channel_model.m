function [hH_matrix,vH_matrix,G_matrix,bar_PHI_matrix,P_phi_matrix,hat_h_matrix,Pr_matrix,f_matrix] = Channel_model(A)

hH_matrix = zeros(1,A.Nt,A.B,A.K);
vH_matrix = zeros(1,A.N,A.R,A.K);
G_matrix = zeros(A.N,A.Nt,A.B,A.R);
bar_PHI_matrix = zeros(A.N,A.N,A.Nphi,A.R);
P_phi_matrix = zeros(A.N,A.N,A.Nphi,A.R);%相移的概率三维矩阵
hat_h_matrix = zeros(1,A.Nt,A.B,A.K);
f_matrix = zeros(A.Nt,1,A.B,A.K);

% Location_B = [1,2,3,4];
% Location_R = [3,4,5,6,7];
% Location_K = [0,1,2,3,4,5];  %BS IRS UE 坐标

PL_BS2UE = zeros(A.B,A.K);
PL_BS2IRS = zeros(A.B,A.R);
PL_IRS2UE = zeros(A.R,A.K);
for i = 1:A.B
    for j = 1:A.K
        PL_BS2UE(i,j) = sqrt((9-1)^2 + (j-1-i)^2);
    end
end
for i = 1:A.B
    for j = 1:A.R
        PL_BS2IRS(i,j) = sqrt((j+2-1)^2 + (5-i)^2);
    end
end
for i = 1:A.R
    for j = 1:A.K
        PL_IRS2UE(i,j) = sqrt((9-i-2)^2 + (j-1-5)^2);
    end
end
% for r = 1:A.R    
%     angle_r(r,:) = 2*pi*rand(1,A.N * A.N);
%     phi_r{r} = diag(exp( j * angle_r(r,:)));
% end

%% 波束成形向量初始化，平均分配
for b = 1:A.B     
    for k = 1:A.K
        a = rand(A.Nt,1);     %中间变量
        f_matrix(:,:,b,k) = a/norm(a,'fro')*sqrt(A.W_max/A.K);
    end
end
clear a;    %清除中间变量

%% 定义各信道  
% BS-UE大尺度*瑞丽衰落
for b = 1:A.B         
    for k= 1:A.K  
        hH_matrix(:,:,b,k) = ((randn(1, A.Nt) + 1i*randn(1, A.Nt))/sqrt(2)*sqrt(PL_BS2UE(b,k)));   %BS-UE实际信道 行向量
    end
end
% 其他大尺度*莱斯衰落
d_ratio_lamda = 1/2;
beta = 3;       %莱斯衰落系数
% BS-IRS
for b = 1:A.B
    for r = 1:A.R
        a_D_r = zeros(A.N,1);
        a_D_t = zeros(A.Nt,1);
        theta_AoA = rand(1)*2*pi;
        theta_AoD = rand(1)*2*pi;
        for nr = 1:A.N
            a_D_r(nr) = exp(i*2*pi*d_ratio_lamda*(nr-1)*sin(theta_AoA));
        end
        for nt = 1:A.Nt
            a_D_t(nt) = exp(i*2*pi*d_ratio_lamda*(nt-1)*sin(theta_AoD));
        end
        H_LOS= a_D_r*a_D_t';
        H_NLOS = ((randn(A.N, A.Nt) + 1i*randn(A.N, A.Nt))/sqrt(2));
        G_matrix(:,:,b,r)= sqrt(PL_BS2IRS(b,r))*(sqrt(beta/(beta+1))*H_LOS+sqrt(1/(beta+1))*H_NLOS);
    end
end
% IRS-UE
for r = 1:A.R
    for k = 1:A.K
        a_D_r = zeros(1,1);
        a_D_t = zeros(A.N,1);
        theta_AoA = rand(1)*2*pi;
        theta_AoD = rand(1)*2*pi;
        for nr = 1:1
            a_D_r(nr) = exp(i*2*pi*d_ratio_lamda*(nr-1)*sin(theta_AoA));
        end
        for n = 1:A.N
            a_D_t(n) = exp(i*2*pi*d_ratio_lamda*(n-1)*sin(theta_AoD));
        end
        H_LOS= a_D_r*a_D_t';
        H_NLOS = ((randn(1, A.N) + 1i*randn(1, A.N))/sqrt(2));
        vH_matrix(:,:,r,k) = sqrt(PL_IRS2UE(r,k))*(sqrt(beta/(beta+1))*H_LOS+sqrt(1/(beta+1))*H_NLOS);
    end
end


%%

% %BS-UE
% %for m = 1:L.B * A.K%L.B    
% for m = 1:A.B * A.K    
%         eval(['hH',num2str(m), ' =rand(1,A.Nt)']);
% end
% %for i = 1:L.B * A.K
%  for i = 1:A.B * A.K
%     hH_matrix(:,:,i) = eval(['hH', num2str(i)]); %hH
%  end

%% 
Pr_matrix = rand(1,A.R);%IRS开关概率

%%

% %IRS-UE
% for m = 1:A.K * A.R    
%         eval(['vH',num2str(m), ' =rand(1,A.N)']);
% end
% 
%  for i = 1:A.K * A.R
%     vH_matrix(:,:,i) = eval(['vH', num2str(i)]); %vH
%  end

 %%
 
 %BS-IRS
% for m = 1:A.B * A.R    
%         eval(['G',num2str(m), ' =rand(A.N,A.Nt)']);
% end
% 
%  for i = 1: A.B * A.R
%     G_matrix(:,:,i) = eval(['G', num2str(i)]); %G
%  end
%  
 %%
 
 %相移矩阵 10种
 for m = 1:A.Nphi * A.R    
        eval(['phi',num2str(m), ' =diag(rand(1,A.N))',';']);  %做50个IRS相移矩阵
end

 for i = 1:A.Nphi * A.R
    bar_PHI_matrix(:,:,i) = eval(['phi', num2str(i)]); %PHI 首先是每个IRS有10个，然后有5个IRS
 end
 
 %概率模型 10个
for m = 1:A.Nphi  
        a = rand(1,A.Nphi);
        a = a./sum(a);
        eval(['P_phi',num2str(m), ' =diag(a)',';']);
end

for j = 1:A.R
     for i = 1:A.Nphi
        P_phi_matrix(:,:,i,j) = eval(['P_phi', num2str(i),';']); %  (10,5)的概率模型
     end
end

%%
%选择相移矩阵部分

%第一种做法
P_Phi_SUM_matrix = zeros(A.N,A.N,A.R);
for i = 1:A.R % 第i个IRS
    aa = eval(['P_phi',num2str(i)],';');
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
        hat_h = hH_matrix(:,:,b,k) + P_V_P_G;
%         eval(['hat_h',num2str((k-1)*A.B+b), ' =hat_h']);
        hat_h_matrix(:,:,(k-1)*A.B+b) = hat_h;        
    end
end

        
