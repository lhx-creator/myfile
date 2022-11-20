function [hat_h,tilde_h,tilde_f] = Channel_model(A,L,F)

hH_matrix = zeros(1,4,4,6);
vH_matrix = zeros(1,10,5,6);
G_matrix = zeros(10,4,4,5);
bar_PHI_matrix = zeros(10,10,10,5);
P_matrix = zeros(10,10,5);
hat_h_matrix = zeros(1,4,4,6);

%%

%BS-UE
%for m = 1:L.B * A.K%L.B    
for m = 1:24    
        eval(['hH',num2str(m), ' =rand(1,4)']);
end
%for i = 1:L.B * A.K
 for i = 1:24
    hH_matrix(:,:,i) = eval(['hH', num2str(i)]); %hH
 end

%% 
P_R = rand(1,5);%IRS开关概率

%%

%IRS-UE
for m = 1:30    
        eval(['vH',num2str(m), ' =rand(1,10)']);
end
%for i = 1:L.B * A.K
 for i = 1:30
    vH_matrix(:,:,i) = eval(['vH', num2str(i)]); %vH
 end

 %%
 
 %BS-IRS
for m = 1:20    
        eval(['G',num2str(m), ' =rand(10,4)']);
end
%for i = 1:L.B * A.K
 for i = 1:20
    G_matrix(:,:,i) = eval(['G', num2str(i)]); %G
 end
 
 %%
 
 %相移矩阵 10种
 for m = 1:50    
        eval(['phi',num2str(m), ' =diag(rand(1,10))']);  %做50个IRS相移矩阵
end
%for i = 1:L.B * A.K
 for i = 1:50
    bar_PHI_matrix(:,:,i) = eval(['phi', num2str(i)]); %PHI 首先是每个IRS有10个，然后有5个IRS
 end
 
 %概率模型 10个
for m = 1:10  
        a = rand(1,10);
        a = a./sum(a);
        eval(['P_phi',num2str(m), ' =diag(a)']);
end
%for i = 1:L.B * A.K
 for i = 1:5
    P_matrix(:,:,i) = eval(['P_phi', num2str(i)]); %P
 end

%%
%选择相移矩阵部分

%第一种做法
P_Phi_SUM_matrix = zeros(10,10,5);
for i = 1:5 % 第i个IRS
    aa = eval(['P_phi',num2str(i)]);
    aaa = zeros(10,10);
    for j = 1:10 %第j个相移矩阵以及概率
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
for b = 1:4
    for k = 1:6      

        P_V_P_G = [];
        for i = 1:5
            P_V_P_G = P_V_P_G + P_R(1,i) * vH_matrix(:,:,i,k) * P_Phi_SUM_matrix(:,:,i) * G_matrix(:,:,b,i);
        end
        hat_h = eval(['hH',num2str(b*k)]) + P_V_P_G;
        eval(['hat_h',num2str(b*k), ' =hat_h']);
        hat_h_matrix(:,:,b*k) = eval(['hat_h',num2str(b*k)]);        
    end
end


