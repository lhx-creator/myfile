function [phi_r,f,h_BS2UE_H,G_BS2IRS,v_IRS2UE_H] = initial(P,Nr,Nt,N,R,B,K,angle_r,phi_r,f,h_BS2UE_H,PL_BS2UE,...
    G_BS2IRS,PL_BS2IRS,v_IRS2UE_H,PL_IRS2UE)
%initial 产生初始化信道、波束成形、IRS随机相位

%% 随机初试化各IRS相位控制矩阵
for r = 1:R    
    angle_r(r,:) = 2*pi*rand(1,N);
    phi_r{r} = diag(exp(i*angle_r(r,:)));
end

%% 波束成形向量初始化，平均分配
for b = 1:B     
    for k = 1:K
        a = rand(Nt,1);     %中间变量
        f{b,k} = a/norm(a,'fro')*sqrt(P/K);
    end
end
clear a;    %清除中间变量

%% 定义各信道  
% BS-UE大尺度*瑞丽衰落
for b = 1:B         
    for k= 1:K  
        h_BS2UE_H{b,k} = ((randn(Nr, Nt) + 1i*randn(Nr, Nt))/sqrt(2)*sqrt(PL_BS2UE(b,k)));   %BS-UE实际信道 行向量
    end
end
% 其他大尺度*莱斯衰落
d_ratio_lamda = 1/2;
beta = 3;       %莱斯衰落系数
% BS-IRS
for b = 1:B
    for r = 1:R
        a_D_r = zeros(N,1);
        a_D_t = zeros(Nt,1);
        theta_AoA = rand(1)*2*pi;
        theta_AoD = rand(1)*2*pi;
        for nr = 1:N
            a_D_r(nr) = exp(i*2*pi*d_ratio_lamda*(nr-1)*sin(theta_AoA));
        end
        for nt = 1:Nt
            a_D_t(nt) = exp(i*2*pi*d_ratio_lamda*(nt-1)*sin(theta_AoD));
        end
        H_LOS= a_D_r*a_D_t';
        H_NLOS = ((randn(N, Nt) + 1i*randn(N, Nt))/sqrt(2));
        G_BS2IRS{b,r} = sqrt(PL_BS2IRS(b,r))*(sqrt(beta/(beta+1))*H_LOS+sqrt(1/(beta+1))*H_NLOS);
    end
end
% IRS-UE
for r = 1:R
    for k = 1:K
        a_D_r = zeros(Nr,1);
        a_D_t = zeros(N,1);
        theta_AoA = rand(1)*2*pi;
        theta_AoD = rand(1)*2*pi;
        for nr = 1:Nr
            a_D_r(nr) = exp(i*2*pi*d_ratio_lamda*(nr-1)*sin(theta_AoA));
        end
        for n = 1:N
            a_D_t(n) = exp(i*2*pi*d_ratio_lamda*(n-1)*sin(theta_AoD));
        end
        H_LOS= a_D_r*a_D_t';
        H_NLOS = ((randn(Nr, N) + 1i*randn(Nr, N))/sqrt(2));
        v_IRS2UE_H{r,k} = sqrt(PL_IRS2UE(r,k))*(sqrt(beta/(beta+1))*H_LOS+sqrt(1/(beta+1))*H_NLOS);
    end
end

end

