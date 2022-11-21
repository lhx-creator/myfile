%代码所涉及参数
function [A,L,F] = Parameter
A.W_b = 10; %基站功耗暂时设置10W
A.W_r = 10; %IRS功耗暂时设置3W
A.W_bh = 1; %回程链路功耗1W
A.W_max = 8; %得看之后波束成型向量设置 暂设
A.Rate_min = 10; %最小传输速率
A.T_min = 2^(A.Rate_min) - 1; %最小信噪比
A.K = 6;%用户的数量
A.sigma = [zeros(6,6);[1,1,1,1,1,1]]; %噪声功率

L.B = 4;%基站的数量
L.Nt = 4;%基站天线数

F.R = 5%IRS的数量
F.N = 10;%IRS的行数
F.Nphi = 10;%相移矩阵可选数

%信道还没有变成复信道
