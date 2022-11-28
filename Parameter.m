%代码所涉及参数
function [A] = Parameter
A.W_b = 10; %基站功耗暂时设置10W
A.W_r = 10; %IRS功耗暂时设置3W
A.W_bh = 1; %回程链路功耗1W
A.W_max = 100; %得看之后波束成型向量设置 暂设
A.Rate_min = 10; %最小传输速率
A.T_min = 2^(A.Rate_min) - 1; %最小信噪比
A.K = 6;%用户的数量
A.sigma = [zeros(6,6);[1,1,1,1,1,1]]; %噪声功率
A.alpha = 0.5;
A.beta = 0.5;

A.B = 4;%基站的数量
A.Nt = 4;%基站天线数

A.R = 5%IRS的数量
A.N = 10;%IRS的行数
A.Nphi = 10;%相移矩阵可选数



%信道还没有变成复信道
