function [fun_leader] = Leader(f_matrix,Pb_matrix,A)

%clc,clear; %

miu = 0.5;
% %Pb_matrix = ones(1,L.B);
% Pb_matrix = rand(1,4);
% %f_matrix = zeros(10,10,L.B,A.K);%L.B A.K
% f_matrix = zeros(4,1,4,6);
% 
% %for m = 1:L.B * A.K%L.B    
% for m = 1:24    
%         eval(['f',num2str(m), ' =rand(4,1)']);
% end
% %for i = 1:L.B * A.K
%  for i = 1:24
%     f_matrix(:,:,i)= eval(['f', num2str(i)]);
%  end
% 
% [A,L,F] = Parameter;
% [hH_matrix,vH_matrix,G_matrix,bar_PHI_matrix,P_matrix,hat_h_matrix] = Channel_model(A,L,F);
% 
% HF = [];
% for i = 1:6 %K
%     for j = 1:4
%         HkFk = 0;
%         HkFk = HkFk + Pb_matrix(j) * abs( hat_h_matrix(:,:,j,i) * f_matrix(:,:,j,i) ); %复信道取模运算
%     end
%     HF(i) = HkFk;
% end
% HF = [HF,0]';
% HF_total = [HF,HF,HF,HF,HF,HF];
% HF_sigma = HF_total + A.sigma;

% [f_matrix] = L1_f( f_matrix, Pb_matrix, miu, hat_h_matrix, HF_sigma);
% [Pb_matrix] = L2_Pb( Pb_matrix, f_matrix, miu, A, hat_h_matrix, HF_sigma );

f_b = zeros(1,6);
fun_leader = 0;
for b = 1:4
    a = 0;
    
    for k = 1:6
        a = a + ( abs( f_matrix(:,:,b,k) ) ).^2
    end
    f_b(b) = a;
    fun_leader = fun_leader + Pb_matrix(1,b) * A.W_b + miu * Pb_matrix(1,b) * f_b(b) + Pb_matrix(1,b) * A.W_bh;
end



