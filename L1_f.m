function [f_matrix] = L1_f( A, f_matrix, Pb_matrix, miu, hat_h_matrix)

square_sum_f = zeros(1,A.B);
cons_square_sum_f = zeros(1,A.B);
PHF = zeros(A.K,1);
HF = zeros(A.K+1,1);
H_matrix = zeros(A.K+1,A.B * A.K * A.Nt);
HF_sigma = zeros(A.K+1,A.K);
%%
% f 优化子问题
cvx_begin 
    variable f_matrix(A.Nt,1,A.B,A.K) complex ;
    expression square_sum_f(A.B);%B
    expression sum_fun;
    expression cons_square_sum_f(A.B);%B
    expression PHF(A.K,1);   
    expression HF(A.K+1);
    expression HF_sigma(A.K+1,A.K);
    
    for b = 1:A.B%B
        a1 = 0;
        for k = 1:A.K%K       
            a1 = a1 + square_pos(norm(f_matrix(:,:,b,k), 'fro'));  %!!!得看看能不能用
        end
        square_sum_f(b) = a1;
    end
    for b = 1:A.B%B
        a2 = 0;
        a2 = a2 + Pb_matrix(1,b) * square_sum_f(b);
    end
    sum_fun = a2;
    minimize( miu * sum_fun )
    subject to
    for i = 1:A.K * A.B
        F_matrix((i-1)*A.Nt+1:A.Nt*i,1) = f_matrix(:,:,i);
    end
    for b = 1:A.B%B
        cons_square_sum_f(b) = 0;
        for k = 1:A.K%K       
            cons_square_sum_f(b) = cons_square_sum_f(b) + square_pos(norm(f_matrix(:,:,b,k), 'fro'));
        end
        cons_square_sum_f(b) <= A.W_max;  
    end
    for k = 1:A.K
        PHF(k) = 0;
        for b = 1:A.B      
            PHF(k) = PHF(k) + Pb_matrix(1,b) *  hat_h_matrix(:,:,b,k) *  f_matrix(:,:,b,k) ;
        end
        imag(PHF(k)) == 0 ;
    end 
    


    for i = 1:A.K%K 行数
        %for m = 1:4% 列数
            for j = 1:A.K %K
                for k = 1:A.B
                    H_matrix(i,(i-1)*4+1:4*i) = hat_h_matrix(:,:,k,j) * Pb_matrix(k);
                    %H_matrix(i,(i-1)*m+1:i*m) = hat_h_matrix(:,:,k,j);
                end
            end
        %end
    end
    HF = H_matrix * F_matrix;

    for k = 1:A.K%   
            HF_sigma(:,k) = HF + A.sigma(:,k);
            norm(HF_sigma(:,k),2) - sqrt((A.T_min+1)/A.T_min) * real(PHF(k)) <= 0; %这里取实数会不会影响结果
    end
cvx_end

