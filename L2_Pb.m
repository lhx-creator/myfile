function [Pb_matrix] = L2_Pb( Pb_matrix, f_matrix, miu, A, hat_h_matrix)
PHF = zeros(A.K,1);
sum_fun = 0;
f_b = zeros(1,A.K);
cons_square_sum_f = zeros(1,A.B);
H_matrix = zeros(A.K+1,A.Nt * A.B * A.K);
HF_sigma = zeros(A.K+1,A.K);

cvx_begin sdp
    variable Pb_matrix(1,A.B) 
    expression cons_square_sum_f(A.B)%B   
    expression square_sum_f(A.B);%B
    expression sum_fun;
    expression PHF(A.K,1);   
    expression HF(A.K+1,1);
    expression HF_sigma(A.K+1,A.K);
    expression H_matrix(A.K+1,A.Nt * A.B * A.K)
    
    
    for b = 1:A.B
        a = 0;
        
        for k = 1:A.K
            a = a + square_pos(norm(f_matrix(:,:,b,k), 'fro'));
        end
        f_b(b) = a;
        sum_fun = sum_fun + Pb_matrix(1,b) * A.W_b + miu * Pb_matrix(1,b) * f_b(b) + Pb_matrix(1,b) * A.W_bh;
    end   
    minimize( sum_fun )
    subject to
    for b = 1:A.B
        Pb_matrix(1,b) <= 1;
        Pb_matrix(1,b) >= 0;
    end

    for i = 1:A.K * A.B
        F_matrix((i-1)*A.Nt+1:A.Nt*i,1) = f_matrix(:,:,i);
    end
    for b = 1:A.B%B
        cons_square_sum_f(b) = 0;
        for k = 1:A.K%K       
            %a3 = a3 + sum_square_abs( f_matrix(:,:,b,k) ) ;
            cons_square_sum_f(b) = cons_square_sum_f(b) + square_pos(norm(f_matrix(:,:,b,k), 'fro'));
        end
        cons_square_sum_f(b) <= A.W_max;  
    end
    for k = 1:A.K
        PHF(k) = 0;
        for b = 1:A.B      
            PHF(k) = PHF(k) + Pb_matrix(1,b) *  hat_h_matrix(:,:,b,k) *  f_matrix(:,:,b,k);
        end
        imag(PHF(k)) == 0 ;
    end 

    for i = 1:A.K%K 行数
        %for m = 1:4% 列数
            for j = 1:A.K %K
                for k = 1:A.B
                    H_matrix(i,(i-1)*A.Nt+1:A.Nt*i) = hat_h_matrix(:,:,k,j) * Pb_matrix(k);
                    %H_matrix(i,(i-1)*m+1:i*m) = hat_h_matrix(:,:,k,j);
                end
            end
        %end
    end
    HF = H_matrix * F_matrix;

    for k = 1:A.K%B   
            HF_sigma(:,k) = HF + A.sigma(:,k);
            norm(HF_sigma(:,k),2) - sqrt((A.T_min+1)/A.T_min) * real(PHF(k)) <= 0; %这里取实数会不会影响结果
    end
    
cvx_end