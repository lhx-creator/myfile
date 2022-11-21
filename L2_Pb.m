function [Pb_matrix] = L2_Pb( Pb_matrix, f_matrix, miu, A, hat_h_matrix, HF_sigma )
PHF = zeros(1,6);
sum_fun = 0;
f_b = zeros(1,6);

cvx_begin
    variable Pb_matrix(1,4) nonnegative
    expression sum_fun
    expression cons_square_sum_f(4)%B   
     
    for b = 1:4
        a = 0;
        
        for k = 1:6
            a = a + ( abs( f_matrix(:,:,b,k) ) ).^2
        end
        f_b(b) = a;
        sum_fun = sum_fun + Pb_matrix(1,b) + miu * Pb_matrix(1,b) * f_b(b) + Pb_matrix(1,b) * A.W_bh;
    end
    
    minimize( sum_fun )
    subject to
    for k = 1:6%B
        a1 = 0;
        for b = 1:4%K       
            a1 = a1 + Pb_matrix(1,b) * abs( hat_h_matrix(:,:,b,k) ) * abs( f_matrix(:,:,b,k) );
        end
        PHF(k) = a1;
    end 
    for k = 1:6%B     
            sqrt( sum( HF_sigma(:,k).^2 ) ) <= PHK(k);
    end
    for b = 1:4
        Pb_matrix(b) <= 1;
    end
    
cvx_end