function [f_matrix] = L1_f( f_matrix, Pb_matrix, miu, hat_h_matrix, HF_sigma)

square_sum_f = zeros(1,4);
cons_square_sum_f(4) = zeros(1,4);
PHF = zeros(1,6);


%%
% f 优化子问题
cvx_begin
    variable f_matrix(4,1,4,6) complex 
    expression square_sum_f(4)%B
    expression sum_fun
    expression cons_square_sum_f(4)%B

    
    for b = 1:4%B
        a1 = 0;
        for k = 1:6%K       
            a1 = a1 + (abs(f_matrix(:,:,b,1))).^2 ;
        end
        square_sum_f(b) = a1;
    end
    for b = 1:4%B
        a2 = 0;
        a2 = a2 + Pb_matrix(1,b) * square_sum_f(b);
    end
    sum_fun = a2;
    minimize( miu * sum_fun )
    subject to
    for b = 1:4%B
        a3 = 0;
        for k = 1:6%K       
            a3 = a3 + (abs(f_matrix(:,:,b,1))).^2 ;
        end
        cons_square_sum_f(b) = a3;
        cons_square_sum_f(b) <= A.W_max;  
    end
    for k = 1:6%B
        a4 = 0;
        for b = 1:4%K       
            a4 = a4 + Pb_matrix(1,b) * abs( hat_h_matrix(:,:,b,k) ) * abs( f_matrix(:,:,b,k) );
        end
        PHF(k) = a4;
    end 
    for k = 1:6%B     
            sqrt( sum( HF_sigma(:,k).^2 ) ) <= PHK(k);
    end
cvx_end

