function [ Pr_matrix ] = F2_Pr( Pr_matrix, A, fq_fun )
a1 = 0;
a3 = 0;

cvx_begin
    variable Pr_matrix(1,5) nonnegative;
    expression sum_PWr(1);
    expression sum_PPWbh(1);
    for r = 1:5%R
        a1 = a1 + Pr_matrix(r) * A.W_r;
    end
    sum_PWr = a1;
    for b = 1:4
        a2 = 0;
        for r = 1:5
            a2 = a2 + Pr_matrix(r) * A.W_bh;
        end
        a3 = a3 + a2;
    end
    sum_PPWbh = a3;
    minimize ( ( -1 ) * fq_fun + A.alpha * sum_PWr + A.beta * sum_PPWbh )
    subject to 
    for r = 1:5
        Pr_matrix(r) <= 1;
    end
    
        