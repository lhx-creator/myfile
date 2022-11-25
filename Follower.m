function [ fun_follower ] = Follower( A, Pr_matrix, fq_fun )

a1 = 0;
a3 = 0;

for r = 1:A.R%R
    a1 = a1 + Pr_matrix(r) * A.W_r;
end
sum_PWr = a1;
for b = 1:A.B
    a2 = 0;
    for r = 1:A.R
        a2 = a2 + Pr_matrix(r) * A.W_bh;
    end
    a3 = a3 + a2;
end
sum_PPWbh = a3;

fun_follower = ( -1 ) * fq_fun + A.alpha * sum_PWr + A.beta * sum_PPWbh




