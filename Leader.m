function [fun_leader] = Leader(f_matrix,Pb_matrix,A)

miu = 0.5;

f_b = zeros(1,A.K);
fun_leader = 0;
for b = 1:A.B
    a = 0;
    
    for k = 1:A.K
        a = a + square_pos(norm(f_matrix(:,:,b,k), 'fro'));
    end
    f_b(b) = a;
    fun_leader = fun_leader + Pb_matrix(1,b) * A.W_b + miu * Pb_matrix(1,b) * f_b(b) + Pb_matrix(1,b) * A.W_bh;
end



