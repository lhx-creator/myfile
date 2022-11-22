function [P_phi_matrix] = F1_P_phi( P_phi_matrix, fq_fun )
P = zeros(1,5);


cvx_begin
    variable P_phi_matrix(10,5) semidefinite;
    minimize ( (-1) * fq_fun )
    subject to
    for j = 1:5%R
        for i = 1:10%I
            P(j) = P(j) + P_phi_matrix(i,j);
        end
        P(j) =1;
    end
    for j = 1:5%R
        for i = 1:10%I
            P_phi_matrix(i,j) <= 1;
        end
    end
cvx_end
    
        