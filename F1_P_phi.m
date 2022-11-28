function [P_phi_matrix] = F1_P_phi(A, P_phi_matrix, fq_fun )
P = zeros(1,A.R);


cvx_begin sdp
    variable P_phi_matrix(A.Nphi,A.R) ;
    expression P(A.R);
    minimize ( (-1) * fq_fun )
    subject to
    for j = 1:A.R%R
        for i = 1:A.Nphi%I
            P(j) = P(j) + P_phi_matrix(i,j);
        end
        P(j) == 1;
    end
    for j = 1:A.R%R
        for i = 1:A.Nphi%I
            P_phi_matrix(i,j) <= 1;
            P_phi_matrix(i,j) >= 0;
        end
    end
cvx_end
    
        