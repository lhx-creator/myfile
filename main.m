clc,clear;
[A] = Parameter;
[hH_matrix,vH_matrix,G_matrix,bar_PHI_matrix,P_phi_matrix,hat_h_matrix,Pr_matrix] = Channel_model(A);
[Pb_matrix,f_matrix,miu,fq_fun,T_initial,gamma,Z ] = Initialization( A,hat_h_matrix );

%第一次计算
[f_matrix] = L1_f( A, f_matrix, Pb_matrix, miu, hat_h_matrix);
[Pb_matrix] = L2_Pb( Pb_matrix, f_matrix, miu, A, hat_h_matrix);
[P_phi_matrix] = F1_P_phi( A, P_phi_matrix, fq_fun );
[ Pr_matrix ] = F2_Pr( Pr_matrix, A, fq_fun );


%第一次的函数值
[fun_leader] = Leader(f_matrix,Pb_matrix,A);
[ fun_follower ] = Follower( A, Pr_matrix, fq_fun );
fun_L = [];
fun_F = [];
fun_L(1) = fun_leader;
fun_F(1) = fun_follower;
for iter = 1:20 
    
    [hat_h_matrix,fq_fun] = Iteration_fun(A,P_phi_matrix,bar_PHI_matrix,Pr_matrix,vH_matrix,G_matrix,hH_matrix,f_matrix,Pb_matrix );
    for i = 2:200
        [f_matrix] = L1_f( A, f_matrix, Pb_matrix, miu, hat_h_matrix);
        [Pb_matrix] = L2_Pb( Pb_matrix, f_matrix, miu, A, hat_h_matrix);
        [fun_L(i)] = Leader(f_matrix,Pb_matrix,A);
        if abs( fun_L(i) - fun_L(i-1) ) <= 1e-2
            break;
        end
    end
    min_L(iter) = fun_L(i);
    for i = 2:200
        [P_phi_matrix] = F1_P_phi( A, P_phi_matrix, fq_fun );
        [ Pr_matrix ] = F2_Pr( Pr_matrix, A, fq_fun );
        [fun_F(i) ] = Follower( A, Pr_matrix, fq_fun );
        if abs( fun_F(i) - fun_F(i-1) ) <= 1e-2
            break;
        end
    end              
    min_F(iter) = fun_F(i);
    
    if iter > 1
        if fun_L(iter) <= fun_L(iter-1)
            x_L = fun_L(iter);
        else
            x_L = fun_L(iter-1);
        end

        if fun_F(iter) <= fun_F(iter-1)
            x_F = fun_F(iter);
        else
            x_F = fun_F(iter-1);
        end

        if ( x_L - fun_L(iter) )/x_L < 10^(-6) && ( x_F - fun_F(iter) )/x_F < 10^(-6)
            break;
        end
    end
end



