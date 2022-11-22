clc,clear;
[A,L,F] = Parameter;
[hH_matrix,vH_matrix,G_matrix,bar_PHI_matrix,P_phi_matrix,hat_h_matrix] = Channel_model( A,L,F );
[Pb_matrix,f_matrix,Pr_matrix,miu,HF_sigma,fq_fun,T_initial,gamma,Z ] = Initialization( A, F, hat_h_matrix );

%第一次计算
[f_matrix] = L1_f( A, f_matrix, Pb_matrix, miu, hat_h_matrix, HF_sigma);
[Pb_matrix] = L2_Pb( Pb_matrix, f_matrix, miu, A, hat_h_matrix, HF_sigma );
[P_phi_matrix] = F1_P_phi( P_phi_matrix, fq_fun );
[ Pr_matrix ] = F2_Pr( Pr_matrix, A, fq_fun );


%第一次的函数值
[fun_leader] = Leader(f_matrix,Pb_matrix,A);
[ fun_follower ] = Follower( A, Pr_matrix, fq_fun );
fun_L = [];
fun_F = [];

for iter = 2:200 
    fun_L(1) = fun_leader;
    fun_F(1) = fun_follower;
    
    [f_matrix] = L1_f( A, f_matrix, Pb_matrix, miu, hat_h_matrix, HF_sigma);
    [Pb_matrix] = L2_Pb( Pb_matrix, f_matrix, miu, A, hat_h_matrix, HF_sigma );
    [P_phi_matrix] = F1_P_phi( P_phi_matrix, fq_fun );
    [ Pr_matrix ] = F2_Pr( Pr_matrix, A, fq_fun );
    [fun_L(iter)] = Leader(f_matrix,Pb_matrix,A);
    [fun_F(iter) ] = Follower( A, Pr_matrix, fq_fun );
    
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
    
    if ( x_L - fun_L(iter) )^2 < 10^(-3) && ( x_F - fun_F(iter) )^2 < 10^(-3)
        break;
    end
end



