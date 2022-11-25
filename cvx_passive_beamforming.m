function [t,THETA] = cvx_passive_beamforming(PSI,p,N,R,K)
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明

% theta_head = [theta.',1].';
% THETA0 = theta_head*theta_head';
% normTHETA0 = norm(THETA0,2);
% [vec_Matrix,value_Matrix] = eig(THETA0);
% [max_eigvalue,max_eig_Index] = max(diag(value_Matrix));
% max_eig_Vec = vec_Matrix(:,max_eig_Index);
cvx_begin sdp quiet
variable t
variable THETA(N*R+1,N*R+1) complex semidefinite hermitian
maximize (t)
subject to
for k =1:K
    t + real(trace(PSI{k}*THETA)) + p{k}<=0;
end
diag(THETA)==ones(N*R+1,1);
% norm_nuc(THETA)-normTHETA0-real(trace(max_eig_Vec*max_eig_Vec'*(THETA-THETA0)))<=0;

cvx_end
end

