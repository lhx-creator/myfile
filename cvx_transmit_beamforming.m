function [f,t] = cvx_transmit_beamforming(B,K,h_head_H,noise_power,f,Nt,y,Pb,alpha)
%利用cvx求解波束成形子问题分时规划凸问题
%Nt=Nt;
cvx_begin quiet
variable t
variable F(Nt,1,B,K) complex
expression sumhf(K,1)
expression sum_abs_sumhf(K,1)
expression power_for_BSb(B,1) 
maximize (t)
subject to 
for k =1:K
    sumhf(k) = 0;
    sum_abs_sumhf(k) = 0;
    for b = 1:B
        sumhf(k) = sumhf(k) + h_head_H{b,k}*F(:,:,b,k);
    end
    for j = 1:K
        if j ~= k
            a = 0;  %中间变量，表示模值内的内积求和
            for b = 1:B
                a = a+h_head_H{b,k}*F(:,:,b,j);
            end
            a = square_abs(a);
            sum_abs_sumhf(k) = sum_abs_sumhf(k)+a;
        end
%         sum_abs_sumhf(k) = sum_abs_sumhf(k) + noise_power * square_pos(norm(y{k}, 'fro'));
    end
end
for b = 1:B
    power_for_BSb(b) = 0;
    for k = 1:K
       power_for_BSb(b) = power_for_BSb(b) + square_pos(norm(F(:,:,b,k), 'fro')); 
    end
end
for k = 1:K
    2*real(y{k}'*sumhf(k))-alpha(k)*(sum_abs_sumhf(k)+noise_power)*square_abs(y{k})>=t;
end
for b = 1:B
    power_for_BSb(b)<=Pb(b);
end
cvx_end
for b = 1:B
    for k=1:K
        f{b,k} = F(:,:,b,k);
    end
end

end

