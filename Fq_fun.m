function [ fq_fun ] = Fq_fun( tilde_f_matrix, tilde_h_matrix, A, gamma, Z )

fq_fun = [];
tilde_fhH = zeros(1,6);
tilde_fh = zeros(1,6);
sum_fhhf = zeros(1,6);

for k = 1:6  
    for b = 1:4%B
        tilde_fhH(k) = tilde_fhH(k) + ctranspose( tilde_f_matrix(:,:,b,k) ) * ctranspose( tilde_h_matrix(:,:,b,k) );
    end
end

for k = 1:6  
    for b = 1:4%B
        tilde_fh(k) = tilde_fh(k) + tilde_h_matrix(:,:,b,k) * tilde_f_matrix(:,:,b,k);
    end
end
for g = 1:6
    for k = 1:6
    sum_fhhf(g) = sum_fhhf(g) + tilde_fh(k) * tilde_fhH(k);
    end
    sum_fhhf(g) = sum_fhhf(g) + A.sigma(g)^2;
end


for k = 1:6%K
    fq_fun  = fq_fun + 2 * sqrt( 1 + gamma(k) ) * real( sum ( tilde_fhH ) * Z(k) ) - ctranspose( Z(k) ) * sum_fhhf(k) * Z(k) + log2( 1 + gamma(k) ) - gamma(k);
       %取实部求和有点问题感觉
end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    