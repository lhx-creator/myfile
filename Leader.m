function [f_matrix,Pb_matrix] = Leader(L,A)

miu = 0.5;
%Pb_matrix = ones(1,L.B);
Pb_matrix = ones(1,4);
%f_matrix = zeros(10,10,L.B,A.K);%L.B A.K
f_matrix = zeros(10,1,4,6);

%for m = 1:L.B * A.K%L.B    
for m = 1:24    
        eval(['f',num2str(m), ' =rand(10,1)']);
end
%for i = 1:L.B * A.K
 for i = 1:24
    f_matrix(:,:,i)= eval(['f', num2str(i)]);
end


