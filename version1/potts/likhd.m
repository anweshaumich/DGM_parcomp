function y = likhd(th, Rg,h_matj,H_mat)
% jth negative conditional loglikelihood

n = size(h_matj,1);
m = length(Rg);
k1 = sum(h_matj,1)* th';
mat_2 = H_mat*th';% mat_2 = mat_2 - max(mat_2); 
mat_2 = reshape(mat_2,m,n)';
k2 = exp(mat_2);
k2 = sum(k2,2);
k3 = sum(log(k2));


y = -k1 + k3 ;


end

