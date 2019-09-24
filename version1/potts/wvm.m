function  [ec, w]  = wvm(th,h_matj,H_mat,Rg)
% th is (theta with indicator 1 and rest 0). it is a row vector
[n,p] = size(h_matj);
m = length(Rg);
mat_1 = H_mat*th'; mat_1 = mat_1 - max(max(mat_1));
mat_2 = exp(mat_1);
mat_3 = repmat(mat_2,1,p);
mat_4 = H_mat.*mat_3;
mat_5 = sm(mat_4,m,n);%blockproc(mat_4,[m p],@(k) sum(k.data))
mat_6 = reshape(mat_2,m,n);
mat_7 = sum(mat_6,1); mat_8 = repmat(mat_7',1,p);
k1 = mat_5./mat_8; k2 = repmat(k1,1,m); ec = reshape(k2',p,m*n)';
mat_9 = repmat(mat_7',1,m); mat_10 = reshape(mat_9',m*n,1);
w = mat_2./mat_10;
%w = sparse(1:m*n,1:m*n,mat_11,m*n,m*n);
end




