function  l  = likehdder(th, Rg,h_matj,H_mat)
%derivative of negative of jth conditional loglikelihood

[n,p] = size(h_matj);
m = length(Rg);
k1 = sum(h_matj,1);
mat_1 = H_mat*th'; mat_1 = mat_1 - max(max(mat_1));
mat_2 = exp(mat_1);
mat_3 = repmat(mat_2,1,p);
mat_4 = H_mat.*mat_3;
mat_5 = sm(mat_4,m,n);%blockproc(mat_4,[m p],@(k) sum(k.data));
mat_6 = reshape(mat_2,m,n);
mat_7 = sum(mat_6,1); mat_8 = repmat(mat_7',1,p);
mat_9 = mat_5./mat_8; k2 = sum(mat_9);

l = -k1 + k2;% returns row vector



end

