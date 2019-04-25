function  l  = likehdder(th, Rg,Z,j)
%derivative of negative of loglikelihood
%   Z is the data matrix, th is the jth row of the parameter matrix
%   associated with jth conditional likelihood, j is the index of the
%   conditional likelihood to be used, m is the number of state spaces

z = Z;  p = size(Z,2);
x = Z(:,j); % identifies the jth variable in n observations
z(:,j) = 1;
X = repmat(x,1,p);
c = X.*z; % creates the derivative of the numerator of the likelihood for each observation
k1 = sum(c);

t = th;
% creating derivative of the proportional constant.
exp_mat = exp((z*t')*Rg);
k2 = sum(exp_mat,2);
k4 = exp_mat(:,2);
k6 = ( (k4./k2)' )*z;


l = -k1 + k6; % returns row vector

end

