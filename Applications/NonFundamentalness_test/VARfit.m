function [B, res, Z, sigmaB]=VARfit(data,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function tries to fit the VAR model with lag p, and returns parameter
%estimates and residual.
% Input: data, with dimension T x d, T is the sample size and d is
%        the dimension of the vector;
%        p, lag of the VAR model;
% Output: theta, the parameter estimates;
%         res,   residual of the fitted model;
%         Z,     design matrix in vector form
%         sigmaB, asymtotic variance of estimator B
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T = size(data,1);
d = size(data,2);
if T<d
    data = data';
    [T,d] = size(data);
end

Y = data((p+1):end,:)';        % response variable
%y = reshape(Y,d*(T-p),1);      % vectorize the response variable
%%%%generate regressors
X = data(1:(T-1),:)';
t = 1;
while t<=(T-p)
    M = X(:,(t+p-1):-1:t);  % Reverse the order of elements
    Z(:,t) = [1; reshape(M,p*d,1)];
%    Z(:,t) = reshape(M,p*d,1);
    t = t+1;
end

%%%%OLS estimation
%beta = kron(inv(Z*Z')*Z,eye(d))*y;   % vectorization of the matrix coefficient
%B = reshape(beta,d,d*p+1);
B = Y*Z'*inv(Z*Z');
res = Y - B*Z;
res = res';
B = B';
