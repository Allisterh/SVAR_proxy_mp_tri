clear all
% warning off
lag = 4;      %% lag used in VAR estimation
h = 1:25;     %% arbitrarily selected autocorrelation order for classic BP test Q(h), chosen by user
alpha = 0.05; %% significance level
q = 2.4;      %% parameter used for deciding which information criteria to use

I = eye(length(h));
Qh = I(h>lag,:)*h';   %% Q(h) can only be used for h > lag

d0 = 12;              %% upper bound of lag selection in the data-driven test for different sample size
d = max(max(h),max(d0));
dd = d0;



res = load('resid.txt');    %%load your data for VAR(lag) estimation

%data = load('dat.txt');
%[B_hat, res, Z] = VARfit(data,lag);  

%res = res'

[m, n] = size(res); % m = observations (T), n = VAR dimension (K)
sigmaE = 1/n*(res*res');
rtsigmaE = inv((eye(m).*sigmaE).^0.5);
%%
X = res(1:(m-1),:)';
t = 1;
while t<=m
    M = X(:,(t+lag-1):-1:t);
    Z(:,t) = [1; reshape(M,lag*n,1)];
    t = t+1;
end



%%

Omga = 1/n*(Z*Z');        
H = [zeros(1,m); sigmaE;zeros(m*(lag-1),m)];
for j=1:d
    e = res(:,(j+1):end);            
    eJ = res(:,1:(n-j));
    G = 1/(n)*e*eJ';
    Q(j,1) = n*sum(diag(G'*inv(sigmaE)*G*inv(sigmaE)));      %% BP 
    QA(j,1) = n/(n-j)*Q(j,1);                                %% LB version
    lamda(j,1) = max(max(abs(rtsigmaE*G*rtsigmaE)));       
end

QAstat = tril(ones(d))*QA;
        
%%%% Q(1) is distributed as weighted sum of chi2 random variables with weights being the eigenvalue of the 
%%%% asymptotic variance of gamma(1) in the paper.
%%%% First estimate the asymptotic variance
e = res(:,2:end);
e1 = res(:,1:(n-1));
ZJ = Z(:,2:end);
sigmaC = 0;
sigmaB = 0;
sigmaCB = 0;
for j = 1:(n-1)
    sigmaC = sigmaC + 1/(n-1)*kron(e1(:,j)*e1(:,j)',e(:,j)*e(:,j)');
    sigmaB = sigmaB + 1/(n-1)*kron(inv(Omga)*ZJ(:,j)*ZJ(:,j)'*inv(Omga),e(:,j)*e(:,j)');
    sigmaCB = sigmaCB + 1/(n-1)*kron(e1(:,j)*ZJ(:,j)'*inv(Omga),e(:,j)*e(:,j)');
end
D = sigmaC + kron(H',eye(m))*sigmaB*kron(H,eye(m)) - sigmaCB*kron(H,eye(m)) - kron(H',eye(m))*sigmaCB';
A = inv(sqrtm(sigmaE));
DD = kron(A,A)*D*kron(A,A);
        
eigDD = eig(DD);
               
%%%% Calculate the data-driven test statistic 
if max(n^0.5*lamda(1:dd)) <= ((q*log(n))^0.5)            
    Lp = QAstat(1:dd) - m^2*log(n)*(1:dd)';            
else
    Lp = QAstat(1:dd) - 2*m^2*(1:dd)';
end
[L,h_hat] = max(Lp);
AQ = QAstat(h_hat);       
              
%%%% Compute the pvalue
Qpv = 1-chi2cdf(QAstat(Qh),m^2*(Qh-lag));    %% p-value of Q(h)
      
global w x                 %% global variable used in the Imhof function
w = eigDD;
x = AQ;
AQpv = 0.5 + quadl(@Imhof,0,10000)/pi;
       
      
    