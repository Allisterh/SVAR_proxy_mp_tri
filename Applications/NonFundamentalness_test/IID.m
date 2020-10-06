clear
clc


%ITRmax is the maximum number of itrations.
ITRmax = 500;


% Range of priliminary bandwidth Pbar.
Pmax = 8;
Pmin = 2;


%P0max is the maximum lag-order to be integrated numerically.
P0max=10;



%Vector of AIC and BIC lag order
aic = ones(ITRmax,1)*nan;
bic = aic;


% Matrix of the Test statistics and P-Values
Test = zeros(8,ITRmax);
PV = zeros(8,ITRmax);



%Parametrization of the bootstrap
alpha=0.36;
beta=0.99;
tau=0.25;
teta=alpha*beta*(1-tau);
kapa= tau*(1-teta)/(1-tau);







for ITR = 1:ITRmax
    
    %Dimension of VAR
    Dim = 2;
    
    
    %Sample size (I throw away 1000 observations to remove the effect of
    %intialization)
    T = 1000 + 250;
    
    
    
    
    %Data Matrices
    Y = zeros(T,Dim);
    h = Y;
    e = Y;
    
    
    
    %  DGP1
    q = lognrnd(0,1,T,1);
    z = lognrnd(0,1,T,1);
    
%     q = q - mean(q);
%     z = z - mean(z);
%     
%     for i = 20:T;
%         Y(i,1) = z(i);
%         Y(i,2) = alpha*Y(i-1,2) + q(i);
%     end;
    
    

    for i = 20:T;
        Y(i,1) = z(i);
        h(i) = 0.001 + 0.8 * h(i-1) + 0.02 * e(i-1)^2;
        e(i) = h(i)^(0.5) * q(i);
        Y(i,2) = alpha*Y(i-1,2) + e(i);
    end;


    x = Y(1001:T,:);
    [P, resid] =  VARfit(x, 8);
    DATA = resid;
    N = size(DATA,1);
    
    ST = std(DATA);
    DATA  = bsxfun(@rdivide, DATA,ST);

    
    %Bartlett Kernel
    WB = zeros(P0max);
    for p = 1:Pmax
        for j = 1:P0max
            WB(j,p) = Bar(j/(p+Pmin))^2;
        end
    end
    
    
    
    
    %Choice of data-driven lags
    M = zeros(Pmax,1);
    MD = M;
    MN = M;
    C = M;
    D = M;
    T = M;
    
    
    E0 = zeros(N);
    for a = 1:N
        for b = 1:N
            E0(a,b) = exp(-norm(DATA(a,:) - DATA(b,:))^2);
        end
    end
    ME0 = mean2(E0);
    
    
    E = zeros(N);
    for a = 1:N
        for b = 1:N
            E(a,b) = exp(-0.5*norm(DATA(a,:) - DATA(b,:))^2);
        end
    end
    ME = mean2(E);
    SRE = sum(E,2);
    SCE = sum(E);
    SRSQ = sum(SRE.^2);
    
    
    
    %Determinant of Phat
    for p = 1:Pmax
        for j = 1:P0max
            MD(p) = MD(p) + (N-j) * (  (N-j)^(-1)*sum(diag(E(1+j:N,1:N-j))) - (N-j)^(-2)*sum(sum(E(1+j:N,1:N-j))) )^2 * WB(j,p);
        end
    end
    
    
    %Numerator of Phat
    SJP = zeros(P0max,1);
    for j = 1:P0max
        S1 = 0;
        S2 = 0;
        for t = j+1:N
            for s = j+1:N
                S1 = S1 + E(t,s)*E(t-j,s-j);
                for r = j+1:N
                    S2 = S2 + E(t,s) * E(t-j,r-j);
                end
            end
        end
        SJP(j) = SJP(j) + (N-j)^(-1) * S1 - 2 * (N-j)^(-2) * S2 + (N-j)^(-3) * sum(sum(E(1+j:N,1+j:N))) * sum(sum(E(1:N-j,1:N-j)));
    end
    
    
    
    for p = 1:Pmax
        for j = 1:P0max
            MN(p) = MN(p) + 2*j^2 * WB(j,p) * SJP(j);
        end
    end
    
    MD0 =  N * (1-ME)^2;
    MD = 2* MD + MD0;
    Phat = 1.3934819 * (N * 2 * MN ./ MD) .^(1/5);
    
    %Bartlett Kernel
    WP = zeros(P0max);
    for p = 1:Pmax
        for j = 1:P0max
            WP(j,p) = Par(j/Phat(p))^2;
        end
    end
    
    %Computing the test statistic
    for i = 1:8
        for j = 1:P0max
            T(i) = T(i) + WP(j,i) * SJP(j);
        end
    end
    
    
    %Mean and Variance of the test statistic
    
    for i = 1:8
        for j = 1:P0max
            C(i) = C(i) + WP(j,i);
            D(i) = D(i) + WP(j,i).^2;
        end
    end
    
    C = C * (1-ME)^2;
    D = D * ( ME0 + ME^2 - 2*N^(-3)* SRSQ )^2;
    
    
    
    %Test Statistic
    for i = 1:8
        Test(i,ITR) = ( T(i)-C(i))./sqrt(2*D(i));
    end
    
    PV(:,ITR) = 1- normcdf(Test(:,ITR));
    
    ITR
    
end

PV1 = [sum((PV < 0.01),2), sum((PV < 0.05),2), sum((PV < 0.1),2)];