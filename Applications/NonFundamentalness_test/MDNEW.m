clear
clc
tic

% Range of priliminary bandwidth Pbar.
Pmax = 8;
Pmin = 2;


%P0max is the maximum lag-order to be integrated numerically.
P0max=Pmin+Pmax;


% Matrix of the Test statistics and P-Values(using Parzen, Daniell, Quadratic
% Kernels)
TP = zeros(Pmax,1);
TD = zeros(Pmax,1);
TQS = zeros(Pmax,1);


x = load('BPADJ.txt');  %%load your data for VAR(p) estimation;
[P, resid] =  VARfit(x, 4);
DATA = resid;
[N Dim]= size(DATA);



ST = std(DATA);
DATA  = bsxfun(@rdivide, DATA,ST);

SMat = zeros(2);
for i = 1:Dim;
    for j = 1:Dim;
        SMat(i,j) = mean(bsxfun(@times, DATA(:,i),DATA(:,j)));
    end;
end;


SS2=0;
for i = 1:Dim;
    for j = 1:Dim;
        SS2 = SS2 + SMat(i,j)^2;
    end;
end;


DS = zeros(N);
for t = 1:N;
    for s = 1:N;
        for d = 1:Dim;
            DS(t,s) = DS(t,s) + DATA(t,d)*DATA(s,d);
        end;
    end;
end;



COVDATA = zeros(1,P0max+1);
Z = zeros(Dim);
for i = 0:P0max
    COVDATA(i+1) = COVX(DATA,i,N,Z);
end





%Choice of data-driven lags
MD = zeros(Pmax,1);


MPar = zeros(Pmax,1);
CP = MPar;
DP = MPar;

MDan = zeros(Pmax,1);
CD = MDan;
DD = MDan;

MQS = zeros(Pmax,1);
CQS = MQS;
DQS = MQS;

E0 = zeros(N);
for a = 1:N
    for b = 1:N
        E0(a,b) = exp(-norm(DATA(a,:) - DATA(b,:))^2);
    end
end
ME0 = mean2(E0);


E = zeros(N);
for a = 1:N;
    for b = 1:N;
        E(a,b) = exp(-0.5*norm(DATA(a,:) - DATA(b,:))^2);
    end;
end;
ME = mean2(E);
SRE = sum(E,2);
SCE = sum(E);
SRSQ = sum(SRE.^2);



%Bartlett Kernel
WB = zeros(P0max);
for i = 1:Pmax;
    for j = 1:(i+Pmin);
        WB(j,i) = Bar(j/(i+Pmax))^2;
    end;
end;



%Determinant of Phat
M0 = N * COVDATA(1) * (1-ME);
M1 = zeros(1,P0max);
for j = 1:P0max;
    M1(j) = M1(j) + 2 * COVDATA(j+1) * ( (N-j)* sum(diag(E(1+j:N,1:N-j))) - sum(sum(E(1+j:N,1:N-j)))) *(N-j)^(-1);
end;


for i=1:Pmax;
    for j = 1:(i+Pmin);
        MD(i) = MD(i) + M1(j)*WB(j,i);
    end;
end;




%Numerator of Phat
SJP = zeros(P0max,1);
SD = zeros(P0max,1);
for j = 1:P0max;
    S1 = 0;
    S2 = 0;
    S3 = 0;
    SD1 = 0;
    SD2 = 0;
    for t = (j+1):N;
        for s = (j+1):N;
            S1 = S1 + DS(t,s) * E(t-j,s-j);
            S2 = S2 + DS(t,s) * sum(E(t-j,1:N-j),2);
            S3 = S3 + DS(t,s) * sum(sum(E(1:N-j,1:N-j)));
            SD1 = SD1 + E(t,s) * E(t-j,s-j);
            for r = (1+j):N;
                SD2 = SD2 + E(t,s) * E(t-j,r-j);
            end;
        end;
    end;
    SJP(j) = (N-j)^(-1) * S1 - 2 * (N-j)^(-2) * S2 + (N-j)^(-3) *S3;
    SD(j) = (N-j)^(-2) * SD1 - 2*(N-j)^(-3) * SD2 + (N-j)^(-4) * sum(sum(E(1+j:N,1+j:N))) * sum(sum(E(1:N-j,1:N-j)));
end;


MN = zeros(Pmax,1);
for i = 1:Pmax;
    for j = 1:(i+Pmin);
        MN(i) = MN(i) + 2*j^2 * WB(j,i) * SJP(j);
    end;
end;

MD = MD + M0;

PhatP = 1.3934819 * (N * 2 * MN ./ MD) .^(1/5);
PhatP = (PhatP > log(N)).*log(N)  + (PhatP <= log(N)).*PhatP;

PhatD = 1.4017333 * (N * 2 * MN ./ MD) .^(1/5);
PhatD = (PhatD > log(N)).*log(N)  + (PhatD <= log(N)).*PhatD;

PhatQS = 1.4223747 * (N * 2 * MN ./ MD) .^(1/5);
PhatQS = (PhatQS > log(N)).*log(N)  + (PhatQS <= log(N)).*PhatQS;

PDim = ceil(max(PhatP,[],1));
PD = min(PDim, P0max);


%Parzen Kernel
WP = zeros(P0max);
WD = zeros(P0max);
WQS = zeros(P0max);
for i = 1:Pmax;
    for j = 1:Pmax;
        WP(j,i) = Par(j/PhatP(i))^2;
        WD(j,i) = Dan(j/PhatD(i))^2;
        WQS(j,i) = QS(j/PhatQS(i))^2;
    end;
end;



%Computing the test statistic
for i = 1:Pmax;
    for j = 1:PD;
        MPar(i) = MPar(i) + WP(j,i) * SJP(j);
        MDan(i) = MDan(i) + WD(j,i) * SJP(j);
        MQS(i) = MQS(i) + WQS(j,i) * SJP(j);
    end;
end;


%Mean and Varianceof the test statistic under Conditional Heteroscedasticity
C1 = zeros(PD,1);

for j = 1:PD
    for t = (j+1):(N-1)
        C1(j) = C1(j) +  norm(DATA(t,:))^2 *( 1 - 2 * mean(exp(-0.5 * (sum(abs(bsxfun(@minus,DATA,DATA(t-j,:))) .^2 , 2)))) + ME );
    end
end


D2P = zeros(Pmax,1);
D2D = zeros(Pmax,1);
D2QS = zeros(Pmax,1);

for i = 1:Pmax;
    for j = 1:PD;
        CP(i) = CP(i) + WP(j,i) * (N-j)^(-1) * C1(j);
        CD(i) = CD(i) + WD(j,i) * (N-j)^(-1) * C1(j);
        CQS(i) = CQS(i) + WQS(j,i) * (N-j)^(-1) * C1(j);
        D2P(i) = D2P(i) +  WP(j,i)^2;
        D2D(i) = D2D(i) +  WD(j,i)^2;
        D2QS(i) = D2QS(i) +  WQS(j,i)^2;
    end;
end;

D20P = D2P * ( ME0 + ME^2 - 2*N^(-3)* SRSQ );
D20D = D2D * ( ME0 + ME^2 - 2*N^(-3)* SRSQ );
D20QS = D2QS * ( ME0 + ME^2 - 2*N^(-3)* SRSQ );




for i = 1:Pmax;
    for j = 1:PD;
        for l = 1:PD;
            if (j~=l);
                DP(i) = DP(i) + WP(j,i)*WP(l,i)*SD(abs(j-l));
                DD(i) = DD(i) + WD(j,i)*WD(l,i)*SD(abs(j-l));
                DQS(i) = DQS(i) + WQS(j,i)*WQS(l,i)*SD(abs(j-l));
            end;
        end;
    end;
end;

DP = SS2 * (DP + D20P);
DD = SS2 * (DD + D20D);
DQS = SS2 * (DQS + D20QS);



%Test Statistic Conditional Heteroscedasticity
for i = 1:Pmax;
    TP(i) = ( MPar(i) - CP(i) )./sqrt( 2*DP(i) );
    TD(i) = ( MDan(i) - CD(i) )./sqrt( 2*DD(i) );
    TQS(i) = ( MQS(i) - CQS(i) )./sqrt( 2*DQS(i) );
end;

PVP = 1- normcdf(TP);
PVD = 1- normcdf(TD);
PVQS= 1- normcdf(TQS);




