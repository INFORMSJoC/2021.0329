
clearvars,clc;

% MM1 Queue
fun = @QueueSim; 
n0 = 10; n2 = 5; 
span1 = 4; span2 = 0.5;
desig = 6 + 0.3*(1:10); desig = desig'; 
nd = length(desig);  %location of design point
Lm = 100;    % number of replications for random covariates point

cu = 0.1; U = 2.5; 

m = [10,15,23,35,53,80]; nm = length(m);    % number of randomly selected covariate
%%% PFS
testa = 40000;    % number of test covariate point for PFS
Del = 0.01;
%%%%%%%
% lose the error tolerance in finding the optimal point
options = optimoptions('fmincon','OptimalityTolerance',0.005,...
    'StepTolerance',0.005,'Display','off');
optseed = 10;  
d = 1;  % context dim.
% %%%%%%%%%%%%%%%%%%%% Sample distribution: Uniform %%%%%%%%%%%%%%%%%%%%%%
thres1 = 0.0002;  % imse target of exponential kernel
thres2 = 0.000075;  % imse target 
[ UPmse1_all,UPmse2_all,UPmse3_all,...
    UPmse4_all,pUPmse1,pUPmse2,pUPmse3,pUPmse4,pmU_mse,m_thresU,...
    UPnum1,UPnum2,UPnum3,UPnum4] ...
    = UniSQ_reg_self( m,desig,span1,span2,testa,Del,n0,cu,U,thres1,thres2 );

% %%%%%%%%%%%%%%%%%%%% Sample distribution: Truncated Normal %%%%%%%%%%%%%%%%%%%%%%
Truva = 3;
[ TPmse1_all,TPmse2_all,TPmse3_all,...
    TPmse4_all,pTPmse1,pTPmse2,pTPmse3,pTPmse4,pmT_mse,m_thresT,...
    TPnum1,TPnum2,TPnum3,TPnum4] ...
    = TrnSQ_reg_self( m,desig,span1,span2,testa,Del,n0,Truva,cu,U,thres1,thres2 );

save testq_reg;

pmT_mse  %coef
m_thresT  % m0
mean(TPnum4), median(TPnum4) 
mean(TPnum3), median(TPnum3)
mean(TPnum2), median(TPnum2)
mean(TPnum1), median(TPnum1)

m_thresU
pmU_mse
mean(UPnum4), median(UPnum4)
mean(UPnum3), median(UPnum3)
mean(UPnum2), median(UPnum2)
mean(UPnum1), median(UPnum1)







