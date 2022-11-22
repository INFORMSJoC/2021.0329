function SKmodel = SKmodelfit(X, Y, F, Vhat, kernind)
% fit a stochastic kriging model for simulation with covariates

numX = size(X,1);
% diagonal intrinsic variance matrix
V = diag(Vhat);
% initialize parameters of tau2
beta = (F'*F)\(F'*Y);
tau2_0 = var(Y-F*beta);

distmat = [];
for i = 1:numX
    distmat = [distmat; sum((X(i,:) - X((i+1):end, :)).^2, 2)];
end
average_distance = mean(distmat);
% initialize parameters of theta
if kernind == 2  % squared exponential
    theta_0 = sqrt(2*average_distance/log(2));
elseif kernind == 1  % exponential
    theta_0 = sqrt(average_distance)/log(2);
elseif kernind == 3  % matern 3/2
    theta_0 = sqrt(3*average_distance)/log(2);
elseif kernind == 4  % matern 5/2
    theta_0 = sqrt(5*average_distance)/log(2);
end

% lower bounds for numerical stability
lbtau = sqrt(0.00000001*tau2_0);    
lbtheta = 0.000001; 
lb = [lbtheta;lbtau];

% optimize the log-likelihood function
options = optimoptions('fmincon','MaxFunEvals',1000000,...
    'MaxIter',500,'Display','notify');
parms = fmincon(@(x) loglikelihood(x,numX,X,F,V,Y,kernind),...
        [theta_0;sqrt(tau2_0)],[],[],[],[],lb,[],[],options); 

SKmodel.tausquared =  (parms(2))^2;
SKmodel.theta = parms(1);
Sigma2  =  Kernval( X,X,parms,kernind ) + V;
Lhat = chol(Sigma2)';
Lhatinv = inv(Lhat);
Sigma2hatinv = Lhatinv'*Lhatinv;
FSigma2inv = F' * Sigma2^(-1);
SKmodel.beta = (FSigma2inv*F) \ (FSigma2inv*Y);
SKmodel.Sigma2inv = Sigma2hatinv;
SKmodel.kernind = kernind;
SKmodel.Znew = Sigma2\(Y-F*beta);
SKmodel.F = F;
SKmodel.X = X;

end



function y = loglikelihood(parms,numX,X,F,V,Y,kernind)

Sigma2  =  Kernval( X,X,parms,kernind ) + V;
eigen_S = eig(Sigma2);
Sigma2inv = Sigma2^(-1);
FSigma2inv = F' * Sigma2inv;
coef = (FSigma2inv*F) \ (FSigma2inv*Y);
y = 0.5*numX*log(2*pi) + 0.5*sum(log(eigen_S)) + 0.5*(Y-F*coef)'*Sigma2inv*(Y-F*coef);
end

