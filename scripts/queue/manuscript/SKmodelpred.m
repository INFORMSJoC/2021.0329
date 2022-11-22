function y = SKmodelpred(SKmodel,Xpred,Fpred)

X = SKmodel.X;
theta = SKmodel.theta;
kernind = SKmodel.kernind;
Znew = SKmodel.Znew;
beta = SKmodel.beta;
tau = sqrt(SKmodel.tausquared);

Spred = Kernval(Xpred, X, [theta,tau], kernind);
y = Fpred*beta + Spred'*Znew;
    
end
