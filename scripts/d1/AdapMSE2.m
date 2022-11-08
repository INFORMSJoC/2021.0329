
clearvars,clc;

% fun: Griewank
fun = @testfun2;
desig = 1:10; desig = desig'; 
nd = length(desig);  %location of design point
n0 = 10;
va = 2; sd = sqrt(va);
span1 = 9; span2 = 1;    
scen = 5.5; sva = 3; ssd = sqrt(sva); ssd2 = ssd; 
Lm = 100;    % number of replications for random covariates point

m = [5,8,12,18,28,42,65,100]; nm = length(m);    % number of randomly selected covariate
%%% PFS
testa = 1000;    % number of test covariate point for PFS
Del = 0.10;
%%%%%%%
options = optimoptions('fmincon','OptimalityTolerance',0.01,...
    'StepTolerance',0.01,'Display','off');
optseed = 10;  
d = 1;  % context dim.
%%%%%%%%%%%%%%%%%%% Sample distribution: Uniform %%%%%%%%%%%%%%%%%%%%%%
[ UPmse1,UPpfs1,UPmse2,UPpfs2,UPmse3,UPpfs3,UPmse4,UPpfs4 ] ...
    = UniS( fun,m,desig,span1,span2,Lm,testa,Del,sd,n0 );

%%%%%%%%%%%%%%%%%%% Sample distribution: Truncated Normal %%%%%%%%%%%%%%%%%%%%%%
Tmu0 = 5.5; Trusd0 = 7; 
[ TTPmse1,TTPpfs1,TTPmse2,TTPpfs2,TTPmse3,TTPpfs3,TTPmse4,TTPpfs4 ] ...
    = TrnS( fun,m,desig,span1,span2,Lm,testa,Del,sd,n0,Tmu0,Trusd0 );

%%%%%%%%%%%%%%%%%%%% Sample distribution: Norm %%%%%%%%%%%%%%%%%%%%%%
[ NPmse1,NPpfs1,NPmse2,NPpfs2,NPmse3,NPpfs3,NPmse4,NPpfs4 ] ...
    = NorS( fun,m,desig,scen,ssd,Lm,testa,Del,sd,n0,ssd2 );

%%%%%%%%%%%%%%%%%%%% Adaptive MSE %%%%%%%%%%%%%%%%%%%%%%
Tmu = 5.5; Trusd = 1;   % Truva = 1; Tmu = 5.5;
pd = makedist('Normal','mu',Tmu,'sigma',Trusd); 
Bd1 = span2; Bd2 = span2+span1;
Tpd = truncate(pd,Bd1,Bd2);
rng(500)
Etestx = random(Tpd,testa,d);   % test covariate generation for PFS
Etesty = fun(desig,Etestx);    % exact testfun value
Eargm = min(Etesty,[],2);    % true optimal design

[ TPmse1,TPpfs1,TPmse2,TPpfs2,TPmse3,TPpfs3,TPmse4,TPpfs4 ] ...
    = TrnS( fun,m,desig,span1,span2,Lm,testa,Del,sd,n0,Tmu,Trusd );

kernind = 1;
[ EPmse1,EPpfs1,Econshortrec1 ] ...
    = AdaS( fun,m,desig,span1,span2,Lm,testa,Del,sd,n0,kernind,Etestx,...
    Etesty,Eargm,options,optseed,d);
kernind = 2;
[ EPmse2,EPpfs2,Econshortrec2 ] ...
    = AdaS( fun,m,desig,span1,span2,Lm,testa,Del,sd,n0,kernind,Etestx,...
    Etesty,Eargm,options,optseed,d);
kernind = 3;
[ EPmse3,EPpfs3,Econshortrec3 ] ...
    = AdaS( fun,m,desig,span1,span2,Lm,testa,Del,sd,n0,kernind,Etestx,...
    Etesty,Eargm,options,optseed,d);
kernind = 4;
[ EPmse4,EPpfs4,Econshortrec4 ] ...
    = AdaS( fun,m,desig,span1,span2,Lm,testa,Del,sd,n0,kernind,Etestx,...
    Etesty,Eargm,options,optseed,d);

%%%%%%%%%%%%%%%%%%%%%% adaptive MSE under uniform distribution %%%%%%%%%%%%%%%%
[ UUPmse1,UUPpfs1,UUPmse2,UUPpfs2,UUPmse3,UUPpfs3,UUPmse4,UUPpfs4 ] ...
    = UniS( fun,m,desig,span1,span2,Lm,testa,Del,sd,n0 );
rng(500)
Etestx = rand(testa,d)*span1+span2;  
Etesty = fun(desig,Etestx);    % exact testfun value
Eargm = min(Etesty,[],2);    % true optimal design
kernind = 1;
[ EEPmse1,EEPpfs1,EEconshortrec1 ] ...
    = AdaS( fun,m,desig,span1,span2,Lm,testa,Del,sd,n0,kernind,Etestx,...
    Etesty,Eargm,options,optseed,d);
kernind = 2;
[ EEPmse2,EEPpfs2,EEconshortrec2 ] ...
    = AdaS( fun,m,desig,span1,span2,Lm,testa,Del,sd,n0,kernind,Etestx,...
    Etesty,Eargm,options,optseed,d);
kernind = 3;
[ EEPmse3,EEPpfs3,EEconshortrec3 ] ...
    = AdaS( fun,m,desig,span1,span2,Lm,testa,Del,sd,n0,kernind,Etestx,...
    Etesty,Eargm,options,optseed,d);
kernind = 4;
[ EEPmse4,EEPpfs4,EEconshortrec4 ] ...
    = AdaS( fun,m,desig,span1,span2,Lm,testa,Del,sd,n0,kernind,Etestx,...
    Etesty,Eargm,options,optseed,d);


save test2;

figure
subplot(1,2,1)
plot(m,TPmse1,'*-')
hold on;
plot(m,TPmse2,'o-')
plot(m,TPmse3,'+-')
plot(m,TPmse4,'d-')
ylabel('Maximal IMSE','Interpreter','latex')
xlabel('m','Interpreter','latex')
title('Griewank (1D): Static Sampling','Interpreter','latex')
legend('Exp','Sq-Exp','Matern3/2','Matern5/2');
set(gca,'XTick', m, 'Xscale','log', 'Yscale','log', 'XMinorTick','off', 'YMinorTick','off')
set(gca,'XMinorGrid','off', 'YMinorGrid','off', 'XGrid','on', 'YGrid','on','GridLineStyle',':')
yt = yticks;
yl = ylim;



subplot(1,2,2)
plot(m,EPmse1,'*-')
hold on;
plot(m,EPmse2,'o-')
plot(m,EPmse3,'+-')
plot(m,EPmse4,'d-')
ylabel('Maximal IMSE','Interpreter','latex')
xlabel('m','Interpreter','latex')
title('Griewank (1D): Adaptive MSE','Interpreter','latex')
legend('Exp','Sq-Exp','Matern3/2','Matern5/2');
set(gca,'Ylim',yl,'YTick',yt)
set(gca,'XTick', m, 'Xscale','log', 'Yscale','log', 'XMinorTick','off', 'YMinorTick','off')
set(gca,'XMinorGrid','off', 'YMinorGrid','off', 'XGrid','on', 'YGrid','on','GridLineStyle',':')





figure
subplot(1,2,1)
plot(m,TPpfs1,'*-')
hold on;
plot(m,TPpfs2,'o-')
plot(m,TPpfs3,'+-')
plot(m,TPpfs4,'d-')
ylabel('IPFS','Interpreter','latex')
xlabel('m','Interpreter','latex')
title('Griewank (1D): Static Sampling','Interpreter','latex')
legend('Exp','Sq-Exp','Matern3/2','Matern5/2');
set(gca,'Ylim',[0.05,0.72])
set(gca,'XTick', m, 'xscale','log', 'yscale', 'log', 'XMinorTick','off', 'YMinorTick','off')
set(gca,'XMinorGrid','off', 'YMinorGrid','off', 'XGrid','on', 'YGrid','on','GridLineStyle',':')
yt = yticks;
yl = ylim;


subplot(1,2,2)
plot(m,EPpfs1,'*-')
hold on;
plot(m,EPpfs2,'o-')
plot(m,EPpfs3,'+-')
plot(m,EPpfs4,'d-')
ylabel('IPFS','Interpreter','latex')
xlabel('m','Interpreter','latex')
title('Griewank (1D): Adaptive MSE','Interpreter','latex')
legend('Exp','Sq-Exp','Matern3/2','Matern5/2');
set(gca,'XTick', m, 'xscale','log', 'yscale', 'log', 'XMinorTick','off', 'YMinorTick','off')
set(gca,'XMinorGrid','off', 'YMinorGrid','off', 'XGrid','on', 'YGrid','on','GridLineStyle',':')
set(gca,'Ylim',yl,'YTick',yt)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
funname = 'Griewank (1D)';
subplot(1,3,3)
plot(m,NPmse1,'*-')
hold on;
plot(m,NPmse2,'o-')
plot(m,NPmse3,'+-')
plot(m,NPmse4,'d-')
ylabel('Maximal IMSE','Interpreter','latex')
xlabel('m','Interpreter','latex')
title(strcat(funname,': Normal Sampling'),'Interpreter','latex')
legend('Exp','Sq-Exp','Matern3/2','Matern5/2');
set(gca,'XTick', m, 'Xscale','log', 'Yscale','log', 'XMinorTick','off', 'YMinorTick','off')
set(gca,'XMinorGrid','off', 'YMinorGrid','off', 'XGrid','on', 'YGrid','on','GridLineStyle',':')
yt1 = yticks;
yl1 = ylim;


subplot(1,3,1)
plot(m,UPmse1,'*-')
hold on;
plot(m,UPmse2,'o-')
plot(m,UPmse3,'+-')
plot(m,UPmse4,'d-')
ylabel('Maximal IMSE','Interpreter','latex')
xlabel('m','Interpreter','latex')
title(strcat(funname,': Uniform Sampling'),'Interpreter','latex')
legend('Exp','Sq-Exp','Matern3/2','Matern5/2');
set(gca,'XTick', m, 'Xscale','log', 'Yscale','log', 'XMinorTick','off', 'YMinorTick','off')
set(gca,'XMinorGrid','off', 'YMinorGrid','off', 'XGrid','on', 'YGrid','on','GridLineStyle',':')
yt2 = yticks;
yl2 = ylim;


subplot(1,3,2)
plot(m,TTPmse1,'*-')
hold on;
plot(m,TTPmse2,'o-')
plot(m,TTPmse3,'+-')
plot(m,TTPmse4,'d-')
ylabel('Maximal IMSE','Interpreter','latex')
xlabel('m','Interpreter','latex')
title(strcat(funname,': Truncated Normal'),'Interpreter','latex')
legend('Exp','Sq-Exp','Matern3/2','Matern5/2');
set(gca,'XTick', m, 'Xscale','log', 'Yscale','log', 'XMinorTick','off', 'YMinorTick','off')
set(gca,'XMinorGrid','off', 'YMinorGrid','off', 'XGrid','on', 'YGrid','on','GridLineStyle',':')
yt3 = yticks;
yl3 = ylim;




figure
subplot(1,3,3)
plot(m,NPpfs1,'*-')
hold on;
plot(m,NPpfs2,'o-')
plot(m,NPpfs3,'+-')
plot(m,NPpfs4,'d-')
ylabel('IPFS','Interpreter','latex')
xlabel('m','Interpreter','latex')
title(strcat(funname,': Normal Sampling'),'Interpreter','latex')
legend('Exp','Sq-Exp','Matern3/2','Matern5/2');
set(gca,'XTick', m, 'xscale','log', 'yscale', 'log', 'XMinorTick','off', 'YMinorTick','off')
set(gca,'XMinorGrid','off', 'YMinorGrid','off', 'XGrid','on', 'YGrid','on','GridLineStyle',':')
yt1 = yticks;
yl1 = ylim;


subplot(1,3,1)
plot(m,UPpfs1,'*-')
hold on;
plot(m,UPpfs2,'o-')
plot(m,UPpfs3,'+-')
plot(m,UPpfs4,'d-')
ylabel('IPFS','Interpreter','latex')
xlabel('m','Interpreter','latex')
title(strcat(funname,': Uniform Sampling'),'Interpreter','latex')
legend('Exp','Sq-Exp','Matern3/2','Matern5/2');
set(gca,'Ylim',yl1,'YTick',yt1)
set(gca,'XTick', m, 'xscale','log', 'yscale', 'log', 'XMinorTick','off', 'YMinorTick','off')
set(gca,'XMinorGrid','off', 'YMinorGrid','off', 'XGrid','on', 'YGrid','on','GridLineStyle',':')


subplot(1,3,2)
plot(m,TTPpfs1,'*-')
hold on;
plot(m,TTPpfs2,'o-')
plot(m,TTPpfs3,'+-')
plot(m,TTPpfs4,'d-')
ylabel('IPFS','Interpreter','latex')
xlabel('m','Interpreter','latex')
title(strcat(funname,': Truncated Normal'),'Interpreter','latex')
legend('Exp','Sq-Exp','Matern3/2','Matern5/2');
set(gca,'Ylim',yl1,'YTick',yt1)
set(gca,'XTick', m, 'xscale','log', 'yscale', 'log', 'XMinorTick','off', 'YMinorTick','off')
set(gca,'XMinorGrid','off', 'YMinorGrid','off', 'XGrid','on', 'YGrid','on','GridLineStyle',':')






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
subplot(1,2,1)
plot(m,UUPmse1,'*-')
hold on;
plot(m,UUPmse2,'o-')
plot(m,UUPmse3,'+-')
plot(m,UUPmse4,'d-')
ylabel('Maximal IMSE','Interpreter','latex')
xlabel('m','Interpreter','latex')
title('Griewank (1D): Static Sampling','Interpreter','latex')
legend('Exp','Sq-Exp','Matern3/2','Matern5/2');
set(gca,'Ylim',[0.027,1])
set(gca,'XTick', m, 'Xscale','log', 'Yscale','log', 'XMinorTick','off', 'YMinorTick','off')
set(gca,'XMinorGrid','off', 'YMinorGrid','off', 'XGrid','on', 'YGrid','on','GridLineStyle',':')
yt = yticks;
yl = ylim;



subplot(1,2,2)
plot(m,EEPmse1,'*-')
hold on;
plot(m,EEPmse2,'o-')
plot(m,EEPmse3,'+-')
plot(m,EEPmse4,'d-')
ylabel('Maximal IMSE','Interpreter','latex')
xlabel('m','Interpreter','latex')
title('Griewank (1D): Adaptive MSE','Interpreter','latex')
legend('Exp','Sq-Exp','Matern3/2','Matern5/2');
set(gca,'Ylim',yl,'YTick',yt)
set(gca,'XTick', m, 'Xscale','log', 'Yscale','log', 'XMinorTick','off', 'YMinorTick','off')
set(gca,'XMinorGrid','off', 'YMinorGrid','off', 'XGrid','on', 'YGrid','on','GridLineStyle',':')





figure
subplot(1,2,1)
plot(m,UUPpfs1,'*-')
hold on;
plot(m,UUPpfs2,'o-')
plot(m,UUPpfs3,'+-')
plot(m,UUPpfs4,'d-')
ylabel('IPFS','Interpreter','latex')
xlabel('m','Interpreter','latex')
title('Griewank (1D): Static Sampling','Interpreter','latex')
legend('Exp','Sq-Exp','Matern3/2','Matern5/2');
set(gca,'Ylim',[0.09,0.71])
set(gca,'XTick', m, 'xscale','log', 'yscale', 'log', 'XMinorTick','off', 'YMinorTick','off')
set(gca,'XMinorGrid','off', 'YMinorGrid','off', 'XGrid','on', 'YGrid','on','GridLineStyle',':')
yt = yticks;
yl = ylim;


subplot(1,2,2)
plot(m,EEPpfs1,'*-')
hold on;
plot(m,EEPpfs2,'o-')
plot(m,EEPpfs3,'+-')
plot(m,EEPpfs4,'d-')
ylabel('IPFS','Interpreter','latex')
xlabel('m','Interpreter','latex')
title('Griewank (1D): Adaptive MSE','Interpreter','latex')
legend('Exp','Sq-Exp','Matern3/2','Matern5/2');
set(gca,'XTick', m, 'xscale','log', 'yscale', 'log', 'XMinorTick','off', 'YMinorTick','off')
set(gca,'XMinorGrid','off', 'YMinorGrid','off', 'XGrid','on', 'YGrid','on','GridLineStyle',':')
set(gca,'Ylim',yl,'YTick',yt)



