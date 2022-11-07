
clearvars,clc;

% De Jong's function: testfun
fun = @testfun1;
desig = 1:10; desig = desig'; 
nd = length(desig);  %location of design point
n0 = 10;
va = 2; sd = sqrt(va);
span1 = 9; span2 = 1;    % test interval [0,4]
scen = 5.5; sva = 3; ssd = sqrt(sva); ssd2 = ssd;    % x^* distribution
Lm = 100;    % number of replications for random covariates point

m = [5,8,12,18,28,42,65,100]; nm = length(m);    % number of randomly selected covariate
%%% PFS
testa = 1000;    % number of test covariate point for PFS
Del = 0.05;
%%%%%%%
options = optimoptions('fmincon','OptimalityTolerance',0.01,...
    'StepTolerance',0.01,'Display','off');
optseed = 10;  
d = 1;  % context dim.
%%%%%%%%%%%%%%%%%%%% Sample distribution: Uniform %%%%%%%%%%%%%%%%%%%%%%
[ UPmse1,UPpfs1,UPmse2,UPpfs2,UPmse3,UPpfs3,UPmse4,UPpfs4 ] ...
    = UniS( fun,m,desig,span1,span2,Lm,testa,Del,sd,n0 );

%%%%%%%%%%%%%%%%%%% Sample distribution: Truncated Normal %%%%%%%%%%%%%%%%%%%%%%
Tmu0 = 5.5; Trusd0 = 7; 
[ TTPmse1,TTPpfs1,TTPmse2,TTPpfs2,TTPmse3,TTPpfs3,TTPmse4,TTPpfs4 ] ...
    = TrnS( fun,m,desig,span1,span2,Lm,testa,Del,sd,n0,Tmu0,Trusd0 );

%%%%%%%%%%%%%%%%%%%% Sample distribution: Norm %%%%%%%%%%%%%%%%%%%%%%
[ NPmse1,NPpfs1,NPmse2,NPpfs2,NPmse3,NPpfs3,NPmse4,NPpfs4 ] ...
    = NorS( fun,m,desig,scen,ssd,Lm,testa,Del,sd,n0,ssd2 );

% %%%%%%%%%%%%%%%%%%%% Adaptive MSE %%%%%%%%%%%%%%%%%%%%%%
Tmu = 5.5; Trusd = 0.25;  
pd = makedist('Normal','mu',Tmu,'sigma',Trusd); 
Bd1 = span2; Bd2 = span2+span1;
Tpd = truncate(pd,Bd1,Bd2);
rng(500)
Etestx = random(Tpd,testa,d);   % test covariate generation for PFS
Etesty = fun(desig,Etestx);    % exact testfun value
Eargm = min(Etesty,[],2);    % true optimal design

[ TPmse1,TPpfs1,TPmse2,TPpfs2,TPmse3,TPpfs3,TPmse4,TPpfs4 ] ...
    = TrnS( fun,m,desig,span1,span2,Lm,testa,Del,sd,n0,Tmu,Trusd);

kernind = 1;
[ EPmse1,EPpfs1,Econshortrec1 ] ...
    = AdaS( fun,m,desig,span1,span2,Lm,testa,Del,sd,n0,kernind,Etestx,...
    Etesty,Eargm,options,optseed,d );
kernind = 2;
[ EPmse2,EPpfs2,Econshortrec2 ] ...
    = AdaS( fun,m,desig,span1,span2,Lm,testa,Del,sd,n0,kernind,Etestx,...
    Etesty,Eargm,options,optseed,d );
kernind = 3;
[ EPmse3,EPpfs3,Econshortrec3 ] ...
    = AdaS( fun,m,desig,span1,span2,Lm,testa,Del,sd,n0,kernind,Etestx,...
    Etesty,Eargm,options,optseed,d );
kernind = 4;
[ EPmse4,EPpfs4,Econshortrec4 ] ...
    = AdaS( fun,m,desig,span1,span2,Lm,testa,Del,sd,n0,kernind,Etestx,...
    Etesty,Eargm,options,optseed,d );

% %%%%%%%%%%%%%%%%%%%% Adaptive MSE va=1 %%%%%%%%%%%%%%%%%%%%%%
T2mu = 5.5; T2rusd = 1;  
p2d = makedist('Normal','mu',T2mu,'sigma',T2rusd); 
Bd1 = span2; Bd2 = span2+span1;
T2pd = truncate(p2d,Bd1,Bd2);
rng(500)
E2testx = random(T2pd,testa,d);   % test covariate generation for PFS
E2testy = fun(desig,E2testx);    % exact testfun value
E2argm = min(E2testy,[],2);    % true optimal design

[ T2Pmse1,T2Ppfs1,T2Pmse2,T2Ppfs2,T2Pmse3,T2Ppfs3,T2Pmse4,T2Ppfs4 ] ...
    = TrnS( fun,m,desig,span1,span2,Lm,testa,Del,sd,n0,T2mu,T2rusd);

kernind = 1;
[ E2Pmse1,E2Ppfs1,E2conshortrec1 ] ...
    = AdaS( fun,m,desig,span1,span2,Lm,testa,Del,sd,n0,kernind,E2testx,...
    E2testy,E2argm,options,optseed,d );
kernind = 2;
[ E2Pmse2,E2Ppfs2,E2conshortrec2 ] ...
    = AdaS( fun,m,desig,span1,span2,Lm,testa,Del,sd,n0,kernind,E2testx,...
    E2testy,E2argm,options,optseed,d );
kernind = 3;
[ E2Pmse3,E2Ppfs3,E2conshortrec3 ] ...
    = AdaS( fun,m,desig,span1,span2,Lm,testa,Del,sd,n0,kernind,E2testx,...
    E2testy,E2argm,options,optseed,d );
kernind = 4;
[ E2Pmse4,E2Ppfs4,E2conshortrec4 ] ...
    = AdaS( fun,m,desig,span1,span2,Lm,testa,Del,sd,n0,kernind,E2testx,...
    E2testy,E2argm,options,optseed,d );

save test1;

% %%%%%%%%%%%%%%%%%%%% Adaptive MSE var=0.25^2%%%%%%%%%%%%%%%%%%%%%%
figure

subplot(1,2,2)
plot(m,EPmse1,'*-')
hold on;
plot(m,EPmse2,'o-')
plot(m,EPmse3,'+-')
plot(m,EPmse4,'d-')
ylabel('Maximal IMSE','Interpreter','latex')
xlabel('m','Interpreter','latex')
title('De Jong (1D): Adaptive MSE','Interpreter','latex')
legend('Exp','Sq-Exp','Matern3/2','Matern5/2');
set(gca,'XTick', m, 'Xscale','log', 'Yscale','log', 'XMinorTick','off', 'YMinorTick','off')
set(gca,'XMinorGrid','off', 'YMinorGrid','off', 'XGrid','on', 'YGrid','on','GridLineStyle',':')
yt = yticks;
yl = ylim;

subplot(1,2,1)
plot(m,TPmse1,'*-')
hold on;
plot(m,TPmse2,'o-')
plot(m,TPmse3,'+-')
plot(m,TPmse4,'d-')
ylabel('Maximal IMSE','Interpreter','latex')
xlabel('m','Interpreter','latex')
title('De Jong (1D): Static Sampling','Interpreter','latex')
legend('Exp','Sq-Exp','Matern3/2','Matern5/2');
set(gca,'Ylim',yl,'YTick',yt)
set(gca,'XTick', m, 'Xscale','log', 'Yscale','log', 'XMinorTick','off', 'YMinorTick','off')
set(gca,'XMinorGrid','off', 'YMinorGrid','off', 'XGrid','on', 'YGrid','on','GridLineStyle',':')


figure
subplot(1,2,2)
plot(m,EPpfs1,'*-')
hold on;
plot(m,EPpfs2,'o-')
plot(m,EPpfs3,'+-')
plot(m,EPpfs4,'d-')
ylabel('IPFS','Interpreter','latex')
xlabel('m','Interpreter','latex')
title('De Jong (1D): Adaptive MSE','Interpreter','latex')
legend('Exp','Sq-Exp','Matern3/2','Matern5/2');
set(gca,'Ylim',[0.02,0.45])
set(gca,'XTick', m, 'xscale','log', 'yscale', 'log', 'XMinorTick','off', 'YMinorTick','off')
set(gca,'XMinorGrid','off', 'YMinorGrid','off', 'XGrid','on', 'YGrid','on','GridLineStyle',':')


subplot(1,2,1)
plot(m,TPpfs1,'*-')
hold on;
plot(m,TPpfs2,'o-')
plot(m,TPpfs3,'+-')
plot(m,TPpfs4,'d-')
ylabel('IPFS','Interpreter','latex')
xlabel('m','Interpreter','latex')
title('De Jong (1D): Static Sampling','Interpreter','latex')
legend('Exp','Sq-Exp','Matern3/2','Matern5/2');
set(gca,'Ylim',[0.02,0.45])
set(gca,'XTick', m, 'xscale','log', 'yscale', 'log', 'XMinorTick','off', 'YMinorTick','off')
set(gca,'XMinorGrid','off', 'YMinorGrid','off', 'XGrid','on', 'YGrid','on','GridLineStyle',':')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
funname = 'De Jong (1D)';
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
set(gca,'XTick', m, 'xscale','log', 'yscale', 'log', 'XMinorTick','off', 'YMinorTick','off')
set(gca,'XMinorGrid','off', 'YMinorGrid','off', 'XGrid','on', 'YGrid','on','GridLineStyle',':')
yt2 = yticks;
yl2 = ylim;


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
set(gca,'XTick', m, 'xscale','log', 'yscale', 'log', 'XMinorTick','off', 'YMinorTick','off')
set(gca,'XMinorGrid','off', 'YMinorGrid','off', 'XGrid','on', 'YGrid','on','GridLineStyle',':')
yt3 = yticks;
yl3 = ylim;


%%%%%%%%%%%%%%%%%%%%%%%% Adaptive MSE va=1 %%%%%%%%%%%%%%%%%%%%%%
figure

subplot(1,2,2)
plot(m,E2Pmse1,'*-')
hold on;
plot(m,E2Pmse2,'o-')
plot(m,E2Pmse3,'+-')
plot(m,E2Pmse4,'d-')
ylabel('Maximal IMSE','Interpreter','latex')
xlabel('m','Interpreter','latex')
title('De Jong (1D): Adaptive MSE','Interpreter','latex')
legend('Exp','Sq-Exp','Matern3/2','Matern5/2');
set(gca,'XTick', m, 'Xscale','log', 'Yscale','log', 'XMinorTick','off', 'YMinorTick','off')
set(gca,'XMinorGrid','off', 'YMinorGrid','off', 'XGrid','on', 'YGrid','on','GridLineStyle',':')
yt = yticks;
yl = ylim;

subplot(1,2,1)
plot(m,T2Pmse1,'*-')
hold on;
plot(m,T2Pmse2,'o-')
plot(m,T2Pmse3,'+-')
plot(m,T2Pmse4,'d-')
ylabel('Maximal IMSE','Interpreter','latex')
xlabel('m','Interpreter','latex')
title('De Jong (1D): Static Sampling','Interpreter','latex')
legend('Exp','Sq-Exp','Matern3/2','Matern5/2');
set(gca,'Ylim',yl,'YTick',yt)
set(gca,'XTick', m, 'Xscale','log', 'Yscale','log', 'XMinorTick','off', 'YMinorTick','off')
set(gca,'XMinorGrid','off', 'YMinorGrid','off', 'XGrid','on', 'YGrid','on','GridLineStyle',':')


figure
subplot(1,2,2)
plot(m,E2Ppfs1,'*-')
hold on;
plot(m,E2Ppfs2,'o-')
plot(m,E2Ppfs3,'+-')
plot(m,E2Ppfs4,'d-')
ylabel('IPFS','Interpreter','latex')
xlabel('m','Interpreter','latex')
title('De Jong (1D): Adaptive MSE','Interpreter','latex')
legend('Exp','Sq-Exp','Matern3/2','Matern5/2');
set(gca,'XTick', m, 'xscale','log', 'yscale', 'log', 'XMinorTick','off', 'YMinorTick','off')
set(gca,'XMinorGrid','off', 'YMinorGrid','off', 'XGrid','on', 'YGrid','on','GridLineStyle',':')
yt = yticks;
yl = ylim;


subplot(1,2,1)
plot(m,T2Ppfs1,'*-')
hold on;
plot(m,T2Ppfs2,'o-')
plot(m,T2Ppfs3,'+-')
plot(m,T2Ppfs4,'d-')
ylabel('IPFS','Interpreter','latex')
xlabel('m','Interpreter','latex')
title('De Jong (1D): Static Sampling','Interpreter','latex')
legend('Exp','Sq-Exp','Matern3/2','Matern5/2');
set(gca,'Ylim',yl,'YTick',yt)
set(gca,'XTick', m, 'xscale','log', 'yscale', 'log', 'XMinorTick','off', 'YMinorTick','off')
set(gca,'XMinorGrid','off', 'YMinorGrid','off', 'XGrid','on', 'YGrid','on','GridLineStyle',':')






