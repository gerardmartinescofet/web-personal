%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% AVERAGE TARGETING IN A NEW KEYNESIAN MODEL %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc;

% PARAMETERS
p.sigma   = 1;
p.rho     = (1.005)^0.25-1;
p.beta    = 1/(1+p.rho);
p.kappa   = 0.3;
p.avpi    = (1.02)^0.25-1;
p.T       = 16+1;
p.tshock  = 1;
p.shock   = -.01;
p.rhou    = sqrt(2)/2;
p.phipi_s = 1.5;
p.phipi_a = 1.5;
p.N       = 4;

% STEADY STATE
p.piSS = p.avpi;
p.iSS  = p.piSS+p.rho;
p.xSS  = ((1-p.beta)/p.kappa)*p.avpi;


% SOLUTION
A       = [1 0 1/p.sigma; -p.kappa 1 0; 0 -p.phipi_s 1];
B       = [1 1/p.sigma 0; 0 p.beta 0; 0 0 0];
C       = [ones(3, p.T).*[p.rho/p.sigma; 0; p.rho+(1-p.phipi_s)*p.avpi]];
C(2, 1) = p.shock;
for i=2:1:p.T
    C(2, i) = p.rhou*C(2, i-1);
end
z = [zeros(3, p.T),[p.xSS; p.piSS; p.iSS]];

for t=0:1:p.T-1
    
    z(:, p.T-t) = A\(B*z(:, p.T+1-t)+C(:, p.T-t));
    
end

% PLOT
% t = [0:1:p.T-1];
% figure(1)
% 
% subplot(3,1,1)
% plot(t, z(1,1:p.T))
% hold on
% drawLine([0,p.xSS],[p.T-1,p.xSS])
% xlabel('t', 'interpreter', 'tex')
% ylabel('x_t', 'interpreter', 'tex')
% 
% subplot(3,1,2)
% plot(t, z(2,1:p.T))
% hold on
% drawLine([0,p.piSS],[p.T-1,p.piSS])
% xlabel('t', 'interpreter', 'tex')
% ylabel('\pi_t', 'interpreter', 'tex')
% 
% subplot(3,1,3)
% plot(t, z(3,1:p.T))
% hold on
% drawLine([0,p.iSS],[p.T-1,p.iSS])
% xlabel('t', 'interpreter', 'tex')
% ylabel('i_t', 'interpreter', 'tex')

% AVERAGE INFLATION

pi_eq      = fsolve(@(x) inflation(x,p)-x, p.piSS*ones(1,p.T));
average_pi = pi_average(pi_eq, p.N, p);

% SOLUTION
A              = [1 0 1/p.sigma; -p.kappa 1 0; 0 -p.phipi_a/(p.N+1) 1];
B              = [1 1/p.sigma 0; 0 p.beta 0; 0 0 0];
C              = [ones(1, p.T).*(p.rho/p.sigma); zeros(1, p.T); zeros(1, p.T)];
C(2, p.tshock) = p.shock;
for t=p.tshock+1:1:p.T
    C(2, t) = p.rhou*C(2, t-1);
end
pi_av = pi_average(pi_eq, p.N, p);
for t=1:1:p.T
    C(3, t) = p.rho+p.phipi_a*(pi_av(t)-(1/(p.N+1))*pi_eq(t))+(1-p.phipi_a)*p.avpi;
end


z1 = [zeros(3, p.T),[p.xSS; p.piSS; p.iSS]];
for t=0:1:p.T-1
    
    z1(:, p.T-t) = A\(B*z1(:, p.T+1-t)+C(:, p.T-t));
    
end

% PLOT
t              = [0:1:p.T-1];
x              = 100*(z(1, 1:p.T));
x1             = 100*(z1(1, 1:p.T));
pi             = 100*((1+z(2, 1:p.T)).^4-1);
pi1            = 100*((1+z1(2, 1:p.T)).^4-1);
average_pi_ann = 100*((1+average_pi).^4-1);
i              = 100*((1+z(3, 1:p.T)).^4-1);
i1             = 100*((1+z1(3, 1:p.T)).^4-1);

figure(2)
subplot(3,1,1)
plot(t, x, 'b-o')
hold on
plot(t, x1, 'R-x')
hold on
drawLine([0,100*p.xSS],[p.T-1,100*p.xSS])
xlabel('t', 'interpreter', 'tex')
xticks([1:1:16])
ylabel('x_t', 'interpreter', 'tex')
legend('x_t^s','x_t^a', 'interpreter', 'tex')

subplot(3,1,2)
plot(t, pi, 'b-o')
hold on
plot(t, pi1, 'R-x')
hold on
plot(t, average_pi_ann, 'g-.*')
hold on
drawLine([0,100*((1+p.piSS)^4-1)],[p.T-1,100*((1+p.piSS)^4-1)])
xlabel('t', 'interpreter', 'tex')
xticks([1:1:16])
ylabel('\pi_t', 'interpreter', 'tex')
legend('\pi_t^s','\pi_t^a','\pi_t^{av}', 'interpreter', 'tex')

subplot(3,1,3)
plot(t, i, 'b-o')
hold on
plot(t, i1, 'r-x')
hold on
drawLine([0,100*((1+p.iSS)^4-1)],[p.T-1,100*((1+p.iSS) ^4-1)])
xlabel('t', 'interpreter', 'tex')
xticks([1:1:16])
ylabel('i_t', 'interpreter', 'tex')
legend('x_t^s','x_t^a', 'interpreter', 'tex')

