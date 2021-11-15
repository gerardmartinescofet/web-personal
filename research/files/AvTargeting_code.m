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
p.shock   = -0.01;
p.rhou    = sqrt(2)/2; %0
p.phipi_s = 1.5;
p.phipi_a = 1.5;
p.N       = 3;

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

% AVERAGE INFLATION

pi_eq      = fsolve(@(x) inflation(x, p.phipi_a, p.N, p.shock, p)-x, p.piSS*ones(1,p.T));
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
legend('$x_t^s$','$x_t^a$', 'interpreter', 'latex', 'FontSize',10)

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
legend('$\pi_t^s$','$\pi_t^a$','$\tilde{\pi_t}$', 'interpreter', 'latex','FontSize',10, 'location', 'east')

subplot(3,1,3)
plot(t, i, 'b-o')
hold on
plot(t, i1, 'r-x')
hold on
drawLine([0,100*((1+p.iSS)^4-1)],[p.T-1,100*((1+p.iSS) ^4-1)])
xlabel('t', 'interpreter', 'tex')
xticks([1:1:16])
ylabel('i_t', 'interpreter', 'tex')
legend('$i_t^s$','$i_t^a$', 'interpreter', 'latex', 'FontSize', 10, 'location', 'east')

%% WELFARE ANALYSIS

weight   = 0.1;
var_x_s  = @(phipi_s)(((p.rhou-phipi_s)/(p.sigma*(1-p.rhou*p.beta)*(1-p.rhou)-(p.rhou-phipi_s)*p.kappa))^2)*(0.01)^2;
var_pi_s = @(phipi_s)((p.sigma*(1-p.rhou)/(p.sigma*(1-p.rhou*p.beta)*(1-p.rhou)-(p.rhou-phipi_s)*p.kappa))^2)*(0.01)^2;
loss_s   = @(phipi_s) weight*var_x_s(phipi_s) + var_pi_s(phipi_s);

phipi = [3.5:0.01:4]; % Change upper and lower bound and determine the minimum visually. This is currently the fastest method, rather than using a global solver.
for i = 1:1:length(phipi)
    l(i) = loss_s(phipi(i));
    r(i) = loss_a([phipi(i), 3], weight, p);
end

figure(1)
plot(phipi, l)
figure(2)
plot(phipi, r)

%% FUNCTIONS (copy-paste each function in separate .m files to run main code)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function inflation = f(guess, phipi_a, N, shock, p) % Computes the equilibrium path of inflation given the guessed equilibrium inflation process by finding a fixed point (see paper).

A              = [1 0 1/p.sigma; -p.kappa 1 0; 0 -phipi_a/(N+1) 1];
B              = [1 1/p.sigma 0; 0 p.beta 0; 0 0 0];
C              = [ones(1, p.T).*(p.rho/p.sigma); zeros(1, p.T); zeros(1, p.T)];
C(2, p.tshock) = shock;
for t=p.tshock+1:1:p.T
    C(2, t) = p.rhou*C(2, t-1);
end
pi_av = pi_average(guess, N, p);
for t=1:1:p.T
    C(3, t) = p.rho+phipi_a*(pi_av(t)-(1/(N+1))*guess(t))+(1-phipi_a)*p.avpi;
end


z = [zeros(3, p.T),[p.xSS; p.piSS; p.iSS]];
for t=0:1:p.T-1
    
    z(:, p.T-t) = A\(B*z(:, p.T+1-t)+C(:, p.T-t));
    
end

inflation = z(2, (1:p.T));

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function takes as input a hypothetical path of inflation "pi", and
% yields as output the path of average inflation in the last N+1 periods 
% (N past periods + today).

function pi_average = f(pi, N, p)

for t=1:1:length(pi)
    if N>=t
        pi_average(t) = (1/(N+1))*sum(pi(1:t))+((N-(t-1))/(N+1))*p.piSS;
    else
        pi_average(t) = (1/(N+1))*sum(pi((t-N):t));
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function takes as input parameters m (vector of phi_a and N), and computes as
% output the variance of inflation in t=0 by performing montecarlo
% simulations in the case of a one-period shock (at t=0) to epsilon_t
% (innovation in the AR(1) for the cost-push shock).


function montecarlo = f(m, p)
for m2=1:1:1000
for m1=1:1:10

shock(m1) = 0.01*randn(1);

% AVERAGE INFLATION SOLUTION

pi_eq      = fsolve(@(x) inflation(x, m(1), round(m(2)), shock(m1), p)-x, p.piSS*ones(1,p.T));
average_pi = pi_average(pi_eq, round(m(2)), p);

% SOLUTION
A              = [1 0 1/p.sigma; -p.kappa 1 0; 0 -m(1)/(round(m(2))+1) 1];
B              = [1 1/p.sigma 0; 0 p.beta 0; 0 0 0];
C              = [ones(1, p.T).*(p.rho/p.sigma); zeros(1, p.T); zeros(1, p.T)];
C(2, p.tshock) = shock(m1);
for t=p.tshock+1:1:p.T
    C(2, t) = p.rhou*C(2, t-1);
end
pi_av = pi_average(pi_eq, round(m(2)), p);
for t=1:1:p.T
    C(3, t) = p.rho+m(1)*(pi_av(t)-(1/(round(m(2))+1))*pi_eq(t))+(1-m(1))*p.avpi;
end


z1 = [zeros(3, p.T),[p.xSS; p.piSS; p.iSS]];
for t=0:1:p.T-1
    
    z1(:, p.T-t) = A\(B*z1(:, p.T+1-t)+C(:, p.T-t));
    
end

x_a(m1)  = z1(1,1);
pi_a(m1) = z1(2,1);
end

var_x_a(m2)  = var(x_a);
var_pi_a(m2) = var(pi_a);
end

montecarlo(1) = mean(var_x_a);
montecarlo(2) = mean(var_pi_a);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function takes as inut the variances of the welfare-relevant output
% gap and of inflation (summarized in vector x - x(1) for the gap, x(2) for
% inflation), and the weight/penalty that the policymaker gives to each
% variance. The output is the welfare loss L.

function loss_a = f(m, weight, p)

R = montecarlo(m, p);
var_x_a  = R(1);
var_pi_a = R(2);

loss_a = weight*var_x_a + var_pi_a;
 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Necessary to draw plots
function [] = drawLine(p1, p2)
theta = atan2( p2(2) - p1(2), p2(1) - p1(1));
r = sqrt( (p2(1) - p1(1))^2 + (p2(2) - p1(2))^2);
line = 0:0.01: r;
x = p1(1) + line*cos(theta);
y = p1(2) + line*sin(theta);
plot(x, y,'LineStyle','--','Color','k')