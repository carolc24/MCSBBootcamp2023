clear all;
close all;

% Modeling phosphorylation-dephosphorylation

% kinase catalyzes I + K <-> IK -> A + K
% phosphatase catalyzes A + P <-> AP -> I + P

% total of both enzymes is constant
P_tot = 1; % uM
K_tot = 1; % uM
% initial conditions
I_int = 1; % uM
A_int = 0; % uM
IK_int = 0; % uM
AP_int = 0; % uM

% binding rates
kon_I = 10; % 1 / (uM * sec)
kon_A = 10; % 1 / (uM * sec)
% unbinding rates
koff_I = 10; % 1 / sec
koff_A = 10; % 1 / sec
% catalysis rates
kcat_I = 10; % 1 / sec
kcat_A = 100; % 1 / sec

% define ODEs
dxdt = @(t,y,K0,P0) [ -kon_A * y(1) * (K0 - y(3)) + koff_I * y(3) + kcat_I * y(4) ; % dI/dt
                      -kon_I * y(2) * (P0 - y(4)) + koff_A * y(4) + kcat_A * y(3) ; % dA/dt
                       kon_A * y(1) * (K0 - y(3)) - koff_I * y(3) - kcat_A * y(3) ; % dIK/dt
                       kon_I * y(2) * (P0 - y(4)) - koff_A * y(4) - kcat_I * y(4)]; % dAP/dt

[T,X] = ode45(@(t,y) dxdt(t,y,K_tot,P_tot) ,[0 5],[I_int;A_int;IK_int;AP_int]);

sum(X(length(X),:)) % this should be 1

figure(1);
plot(T,X);
xlabel("time");
ylabel("concs");
legend("I","A","IK","AP");

figure(2);
plot(T,(X(:,2))./(X(:,1) + X(:,2))); % K/P ratio
xlabel("time");
ylabel("[A]/([I]+[A])");

%%

% parameter sweep on K_tot

Af_1 = zeros(21,1); % when I_tot = 1 uM
Af_100 = zeros(21,1); % when I_tot = 100 uM

for i=1:21
   K_tot = 10^(i*0.25-3.25); % vary K_tot over range 0.001 - 100
   [T1,X1] = ode45(@(t,y) dxdt(t,y,K_tot,P_tot),[0 5],[1,0,0,0]); % low I simulation
   [T2,X2] = ode45(@(t,y) dxdt(t,y,K_tot,P_tot),[0 5],[100,0,0,0]); % high I simulation
   Af_1(i) = X1(length(T1),2) / sum(X1(length(T1),:)); % normalized by total protein count
   Af_100(i) = X2(length(T2),2) / sum(X2(length(T2),:)); 
end

logKtot = -3:0.25:2;
plot(logKtot, Af_1, logKtot, Af_100);
xlabel("log_{10} K_{tot}")
ylabel("steady state [A]/total protein")
legend("I_{tot} = 1","I_{tot} = 100")