% Curve Fitting - Modified Logistic Growth
% UCI COSMOS 2022:  Tissue and Tumor Modeling (Cluster 3)

% clear all

% For the user to run this curve-fitting routine, the workspace must
% contain the vectors Time and Population. These vectors can be either data
% loaded from GroupX_ModifiedLogisticData.mat from week 1, or experiemental
% data manually entered with the appropriate names.

Time = [0,30,65,95,125,155]; % minutes
Population = [-.021,.009,-.012,-.014,.031,.003; % no bacteria 0 glucose
              -.004,.005,-.009,-.003,.008,.004; % no bacteria 250uM glucose
              .004,.013,-.026,-.004,.023,-.002; % no bacteria 500 uM glucose
              .035,.074,.099,.208,.407,.531; % bacteria 0 glucose
              .040,.073,.102,.191,.365,.571; % bacteria 250uM glucose
              .044,.054,.109,.173,.371,.549]; % bacteria 500uM glucose

figure
plot(Time,Population');
xlabel("time (min)")
ylabel("OD (600 nm)")
legend("B-0","B-250","B-500","B+0","B+250","B+500")

%% User-defined inputs

% Initial guess in the order of:
% [lambdam; theta; InitialPopulation; alpha]
% Note: Order here is important.
initialguess = [0.5; 2; 0.01; 0.2];

%% Begin Curve Fit

FinalTime = Time(end);

% Transpose population vector and time vectors if they are rows
if size(Population,1) == 1
    Population = Population';
end
if size(Time,1) == 1
    Time = Time';
end

% Find parameter values minimizing error.
min_0 = fminsearch(@(x) CF_Error_ModLog2022(x,Time,Population(4,:)'),initialguess);

lambdam = min_0(1);            % Growth rate
theta = min_0(2);              % Carrying Capacity
InitialPopulation = min_0(3);  % Initial Population
alpha = min_0(4);              % Exponent

[h,N,times] = CF_Sim_ModLog2022(lambdam,theta,InitialPopulation,alpha,Time,true);

%% Plotting

figure,plot(Time,Population(4,:)','ro',times,N,'r')
err_0 = sum((Population(4,:)' - h).^2);
text(Time(4)+5,Population(4,4),['SSE = ', num2str(err_0)],'Color','r');
xlabel('Time','FontWeight','bold')
ylabel('Population','FontWeight','bold')
grid on
hold on

% 250 uM glucose
min_250 = fminsearch(@(x) CF_Error_ModLog2022(x,Time,Population(5,:)'),initialguess);
[h,N,times] = CF_Sim_ModLog2022(min_250(1), min_250(2), min_250(3), min_250(4),Time,true);
err_250 = sum((Population(5,:)' - h).^2);
plot(Time,Population(5,:)','bx',times,N,'b')
text(Time(4)+5,Population(5,4),['SSE = ', num2str(err_250)],'Color','b');

% 500 uM glucose
min_500 = fminsearch(@(x) CF_Error_ModLog2022(x,Time,Population(6,:)'),initialguess);
[h,N,times] = CF_Sim_ModLog2022(min_500(1), min_500(2), min_500(3), min_500(4),Time,true);
err_500 = sum((Population(6,:)' - h).^2);
plot(Time,Population(6,:)','m+',times,N,'m')
text(Time(4)+5,Population(6,4),['SSE = ', num2str(err_500)],'Color','m');

title('Modified Logistic Growth Curve Fit','FontSize',12,'FontWeight','bold')
set(gca,'FontWeight','bold')
legend('Data (no glucose)','Curve fit (no glucose)', 'Data (250uM glucose)','Curve fit (250uM glucose)', ...
    'Data (500uM glucose)', 'Curve fit (500uM glucose)')

