function error = derivative_estimator(noise)

% ------------------------------------------------------
A     = 1.1; % fluorescence intensity units
omega = 2.6; % rad/s
A_0   = 0.01;

u=@(t) A*sin(omega*t)+A_0;

tArray = linspace(0,1.6,200);
uArray = u(tArray); % an array of samples of u
% ------------------------------------------------------

% analytical solutions (in real life, we might not know these)
dudtExact      =  A*omega*cos(omega*tArray);
du2dt2Exact    = -A*omega^2*sin(omega*tArray);
du3dt3Exact    = -A*omega^3*cos(omega*tArray);

% Take the sample and add a bit of noise
uObserved = uArray + (noise)*randn(size(tArray));

%%

% create finite-difference derivatives for first derivative, second derivative and third derivative
dudt   = diff(uObserved)./diff(tArray);
du2dt2 = diff(dudt)./diff(tArray(1:end-1));
du3dt3 = diff(du2dt2)./diff(tArray(1:end-2));

%%

% function to calculate accuracy

err = @(x_pred, x_actual) sum(abs(x_pred - x_actual(1:length(x_pred)))) / length(x_pred);

e1 = err(dudt,dudtExact);
e2 = err(du2dt2,du2dt2Exact);
e3 = err(du3dt3,du3dt3Exact);

error = [e1 e2 e3];