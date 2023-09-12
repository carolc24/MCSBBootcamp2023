% FitzHugh-Nagumo model

% Note: I did some of this last week...
% but I'm going through it again step-by-step today.

% v = x(1), w = x(2) 
dxdt = @(t,x,e,a,b) [x(1) - 1/3 * x(1).^3 - x(2);
                     e*(x(1) + a - b*x(2))     ];

tf = 500;

% oscillations

e = 0.08;
a = 0.5;
b = 0.2;

[T,X] = ode45(@(t,y) dxdt(t,y,e,a,b), [0,tf], [3.14,1.59]);

figure(1)
plot(T,X)
xlabel("time")
legend("v","w")
%%
% excitability
a = 1.0;

[T1,X1] = ode45(@(t,y) dxdt(t,y,e,a,b), [0,tf],[-1.5,-0.5]);
[T2,X2] = ode45(@(t,y) dxdt(t,y,e,a,b), [0,tf],[0.0,-0.5]);

figure(2)
plot(T1,X1,T2,X2)
xlabel("time")
legend("v_1","w_1","v_2","w_2")

x_ss = X1(length(X1),:);

%%
% perturbation
% different impulse durations

I0 = 1.0;
tStart = 40;
tEnd = 47;
Ifunc = @(t,t0,t1) I0*(t>t0)*(t<t1);

dxdt_p = @(t,x,e,a,b,t0,t1) [x(1) - 1/3 * x(1).^3 - x(2) + Ifunc(t,t0,t1);
                       e*(x(1) + a - b*x(2))                ];

[T1,X1] = ode45(@(t,y) dxdt_p(t,y,e,a,b,tStart,tEnd-5),[0,tf],x_ss);
[T2,X2] = ode45(@(t,y) dxdt_p(t,y,e,a,b,tStart,tEnd),[0,tf],x_ss);
[T3,X3] = ode45(@(t,y) dxdt_p(t,y,e,a,b,tStart,tEnd+200),[0,tf],x_ss);

figure(3);
plot(T1,X1,T2,X2,T3,X3);
xlabel("time")
legend("v (short)","w (short)", "v (med)", "w (med)","v (long)","w (long)")

%%
% neural network
% matrix of interactions
N = 10; % number of neurons
tf = 100;

I0 = zeros(N,1); % N elements
I0(4) = 1;
tStart = 40;
tEnd = tf;
Ifunc_v = @(t) I0*(t>tStart)*(t<tEnd);

dvdt_large = @(t,v,w,D) v - 1/3 * v.^3 - w + Ifunc_v(t) + D*(circshift(v,1) - 2*v + circshift(v,-1));
dwdt_large = @(t,v,w,e,a,b) e*(v + a - b*w);
dxdt_large = @(t,x,N,D,e,a,b) [dvdt_large(t,x(1:N),x(N+1:2*N),D); 
                               dwdt_large(t,x(1:N),x(N+1:2*N),e,a,b)];

D = 0.1; % connectivity
v_int = x_ss(1) * ones(1,N);
w_int = x_ss(2) * ones(1,N);

[T,X] = ode45(@(t,y) dxdt_large(t,y,N,D,e,a,b), [0,tf],[v_int';w_int']);

figure(4);
plot(T,X(:,1:N));
xlabel("time")
ylabel("membrane potential")

% movie
for nt=1:numel(T)
figure(5); clf; hold on; box on;
plot(X(nt,1:N));
set(gca,'ylim', [-2.5,2.5])
xlabel('Cell');
ylabel('Voltage')
drawnow;
end