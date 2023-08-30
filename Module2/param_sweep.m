% parameter sweep of r in logistic equation
xdot = @(x,r,K) x + r*(1 - x/K)*x; 

nMax=1000;
K = 0.6;
x = zeros(1,nMax);
x(1)=0.2;

% Run lots of simulations for different values of r
figure(1)
xlabel('r')
ylabel('steady state x')
for r=2:0.01:3
    % run a simulation
    for n=2:nMax
        x(n) = xdot(x(n-1),r,K);
    end
    plot(r,x(nMax-64:nMax),'o')
    hold on
end