% rabbit population growth
% bound by carrying capacity
nMax = 1000; % length of sim
r = sqrt(8); % max growth rate
K = 0.6; % carrying capacity

x = zeros(nMax,1); % thousands of rabbits
x(1,:) = 0.2; % initial population

xdot = @(x,r,K) x + r*(1 - x./K).*x; 

for n=2:nMax
    
    x(n,:) = xdot(x(n-1,:),r,K);
    
end % finished loop

% THE MODEL ^
% ------------------------------------------
% THE BEHAVIOR / THE OUTPUT ? 

figure(1); 
plot(x(nMax-100:nMax),'-ok');
ylabel('pop. size (thousands)')
xlabel('months')