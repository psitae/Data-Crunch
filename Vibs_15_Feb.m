m = 4; %kg
k = 6; %N/m
c = [1 5 10]; %N*s/m

time_interval = 10; %s
t = linspace(0,time_interval,1000); %s
dt = time_interval/1000; %s

v0 = 4; %m/s right is positive
x0 = 2; %m

history = zeros(3,1000);

for i = 1:3 
    x = x0;  
    v = v0;  
    for j = 1:1000   % this is a numerical differentiation
        history(i,j) = x;
        a = -x*k/m - c(i)*v/m;
        v = a*dt + v;
        x = v*dt + x;
    end
end

plot( t, history(1,:), t, history(2,:), t, history(3,:), 'green' )
title('６ N/m spring with varing damping constants')
xlabel('Time (s)')
ylabel('Position (m)')
legend('C = 1 N*s/m','C = 5 N*s/m','C = 10 N*s/m')
grid on
        
        
