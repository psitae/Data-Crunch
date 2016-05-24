m = 4; %kg
k = 6; %N/m
c = 2*sqrt(24); %N*s/m

dt = .001; %s

v0 = -4; %m/s right is positive
x0 = 2; %m

history = zeros(1,10000);

x = x0;  
v = v0; 
j = 1;

while abs(v) > 1e-2 % this is a numerical differentiation
    history(j) = x;
    a = -x*k/m - c*v/m;
    v = a*dt + v;
    x = v*dt + x;
    j = j + 1;
end

history(find(history==0)) = NaN;
plot( history(:) )
title('６ N/m spring with overshoot')
xlabel('Time (ms)')
ylabel('Position (m)')
grid on
        
fprintf('The overshoot is %1.5f meters\n',x)
