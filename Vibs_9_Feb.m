clc
close all
clear 

init = linspace(.001,.9*pi,100);
period = zeros(1,100);
dt = .001;


%assume L = 1

for i = 1:100
    theta = init(i);
    deriv = 0;
    deriv2 = 0;
    counter = 0;
    while counter < 10 || deriv < 1e-6
        deriv2 = -4.9*sin(theta);
        deriv = deriv + dt*deriv2;
        theta = theta + dt*deriv; 
        counter = counter + 1;
    end
    period(i) = counter * 1000 * 2;
    
end

figure(2)
plot(period)
  
