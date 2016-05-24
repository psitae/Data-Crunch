close all
clear
clc

mu = 3.986e5;
re = 6397; 
h1 = 4669;
r1 = h1 + re;
h2 = 6368;
r2 = h2 + re;
theta = 86*pi/180; 
d = sqrt( r1^2 + r2^2 - 2*r1*r2*cos(theta) );
true_time = 48.5*60;
a = (r1 + r2)/2;

TOF = zeros(2, 25);
for i = 1:25
	alpha = acos( 1 - (r1 + r2 + d)/(2*a) );
	beta = acos( 1 - (r1 + r2 - d)/(2*a) );
	TOF(1,i) = a;
	TOF(2,i) = sqrt(a^3/mu)*( (alpha - sin(alpha)) - (beta - sin(beta)) ); 
	a = a * 1.045;
end

plot( TOF(1,:) , TOF(2,:) )
title 'Time of flight vs. Semi-major Axis'
xlabel 'Semi-major Axis (km)'
ylabel 'Time of flight (s)'
grid on
hold on

TOF_exact= 1e6;
a = (r1 + r2)/2;

while abs( TOF_exact - true_time ) > 1
	alpha = acos( 1 - (r1 + r2 + d)/(2*a) );
	beta = acos( 1 - (r1 + r2 - d)/(2*a) );
	TOF_exact= sqrt(a^3/mu)*( (alpha - sin(alpha)) - (beta - sin(beta)) ); 
        a = a * 1.004; % changed from .45% increments
end


plot( a, TOF_exact, 'r*')
legend('Approximate', sprintf('a = %4.0f km, accurate to one second', a))

Energy = -mu/(2*a);
phi = alpha - beta
ua = atan( sin(phi)^-1*(cos(phi) - (a-r2)/(a-r1)) )
ecc = (a-r1)/(a*cos(ua));
p = a*(1 - ecc^2);
H = sqrt(p*mu);
%H = sqrt(mu^2*(ecc^2-1)/(2*Energy))
Vp = H/(a*(1-ecc))
Va = H/(a*(1+ecc))

% check work using mean anomaly

ub = cos(ua)*(a-r2)/(a-r1);
TOF2 = sqrt(a^3/mu)*( (ua - ecc*sin(ua)) - (ub - ecc*sin(ub)) )
errorp = (TOF_exact - TOF2)/TOF_exact;
fprintf('First flight time: %4.2f\n', TOF_exact)
fprintf('Error is %4.2f%%\n', 100*errorp)

nu1 = acos( ecc^-1*(p/r1 - 1) );
nu2 = acos( ecc^-1*(p/r2 - 1) );
rp = a*(1 - ecc);

ub = acos( (ecc+cos(nu2))/(1 + ecc*cos(nu2)) );
Ma = ua - ecc*sin(ua);
Mb = ub - ecc*sin(ub);
TOF3 = sqrt(a^3/mu)*(Mb-Ma)
t = linspace(0, 2*pi, 1000);
earthx = re*cos(t);
earthy = re*sin(t);
c = ecc*a;
pathx = a*cos(t) - c;
b = sqrt( a^2 - c^2 );
pathy = b*sin(t);

point1x = r1 * cos(nu1);
point1y = r1 * sin(nu1);
point2x = r2 * cos(nu2);
point2y = r2 * sin(nu2);

figure (2)
plot(earthx, earthy, 'green', pathx, pathy, ':', ...
    point1x, point1y, 'r*', ...
    point2x, point2y, 'r*')
grid on
title 'Two Observations'
xlabel 'Distance (km)'
legend 'Earth' 'Calculated Orbit Path' 'Observation 1' ...
    'Observation 2' 'Location' 'Southwest'
ylabel 'Distance (km)'

