close all
clear
clc

r0 = [6713e3 0
      7343e3 0
     -7903e3 0];  % initial position
 
v0 = [0  8011
      0  8614 
      0 -6970];   % initial velocity
u = 3.986e14;
dt = 1;

% The end goal of the next section is to find E and H
% so that the nature of the orbits can be determined

num = 2;
nu = atan(r0(num,2)/r0(num,1))
flight_angle = arctan( v0(num,1) , v0(num,2) )
elevation_angle = (pi/2 + nu - flight_angle)
E = mag(v0(num,:))^2/2 - u/mag(r0(num,:))
H = mag(r0(num,:))*mag(v0(num,:))* cos (elevation_angle)

% this is a basic numerical differentiation
history = zeros(2,20e4,3);

for j = 1:3
    r = r0(j,:);
    v = v0(j,:);
    for i = 1:15000
        history(1,i,j) = r(1);
        history(2,i,j) = r(2);
        acc = ( -u/mag(r)^3 ) * r;
        v = acc*dt + v;
        r = v*dt + r;
        if j == 3 && abs( v(1) ) < 1
            disp('break!')
            v
            arctan( r(1), r(2) )
            break
        end
    end
end

% display the paths saved in history using plot

theta = linspace(0,2*pi,1000);
earthx = 6.378e6 * cos(theta);
earthy = 6.378e6 * sin(theta);

history(find(history==0))= NaN; %delete remaining zeros

f = figure(1);
plot(history(1,:,1) , history(2,:,1), ...
     history(1,:,2) , history(2,:,2), 'red', ...
     history(1,:,3) , history(2,:,3), 'magenta', ...
     earthx, earthy, 'green',...
     'MarkerSize', 10)

range = 20e6;
axis([-range range -range range])
f.Position = ([500 300 700 620]);
legend 'Initial Orbit' 'Final Orbit' 'Transfer' 'Earth'
grid on

