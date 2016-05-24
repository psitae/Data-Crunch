clc; close all

% raw data
init_wgt = [2001.8 1989.2 1976.2 1961.8 1944.2 1923.7 1906.6 1896.6 1892.4 1885.3 1879.5 1866.1];
final_wgt = [1994.8 1983.4 1967.6 1953.1 1934.9 1916.6 1901.5 1893.4 1889.2 1881.7 1875.8 1862.6];
oil_pressure = [ 1510 1520 1510 1490 1490 1320 890 490 405 320 240 180];
oil_flow = [2.9 3.1 3.35 3.75 4.2 3.3 3.4 2.8 3.2 3.75 4.1 4.6];
cyl_temp = [268 335 370 393 410 410 383 351 334 320 313 307];
oil_temp = [89 92 95 98 102 110 113 115 115 115 114 114];
RPM      = [2000:300:3200, 2600, 2600, 2000:300:3200];
%{
weights are in grams
pressure is measured in psi
flow rate is measure in gpm
temperatures are in farenheit
entries 1 - 5 are the first data set, then 6 and 7 are runs
two are three, 8 - 12 are the last set
%}

% make appropriate unit conversrions
% underscores are usefully read as "to"

RPS = RPM/60;
displacement = .206; % Liters
heating_value = 47.3e6; % J/kg
oil_density = 800; % kg/m^3

psi_pa = 6.894757e3; 
watt_hp = 1/745.7;
nm_ftlb = .7376;
g_lbf = 1/453.59;
lbf_lbm = 32.2;
difference = (init_wgt - final_wgt)/30; % g/s
pressure = psi_pa * oil_pressure; %Pa
mass_flow  = oil_flow / (60); % kg/s
vol_flow  = mass_flow / oil_density;


cyl_temp = (cyl_temp - 32)/1.8 + 273; % kelvin
oil_temp = (oil_temp - 32)/1.8 + 273;

% begin calculations 
eff = .8; % pump efficiency

brake_power_SI = pressure .* vol_flow / eff; % watts
brake_power = brake_power_SI * watt_hp; % horse power
brake_torque_SI = brake_power_SI ./ (2*pi*RPS); % N*m
brake_torque = brake_torque_SI * nm_ftlb; % ft * lbf
fuel_cons = difference * g_lbf * 3600; % lbf / hour
spc_fl_cns = fuel_cons * lbf_lbm ./ brake_power; % lbm / hour / horsepower
bmep_SI = brake_power_SI * 2 ./ (displacement * RPS); % Pa
bmep = bmep_SI * 1/psi_pa * 1000; % psi
mass_flow_SI = difference / 1000; % kg / s
conv_eff = brake_power_SI ./(mass_flow_SI * heating_value)*100; % percent

fst = RPM(1:5);
y = 1:5;
rd = RPM(8:12);
z = 8:12;

% this part creates plots of the calculated data

figure (3) % power
plot(fst, brake_power(y), 'redx-', ...
    2600, brake_power(6), 'kx', 2600, brake_power(7),'mx',...
    rd, brake_power(z), 'bluex-')
title({'Figure 3', 'Power vs. RPM with Different Amounts of Throttle'})
legend 'Full Throttle' '3/4 Throttle' 'Half Throttle' '1/4 Throttle' ...
    'Location' 'Northwest'
xlabel 'RPM'; ylabel 'Horsepower (550 ft*lbf/s)'
grid on

figure (4) % torque
plot(fst, brake_torque(y), 'redx-', ...
    2600, brake_torque(6), 'kx', 2600, brake_torque(7),'mx',...
    rd, brake_torque(z), 'bluex-')
title({'Figure 4', 'Torque vs. RPM with Different Amounts of Throttle'})
legend 'Full Throttle' '3/4 Throttle' 'Half Throttle' '1/4 Throttle' ...
    'Location' 'best'
axis( [ 2000 3200 0 3] )
xlabel 'RPM'; ylabel 'Torque (ft*lbf)'
grid on

figure (5) % fuel consumption
plot([2000,2600:300:3200], fuel_cons([1,3:5]), 'redx-', ...
    2600, fuel_cons(6), 'kx', 2600, fuel_cons(7),'mx',...
    rd, fuel_cons(z), 'bluex-', 2300, fuel_cons(2), 'redx')
title({'Figure 5', 'Fuel Consumption vs. RPM with Different Amounts of Throttle'})
legend 'Full Throttle' '3/4 Throttle' 'Half Throttle' ...
    '1/4 Throttle' 'Outlier' 'Location' 'Northwest'
xlabel 'RPM'; ylabel 'Fuel Consumption (lbf/hr)'
grid on

figure (6) % Breake Specific Fuel Consumption
plot([2000,2600:300:3200], spc_fl_cns([1,3:5]), 'redx-', ...
    2600, spc_fl_cns(6), 'kx', 2600, spc_fl_cns(7),'mx',...
    rd, spc_fl_cns(z), 'bluex-', 2300, spc_fl_cns(2), 'redx')
title({'Figure 6', 'Break Specific Fuel Consumption vs. RPM with Different Amounts of Throttle'})
legend 'Full Throttle' '3/4 Throttle' 'Half Throttle' '1/4 Throttle' ...
    'Outlier' 'Location' 'Northwest'
xlabel 'RPM'; ylabel 'Specific Fuel Consumpion (lbm/hr/hp)'
grid on

figure (7) % efficiency
plot([2000,2600:300:3200], conv_eff([1,3:5]), 'redx-', ...
    2600, conv_eff(6), 'kx', 2600, conv_eff(7),'mx',...
    rd, conv_eff(z), 'bluex-', 2300, conv_eff(2), 'redx')
title({'Figure 7', 'Efficiency vs. RPM with Different Amounts of Throttle'})
legend 'Full Throttle' '3/4 Throttle' 'Half Throttle' '1/4 Throttle' ...
    'Outlier' 'Location' 'Northeast'
axis( [ 2000 3200 0 20] )
xlabel 'RPM'; ylabel 'Efficiency (%)'
grid on

figure (8) % Brake Mean Effective Pressure
plot(fst, bmep(y), 'redx-', ...
    2600, bmep(6), 'kx', 2600, bmep(7),'mx',...
    rd, bmep(z), 'bluex-')
title({'Figure 8', 'Brake Mean Effective Pressure vs.' ...
    'RPM with Different Amounts of Throttle'})
legend 'Full Throttle' '3/4 Throttle' 'Half Throttle' '1/4 Throttle' ...
    'Location' 'best'
%axis( [2000 3900 0 300] )
xlabel 'RPM'; ylabel 'BMEP (psi)'
grid on

% further required calculations

eff_theory = 1 - 1/6^(.3);
stroke = 0.0619; % m
bore = 0.0651; % m
area_rod = pi * (bore/5)^2 / 4;
w = 2900 * 2 * pi / 60; % rad/s
sparks = 2900 / 2 / 60;
max_v = stroke/2*w;
stress = 0.3 * w^2 * stroke / 2 /  area_rod;
fprintf('The theoretical efficiency is %5.2f.\n', eff_theory)
fprintf('The maximum piston speed is %4.2f m/s \n', max_v)
fprintf('The number of sparks per second is %4.1f.\n', sparks)
fprintf('The approximate tensile stress in the connecting rod is %4.2f MPa.', stress/1e6)








