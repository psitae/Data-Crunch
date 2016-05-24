tic
clear
close all
clc

% reading the directory into a cell array is platform specific
if ispc % WINDOWS
    [aoeu, files0] = system( 'dir /B /S *.lvm');
    files1 = textscan( files0, '%s', 'delimiter', '\');
    files2 = files1{1};
    files = files2(4:4:end);
elseif isunix % Linux
    [aoeu, files0] = system( 'ls -l *.lvm');
    files1 = textscan( files0, '%s', 'delimiter', ' ');
    files = sort( files1{1}(9:9:end) )
end

psi_pa = 6.894757e3;
d = 2.195 - 1.647; % m
gam = 7/5;
R = 287;         % J/( kg * K )  
gam_he = 5/3;    
R_he = 2078;     % J/( kg * K )
atm = 101.3e3;   % Pa
T1 = 290;        % K
cp1 = 132890;    % Pa/V
cp2 = 243390;    % Pa/V
long = length(files);
[ press, data ] = deal( cell( 1, long ) );
[p2_rat, p4_rat, mach] = deal( zeros( 3, long) );
count = 1;

for i=1:long
    % read data from file
    fid = fopen( files{i} );
    for j = 1:23, fgetl(fid); end
    data{i} = fscanf(fid, '%f %f', [2,inf])';
    
    % characterize data
    p0 = textscan( files{i}, '%f' ); 
    p4 = p0{1}*psi_pa + atm;
    type = textscan( files{i}, '%s', 'delimiter', '_');
    type = strrep( type{1}{2}, '.lvm', '');
    if strcmp( type, 'open' )
        fig_point(count) = i;
        count = count + 1;
    end
    if strcmp( type, 'helium' )
        R4 = R_he;
        gam4 = gam_he;
    else
        R4 = R;
        gam4 = gam;
    end
  
    % calculate mach number of shock wave
    del_t = (    find( data{i}(:,2) > 0.05, 1)  ...
               - find( data{i}(:,1) > 0.05, 1)    )/1e6;
    c1 = sqrt(gam*R*T1);
    w = d/del_t;
    mach(3,i) = w/c1;
    
    % calculate p2 from mach
    p2_p1 = 1 + (mach(3,i)^2-1)*2*gam/(1+gam);                       %% eq 4,10
    p2_rat(1,i) = p2_p1;
    p2  = p2_p1 * atm;
    
    % process transducer data
    max_v(1) = max( data{i}(:,1) );
    max_v(2) = max( data{i}(:,2) );
    p2_p1_meas = mean( max_v .* [cp1 cp2] + atm );
    if strcmp(type, 'open') && ( p0{1} == 6 || p0{1} == 8 )
        cp(1) = ( p2 - atm ) / max_v(1);
        cp(2) = ( p2 - atm ) / max_v(2);
%         fprintf('For %20s: cp1 is %7.0f and cp2 in %7.0f\n', files{i}, cp(1), cp(2) )
    end
    
    %  calculate mach number using numerical iteration
    c4 = sqrt( gam4 * R4 * T1 );                                  %% eq 2
    mach_g = 0.01;
    guess = 0;        
    term2 = ( gam4 - 1 ) * c4 / ( (gam + 1) * c1 );     
    raise = -2*gam4/(gam4 - 1);
    while abs( guess*atm - p4 ) > 1000
        term1 = ( 2*gam*mach_g^2 - gam + 1 ) / ( gam + 1 );
        term3 = mach_g - mach_g^-1;
        guess = term1 * ( 1 - term2 * term3 )^raise;              %% eq 3
        mach_g = mach_g + 0.001;
        if mach_g > 2
            disp('A Loop didn''t work')
            break 
        end   
    end
    mach(1,i) = mach_g;
           
    % calculate theoretical p2_p1 using the shock tube equation
    guess = 0;
    ratio_g = 0;
    while abs( guess*atm - p4 ) > 500
        top = (gam4 - 1) * c1/c4 * (ratio_g - 1);
        bottom = sqrt( 2 * gam * (2*gam+(gam+1)*(ratio_g-1)) ); 
        guess = ratio_g * ( 1 - top/bottom )^raise;               %% eq 1
        ratio_g = ratio_g + 0.001;
        if ratio_g > 5
            disp('Loop didn''t work')
            break 
        end
    end
    
    % prepare data for figures
    mach(2,i) = sqrt( (gam+1)/(2*gam)*(p2_p1_meas - 1) + 1 )/c1;
    p4_rat(1,i) = p4/atm;
    p2_rat(2,i) = ratio_g;
    p2_rat(3,i) = abs( (p2_rat(2,i) - p2_rat(1,i)) / p2_rat(1,i) ); 
    string = strrep( files{i}, '_', ' ');
    string = strrep( string, '.lvm', '');
    string = strrep( string, 'psi', ' psi');
    psi = p0{1};
    name = sprintf('Figure %d', psi/4 + 0.5);    

    % in loop graphing   
    figure(i+6)
    press{i}(:,1) = smooth( data{i}(:,1)*cp1, 60 )/1000;
    press{i}(:,2) = smooth( data{i}(:,2)*cp2, 60 )/1000;
    x = ( 1:length(data{i}) )/1000;
    plot( x, press{i}(:,1), x, press{i}(:,2))
    title( {sprintf('Figure %d', i+6), 'Corrected Pressure Trace', string} )
    xlabel 'Milliseconds'
    ylabel 'Guage Pressure (kPa)'
    legend 'dist = 1.647 m' 'dist = 2.195 m' 'Location' 'SE'
    grid on
    
    if strcmp( type, 'open') && mod( p0{1} + 2, 4) == 0 
        abs_press(:,1) = ( data{i}(:,1)*cp1 + atm ) / 1000;
        abs_press(:,2) = ( data{i}(:,2)*cp1 + atm ) / 1000;
        figure(psi/4 + 0.5)
        plot(x, smooth(abs_press(:,1), 60), ...
             x, smooth(abs_press(:,2), 60) )
        title( {name 'Absolute Pressure vs. Time for Two Locations' string} )
        xlabel 'Milliseconds'
        ylabel 'Absolute Pressure (kPa)'
        axis([0 25 60 150 ])
        grid on
        legend 'dist = 1.647 m' 'dist = 2.195 m' 'Location' 'SE'
    end
    

    

end

% graph pressure ratios
figure(1)
less = [1:4 7];
plot( p4_rat(1,less), p2_rat(1,less), '*', ...
      p4_rat(1,less), p2_rat(2,less), '*')
grid on
title( {'Figure 1' 'Shock Strength Pressure Ratio vs. Initial Pressure Ratio'} )
legend 'Measured Values' 'Theoretical Values' 'Location' 'SE'
A = xlabel('$\frac{P4}{P1}$');
B = ylabel('$\frac{P2}{P1}$');
B.Rotation = 0;
latex( A, 1 )
latex( B, 1 )

cols1 = { 'P4_P1' 'Measured' 'Theoretical' 'Percent_Error' };
table1 = table( p4_rat(1,less)', p2_rat(1,less)', p2_rat(2,less)', ...
    100*p2_rat(3,less)', 'VariableNames', cols1);
fprintf('\nTable 1 - P2_P4 values as a function of P4_P1\n')
disp(table1)

% graph mach numbers
figure(5)
plot( p4_rat(1,fig_point), mach(1,fig_point), '*', ...
    p4_rat(1,fig_point), mach(2,fig_point), '*', ...
    p4_rat(1,fig_point), mach(3,fig_point), '*')
hold on
grid on
title({'Figure 5' 'Mach Numbers of Different Methods'})
C = xlabel('$\frac{P4}{P1}$');
latex(C, 1)
ylabel 'Mach Number'
legend 'Theoretical' 'Measured Pressure' 'Measured Time' ...
    'Location' 'SE'
 
cols2 = { 'P4_P1' 'Theoretical' 'Measured_Pressure' 'Percent_Error1' ...
    'Measured_Time' 'Percent_Error2'};
table2 = table( p4_rat(1,fig_point)', mach(1,fig_point)', mach(2,fig_point)', ...
    100 * abs( (p4_rat(1,fig_point) - mach(1,fig_point))./p4_rat(1,fig_point) )', ...
    mach(3,fig_point)', 100 * abs( (p4_rat(1,fig_point) - mach(3,fig_point))./p4_rat(1,fig_point) )', ...
    'VariableNames', cols2);
fprintf('\nTable 2 - Mach Numbers of Different Methods\n')
disp(table2)

% calculations for reflected shock using theory
m = mach(1,5);
% from temperature of state 2
t2_t1 = 1 + 2*(gam-1)*(gam*m^2 + 1)*(m^2-1) / ( (gam+1)^2 * m^2 );
t2 = t2_t1 * 290;
c2a = sqrt( R * gam * t2 );
% from ratio of speeds of sounnd
drift = 2*c1/(gam + 1)*(m - 1/m);
top1 = 2*gam*m^2 - gam + 1;
top2 = (gam - 1)*m^2 + 2;
bottom = (gam + 1)^2*m;
c2b = c1 * sqrt( top1 * top2 / bottom );
% close results, so average them
c2 = mean([c2a c2b]);
m_refl = roots([ 1 -drift*(gam+1)/(2*c2) -1  ]);
m_refl =  m_refl( m_refl > 0 );
p2 =  p2_rat(2,5) * atm;
p5 = p2 * (2*gam*m_refl^2 - gam + 1) / (gam + 1);

% compare with inspection of data 
p5_meas = mean( press{5}(press{5}>50) ) * 1000 + atm;

% report results
fprintf( 'The theoretical value for P5 is %5.2f kPa.\n', p5/1000 )
fprintf( 'The measured pressure is %5.2f kPa.\n', p5_meas/1000 )
fprintf( 'That''s %3.2f%% error!\n\n',  abs( (p5 - p5_meas )/p5 )*100 ) 


toc
