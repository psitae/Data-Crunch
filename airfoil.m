close all
clc

wings = {wing2412, wing1408, wing4412a, wing4412b, none};
in_m = 0.0254;
c = 2.5*in_m;
S = 2.5*5*in_m^2;
rho_w = 1000; % m^3/kg water
rho_a = 1.225; % air
atm = 101.3e3; % Pa
nu = 1.5e-5; % m^2/s


% this next loop performs the necessary calculations and
% makes a table for each wing
for i = 1:4
    turn = zeros(wings{i}.length,2);
    non_dims = zeros(wings{i}.length,2);
    wings{i}.forceOpen = wings{i}.forceOpen - wings{i}.forceClosed;
    wings{i}.pressOpen = wings{i}.pressOpen * in_m * 9.81 * rho_w;    
    wings{i}.('v') = sqrt( 2*wings{i}.pressOpen/rho_a );
    wings{i}.('Re') = wings{i}.v * c / nu;
    for j = 1:wings{i}.length
        turn(j,1) = wings{i}.forceOpen(j,1) * cos(wings{i}.sting(j)) ...
            - wings{i}.forceOpen(j,2) * sin(wings{i}.sting(j)); % lift
        turn(j,2) = wings{i}.forceOpen(j,2) * cos(wings{i}.sting(j)) ...
            + wings{i}.forceOpen(j,1) * sin(wings{i}.sting(j)); % drag
        non_dims(j,1) = turn(j,1) / (wings{i}.pressOpen * S); % Cl
        non_dims(j,2) = turn(j,2) / (wings{i}.pressOpen * S); % Cd
    end
    wings{i}.('lift') = turn(:,1);
    wings{i}.('drag') = turn(:,2);
    wings{i}.('cl') = non_dims(:,1);
    wings{i}.('cd') = non_dims(:,2);    
    % now to make a table
    fprintf('\n\n------------------------------------------------------------------')
    fprintf('\n\t\tThe Reynold''s Number for %s is %1.2e\n', names{i}, wings{i}.Re)
    cols = {'Sting' 'Attack' 'Lift' 'Drag' 'Cl' 'Cd'};
    T = table(wings{i}.sting*180/pi, wings{i}.attack, ...
        wings{i}.lift, wings{i}.drag, wings{i}.cl, wings{i}.cd, ...
        'VariableNames', cols);
    eval([names{i} '_table_of_values = T'])
end
% remove an outlier
wings{1}.cd(16) = NaN;
wings{1}.cl(16) = NaN;
wings{1}.attack(16) = NaN;


ttl = {'Lift' 'Drag'};
grp = {'cl' 'cd'};
for i = 1:2
    figure(i)
    for j = 1:4
        y = eval(['wings{j}.' grp{i}]);
        plot(wings{j}.attack, y, 'o-')
        hold on
        grid on
        title( sprintf('%s Coefficient vs. Angle of Attack', ttl{i}) )
        xlabel 'Angle of Attack (Degrees)'
        ylabel( sprintf('%s Coefficient', ttl{i}) )
    end
    if i==1
        theory_x = linspace(-2.5,22,1000);
        theory_y = (theory_x*pi/180)*2*pi;
        plot(theory_x, theory_y, '--')
        legend(names{1:4}, 'Theoretical Thin Airfoil', 'Location', 'Se')
    else
        legend(names{1:4}, 'Location', 'Se')
    end
    
end

figure (3)
for i=1:4
    plot(wings{i}.cl, wings{i}.cd, 'o-')
    hold on
    grid on
    title 'Drag Polar'
    xlabel 'Lift Coefficient'
    ylabel 'Drag Coefficient'
end
legend(names{1:4}, 'Location', 'nw')
    



        