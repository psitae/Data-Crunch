clear
clc
%------------- THE FIRST TASK IS TO ORGANIZE THE DATA ----------------
%---------------------------------------------------------------------
% the following algorithm searches the raw data for blocks of data
% of length = 4 and saves them into a cell array 

files = {'many.xlsx' 'NACA_2412.xlsx' 'NACA_4412.xlsx'};
counter = 0;
data = cell(1,20);

for i = 1:3
    f = files{i};
    [nothing, sheets] = xlsfinfo(f);
    for j = 1:2
        tic
        temp = xlsread(f, sheets{j});
        toc
        s = size(temp);
        go_down = 1;
        while 1
            if go_down > s(1), break, end
            test = zeros(1,4);
            for k = 1:s(2)-3
%               fprintf('testing line %d, col %d', go_down, k)
                for m = 0:3, test(m+1) = ~isnan( temp(go_down,k+m) ); end
%               disp(test)
                if all( test ) % four numbers in a row?
                    status = 'entered';
%                   fprintf('entered --------------------\n\n\n')
                    counter = counter + 1;   
                    number_down = 1;
                    while 1    % how far down does in go?
%                       fprintf('Current box (%d, %d) %d\n',number_down+go_down,k,temp(number_down+go_down,k))
                        if number_down + go_down > s(1), break
                        elseif ~isnan( temp(number_down+go_down,k) ) == 1 
                            number_down = number_down + 1;
                        else
                            break
                        end
                    end
                    data{counter} = temp(go_down:go_down+number_down-1,k:k+3);
%                   fprintf('length of section is %d\n', number_down)
                end
            end
            if strcmp(status,'entered')
                go_down = go_down + number_down;
                status = 'renewed';
            else
                go_down = go_down + 1;
            end
        end 
    end
end


% the data from data matrix is then converted to structures for each wing

fields = {'pressOpen' 'forceOpen' 'pressClosed' 'forceClosed' ...
    'attack' 'sting' 'length'};
names = {'wing2412' 'wing1408' 'wing4412a' 'wing4412b' 'none' };
location = [ 12 11 8 7 10 9 14 13 2 1];

for i = 1:5
    a = location(2*i-1);
    b = location(2*i);
    press = data{a}(:,4);
    stats = { mean( press(~isnan(press)) ), ...
        [ data{a}(:,2), data{a}(:,3) ], ...
        mean(data{b}(:,4)), ...
        [ data{b}(:,2), data{b}(:,3) ], ...
        data{b}(:,1), ...
        (data{b}(:,1)-9)*pi/180, ...
        length(data{b}) };
    stats = cell2struct(stats, fields, 2);
    eval([names{i} '=stats'])
    clc
end

% one last fix needs to be made
wing4412b.forceOpen = wing4412b.forceOpen(16:-1:1,:);

disp('The data is now organized!')
clearvars -except wing2412 wing1408 wing4412a wing4412b none names


    
