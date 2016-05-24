clear
clc

string = {'sine wave' 'square wave' 'resonance' ...
    '1st location mode 2' '2nd location end of beam'};

tdata = cell(1,5);              
fdata = cell(1,5);
for j = 1:5                     % extract data from files w/ fscanf
    a = fopen(string{1,j});
    for i = 1:5, skip = fgetl(a); end
    temp_tdata = fscanf(a, '%f %*c %f %*c %f', [3, inf]);
    for i = 1:2, skip2 = fgetl(a); end
    temp_fdata = fscanf(a, '%f %*c %f %*c %f', [3, inf]);
    tdata{1,j} = temp_tdata';
    fdata{1,j} = temp_fdata';
    fclose(a);
end
                              
for i = 1:5                     % convert to millivolts
    tdata{1,i}(:,2:3) = 1000 * tdata{1,i}(:,2:3);
end                          
                                
i=1;                            % remove zeros
counter = 1;
removed = zeros(5,1e3); % the purpose of this matrix is the check if 
for j=1:5                           % the data has been removed from the 
    while i<length(tdata{1,j})                  % right place
        if tdata{1,j}(i,2)==0 && tdata{1,j}(i,3)==0
            tdata{1,j}(i,:) = [];
            removed(counter)=i;
            counter = counter + 1;
        end
        i = i + 1;
    end
end
%%                                create sine wave and square wave graphs
                                                       
for i = 1:2
    figure(i)
    plot(fdata{1,i}(:,1),fdata{1,i}(:,2))
    title(string{1,i})
    xlabel('Frequency (Hz)')
    ylabel('Amplitude')
    axis([0 40 0 .06])
    grid on
end

[Amax, inx] = max(fdata{1,1}(:,2)); % to create sinusoid wave from
freq_max = fdata{1,1}(inx,1);      % laplace transform information
domain = linspace(0,2,1e4);
range = Amax * 1e3 * sin(freq_max*domain*2*pi);

figure(3)
plot(tdata{1,1}(:,1), tdata{1,1}(:,2),':', domain, range,'--')
title('Sine Wave with Spectrum Analysis')
legend('Measured Data', 'Laplace Transform Results')
xlabel('Time (s)')
ylabel('Voltage (mV)')
grid on

[peaks, sq_max] = findpeaks(fdata{1,2}(:,2));
[peaks, sq_max] = remove_zeros(peaks,sq_max);

wave1 = 1e3 * peaks(1) * sin(domain*2*pi*2);
wave2 = zeros(1,1e4);
wave3 = zeros(1,1e4);

for i=1:5
    wave2 = wave2 + 1e3 * peaks(i) * sin(domain*2*pi*(4*i-2));
end    

i=1;
while i<length(peaks)
    wave3 = wave3 + 1e3 * peaks(i) * sin(domain*2*pi*(4*i-2));
    i = i + 1;
end

figure(4)
subplot(3,1,1)
plot(domain, wave1, 'g', tdata{1,2}(:,1), tdata{1,2}(:,2),':')
title('Square Wave and N=1')
legend('N=1','Measured Data')
xlabel('Time (s)')
ylabel('Voltage (mV)')
grid on

subplot(3,1,2)
plot(domain, wave2, 'g', tdata{1,2}(:,1), tdata{1,2}(:,2),':')
title('Square Wave and N=5')
legend('N=5','Measured Data')
xlabel('Time (s)')
ylabel('Voltage (mV)')
grid on

subplot(3,1,3)
plot(domain, wave3, 'g', tdata{1,2}(:,1), tdata{1,2}(:,2),':')
title('Square Wave and IFT')
legend('IFT','Measured Data')
xlabel('Time (s)')
ylabel('Voltage (mV)')
grid on

%%                               Free Vibration

% constants
rho = 2700;    % kg/m^3
E = 68.9e9;    % Pa
I = 4.336e-9;  % m^-4
L = 32*.0254;  % m
p = 8.710e-1;  % kg/m

eig = [1.875 4.694 7.855]'; % not correct
nfreq = eig.^2 * sqrt( E*I / (p*L^4) )/(2*pi) % Hz

counter = 5;
for i = 5:-1:4     % plots of time and frequency
    figure(counter)
    plot(fdata{1,i}(:,1),fdata{1,i}(:,2),fdata{1,i}(:,1),fdata{1,i}(:,3))
    title(string{1,i})
    xlabel('Frequency (Hz)')
    ylabel('Amplitude')
%     axis([0 40 0 .06])
    grid on
    
    figure(counter + 1)
    plot(tdata{1,i}(:,1),tdata{1,i}(:,2),tdata{1,i}(:,1),tdata{1,i}(:,3))   
    title(string{1,i})
    xlabel('Time (s)')
    ylabel('Voltage (mV)')
%     axis([.3 .7 -700 700])
    grid on
    counter = counter + 2;
end

% averaging the distance
% between all peaks (better than just 7). The results
% are graphed and analyzed

[maxes, inx] = findpeaks (tdata{1,4}(:,3),'MinPeakDistance',142);
t_crit = tdata{1,4}(inx,1);

% fit data was fitted and graphed
domain = linspace(0,5,1e4);
fit = 706.3*exp(-.349*domain);

figure (9)
plot(t_crit, maxes, domain, fit)
title('Peaks of Impulse Response w/ Exponential Envelope Fit Curve')
xlabel('Time (s)')
ylabel('Voltage (mV)')
legend('Relative Maxima', 'f(x) = 706.3*exp(-.349x)')
grid on

el = length(t_crit);
delta_t = zeros(1,el-1);

for i = 1:el-1
    delta_t(i) = t_crit(i+1) - t_crit(i);
end
delta_t(end) = [];
av_t = mean( delta_t );
nfreq_data = 1/av_t; %Hz

for i=1:el-1
    d(i) = log ( maxes(i)/maxes(i+1) );
    xi(i) = d(i)/(2*pi);
end
nextmax = maxes; %rearrange sizes so they all fit in the table
nextmax(1) = [];
maxes(end) = [];

result = table([1:48]', maxes, nextmax, nextmax-maxes, d', xi', ...
    'VariableNames', {'Peak_Number' 'Peak' 'Next_Peak' 'Difference' ...
    'Logarithmic_Decent' 'Xi'} );

%%                          Resonance

a = fopen('forced vibration');
temp = textscan(a, '%s', 'Delimiter','\n');
data_string = cell(1,41);
for i=6:4:167
    data_string((i-2)/4) = temp{1,1}(i,1);
end

adata = zeros(41,5);
for i = 1:41    
    adata(i,:) = str2num(data_string{1,i});
end
   
result2 = table(adata(:,2),adata(:,3),adata(:,5), ...
    'VariableNames', {'Frequency' 'Input_Amplitude' 'Output_Amplitude'} );
    
amp_r = zeros (1,41);
for i=1:41
    amp_r(i) = ( adata(i,5)*adata(41,3) ) / ( adata(i,3)*adata(41,5) );
end

freq=12.5;
freq_r = adata(:,2)/(freq);
domain = linspace(0,50,1e4);
ksi = mean(xi);
magnitude = ( ( 1-(domain/freq).^2 ).^2 + (2*ksi*domain/freq).^2 ).^-0.5;

figure(10)
plot(freq_r,amp_r, domain*4/50, magnitude)
axis( [0 1.5 0 40] )

figure(11)
subplot(2,1,1)
plot(tdata{1,3}(:,1),tdata{1,3}(:,2),tdata{1,3}(:,1),tdata{1,3}(:,3))
legend('Excitation', 'Response')
title('Time Response 0-2s')
grid on

subplot(2,1,2)
plot(tdata{1,3}(:,1),tdata{1,3}(:,2),tdata{1,3}(:,1),tdata{1,3}(:,3))
legend('Excitation', 'Response')
title('Resonance Time Response 0-0.5s')
axis( [0 .5 -1000 1000] )
grid on


