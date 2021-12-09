%%  Term_Project_Part_1.m
%   Author: Edward Hackett

close all;
clear all;

r1 = 8;
v1 = 22;
r2 = 90;
v2 = 17;
%% Define the system parameters

% Using 77 GHz for typical automotive radar

fc = 77e9;
c = 3e8;            % speed of light
lambda = c/fc;      % wavelength

% System Goals
% range_max: maximum desired range;
% tm: sweep time based on maximum unambiguous range (range_max) and a
% factor multiplier of 5.5

range_max = 200;
tm = 5*range2time(range_max,c);
% tm = 2e-3;
% range_res: 1 meter range resolution
% bw: define the required bandwidth for each frequency sweep
% sweep_slope: the sweep slope

range_res = 1;
bw = rangeres2bw(range_res,c);
bw = 200e6;
% sweep_slope = bw/tm;
sweep_slope = bw/tm;

% fr_max: the beat frequency corresponding to the maximum range
% fr_max = range2beat(range_max,sweep_slope,c);
fr_max = (2*sweep_slope*range_max)/c;

%% Define the target
% v_max: the max velocity of the target
% fd_max: the max Doppler frequency of the target
% fb_max: the max beat frequency of the target
% fs: the sample frequency defined as twice the beat frequency and
% bandwidth

v_max = 230*1000/3600;
fd_max = speed2dop(2*v_max,lambda);

fb_max = fr_max+fd_max; % beat signal after mixing range and doppler freq

fs = max(2*fb_max,bw); % sample frequency

fif_val = fr_max; %fb_max + fd_max;
T = c/(4*v_max*fc);
Ns = ceil((4*bw*range_max)/c);
dV = 0.8;
L = ceil(c/(2*fc*dV*T)) ; % No.of Chirps


Ts = tm; %T/Ns; % time sample
Fs = fs; %1/Ts; % frequency sample
%Therefore
t=0:Ts:T-Ts;
%% Big time scale
time_scale = zeros(1,L*Ns);
time_scale(1:length(t)) = t(1:end);

%% For No.of chirps = L
 for i=1:L-1
    time_scale((i*length(t))+1:(i+1)*length(t)) = t + (T*i);
 end
%  time_scale=0:Ts:T*L-Ts;
 
% td=1e-6;
td1=2*(r1+v1.*t)/c;
% R= c*td/2
f_t = fc + sweep_slope*t;
% f_r = f0 + m*(t-td)/2;

%% For L chirps 
t1 = 2*r1/c; 
t2 = 2*r2/c;

phi0 = 2*pi*fc*t1 - pi*sweep_slope*(t1^2);
phi2 = 2*pi*fc*t2 - pi*sweep_slope*(t2^2);

t=time_scale;
% time delay of signal
% tau = (2*(R+v*t))/c

td1=2.*(r1+(v1.*t))./c;
td2=2.*(r2+(v2.*t))./c;




f_t = repmat(f_t,1,L);
% New -----------
f_r = zeros(size(f_t));
% f_r = (bw/2).*sawtooth((t-td1)*2*pi)+2*v1*fc/c;
n = ceil(t1/Ts);
f_r(n+1:end) = f_t(1:end-n);
f_r = f_r + fd_max;
%----------
% f_r = repmat(f_r,1,L);
f_if = f_t-f_r;
% f_if(1:n) = 0;
st = (bw/2)*sawtooth(2*pi.*f_t.*t);
rt = (bw/2)*sawtooth(2*pi.*f_r.*t);
% rt = cos(2*pi.*f_r.*(t+td));  %%%%%%%%
% rt = cos(2*pi*(f0(t-td) + m*((t-td)^2)/2));

                                                                            % t = time_scale;
% f_r2 = (bw/2).*sawtooth((t-td2)*2*pi)+2*v2*fc/c;                                                                            % st = repmat(st,1,L);
                                                                            % rt = repmat(rt,1,L);
                                                                            % f_t = repmat(f_t,1,L);
                                                                            % f_r = repmat(f_r,1,L);
f_r2 = zeros(size(f_t));
n2 = ceil(t2/Ts);
f_r2(n2+1:end) = f_t(1:end-n2);
f_r2 = f_r2 + fd_max;
%----------
% f_r = repmat(f_r,1,L);
f_if2 = f_t-f_r2;
% f_if(1:n) = 0;
st2 = (bw/2)*sawtooth(2*pi.*f_t.*t-td2);
rt2 = (bw/2)*sawtooth(2*pi.*f_r2.*t-td2);

fif = (rt.*st +rt2.*st)/4;
fif_lpf = lowpass(fif,max(f_if),2*bw,'Steepness',0.8);
% fif_lpf_up = lowpass(fif(:,1:2:end),max(f_if),2*bw,'Steepness',0.8);
% fif_lpf_down = lowpass(fif(:,2:2:end),max(f_if),2*bw,'Steepness',0.8);
% fif_lpf = (fif_lpf_up+fif_lpf_down)./2;

%% Final IF signal
% x = sawtooth(2*pi*50*t,1/2);
% fif_the = 0.5*sawtooth(phi0 + phi2 + 2*pi*fif_val.*t);
% fif_the = 0.5*cos(phi0 + phi2 + 2*pi*fif_val.*t-td);
fif_the = 0.5*sawtooth(phi0 + + phi2 + 2*pi.*f_if.*t);
% fif_the = 0.5*cos(phi0 + 2*pi.*f_if.*t);

%% Post processing
% radar_mat = reshape(fif_the,Ns,L); %% Using the Theoretical IF Signal
radar_mat = reshape(fif_lpf,Ns,L); %% Using the Mixed and LPF IF Signal
%% Window function
window_1D = hann(size(radar_mat,1));
window_2D = hann(size(radar_mat,2));
%% FFT
rfft = (fft(radar_mat.*window_1D,[],1));
rfft = rfft./max(max(rfft)); %Normalization
rfft = rfft(1:size(rfft)/2,:);
% zeroPadding = zeros(size(rfft));
% rfft = vertcat(rfft,zeroPadding);
vfft = fft(rfft.*window_2D',[],2);
%% Normalization
% normalize = vfft./max(max(vfft));
%vfft = fftshift(vfft,2)
vfft = vfft./max(max(vfft));

% vfft1 = postProcess(fif_lpf_up, f_if,bw, Ns, L, range_res, v_max, range_max);
% vfft2 = postProcess(fif_lpf_down, f_if,bw, Ns, L, range_res, v_max, range_max);
% 
% vfft = [vfft1 vfft2];

%% Range and Velocity vectors
R = 0:range_res:range_max-range_res;
V = linspace(-v_max, v_max, L);
figure(4);
h=imagesc(V,R,20*log10(abs(fftshift(vfft,2))),[-60 0]);
cb = colorbar;
set(gca,'YDir','normal')
xlabel('Velocity (m/s)');
ylabel('Range (m)');

deltaR = rdcoupling(fd_max,sweep_slope,c)

%% Plots
% xlimit = 2*T;
xlimit = T/2;
f0 = fc;
%-------Fig 3 For Big time scale----------%
figure
subplot(511)
plot(t,st);
xlim([0 xlimit])
title("Received signal as st = cos(2*pi.*f_t.*t);")
subplot(512)
plot(t,rt);
xlim([0 xlimit])
title("Received signal as rt = cos(2*pi.*f_r.*t);")
subplot(513)
plot(t,fif);
xlim([0 xlimit])
title("IF after Mixing")
subplot(514)
plot(t,fif_lpf);
xlim([0 xlimit])
title("IF after LPF")
subplot(515)
plot(t,fif_the);
xlim([0 xlimit])
title("IF from fif_the = 0.5*cos(phi0 + 2*pi*fif_val.*t);")
figure
subplot(211)
plot(t,f_t);
% xlim([0 T/10])
hold on;
grid on;
plot(t,f_r);
ylim([f0-bw f0+(2*bw)]);
xlim([0 T*5])
legend('f_t','f_r')
subplot(212)
plot(t,f_if);
grid on;
xlim([0 T*5])
%% Generate a table of the system parameters
%%
% The following table summarizes the radar parameters.
% 
%  System parameters            Value
%  ----------------------------------

row = cell(9,1);
% row = cell(8,1);
row(1,1) = {'System parameters & Value \\ \hline \hline'};
row(2,1) =  {['Operating frequency (GHz) & ', num2str(fc/10e8), ' \\ \hline']} ;
row(3,1) = {['Maximum target range (m) & ', num2str(range_max), ' \\ \hline']};
row(4,1) = {['Range resolution (m)  & ', num2str(range_res), ' \\ \hline']};    
row(5,1) = {['Maximum target speed (km/h)  & ', num2str(v_max), ' \\ \hline']};
row(6,1) = {['Sweep time (microseconds) &  ', num2str(tm*10e5), '\\ \hline']};
row(7,1) = {['Sweep bandwidth (MHz) & ', num2str(bw), '\\ \hline']}   ;
row(8,1) = {['Maximum beat frequency (MHz) & ', num2str(fb_max), '\\ \hline']};
row(9,1) = {['Sample rate (MHz)  & ', num2str(fs), '\\ \hline']} ;
r = cell2table(row);

writecell(row,'./system_parameters.tex', 'FileType', 'text')
% fileID = fopen('./system_parameters.tex','w');
% for i = 1:9
%     fprintf(fileID, row(i,1));
% end
% fprintf(fileID, r);

%% Generate the FMCW waveform

% With all the information above, one can set up the FMCW waveform used
% in the radar system.

waveform = phased.FMCWWaveform('SweepTime',tm,'SweepBandwidth',bw,...
    'SampleRate',fs);
phased.LinearFMWaveform()



%%
% This is a up-sweep linear FMCW signal, often referred to as sawtooth
% shape. One can examine the time-frequency plot of the generated signal.

% sig = waveform();
% subplot(211); plot(0:1/fs:tm-1/fs,real(sig));
% xlabel('Time (s)'); ylabel('Amplitude (v)');
% title('FMCW signal'); axis tight;
% subplot(212); spectrogram(sig,32,16,32,fs,'yaxis');
% title('FMCW signal spectrogram');

%%

%%
T = tm;
t = 0:1/fs:T-1/fs;
f_t = bw/2*sawtooth(t*fb_max*pi);
r1 = 8;
v1 = 22;
r2 = 113;
v2 = 17;
delay1 = 2*(r1+v1*t)/c;
delay2 = 2*(r2+v2*t)/c;
f_r1 = (bw/2).*sawtooth((t-delay1)*2000*pi)+2*v1*fc/c;
f_r2 = (bw/2).*sawtooth((t-delay2)*2000*pi)+2*v2*fc/c;

offset = 0;
signal_t = cos(2*pi*(fc).*t+2*pi*cumsum(f_t)*dt+offset);
signal_r1 = cos(2*pi*fc.*(t-delay1)+2*pi*cumsum(f_r1)*dt+offset);
signal_r2 = cos(2*pi*fc.*(t-delay2)+2*pi*cumsum(f_r2)*dt+offset);
signal_r = (signal_r1+signal_r2)/2;

%mixing - normalized sgnal
mixed = (signal_t.* (signal_r)/2)/max(signal_t.* (signal_r)/2)+noise;
%%% time domain signals
figure;
subplot(3,1,1)
plot(t*10^3,f_t/1e6,t*10^3,real(mixed),'+-')
legend('Tx','Rx')
xlabel('time [ms]');
ylabel('frequency [MHz]');
xlim([0,3])
subplot(3,1,2)
plot(t*10^3,f_t/1e6,t*10^3,(f_r1+f_r2)/2e6,'+-')
%axis([-100e3,100e3,0,1e-3)]
legend('Tx','Rx')
xlabel('time [ms]');
ylabel('frequency [MHz]');
%xlim([0.992,1.004])
%ylim([99,100.2])
%legend('Tx','Rx')
subplot(3,1,3)
plot(t*10^3,f_t/1e6+(f_r1+f_r2)/2e6)

function vfft = postProcess(fif, f_if,bw, Ns, L, range_res, v_max, range_max)

    fif_lpf = lowpass(fif,max(f_if),2*bw,'Steepness',0.8);
    %% Final IF signal
    
    
    %% Post processing
    % radar_mat = reshape(fif_the,Ns,L); %% Using the Theoretical IF Signal
    radar_mat = reshape(fif_lpf,Ns,L/2); %% Using the Mixed and LPF IF Signal
    %% Window function
    window_1D = hann(size(radar_mat,1));
    window_2D = hann(size(radar_mat,2));
    %% FFT
    rfft = (fft(radar_mat.*window_1D,[],1));
    rfft = rfft./max(max(rfft)); %Normalization
    rfft = rfft(1:size(rfft)/2,:);
    % zeroPadding = zeros(size(rfft));
    % rfft = vertcat(rfft,zeroPadding);
    vfft = fft(rfft.*window_2D',[],2);
    %% Normalization
    % normalize = vfft./max(max(vfft));
    %vfft = fftshift(vfft,2)
    vfft = vfft./max(max(vfft));
    %% Range and Velocity vectors
%     R = 0:range_res:range_max-range_res;
%     V = linspace(-v_max, v_max, L);
%     figure(4);
%     h=imagesc(V,R,20*log10(abs(fftshift(vfft,2))),[-60 0]);
%     cb = colorbar;
%     set(gca,'YDir','normal')
%     xlabel('Velocity (m/s)');
%     ylabel('Range (m)');

end
