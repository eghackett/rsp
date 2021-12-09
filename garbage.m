close all;
clear all;
clc;
%% Params
c=3e8;
%%% Transmit side params
f0 = 10e9;
% dR = 15e-2;
% Rmax = 7.5e3;
% dV = 0.94;
% Vmax = 7.5;
B = 1e9;
T = 1e-5;
Ns = 2048;
chirps = 64;

% fdmax=1/(1*T);
%%% Receive side params
R1 = 70;
v1 = 50; % Corresponds to 43.2kmph
R2 = 30;
v2 = 20; % Co
%% Derived Params
%----------------
% If Vmax, Rmax, dR, dV specified
% T = c/(4*Vmax*f0);
% B=c/(2*dR);
% Ns = (4*B*Rmax)/c;
% L = ceil(c/(2*f0*dV*T)) ; % No.of Chirps
%-------------------------------
%If B, T, Ns, L are specified
Vmax =  240*1000/3600;
dR=c/(2*B);
m = B/T;
Rmax = Ns*c/(4*B);
dV = c/(2*f0*chirps*T);
%-------------------------------
% n=ceil(log10(Vmax));
% factor = roundn(Vmax,n)
% v1 = factor-v1;
%-------------------------------
t0 = 2*R1/c;
phi0 = 2*pi*f0*t0 - pi*m*(t0^2);
fb = 2*R1*m/c;
fd = -2*v1*f0/c;
fd2 = -2*v2*f0/c;
% fif_val1 = fb + fd;           % For comparison purpose
% fif_val2 = m*t0 +f0*2*v/c;
% v1=v1-T*1e8/2;
fif_val= fb + fd;
Ts = T/Ns;
Fs = 1/Ts;
%Therefore
t=0:Ts:T-Ts;
% 
%% Big time scale
time_scale = zeros(1,chirps*Ns);
time_scale(1:length(t)) = t(1:end);
%% For No.of chirps = L
 for i=1:chirps-1
    time_scale((i*length(t))+1:(i+1)*length(t)) = t + (T*i);
 end
%  time_scale=0:Ts:T*L-Ts;
 
% td=1e-6;
td=2*(R1+v1.*t)/c;
td2=2*(R2+v2.*t)/c;
% R= c*td/2
f_t = f0 + m*t;
% f_r = f0 + m*(t-td)/2;
%% For L chirps 
t=time_scale;
td=2*(R1+v1.*t)/c;

f_t = repmat(f_t,1,chirps);
% New -----------
f_r = zeros(size(f_t));
n = ceil(t0/Ts);
f_r(n+1:end) = f_t(1:end-n);
f_r = f_r + fd;


fr2 = zeros(size(f_t));
n = ceil(t0/Ts);
fr2(n+1:end) = f_t(1:end-n);
fd = (2*v2/c)*f0;
fr2 = fr2+fd2;



%----------
% f_r = repmat(f_r,1,L);
f_if = f_t-f_r;
% f_if(1:n) = 0;
st = cos(2*pi.*f_t.*t);
rt = cos(2*pi.*f_r.*t);
rt2 = cos(2*pi.*fr2.*t);
% rt = cos(2*pi.*f_r.*(t+td));  %%%%%%%%
% rt = cos(2*pi*(f0(t-td) + m*((t-td)^2)/2));
                                                                            % t = time_scale;
                                                                            % st = repmat(st,1,L);
% rt = rt+rt2;                                                                            % rt = repmat(rt,1,L);
                                                                            % f_t = repmat(f_t,1,L);
                                                                            % f_r = repmat(f_r,1,L);
fif = rt.*rt.*st;
fif_lpf = lowpass(fif,max(f_if),2*B,'Steepness',0.8);
%% Final IF signal
fif_the = 0.5*cos(phi0 + 2*pi*fif_val.*t);
% fif_the = 0.5*cos(phi0 + 2*pi.*f_if.*t);
%% Plots
% xlimit = 2*T;
xlimit = T/2;
%-------Fig 3 For Big time scale----------%
figure(3)
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
figure(5)
subplot(211)
plot(t,f_t);
% xlim([0 T/10])
hold on;
grid on;
plot(t,f_r);
ylim([f0-B f0+(2*B)]);
xlim([0 T*5])
legend('f_t','f_r')
subplot(212)
plot(t,f_if);
grid on;
xlim([0 T*5])
%% Post processing
radar_mat = reshape(fif_the,Ns,chirps); %% Using the Theoretical IF Signal
% radar_mat = reshape(fif_lpf,Ns,chirps); %% Using the Mixed and LPF IF Signal
%% Window function
window_1D = hann(size(radar_mat,1));
window_2D = hann(size(radar_mat,2));
%% FFT
rfft = (fft(radar_mat.*window_1D,[],1));
rfft = rfft./max(max(rfft)); %Normalization
rfft = rfft(1:size(rfft)/2,:);
% zeroPadding = zeros(size(rfft));
% rfft = vertcat(rfft,zeroPadding);
vfft = fft((rfft.*window_2D'),[],2);
%% Normalization
% normalize = vfft./max(max(vfft));
%vfft = fftshift(vfft,2)
vfft = vfft./max(max(vfft));
% vfft = fftshift(fft2(radar_mat), 2);
%% Range and Velocity vectors
R = 0:dR:Rmax-dR;
V = linspace(-Vmax, Vmax, chirps);
figure(4);
h=imagesc(V,R,20*log10(abs(fftshift(vfft,2))),[-60 0]);
cb = colorbar;
set(gca,'YDir','normal')
xlabel('Velocity (m/s)');
ylabel('Range (m)');