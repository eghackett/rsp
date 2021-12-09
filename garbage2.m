%%  Term Project Part 1
%   Author: Edward Hackett
%   Date: 10/26/2021

%% Clear the workspace

close all;
clear all;

%% System Parameters

c = 3e8;            %speed of light
range_max = 180;    %max detection range
range_res = 1;
tm = 5.5*(2*range_max/c); %sweep time
%tm is 7.2e-6 s
% bw = 200e6;         %sweep bandwidth
% bw = rangeres2bw(range_res,c);
bw = c/(2*range_res);

sweep_slope = bw/tm;
v_max = 180*1000/3600;  %target max velocity
fc = 77e9;          %radar frequency
lambda = c/fc;      %radar wavelength
fs = 72e6;          %sampling rate
ts = 1/fs;
t = 0:ts:tm;
%sampling rate based on ADC datasheet
chirps = 64;        %frame size
samples = ceil(tm*fs);  %samples in one chirp

%% target
% s = rng('hackett');
% myStream = RandStream('hackett');
% s = rand(myStream,1,5)

% a = 10;
% b = range_max;
% R1 = (range_max-a)*rand + a
% V1 = (v_max-a)*rand + a

R0 = 110; %range in meters
V = 36; %radial velocity, m/s
R2 = 50; %range in meters
V2 = 45; %radial velocity, m/s


[mix1, fb1, txChirp1] = sigGen(samples, chirps, R0, V, fc, bw, tm, sweep_slope, fs);
[mix2, fb2, txChirp2] = sigGen(samples, chirps, R2, V2, fc, bw, tm, sweep_slope, fs);
% mix = sigGen(samples, chirps, [R0;R2], [V;V2], fc, bw, tm, sweep_slope, fs);
mix = (mix2+mix1)/2;

slow_time = zeros(1,chirps*samples);
slow_time(1:length(t)) = t(1:end);
 for i=1:chirps-1
    slow_time((i*length(t))+1:(i+1)*length(t)) = t + (tm*i);
 end

fd1 = (2*V/c)*fc;
fd2 = (2*V2/c)*fc;
fb1 = reshape(fb1,1,samples*chirps);
txChirp1 = reshape(txChirp1,1,samples*chirps);
fb2 = reshape(fb2,1,samples*chirps);
% fb1 = fb1+fd1;
% fb2 = fb2+fd2;
rt1 = cos(2*pi.*fb1.*slow_time);
rt2 = cos(2*pi.*fb2.*slow_time);
figure('name','fbFig');
plot(slow_time(1,1:1500), txChirp1(1,1:1500))
% mix = awgn(mix,5, 'measured');

a = 2*pi*fc.*slow_time ...
    + pi*(bw/tm)*slow_time.^2;

xt = cos(a); % transmitted FMCW
figure
plot(slow_time(1:1500),xt(1:1500));

finst = fc+(bw/ts).*t; % instantaneous frequency
finst = repmat(finst,1,chirps);
figure
plot(slow_time(1:1500),finst(1:1500));

td = 2 * R0 / c;
td2 = 2 * R2 / c;

zp = zeros(1,ceil(td/ts));
zp2 = zeros(1,ceil(td2/ts));

% The instantaneous frequency of the receved signal
frx = (2*V*fc/c)+(bw/ts).*(t-td);
frx = repmat(frx,1,chirps);
% Zero padding corresponds to delay
frx = [zp, frx(1,length(zp):end-1)];

% The instantaneous frequency of the receved signal
frx2 = fc+(bw/ts).*(t-td2);
frx2 = repmat(frx2,1,chirps);
% Zero padding corresponds to delay
frx2 = [zp2, frx2(1,length(zp):end)];


figure
plot(slow_time(1:1500),finst(1:1500),slow_time(1:1500),frx(1:1500),slow_time(1:1500),frx2(1:1500));
legend('finst', 'frx')


a1 = 2*pi*fc.*(slow_time-td) ...
    + pi*(bw/tm)*(slow_time-td).^2;

xr = cos(a1) ;%xt.*(slow_time-td)/2;

figure
plot(slow_time(1:1500),xr(1:1500))


f_if = finst-frx;

% The Chirped Transmitted Signal
% st = cos(2*pi.*finst.*slow_time);
st = cos(2*pi*fc.*slow_time+pi*(bw/ts).*slow_time.^2);
% st = repmat(st', 1, chirps);
st = reshape(st, samples, chirps);
% The Received Signal
% rt = cos(2*pi.*frx.*slow_time);
% rt = cos(2*pi*fc.*(t-td)+pi*(bw/ts).*(t-td).^2);



rt = zeros(samples, chirps); %mixer output
for i=1:1:chirps
        td = 2 .* R0 / c; %round trip delay
        phi0 = 4*pi*fc.*R0/c; %inital phase
        t = 0; % Reset
    
        for j=1:1:samples
%             a = (-2*pi*fc*2.*V*i*tm/c ...    %phase shift
%                  -2*pi*(2.*V*(fc+i*bw)/c + sweep_slope.*td)*t);   %frequency
%             a = -2*pi*fc*2*V*i*tm/c ...    %phase shift
%                  -2*pi*((2*V*(fc+i*bw)/c + sweep_slope*td)*t);   %frequency
            a = -2*pi*fc*2*V*i*tm/c ...    %phase shift
                 -2*pi*((sweep_slope*td)*t);   %frequency

%             rt(j,i) = 0.5*cos(2*pi*V*2*i*fc*tm/c+pi*(bw/ts).*(t-td).^2); %% ALMOST1
            frx = ((2*V*fc/c)-(bw/ts)*sweep_slope*td)*t;
            rt(j,i) = 0.5*cos(a);
%             rt(j,i) = 0.5*cos(2*pi*V*2*i*fc*tm/c+pi*V*(bw/ts).*(t-td).^2); %% ALMOST
%             rt(j,i) = 0.5*cos(-2*pi*V*2*i*fc*tm/c-(2*pi*((2*V/c)*(fc)+bw/ts)*td)*t);
%             0.5*cos(2*pi*V*2*i*fc.*(t-td)+pi*i*(bw/ts).*(t-td).^2);
%                     2*V*fc/c+(bw/ts).*(t-td);
%             rt(j,i) = 0.5*cos(a);
            
            t = t + 1/fs;
        end

end

rt2 = zeros(samples, chirps); %mixer output
for i=1:1:chirps
        td2 = 2 .* R2 / c; %round trip delay
        
        t = 0; % Reset
    
        for j=1:1:samples

            a = -2*pi*fc*2*V2*i*tm/c ...    %phase shift
                 -2*pi*((sweep_slope*td2)*t);   %frequency

%             
            
            rt2(j,i) = 0.5*cos(a);

            
            t = t + 1/fs;
        end

 end





% Mixing the transmitted and received signal to produce the IF or Beat
% Frequency


finst = fc+(bw/ts).*slow_time;

fi = reshape(finst, samples, chirps);
fif = ((rt+rt2)./2).*fi;
% fif = fi.*conj(rt);
% Pass the received signal through a lowpass filter 
fif_lpf = lowpass(fif,max(f_if),2*bw,'Steepness',0.8);

yt = fif_lpf;
% yt = (xr.*xt);
% yt = (xt.* (xr)/2)/max(xt.* (xr)/2);

% figure; plot(slow_time(1:1500),yt(1:1500))

% fbeat = ((bw/tm)*td);
% fbeat = finst(8)-frx(1);%mean(finst)-mean(frx);
fbmax = (2*range_max/tm) + (2*v_max*fc/c);
fbmax = range2beat(range_max,sweep_slope);
% bp1 = fbeat-10e1;%(fbeat/10);
% bp2 = fbeat+10e1;%(fbeat/10);
min(yt)
max(yt)
% bandpass(yt,[1.5e6 3e6],fs)
% yt = bandpass(yt,[34.5463e6 34.9463e6],fs);
% yt = bandpass(yt,[1e6 5e6],fs);
% lowpass(yt,fbmax,fs);
yt = lowpass(yt,fbmax,fs);
% yt = lowpass(yt,range2beat(range_max,sweep_slope,c),fs);



% figure; plot(slow_time(1:1500),yt(1:1500))

yt = reshape(yt,samples,chirps);

ytfft = fftshift(fft2(yt),2);

% window_1D = hann(size(yt,1));
% window_2D = hann(size(yt,2));
% 
% ytfft = (fft(yt.*window_1D, samples, 1));
% ytfft = ytfft./max(max(ytfft)); %Normalization
% % ytfft = ytfft(1:size(ytfft)/2,:);
% % figure; plot(slow_time,ytfft)
% zeroPadding = zeros(size(ytfft));
% rfft = vertcat(ytfft,zeroPadding);
% 
% ytfft = fftshift(fft(ytfft.*window_2D', chirps, 2));

% mix = [mix zeros(1,9*length(mix))];
%% Form the range-Doppler map (RDM)


% RDM axes
rangeBinAxis = (0:samples-1).*c/(2*bw);
dopplerBinSize = (1/tm)/chirps;
velocityBinAxis = (-chirps/2:chirps/2-1).*dopplerBinSize*lambda/2;

% 2D FFT to perform range and Doppler compression (i.e. form the RDM)
% rdm = fftshift(fft2(mix), 2);
fbs = finst-frx;%finst(8:end)-frx(1:end-7);
% fbs = [fbs,zeros(1,7)];
fbs = reshape(fbs,samples,chirps);



rdm = fftshift(fft2(fbs), 2);


% Plot the RDM for the valid ranges of interest - targets ahead of you
figure;
h = chirps*1.5;
surf(velocityBinAxis, rangeBinAxis(1:ceil(samples/2)), 20*log10(abs(ytfft(1:ceil(samples/2), :))));
% surf(velocityBinAxis, rangeBinAxis, 20*log10(abs(rdm)));  % See the entire spectrum

xlabel("Velocity (m/s)");
ylabel("Range (m)");
axis tight;
shading flat;
view(0, 90);
colorbar;


sig_mix = reshape(mix,1,chirps*samples);
figure; plot(1:samples,sig_mix(1:samples))
title('signal to dechirp')

lambda = c/fc;  
v_max = 150*1000/3600; 
fr_max = (2*sweep_slope*range_max)/c;
% fb_max = fr_max+fd_max; % beat signal after mixing range and doppler freq
% fs = max(2*fb_max,bw);
% fd_max = speed2dop(2*v_max,lambda);
% 
% a = (-2*pi*fc*2.*V*i*tm/c -2*pi*(2.*V*(fc+i*bw)/c + sweep_slope.*td)*t);   %frequency
% sig_og = 0.5*cos(a);



t0 = 2*R0/c;
ft = fc+sweep_slope*t;
ft = repmat(ft,1,chirps);

fr = zeros(size(ft));
n = ceil(t0/ts);
fr(n+1:end) = ft(1:end-n);
fd = (2*V/c)*fc;
fr = fr+fd;

t0 = 2*R2/c;
ft = fc+sweep_slope*t;
ft = repmat(ft,1,chirps);

fr2 = zeros(size(ft));
n = ceil(t0/ts);
fr2(n+1:end) = ft(1:end-n);
fd = (2*V2/c)*fc;
fr2 = fr2+fd;

t = slow_time;
figure; 
plot(t,ft,t,fr, t, fr2)
% hold on;
grid on;
legend('Tx Signal','Target 1', 'Target 2')
% plot();
ylim([fc fc+(bw)]);
xlim([0 tm*2])

%%
t = 0:ts:tm;
t0 = 2*R2/c;
ft = fc+sweep_slope*t;
ft = repmat(ft,1,chirps);

fr2 = zeros(size(ft));
n = ceil(t0/ts);
fr2(n+1:end) = ft(1:end-n);
fd = (2*V2/c)*fc;
fr2 = fr2+fd;

% t = slow_time;
% figure; 
% plot(t,ft,t,fr)
% % hold on;
% grid on;
% legend('ft','fr')
% % plot();
% ylim([fc fc+(bw)]);
% xlim([0 tm*2])

st = cos(2*pi.*ft.*slow_time);
rt1 = cos(2*pi.*fr.*slow_time);
rt2 = cos(2*pi.*fr2.*slow_time);

Srx = (st+rt1+rt2);

figure;
plot(slow_time(1,1:250),rt1(1,1:250))

%% Sig Gen function
function [mix, fb, txChirp] = sigGen(samples, chirps, R0, V, fc, bw, tm, sweep_slope, fs)
    t = 0; %time
    c = 3e8;
    mix = zeros(samples, chirps); %mixer output
    fb = zeros(samples, chirps);
    txChirp = zeros(samples, chirps);
    for i=1:1:chirps
        td = 2 .* R0 / c; %round trip delay
        phi0 = 4*pi*fc.*R0/c; %inital phase
        t = 0; % Reset
    
        for j=1:1:samples
%             a = (-2*pi*fc*2.*V*i*tm/c ...    %phase shift
%                  -2*pi*(2.*V*(fc+i*bw)/c + sweep_slope.*td)*t);   %frequency
            a = -2*pi*fc*2*V*i*tm/c ...    %phase shift
                 -2*pi*((2*V*(fc+i*bw)/c + sweep_slope*td)*t);   %frequency
            mix(j,i) = 0.5*cos(a);
            fb(j,i) = ((2*V*(fc+i*bw)/c + sweep_slope*td)*t);
%             txChirp(j,i) = -2*pi*fc*2*V*i*tm/c ...    %phase shift
%                  -2*pi*((2*V*(fc+i*bw)/c + sweep_slope*1)*t);
            txChirp(j,i) = -2*pi*fc*2*V*i*tm/c ...    %phase shift
                 -2*pi*((2*V*(fc+i*bw)/c + sweep_slope*1)*t);
            t = t + 1/fs;
        end

    end
end

% function sigplot(R2, c, fc, sweep_slope, t, )
%     t0 = 2*R2/c;
%     ft = fc+sweep_slope*t;
%     ft = repmat(ft,1,chirps);
%     
%     fr = zeros(size(ft));
%     n = ceil(t0/ts);
%     fr(n+1:end) = ft(1:end-n);
%     fd = (2*V2/c)*fc;
%     fr = fr+fd;
%     
%     t = slow_time;
%     figure; 
%     plot(t,ft,t,fr)
%     % hold on;
%     grid on;
%     legend('ft','fr')
%     % plot();
%     ylim([fc fc+(bw)]);
%     xlim([0 tm*2])
% end