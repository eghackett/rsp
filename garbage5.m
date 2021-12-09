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

fbmx = (2*range_max*bw/tm) + ((2*v_max/c)*fc);

ts = 1/fs;
t = 0:ts:tm;
%sampling rate based on ADC datasheet
chirps = 64;        %frame size
samples = ceil(tm*fs);  %samples in one chirp

%% target


R0 = 110; %range in meters
V = 36; %radial velocity, m/s
R2 = 50; %range in meters
V2 = 45; %radial velocity, m/s


mix1 = theoreticalMix(samples, chirps, R0, V, fc, bw, tm, sweep_slope, fs);
mix2 = theoreticalMix(samples, chirps, R2, V2, fc, bw, tm, sweep_slope, fs);
% mix = sigGen(samples, chirps, [R0;R2], [V;V2], fc, bw, tm, sweep_slope, fs);
mix = (mix2+mix1)/2;

slow_time = zeros(1,chirps*samples);
slow_time(1:length(t)) = t(1:end);
 for i=1:chirps-1
    slow_time((i*length(t))+1:(i+1)*length(t)) = t + (tm*i);
 end

fd1 = (2*V/c)*fc;
fd2 = (2*V2/c)*fc;


figure
% plot(slow_time(1:1500),xt(1:1500));

finst = fc+(bw/ts).*t; % instantaneous frequency
finst = repmat(finst,1,chirps);
figure
plot(slow_time(1:1500),finst(1:1500));

td = 2 * R0 / c;
td2 = 2 * R2 / c;

zp = zeros(1,ceil(td/ts));
zp2 = zeros(1,ceil(td2/ts));

% The instantaneous frequency of the receved signal
frx = fc+(bw/ts).*(t-td);

frx = repmat(frx,1,chirps);

% The instantaneous frequency of the receved signal
frx2 = fc+(bw/ts).*(t-td2);
frx2 = repmat(frx2,1,chirps);

figure
plot(slow_time(1:1500),finst(1:1500),slow_time(1:1500)+td,frx(1:1500),slow_time(1:1500)+td2,frx2(1:1500));
legend('$$f_0$$', 'Target 1','Target 2','interpreter','latex')
ylabel('$$f_T(t)$$','interpreter','latex')
xlabel('t','interpreter','latex')

a1 = 2*pi*fc.*(slow_time-td) ...
    + pi*(bw/tm)*(slow_time-td).^2;

xr = cos(a1) ;%xt.*(slow_time-td)/2;

figure
plot(slow_time(1:1500),xr(1:1500))


f_if = finst-frx;

% The Chirped Transmitted Signal

st = cos(2*pi*fc.*slow_time+pi*(bw/ts).*slow_time.^2);

st = reshape(st, samples, chirps);

rt = target(samples, chirps, R0, V, fc, tm, ts, sweep_slope);
rt2 = target(samples, chirps, R2, V2, fc, tm, ts, sweep_slope);



% Mixing the transmitted and received signal to produce the IF or Beat
% Frequency


finst = fc+(bw/ts).*slow_time;
myMix = ((rt+rt2));
myMix = awgn(myMix,0);
fi = reshape(finst, samples, chirps);
% fif = ((rt+rt2)./2).*fi;
fif = myMix.*fi;
fbmax = range2beat(range_max,sweep_slope);


yt = lowpass(fif,fbmax,fs);



yt = reshape(yt,samples,chirps);

ytfft = fftshift(fft2(yt),2);

mixFFT = fftshift(fft2(awgn(mix,0)),2);


%% Form the range-Doppler map (RDM)




rangeBinAxis = (0:samples-1).*c/(2*bw);
dopplerBinSize = (1/tm)/chirps;
velocityBinAxis = (-chirps/2:chirps/2-1).*dopplerBinSize*lambda/2;



figure;
h = chirps*1.5;
surf(velocityBinAxis, rangeBinAxis(1:ceil(samples/2)), 20*log10(abs(ytfft(1:ceil(samples/2), :))));
% surf(velocityBinAxis, rangeBinAxis, 20*log10(abs(ytfft)));  % See the entire spectrum

xlabel("Velocity (m/s)");
ylabel("Range (m)");
axis tight;
shading flat;
view(0, 90);
colorbar;




%% Sig Gen function
function mix = theoreticalMix(samples, chirps, R0, V, fc, bw, tm, sweep_slope, fs)
    t = 0; %time
    c = 3e8;
    ts = 1/fs;
    td = 2 .* R0 / c; %round trip delay
    t = repmat([0:ts:tm]',1,chirps);
    chirpMat = repmat(1:chirps,samples,1);
    a = -2*pi*fc*2*V*chirpMat*tm/c -2*pi*((2*V.*(fc+chirpMat.*bw)/c + sweep_slope*td).*t); 
    mix = 0.5*cos(a);
end

function rt = target(samples, chirps, R2, V2, fc, tm, ts, sweep_slope)
    c = 3e8;
    td2 = 2 .* R2 / c; %round trip delay
    t = repmat([0:ts:tm]',1,chirps);
    chirpMat = repmat(1:chirps,samples,1);
    a = -2*pi*fc*2*V2.*chirpMat*tm/c -2*pi*((sweep_slope*td2).*t);   %frequency
    rt = 0.5*cos(a);
 end

