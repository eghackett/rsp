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

rt1 = cos(2*pi.*fb1.*slow_time);
rt2 = cos(2*pi.*fb2.*slow_time);
figure('name','fbFig');
plot(slow_time(1,1:1500), txChirp1(1,1:1500))



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




% rt = zeros(samples, chirps); %mixer output
% for i=1:1:chirps
%         td = 2 .* R0 / c; %round trip delay
%         phi0 = 4*pi*fc.*R0/c; %inital phase
%         t = 0; % Reset
%     
%         for j=1:1:samples
% 
%             a = -2*pi*fc*2*V*i*tm/c ...    %phase shift
%                  -2*pi*((sweep_slope*td)*t);   %frequency
% 
%             frx = ((2*V*fc/c)-(bw/ts)*sweep_slope*td)*t;
%             rt(j,i) = 0.5*cos(a);
% 
%             t = t + 1/fs;
%         end
% 
% end

% rt2 = zeros(samples, chirps); %mixer output
% for i=1:1:chirps
%         td2 = 2 .* R2 / c; %round trip delay
%         
%         t = 0; % Reset
%     
%         for j=1:1:samples
% 
%             a = -2*pi*fc*2*V2*i*tm/c ...    %phase shift
%                  -2*pi*((sweep_slope*td2)*t);   %frequency
% 
%             rt2(j,i) = 0.5*cos(a);
% 
%             
%             t = t + 1/fs;
%         end
% 
%  end
rt = target2(samples, chirps, R0, V, fc, tm, ts, sweep_slope);
rt2 = target2(samples, chirps, R2, V2, fc, tm, ts, sweep_slope);



% Mixing the transmitted and received signal to produce the IF or Beat
% Frequency


finst = fc+(bw/ts).*slow_time;
myMix = ((rt+rt2)./2);
myMix = awgn(myMix,0);
fi = reshape(finst, samples, chirps);
% fif = ((rt+rt2)./2).*fi;
fif = myMix.*fi;
fbmax = range2beat(range_max,sweep_slope);


yt = lowpass(fif,fbmax,fs);



yt = reshape(yt,samples,chirps);

ytfft = fftshift(fft2(yt),2);


%% Form the range-Doppler map (RDM)



rangeBinAxis = (0:samples-1).*c/(2*bw);
dopplerBinSize = (1/tm)/chirps;
velocityBinAxis = (-chirps/2:chirps/2-1).*dopplerBinSize*lambda/2;



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




%% Sig Gen function
function [mix, fb, txChirp] = sigGen(samples, chirps, R0, V, fc, bw, tm, sweep_slope, fs)
    t = 0; %time
    c = 3e8;
    mix = zeros(samples, chirps); %mixer output
    fb = zeros(samples, chirps);
    txChirp = zeros(samples, chirps);
    for i=1:1:chirps
        td = 2 .* R0 / c; %round trip delay
        
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

function rt2 = target2(samples, chirps, R2, V2, fc, tm, ts, sweep_slope)


    c = 3e8;
    
    td2 = 2 .* R2 / c; %round trip delay
    
    t = repmat([0:ts:tm]',1,chirps);
    
    
    chirpMat = repmat(1:chirps,samples,1);
    
    a = -2*pi*fc*2*V2.*chirpMat*tm/c ...    %phase shift
                     -2*pi*((sweep_slope*td2).*t);   %frequency
    rt2 = 0.5*cos(a);


 end

