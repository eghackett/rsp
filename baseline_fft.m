function [mix, mixFFT] = baseline_fft(mk_heatmap)

    %%  Term Project Part 1
    %   Author: Edward Hackett
    %   Date: 10/26/2021
    
    %% Clear the workspace
    
%     close all;
%     clear all;
    if nargin < 1
        mk_heatmap = 1;
    end
    %% System Parameters
    
    c = 3e8;            %speed of light
    range_max = 180;    %max detection range
    range_res = 1;
    tm = 5.5*(2*range_max/c); %sweep time
    
    bw = c/(2*range_res);
    
    sweep_slope = bw/tm;
    v_max = 180*1000/3600;  %target max velocity
    fc = 77e9;          %radar frequency
    lambda = c/fc;      %radar wavelength
    fs = 72e6;          %sampling rate
    
    fbmx = (2*range_max*bw/tm) + ((2*v_max/c)*fc);
    
    ts = 1/fs;
    t = 0:ts:tm;
    
    chirps = 64;        
    
    samples = ceil(tm*fs);  %samples in one chirp
    
    %% target
    
    R0 = 110; %range in meters
    V = 36; %radial velocity, m/s
    R2 = 50; %range in meters
    V2 = 45; %radial velocity, m/s
    
    s_time = zeros(1,chirps*samples);
    s_time(1:length(t)) = t(1:end);
    chirpMat = tm.*repmat(0:chirps-1,samples,1);
    s_time = repmat([0:ts:tm]',1,chirps);
    s_time = s_time+chirpMat;
    s_time = reshape(s_time, 1, chirps*samples);
    
    
    fd1 = (2*V/c)*fc;
    fd2 = (2*V2/c)*fc;
    
    finst = fc+(bw/ts).*t; % instantaneous frequency
    finst = repmat(finst,1,chirps);
    figure
    plot(s_time(1:1500),finst(1:1500));
    
    td = 2 * R0 / c;
    td2 = 2 * R2 / c;
    
    %% Instantaneous Frequencies
    
    % The instantaneous frequency of the receved signal
    frx = fc+(bw/ts).*(t-td);
    
    frx = repmat(frx,1,chirps);
    
    % The instantaneous frequency of the receved signal
    frx2 = fc+(bw/ts).*(t-td2);
    frx2 = repmat(frx2,1,chirps);
    
    figure
    plot(s_time(1:1500),finst(1:1500),s_time(1:1500)+td,frx(1:1500),s_time(1:1500)+td2,frx2(1:1500));
    legend('$$f_0$$', 'Target 1','Target 2','interpreter','latex')
    ylabel('$$f_T(t)$$','interpreter','latex')
    xlabel('t','interpreter','latex')
    
    xr = cos(2*pi*fc.*(s_time-td) + pi*(bw/tm)*(s_time-td).^2) ;
    
    figure
    plot(s_time(1:1500),xr(1:1500))
    
    f_if = finst-frx;
    
    %% Signal Generation
    
    st = cos(2*pi*fc.*s_time+pi*(bw/ts).*s_time.^2);
    
    st = reshape(st, samples, chirps);
    
    rt = target(samples, chirps, R0, V, fc, tm, ts, sweep_slope);
    rt2 = target(samples, chirps, R2, V2, fc, tm, ts, sweep_slope);
    
    % Mixing the transmitted and received signal to produce the IF or Beat
    % Frequency
    
    %% Theoretical Received Signal
    
    mix1 = theoreticalMix(samples, chirps, R0, V, fc, bw, tm, sweep_slope, fs);
    mix2 = theoreticalMix(samples, chirps, R2, V2, fc, bw, tm, sweep_slope, fs);
    
    mix = (mix2+mix1)/2;
    
    %% Signal Processing
    
    finst = fc+(bw/ts).*s_time;
    myMix = ((rt+rt2));
    myMix = awgn(myMix,0);
    fi = reshape(finst, samples, chirps);
    % fif = ((rt+rt2)./2).*fi;
    fif = myMix.*fi;
    fbmax = range2beat(range_max,sweep_slope);
    
    
    %% Signal Filtering
    
    yt = lowpass(fif,fbmax,fs);
    
    yt = reshape(yt,samples,chirps);
    
    %% 2D-FFT
    ytfft = fftshift(fft2(yt),2);
    
    mixFFT = fftshift(fft2(awgn(mix,0)),2); % Theoretical
    
    %% Heatmap
    if mk_heatmap == 1
        make_heatmap(samples, bw, tm, chirps, mixFFT, lambda)
    end
    
end




