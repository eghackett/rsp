function mix = theoreticalMix(samples, chirps, R, V, fc, bw, tm, sweep_slope, fs)
    t = 0; %time
    c = 3e8;
    ts = 1/fs;
    td = 2 .* R / c; %round trip delay
    t = repmat([0:ts:tm]',1,chirps);
    chirpMat = repmat(1:chirps,samples,1);
    a = -2*pi*fc*2*V*chirpMat*tm/c -2*pi*((2*V.*(fc+chirpMat.*bw)/c + sweep_slope*td).*t); 
    mix = 0.5*cos(a);
end
