function rt = target(samples, chirps, R, V2, fc, tm, ts, sweep_slope)
    c = 3e8;
    td = 2 .* R / c; 
    t = repmat([0:ts:tm]',1,chirps);
    chirpMat = repmat(1:chirps,samples,1);
    a = -2*pi*fc*2*V2.*chirpMat*tm/c -2*pi*((sweep_slope*td).*t);   %frequency
    rt = 0.5*cos(a);
 end
