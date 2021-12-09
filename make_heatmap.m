function make_heatmap(samples, bw, tm, chirps, ytfft, lambda)

 %% Heatmap
    c = 3e8;
    
    rangeBin = (0:samples-1).*c/(2*bw);
    dopplerBin = (1/tm)/chirps;
    velocityBin = (-chirps/2:chirps/2-1).*dopplerBin*lambda/2;
    
    figure;
    surf(velocityBin, rangeBin(1:ceil(samples/2)), 20*log10(abs(ytfft(1:ceil(samples/2), :))));
    xlabel("Velocity (m/s)");
    ylabel("Range (m)");
    axis tight;
    shading flat;
    view(0, 90);
    colorbar;

end