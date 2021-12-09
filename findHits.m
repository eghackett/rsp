function [hits, r,c] = findHits(mix,samples)
    mixFFT = fftshift(fft2(mix),2);
    rdb = 20*log10(abs(mixFFT(1:ceil(samples/2), :)));
    rdbComplex = mixFFT(1:ceil(samples/2), :);
    hits = rdbComplex(rdb>65);
    [r,c] = find(rdb>65);

end