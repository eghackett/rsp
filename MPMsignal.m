function mix = MPMsignal(A, R,d,theta, samples, chirps, V, fc, bw,tm,sweep_slope,fs)
    n = 0:19;
    R = R+d*sind(theta(1))*n;
    
    R = R';
    mix = zeros(samples, chirps);
%     mix = arrayfun(@(x) theoreticalMix(samples, chirps,x, V, fc, bw, tm, sweep_slope, fs), R,'UniformOutput',false);
    for i=1:20
        mix1 = theoreticalMix(samples, chirps, R(i,1), V, fc, bw, tm, sweep_slope, fs);
        
        mix(:,i) = mix1(:,1);
        
    end
%     mix = cell2mat(mix);
    mix = A(:,1)*mix(:,1)';

end