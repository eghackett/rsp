function y = MPM(d, theta, N, s)

    if nargin < 1
        c = 3e8;
        lambda = (c/77e9);
        d = lambda/2;
        N = 20;
        theta = [15 63];

    end

    % samples/3 < P < samples/2
    p = 8; 

    configuration; %load configuration profile

%     samples = 20;

    theta = deg2rad(theta);
    n = [0:N-1]';
    A = exp(1j.*n*2*pi*d.*sin(theta)/lambda); %Steering vector; no of rows should equal N and col = no of theta
    R0d = R0+d*sind(theta(1))*n;
    R2d = R2+d*sind(theta(1))*n;
    
    mix11 = zeros(samples, chirps);
    mix22 = zeros(samples, chirps);
    for i=1:N
        mix11 = theoreticalMix(samples, chirps, R0d(i,1), V, fc, bw, tm, sweep_slope, fs);
        mix22 = theoreticalMix(samples, chirps, R2d(i,1), V2, fc, bw, tm, sweep_slope, fs);
        mix1(:,i) = mix11(:,1);
        mix2(:,i) = mix22(:,1);
    end

    hit1 = zeros(1, N);
    hit2 = zeros(1, N);
    hits = zeros(2,N);
    for i=1:N
        mix3 = theoreticalMix(samples, chirps, R0d(i,1), V, fc, bw, tm, sweep_slope, fs);
        mix4 = theoreticalMix(samples, chirps, R2d(i,1), V2, fc, bw, tm, sweep_slope, fs);
        mix = (mix3+mix4)/2;
%         mix = awgn(mix,0);
        hits(:,i) = findHits(mix,samples);
        hit1(:,i) = hits(1,i);
        hit2(:,i) = hits(2,i);
    end

    mixTest1 = MPMsignal(A(:,1),R0,d,theta, samples, chirps, V, fc, bw,tm,sweep_slope,fs);%theoreticalMix(samples, chirps, R0, V, fc, bw, tm, sweep_slope, fs);
    mixTest2 = MPMsignal(A(:,2),R2,d,theta, samples, chirps, V, fc, bw,tm,sweep_slope,fs);
%     mixTest2 = theoreticalMix(samples, chirps, R2, V, fc, bw, tm, sweep_slope, fs);
    mixTest = (mixTest1+mixTest2)./2;
%     mix = [mix1 mix2]

%     mix = A(:,1)*downsample(mix1(:,1)',25);
    mix = A(:,1)*mix1(:,1)'+A(:,2)*mix2(:,1)';
    
    mix = awgn(mixTest,20); % TODO does this suggest that I did not need to factor in the time delay in the range variable when producing each signal
    
%     fftshift(fft2(mix),2);
%     Y = zeros(N-p,p+1);
%     Y1 = Y(:,1:end-1);
%     Y2 = Y(:,2:end);
    Hcol = p;
    Hrow = N-p+1;
%     Y = hankel(mix(1:end-p),mix(end-p:end));
%     Y = hankel(hit1(1:end-p),hit1(end-p:end));

%     hit2=hit1;
%     hit1=hit2;
    for i=1:N-p+1
        Y(i,:) = hit1(i:p+i-1);
    end

    for i=1:N-p+1
        YY(i,:) = hit1(i:p+i-1);
    end
    
    Y = hankel(hit1(1:end-p),hit1(end-p:end));
    Y1 = Y(:,1:end-1);
    Y2 = Y(:,2:end);
    l = eig(pinv(Y1)*Y2);
    Z = zeros(length(hit1),p);
    for i=1:length(l)
        Z(:,i) = transpose(l(i).^(0:length(hit1)-1));
    end
    rZ = real(Z);
    iZ = imag(Z);

    rZ(isinf(rZ))=realmax*sign(rZ(isinf(rZ)));
    iZ(isinf(iZ))=realmax*sign(iZ(isinf(iZ)));

    Z = rZ+1i*iZ;
    h=Z\hit1';
    Amp = abs(h);
    theta = atan2d(imag(max(h)),real(max(h)))
    theta = asind(imag(log(max(h)))/(2*pi*d/lambda))
%     l = max(l);
    
    phi = asin(atan2(imag(l),real(l))*lambda/(2*pi*d))
    abs(phi)
    rad2deg(abs(phi))
    MPMangle = (rad2deg(abs(phi)))

    [U,S,Vv] = svd(Y);
    U0 = U(1:end-1,:);
    U1 = U(2:end,:);
    Y1 = U0*S*Vv';
    Y2 = U1*S*Vv';
%     Y = hankel(mix(1:p),mix(1:N-p+1));

%     Y = zeros(p,N-p+1);
%     Y1 = Y(1:end-1,:);
%     Y2 = Y(2:end,:);
%     Y1 = Y(:,1:end-1);
%     Y2 = Y(:,2:end);
%     Y1 = Y(:,1:end-1);
%     Y2 = Y(:,2:end);
%     Y1 = Y(1,:);
%     Y2 = Y(N-p+1,:);
%     
% 
%     l_eig = eig(Y1\Y2);
    Yt = (U0'*U0/(U1'*U0));%/eye(Hrow);
%     Yt = Y2*pinv(Y1);
    [v, l_eig] = eig(Yt);
%     [l_eig,ind] = sort(diag(l_eig));
%     Ds = l_eig(ind,ind);

% svd
    l_eig = mean(l_eig');%l_eig(l_eig>0);
%     [l_eig,ind] = sort(diag(l_eig));

    %% This works
    phi = asin(atan2(imag((l_eig)),real((l_eig)))*lambda/(2*pi*d))
    abs(phi)
    rad2deg(abs(phi))
    MPMangle = max(rad2deg(abs(phi)))

    %%
    mix1 = theoreticalMix(samples, chirps, R0, V, fc, bw, tm, sweep_slope, fs);
    mix2 = theoreticalMix(samples, chirps, R2, V2, fc, bw, tm, sweep_slope, fs);
    
%     mix = (mix2+mix1)/2;
    
    ogSignal = theoreticalMix(samples, chirps, 0, V, fc, bw, tm, sweep_slope, fs);
    ogSignal = ogSignal;
    hk = ogSignal\mix';
    thetaProny = atan(imag(hk)./real(hk));
    Prony = rad2deg(thetaProny);

    fk = sin(theta)./lambda;
%     mix2 = A*mix2;

%     mix = (mix2+mix1)/2;

%     mix = awgn(mix,0);

end