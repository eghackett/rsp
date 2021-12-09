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
        mymix(:,i)=mix(i,:)';
        mymix2(:,i)=mix(:,i);
%         mix = awgn(mix,0);
        hits(:,i) = findHits(mix,samples);
        hit1(:,i) = hits(1,i);
        hit2(:,i) = hits(2,i);
    end



    x = mymix(1,:);

    [a,b,c] = mpmfun(x,p,d,lambda);

    target2 = 90-max(c(1:p/2))
    target1 = max(c(1+p/2:end))

    target1Error = abs(15-target1)/15
    target2Error = abs(63-target2)/63



    %%
    x = mymix2(1,:);
    pronyTheta = pronyfun(x,p,ts)
    target1_prony = pronyTheta(1)
    target2_prony = 90-pronyTheta(end)

    hk = ogSignal\mix';
    thetaProny = atan(imag(hk)./real(hk));
    Prony = rad2deg(thetaProny);

    fk = sin(theta)./lambda;
%     mix2 = A*mix2;

%     mix = (mix2+mix1)/2;

%     mix = awgn(mix,0);

end