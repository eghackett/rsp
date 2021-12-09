function y = MPMv3(d, theta, N, s)

    if nargin < 1
        c = 3e8;
        lambda = (c/77e9);
        d = lambda/1.99;
        N = 20;
        theta = [15 63];

    end

    % samples/3 < P < samples/2
    p = 10; 

    configuration; %load configuration profile

%     samples = 20;

    theta = deg2rad(theta);
    n = [0:N-1]';
    A = exp(1j.*n*2*pi*d.*sin(theta)/lambda); %Steering vector; no of rows should equal N and col = no of theta
    R0d = R0+d*sin(theta(1))*n;
    R2d = R2+d*sin(theta(2))*n;
    
    mix11 = zeros(samples, chirps);
    mix22 = zeros(samples, chirps);
    for i=1:N
        mix11 = theoreticalMix(samples, chirps, R0d(i,1), V, fc, bw, tm, sweep_slope, fs);
        mix22 = theoreticalMix(samples, chirps, R2d(i,1), V2, fc, bw, tm, sweep_slope, fs);
        mix1(:,i) = mix11(:,1);
        mix2(:,i) = mix22(:,1);
    end

    r = zeros(2, N);
    c = zeros(2, N);
    hit1 = zeros(1, N);
    hit2 = zeros(1, N);
    hits = zeros(2,N);
    for i=1:N
        mix3 = theoreticalMix(samples, chirps, R0d(i,1), V, fc, bw, tm, sweep_slope, fs);
        mix4 = theoreticalMix(samples, chirps, R2d(i,1), V2, fc, bw, tm, sweep_slope, fs);
        mix = (mix3+mix4)/2;
%         mix = A*mix;
        mymix(:,i)=mix(i,:)';
        mymix2(:,i)=mix(:,i);
%         mix = awgn(mix,0);
        [hits(:,i), r(:,i),c(:,i)] = findHits(mix,samples);
        hit1(:,i) = A(i,1)*hits(1,i);
        hit2(:,i) = A(i,2)*hits(2,i);
    end

%     hit1 = A(:,1)*hit1;
%     hit1 = A(:,2)*hit2;


    x = hit1(1,:);
%     x = mean(x);
    myeig = exp(1j*(2*pi.*d*sin(deg2rad(15)))/lambda);
    phi = asin(atan2(imag(myeig),real(myeig)).*lambda./(2*pi.*d));
    MPMangle = (rad2deg((phi)));

    [a,b,c] = mpmfun(x,p,d,lambda);
    [a2,b2,c2] = mpmfun(hit2,p,d,lambda);

%     target2 = 90-max(c(1:p/2))
    target1 = c%max(c(1+p/2:end))
    target2 = c2%max(c2(1+p/2:end))

%     target1Error = abs(15-target1)/15
%     target2Error = abs(63-target2)/63



    %%
    x = hit1;%mymix2(1,:);
    pronyTheta = pronyfun(x,p,ts)
%     target1_prony = pronyTheta(1)
%     target2_prony = pronyTheta(end)

%     hk = x\mix';
%     thetaProny = atan(imag(hk)./real(hk));
%     Prony = rad2deg(thetaProny);

%     fk = sin(theta)./lambda;

% x = sin(2*pi*f1*t) + sin(2*pi*f2*t) + sin(2*pi*f3*t);
mixFFT = fftshift(fft2(mix),2);
    rdb = 20*log10(abs(mixFFT(1:ceil(samples/2), :)));
    x = mixFFT(1:ceil(samples/2), :);
%Generate random samples
M = 20; %Approximately 10% of total signal
k = randperm(N);
m = k(1:M);
b = x(sort(m));

end