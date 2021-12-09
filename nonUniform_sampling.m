%nonUniform_Sampling.m - Non-Uniform Sampling
clear; clc; close all
f1 = 1e3; f2 = 2e3; f3 = 4e3; %Hz
fs = 10e3; %Hz
T = 1/fs; %Period in seconds
%Number of samples
N = 1028; %For .1s of data
t = [0:N-1]'*T;
%Full-length and sampled signals
x = sin(2*pi*f1*t) + sin(2*pi*f2*t) + sin(2*pi*f3*t);
%Generate random samples
M = 151; %Approximately 10% of total signal
k = randperm(N);
m = k(1:M);
b = x(sort(m));
%Plot time-domain signal
figure;subplot(3,2,1);plot(t,x,'r-',t(m),b,'bo');
xlabel('Time (s)');ylabel('x(t)');title('Sum of sinusoids');
%Plot (sparse) Fourier domain signal
X = fft(x);
subplot(3,2,2);plot(linspace(0,fs,length(X)),abs(X));
xlabel('Frequency (Hz)');ylabel('|X(f)|');
%Generate linear transform matrix
PHI = fft(eye(N));
S = zeros(M,N);
S(sub2ind(size(S),1:M,sort(m))) = 1;
A = S*(1/N)*PHI';

%Naive l2 minimization solution
l2 = pinv(A)*b;
subplot(3,2,4);plot(linspace(0,fs,length(X)),abs(l2));
xlabel('Frequency (Hz)');ylabel('|X_l_2(f)|');
%Reconstruct signal
xl2 = (1/N)*real(PHI'*l2);
subplot(3,2,3);plot(t,xl2);
xlabel('Time (s)');ylabel('x_l_2(t)');
%Approximate l1 solution using OMP
l1 = OMPnorm(A,b,floor(M/4),0);
subplot(3,2,6);plot(linspace(0,fs,length(X)),abs(l1));
xlabel('Frequency (Hz)');ylabel('|X_l_1(f)|');