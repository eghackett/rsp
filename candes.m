%Emmanuel Candes, California Institute of Technology, June 6 2007, IMA Summerschool.
%Transcribed by Jon Dattorro.
%Fails using SDP solver SDPT3  on 7th consecutive run after Matlab R2007b startup.  CVX version 1.2 (build 656).
%Fails using SDP solver Sedumi on 4th consecutive run after Matlab R2007b startup.  CVX version 1.2 (build 656).
clear all, close all                
n = 512;                            % Size of signal
m = 64;                             % Number of samples (undersample by a factor 8)

k = 0:n-1;  t = 0:n-1;
F = exp(-i*2*pi*k'*t/n)/sqrt(n);    % Fourier matrix
freq = randsample(n,m);
A = [real(F(freq,:)); 
     imag(F(freq,:))];              % Incomplete Fourier matrix

S = 28;
support = randsample(n,S);
x0 = zeros(n,1); 
x0(support) = randn(S,1);
b = A*x0;

% Solve l1 using CVX
cvx_quiet(true);
%cvx_solver('sedumi'); 
cvx_begin
    variable x(n);
    minimize(norm(x,1));
    A*x == b;
cvx_end

norm(x - x0)/norm(x0)
figure, plot(1:n,x0,'b*',1:n,x,'ro'), legend('original','decoded')