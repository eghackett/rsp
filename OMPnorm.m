function [A]=OMPnorm(D,X,L,eps)
%OMPNORM Sparse coding of a group of signals based on a 
% dictionary and specified number of atoms 
%
% [USAGE]
% A = OMPnorm(D,X,L,DEBUG)
%
% [INPUTS]
% D M x K overcomplete dictionary
% X Data vector or matrix. Must be of size M x P, where P is 
% number of signals to code 
% L Maximum number of coefficients for each signal
% eps Stopping criterion (stop when l2 norm of residual is less 
% than eps)
%
% [OUTPUTS]
% A Sparse coefficient matrix of size K x P
%
P=size(X,2);
[M K]=size(D);
%Ensure that the dictionary elements have been normalized
Dnorm = D./repmat(sqrt(sum(D.^2)),[M 1]);
A = zeros(K,P);
for k=1:1:P,
 x = X(:,k);
 residual = x;
 indx = zeros(L,1);
 resids = zeros(K+1,1);
 for j = 1:1:L
 proj = Dnorm'*residual;
 [maxVal,pos] = max(abs(proj));
 pos = pos(1);
 indx(j) = pos;
 a = pinv(D(:,indx(1:j)))*x;
 residual = x-D(:,indx(1:j))*a;
 resids(j+1) = sum(residual.^2);
 %Break if error has reach sufficient value
 if sum(residual.^2) < eps
 break;
 end
 end;
 A(indx(1:j),k)=a;
end