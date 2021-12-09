%lpNorm.m 
x = -1.25:.01:1.25;
m = -.66; b = 1;
yl = m*x+b;
p = [.4 .9 1 1.5 2];
figure;
for ii = 1:length(p)
 [Np ind] = min((abs(yl).^p(ii)+abs(x).^p(ii)).^(1/p(ii)));
 xp = linspace(-Np,Np,1000);
 yp = (Np^p(ii)-abs(xp).^p(ii)).^(1/p(ii));
 subplot(1,length(p),ii);
 plot(xp,yp,'-k',xp,-yp,'-k');hold on;
 plot(x,yl,'-r'); xlim([-1.25 1.25]);ylim([-1.5 1.5]);axis square;
 title(sprintf('P = %.1f',p(ii)));xlabel('x_1');ylabel('x_2');
end