function [theta1, theta2, MPMangle] = mpmfun(x,p,d,lambda)

    N = length(x);
    Y = hankel(x(1:end-p+1),x(end-p+1:end));
%     Y=Y';
%     Y1 = Y(:,1:end-1);
%     Y2 = Y(:,2:end);
    Y1 = Y(1:end-1,:); %delete last row
    Y2 = Y(2:end,:); %delete first row
    mon = pinv(Y1)*Y2;
    [a,b,c] = svd(mon);
    l = eig(mon);
    Z = zeros(N,p);
    for i=1:length(l)
        Z(:,i) = transpose(l(i).^(0:length(x)-1));
    end
    rZ = real(Z);
    iZ = imag(Z);

    rZ(isinf(rZ))=realmax*sign(rZ(isinf(rZ)));
    iZ(isinf(iZ))=realmax*sign(iZ(isinf(iZ)));

%     Z = rZ+1i*iZ;
    h=Z\x';
    h = sum(h);
    Amp = abs(h);
    theta1 = atan2(imag((h)),real((h)));
    theta2 = asin(imag(log((h)))/(2*pi*d/lambda));
    theta1 = rad2deg(theta1);
    theta2 = rad2deg(theta2);
%     l = max(l);
    oldL = l;
    l = mean(l);
%     [a,b,c] = svd(l)
    phi = asin(atan2(imag(l),real(l))*lambda/(2*pi*d));
    abs(phi);
    rad2deg(abs(phi));
    MPMangle = (rad2deg(abs(phi)));

%     myeig = exp(1j*(2*pi.*d*sin(deg2rad(13)))/lambda);
%     phi = asin(atan2(imag(myeig),real(myeig)).*lambda./(2*pi.*d));
%     MPMangle = (rad2deg((phi)))

end