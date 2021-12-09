

alfa = log(abs(l_eig))/ts;
    freq = atan2(imag(l_eig),real(l_eig))/(2*pi*ts);

    Z = zeros(N,p);
    for i=1:length(l_eig)
        Z(:,i) = transpose(l_eig(i).^(0:N-1));
    end

    rZ = real(Z);
    iZ = imag(Z);

    rZ(isinf(rZ))=realmax*sign(rZ(isinf(rZ)));
    iZ(isinf(iZ))=realmax*sign(iZ(isinf(iZ)));

    Z = rZ+1i*iZ;
    h = Z\mix;
    Amp = abs(h);
    theta = atand(imag(h)/real(h));

    a = lambda/(2*pi*d);
    asind(a.*theta);

    theta2 = asind(imag(log(diag(Z)))/(2*pi*d));
    theta = asind((lambda/(2*pi*d))*(iZ./rZ));

    mean(abs(theta2))

    sum(theta)

    phi = acosd(lambda*log(l_eig)/(1j*2*pi*d))
    abs(phi)