function theta = pronyfun(x,p,Ts)
    N = length(x);
    T = toeplitz(x(p:N-1),x(p:-1:1));
    a = -T\x(p+1:N)';
    indeterminate_form = sum(isnan(a) | isinf(a));
    if (indeterminate_form)
        Amp = []; alfa = []; freq = []; theta = [];
        return;
    end
    c = transpose([1;a]);
    r = roots(c);
    alfa = log(abs(r))/Ts;
    freq = atan2(imag(r),real(r))/(2*pi*Ts);
    alfa(isinf(alfa))=realmax*sign(alfa(isinf(alfa)));

    len_vandermonde=p;
    Z = zeros(len_vandermonde,p);
    for i=1:length(r)
        Z(:,i)=transpose(r(i).^(0:len_vandermonde-1));
    end
    rZ = real(Z);
    iZ = imag(Z);
    rZ(isinf(rZ))=realmax*sign(rZ(isinf(rZ)));
    iZ(isinf(iZ))=realmax*sign(iZ(isinf(iZ)));
    Z = rZ+1i*iZ;

    h = Z\x(1:len_vandermonde)';
%     indeterminate_form = sum(sum(isnan(Z) | isinf(Z)));
%     if (indeterminate_form)
%         Amp = []; alfa = []; freq = []; theta = [];
%         return;
%     else
%         
%     end
    Amp = abs(h);
    theta = atan2d(abs(imag(h)),abs(real(h)));

end