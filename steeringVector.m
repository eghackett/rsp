function A = steeringVector(theta, N)
    theta = deg2rad(theta);
    n = [0:N-1]';
    A = exp(1j.*n*2*pi*d.*sin(theta)/lambda); %Steering vector; no of rows should equal N and col = no of theta

end