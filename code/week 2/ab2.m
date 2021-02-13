% adam-bashforth method
% x2 = x1 + (1/2)*(3f(x1) - f(x0))

% Crank-Nicolson (order 2)
% vec = cn2(x0,h,k)

function x1 = ab2(x1,x0, h, k)
% f(y,u), h = timestep, k=kappa

    Y0 = x0(1); U0 = x0(2);
    Y1 = x1(1); U1 = x1(2);

    x1(1) = Y1 + (h/2)*(3*U1 - U0);
    x1(2) = U1 + h*(3*U1*(1-U1-Y1)-U0*(1-U0-Y0))/(2*k);

end
