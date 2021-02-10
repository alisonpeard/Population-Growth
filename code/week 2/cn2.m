% Crank-Nicolson (order 2)
% vec = cn2(x0,h,k)

function x0 = cn2(x0, delta_t, kappa)
% fix order of y,u to match

    k = kappa;
    h = delta_t;
    
    Y0 = x0(1); % mixed these up a little
    U0 = x0(2);
    
    f = @(y,u) [y - Y0 - (h/2)*(u+U0); u - U0 - (h/2)*( (u/k)*(1-u-y) + (U0/k)*(1-U0-Y0))];
    
    fdy = @(y,u) [1; (h/2)*(u/k)];   
    fdu = @(y,u) [-h/2; 1 - (h/2)*(1 - 2*u - y)/k];

    J = @(y,u) [fdy(y,u) fdu(y,u)];                 

    y = x0(1); u = x0(2);                      

    TOL = 1e-12; counter = 0;
    
    while ( norm(f(y,u)) > TOL && counter < 500  )   
        y = x0(1); u = x0(2);                       
        x0 = x0 - J(y,u)\f(y,u);                  
        counter = counter + 1;
    end
  

end