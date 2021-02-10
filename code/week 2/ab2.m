% adam-bashforth method

% Crank-Nicolson (order 2)
% vec = cn2(x0,h,k)

function x1 = ab2(x0,x1, delta_t, kappa)
% fix order of y,u to match

    k = kappa;
    h = delta_t;
    
    
    Y0 = x0(1); U0 = x0(2);
    Y1 = x1(1); U1 = x1(2);

    f = @(y,u) [y - Y1 - (h/2)*(3*U1 - U0); u - U0 - (h/2)*( (3*U1/k)*(1-U1-Y1) - (U0/k)*(1-U0-Y0))];
    
    fdy = @(y,u) [1; 0];   
    fdu = @(y,u) [0; 1];

    J = @(y,u) [fdy(y,u) fdu(y,u)]; 


    y = x1(1); u = x1(2);                      

    TOL = 1e-12; counter = 0;
    
    while ( norm(f(y,u)) > TOL && counter < 500  )   
        y = x1(1); u = x1(2);                       
        x1 = x1 - J(y,u)\f(y,u);                  
        counter = counter + 1;
    end
  

end