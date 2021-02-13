% Newton solver for theta method
% f = f(y,u), h = timestep, k = kappa
% theta = 1 => fully implicit
% test: x0 = theta_method([0; 0.5], 0.1, 1, 1/2)

function x0 = theta_method(x0, h, k, theta)

    Yn = x0(1);
    Un = x0(2);
    
    f = @(y,u) [y - Yn - h*((1-theta)*Un + theta*u); u - Un - h*( (1-theta)*Un*(1-Un-Yn) + theta*u*(1-u-y) )/k];
    
	fdy = @(y,u) [1; h*theta*u/k];   
    fdu = @(y,u) [-h*theta*u; 1 - h*theta*(1-2*u-y)/k];
    
    J = @(y,u) [fdy(y,u) fdu(y,u)];                 

    y = x0(1); u = x0(2);                      

    TOL = 1e-12; counter = 0;
    
    while ( norm(f(y,u)) > TOL && counter < 500  )   
        y = x0(1); u = x0(2);                       
        x0 = x0 - J(y,u)\f(y,u);                  
        counter = counter + 1;
    end
 
end