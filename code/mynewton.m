
function x0 = mynewton(x0, delta_t, kappa)
% x0 a column vector
% [counter,x0] = mynewton([1 0])

    theta = 0.5;
    k = kappa;
    h = delta_t;
    
    Un = x0(1);
    Yn = x0(2);
    
    f = @(x,y) [x - Un - h*((1-theta)*Un*(1-Un-Yn)+theta*(1-x-y))/k; y - Yn - h*((1-theta)*Un + theta*x)];
    
    fdx = @(x,y) [1 - h*theta*(1-2*x-y)/k; -h*theta];       % define d/dx column of J
    fdy = @(x,y) [h*theta*x/k; 1];    % define d/dy column of J
    J = @(x,y) [fdx(x,y) fdy(x,y)];                  % define J

    x = x0(1); y = x0(2);                            % define (x,y) to reduce # of computations

    TOL = 1e-12; counter = 0;
    
    % solve for successive roots using xk+1 = xk - J\f
    while ( norm(f(x,y)) > TOL && counter < 500  )    % I think there should also be a condition on J but don't know what
        x = x0(1); y = x0(2);                       
        x0 = x0 - J(x,y)\f(x,y);                     % Newton-Rhapson iteration
        counter = counter + 1;
    end
    
    %x0 = round(complex(x0(1),x0(2)),4);                             % imaginary format

end