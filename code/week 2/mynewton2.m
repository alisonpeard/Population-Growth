% fixing bugs

function x0 = mynewton2(x0, delta_t, kappa, theta)
% x0 a column vector
% theta = 1 => fully implicit
% [counter,x0] = mynewton([1 0])

    k = kappa;
    h = delta_t;
    
    Un = x0(1);
    Yn = x0(2);
    
    f = @(x,y) [x - Un - h*((1-theta)*Un*(1-Un-Yn)+theta*x*(1-x-y))/k; y - Yn - h*((1-theta)*Un + theta*x)];
    
    fdx = @(x,y) [1 - h*theta*(1-2*x-y)/k; -h*theta*x];
    fdy = @(x,y) [h*theta*x/k; 1];   
    J = @(x,y) [fdx(x,y) fdy(x,y)];                 

    x = x0(1); y = x0(2);                      

    TOL = 1e-12; counter = 0;
    
    while ( norm(f(x,y)) > TOL && counter < 500  )   
        x = x0(1); y = x0(2);                       
        x0 = x0 - J(x,y)\f(x,y);                  
        counter = counter + 1;
    end
     

end