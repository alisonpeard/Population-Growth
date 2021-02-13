% adam-bashforth method

% Crank-Nicolson (order 2)
% vec = cn2(x0,h,k)

function x1 = ab2(x1,x0, h, k)
% fix order of y,u to match
    
    
    func = @(v)[v(1);(v(2)/k)*(1-v(1)-v(2))];
        
    x1 = x0 + (h/2)*((3*func(x1))-(func(x0)));
    

    
  
  

end