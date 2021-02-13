clear all;
u0 = 0.1;
y0 = 0;
k = 0.2;
h0 = 0.1;
tol = 10^-5;
Tmax = 10;
c = -1/12;
chat = 5/12;

theta = 1/2;

t = h0;
h = h0;
v0 = [u0;y0];
v1 = theta_method(v0,h0,k,theta)';
Uvec(1) = v0(1);
Uvec(2) = v1(1);
Tvec(1) = 0;
Tvec(2) = h0;
Hvec(1) = h0;
func = @(v)rhs(v,k);


while t<Tmax
    E = 10*h*tol;
    while E > h*tol
        % Theta
        v2 = theta_method(v1,h,k,theta)';
        % AB
        v2hat = ab_method(v0,v1,h,func);
        v0 = v1;
        v1 = v2;
        E = abs(c/(chat-c))*norm(v2-v2hat);
        
        if E > h*tol
            h = h/2;
            v0 = (v0+v1)/2;
        end  
    end
    t = t+h;
    Tvec(end+1) = t;
    Hvec(end+1) = h;
    Uvec(end+1) = v1(1);
    
        if E <= h*tol/10
            h = h*2;
            v0 = (2*v0)-v1;
        end
end
    
