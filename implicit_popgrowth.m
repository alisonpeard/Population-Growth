%SC Week 2- Implicit euler with trapezium rule to solve the ODE
function [time,pop] = implicit_popgrowth(h, tf, u0, T0, p)

%set up step size, initial conditions, parameter k, length of simulation

k = p;
N = tf/h; 
t = 0:h:tf;

%pre-allocate arrays
u = zeros(size(t)); %population size
T = zeros(size(t)); %composite trapezium rule integral approx

%assign initial values
u(1) = u0;
T(1) = T0;

% u' = (u\k)(1-u-int_0_t(u(s))ds) = f(t,u)
%Implicit Euler method:
%u_n+1 = u_n + h(u_n+1/k)f(t_n+1, u_n+1)
%T_n+1 = T_n + (h/2)(u_n + u_n+1)
%Below, we have solved a quadratic for u_n+1 explicitly

for j = 1:N
   
    a = 1 + (h/2);
    b = - (1 - T(j) - (h/2)*u(j) - (k/h));
    c = -(k/h)*u(j);
   
    %implicit Euler
    u(j+1) = (- b + sqrt(b^2 - 4*a*c))/(2*a);
    T(j+1) = T(j) + (h/2)*(u(j) + u(j+1));
   
end

time = t;
pop = u;

end


