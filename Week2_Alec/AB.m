clear all;

h = 1e-4;
Tmax = 5;
N = ceil(Tmax/h);
time = 0:h:Tmax;

k = 0.2;
func = @(v)rhs(v,k);
theta = 0.5;

% ICs
u0 = 0.8;
y0 = 0;

% midpoint method for first step
u1 = theta*exp_u(u0,y0,h,k) + (1-theta)*imp_u(u0,y0,h,k);
y1 = theta*exp_y(u0,y0,h) + (1-theta)*imp_y(u1,y0,h);

% AB Scheme
v0 = [u0;y0];
v1 = [u1;y1];

uvec = zeros(1,N+1);
uvec(1) = v0(1);
uvec(2) = v1(1);

for i = 3:N+1
    
    v2 = ((h/2)*((3*func(v1))-func(v0)))+v1;
    uvec(i) = v2(1);
    
    v0 = v1;
    v1 = v2;
    
end

plot(time,uvec)