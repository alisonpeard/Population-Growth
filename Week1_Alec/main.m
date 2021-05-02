clear all;
gcf;
hold on
theta = 0.5;

Tmax = 10;
h = 1e-2;
N = ceil(Tmax/h);
time = 0:h:Tmax;


uvec = zeros(N,1);
Tvec = zeros(N,1);


kvec = [0.5, 1, 1.5];
for j = 1:length(kvec)
    k = kvec(j);
    u0 = 1.5;
    T0 = 0;
    uvec(1) = u0;
    Tvec(1) = T0;
for i = 2:N+1
   u1 = theta*exp_u(u0,T0,h,k) + (1-theta)*imp_u(u0,T0,h,k);
   T1 = T0 + ((h/2)*(u1+u0));
    
   u0 = u1;
   T0 = T1;
   uvec(i) = u0;
   Tvec(i) = T0;
end
    plot(time,uvec)
    hold all
end