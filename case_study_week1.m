%explicit euler scheme with trapezium rule
clc,clear,close all


% define timestep
h=0.001;

%define parameters k u0
k=1; %1 0.1 0.01
u0=1.5; %0.1, 1, 2



% set up discrete timesteps
tsteps=(0:h:10)';
N=length(tsteps);

% initialise arrays to store numerical solutions and apply initial
% condition
uexp=zeros(N,1); uexp(1)=u0; trap(1)=0;

for n=1:N-1
    uexp(n+1)=uexp(n)+(h*uexp(n)/k)*(1-uexp(n)-trap(n));
    trap(n+1)=trap(n)+(h/2)*(uexp(n)+uexp(n+1));
end

figure
plot(tsteps,uexp,'r')
xlabel('t')
ylabel('U_t')

hold on

%%
u0=1;
uimp=zeros(N,1); uimp(1)=u0; trap(1)=0;
for n=1:N-1
    uimp(n+1)=qua((h/k)+((h^2)/(2*k)),1-(h/k)+(h*trap(n)/k)+((h^2)*uimp(n)/(2*k)),-1*uimp(n));
    trap(n+1)=trap(n)+(h/2)*(uimp(n)+uimp(n+1));
end

plot(tsteps,uimp,'b')
xlabel('t')
ylabel('U_t')


%%
%use newton's method
u0=2;
uimp2=zeros(N,1); trap=uimp2;
uimp2(1)=u0;
trap(1)=0;

for n=1:N-1
    [uimp2(n+1),trap(n+1)]=get_new(uimp2(n),trap(n),h,k);
end

plot(tsteps,uimp2,'k')
xlabel('t')
ylabel('U_t')
legend('Explicit Euler, u_0=0.1','Implicit Euler, u_0=1','Implicit2 Euler, u_0=2')


%%
function x= qua(a,b,c)
x = (-b + sqrt(b^2 - 4 * a * c))/(2*a);
%x = (-b - sqrt(b^2 - 4 * a * c))/(2*a);
end

function [u,v]=get_new(un,vn,h,k)

TOL=1e-10;
u=un; v=vn;
f=[u-un-(((h*u)*(1-u-v))/k); v-vn-((h*(un+u))/2)];
while norm(f)>TOL
    J=[1+(h/k)*v+((2*h*u)/k) (h*u)/k;(-1*h)/2  1];
    uk=[u;v]-J\f;
    u=uk(1); v=uk(2);
    f=[u-un-(((h*u)*(1-u-v))/k); v-vn-((h*(un+u))/2)];
end
    
end
