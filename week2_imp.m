%implicit
% set up discrete timesteps
%explicit euler scheme with trapezium rule
clc,clear,close all
% define timestep
u0=0.1;
h=0.001;
tsteps=(0:h:10)';
N=length(tsteps);
kvec=[0.1 1 2];
for i=1:length(kvec)
k=kvec(i);
uimp=zeros(N,1); 
uimp(1)=u0;
y(1)=0;

for n=1:N-1
    [y(n+1),uimp(n+1)]=get_new(y(n),uimp(n),h,k);
end

 plot(tsteps,uimp);
 hold on;
end
title('Implicit solution with composite trapezoidal, U_0= '+ string(u0))
xlabel('t')
ylabel('U_t')
legendStrings= "k=" +string(kvec);
legend(legendStrings)

function [u,v]=get_new(un,vn,h,k)

TOL=1e-10;
u=un; v=vn;
f=[u-(h*v)-un; v-vn-(h/k)*v*(1-v-u)];
while norm(f)>TOL
    J=[1 -h;(v*h)/k  1-(h/k)+(h*u)/k+(2*h*v)/k];
    uk=[u;v]-J\f;
    u=uk(1); v=uk(2);
    f=[u-(h*v)-un; v-vn-(h/k)*v*(1-v-u)];
end
    
end


