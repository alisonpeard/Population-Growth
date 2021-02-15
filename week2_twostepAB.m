%the AB2 method with initialization with RK2 method.
%y'+5y = 0

%Where, the initial conditions are,

%y(0) = 0.5
close all;

clc;
h = 0.02;     % Step size
Tmax = 0.5;    % Maximum time
N = Tmax / h;  % Maximum number of steps

t = linspace(0,0.5,N+1);  % Time range

% Analytical solution of the differential equation
y_real = 0.5*exp(-5*t);
plot(t,y_real);
hold on

%Numerical solution
f=@(t,y)  -5*y; 
% Initial Conditions
Y=zeros(N,1);
Y(1) = 0.5;
% Initialization with Runge-Kutta method 
%(modified euler which is 2nd oder accurate)
    Y(2) = Y(1) + h*f(t(1)+0.5*h,Y(1)+h*f(t(1),Y(1)));

% Second Order Adams-Bashforth method steps
for n=1:N
    Y(n+1) = Y(n) + 3/2*h*f(t(n),Y(n)) - h/2*f(t(n),Y(n));
end

plot(t,Y,'o:')
legend('Exact Solution','Adams-Bashforth Solution','Location','NorthEast')
title('When h = 0.02')
xlabel('t')
ylabel('y')




%%
%Let's consider the following differential equation:

%y'' + 10 y' + 500y = 0

%Where, the initial conditions are,

%y(0) = -0.025, and y'(0) = -1
close all;

clc;
u0=0.5;
h=0.1;
%t=(0:h:10)';
%N=length(t);
k=1;
  
Tmax = 10;    % Maximum time
N = Tmax / h;  % Maximum number of steps
t = linspace(0,10,N+1);  % Time range

% Analytical solution of the differential equation
%y_real = -(((9.*sqrt(19))/760).*exp(-5.*t).*sin(5.*sqrt(19).*t)) - ((1/40).*exp(-5.*t).*cos(5.*sqrt(19).*t));
%plot(t,y_real);
%hold on

%Numerical solution
f=@(t,y) [y(2); (y(2)/k)*(1-y(2)-y(1))]; % Governing system of equations

% Initial Conditions
alpha=1;
Y=zeros(2,N);
Y(:,1) = [0; u0];
% Initialization with second order Runge-Kutta method
k1 = h.*f(t(1),Y(:,1));
    k2 = h.*f(t(1)+alpha.*h, Y(:,1)+alpha.*k1);
    Y(:,2) = Y(:,1) + (1-1/2/alpha).*k1 + k2/2/alpha;

% Second Order Adams-Bashforth method steps
for i=2:N
    Y(:,i+1) = Y(:,i) + 3/2*h*f(t(i),Y(:,i))- h/2*f(t(i-1),Y(:,i-1));
end

plot(t,Y(2,:),'o:')
legend('Adams-Bashforth Solution','Location','NorthEast')
title('When k=1, u_0=0.5')
xlabel('t')
ylabel('u')
