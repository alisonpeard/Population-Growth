clc;clear;
%%initialise parameters
c=-1/12; %error constant for trap rule
c2=5/12; %error constant of AB2
u(1)=0.5;
y(1)=0;
k=1;
h=0.01;
TOL=10^(-5);
%TOL=1e-7;
Tmax=10;

%%
tsteps(1)=0; % tsteps stores the values of the timesteps
ustore(1)=u(1);
n=1;
while tsteps(n)< Tmax
    E=2*h*TOL;
    while E>h*TOL
         
         [y(n+1),u(n+1)]=trap(y(n),u(n),h,k);
         if n~=1
         [ytil(n+1),util(n+1)]=ab(y(n-1),u(n-1),y(n),u(n),h,k,2);
         else
             [ytil(n+1),util(n+1)]=ab(0,0,y(n),u(n),h,k,1);
         end
         
         E=abs(c/(c+c2))*norm([y(n+1),u(n+1)]-[ytil(n+1),util(n+1)]);
         
         if E>h*TOL
             h=h/2;
             %remesh
             if n~=1
                 yhat(n-1)=0.5*(y(n-1)+y(n));
                 uhat(n-1)=0.5*(u(n-1)+u(n));
                 [y(n+1),u(n+1)]=trap(y(n),u(n),h,k);
                
                 [ytil(n+1),util(n+1)]=ab(yhat(n-1),uhat(n-1),y(n),u(n),h,k,2);
                 
             else
                 [y(n+1),u(n+1)]=trap(y(n),u(n),h,k);
                 [ytil(n+1),util(n+1)]=ab(0,0,y(n),u(n),h,k,1);
             end
         end
    end
        tsteps(n+1)=tsteps(n)+h;
        %ustore(n+1)=u(n+1);
        
        if E<= (1/10)*h*TOL
            h=2*h;
            %remesh
            if n~=1
                 yhat(n-1)=2*y(n-1)-y(n);
                 uhat(n-1)=2*u(n-1)-u(n);
                  [y(n+1),u(n+1)]=trap(y(n),u(n),h,k);
                
                  [ytil(n+1),util(n+1)]=ab(yhat(n-1),uhat(n-1),y(n),u(n),h,k,n);
                 
             else
                 [y(n+1),u(n+1)]=trap(y(n),u(n),h,k);
                 [ytil(n+1),util(n+1)]=ab(0,0,y(n),u(n),h,k,n);
             end
        end
        
        ustore(n+1)=u(n+1);
        n=n+1;
end

%%
%plot
figure(1)
plot(tsteps,ustore,'r')
%axis([0,5,-1,1])
ylabel('Solution')
figure(2)

semilogy(tsteps(1:end-1),tsteps(2:end)-tsteps(1:end-1))
ylabel('Timestep length')

%%   
function [y2,u2]=trap(y,u,h,k)

    
    f=[-h*u; -(h/k)*u*(1-u-y)];
    J=[1, -h;u*h/k, 1-(h/k)+(h/k)*y+2*(h/k)*u];
    
    uk=[y;u]-J\f;
    y2=0.5*uk(1)+(1-0.5)*(y+h*u); 
    u2=0.5*uk(2)+(1-0.5)*(u+h*(u/k)*(1-u-y));
    
    

end
function [ytil2,util2]=ab(y0,u0,y1,u1,h,k,n)
f=@(y) [y(2); (y(2)/k)*(1-y(2)-y(1))]; % Governing system of equations

% Initial Conditions
alpha=1;
Y(:,1)=[y0,u0];
Y(:,2)=[y1,u1];
%Y(:,1) = [0; u0];
  if n==1
    % Initialization with second order Runge-Kutta method
    k1 = h.*f(Y(:,2));
    k2 = h.*f(Y(:,2)+alpha.*k1);
    Y(:,3) = Y(:,2) + (1-1/2/alpha).*k1 + k2/2/alpha;
    ytil2=Y(1,3);
    util2=Y(2,3);
  else

% Second Order Adams-Bashforth method steps

    Y(:,3) = Y(:,2) + 3/2*h*f(Y(:,2))- h/2*f(Y(:,1));
    ytil2=Y(1,3);
    util2=Y(2,3);
 end
end

