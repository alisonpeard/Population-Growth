%Milne device - adaptive step size 

clear all

k = 6;
tmax = 10;
h = 0.1;
u0 = 0.5;
TOL = 10^-5;
c_ab2 = 5/12;
c_trap = -1/12;
theta = 0.5;

t(1) = 0; %initial time
u(1) = u0; %initial system conditions
y(1) = 0; 
h_list(1) = h;

%need 2 initial conditions to use adams-bashford
 t(2) = h;
 u(2) = theta*exp_u(u(1),y(1),h,k) + (1-theta)*imp_u(u(1),y(1),h,k);
 y(2) = theta*exp_y(u(1),y(1),h) + (1-theta)*imp_y(u(2),y(1),h);
 h_list(2) = h;

while t(end) < tmax

    E = 10*h*TOL;
    
    while E > h*TOL

        %compute u_n+s with TR
        u1 = theta*exp_u(u(end),y(end),h,k) + (1-theta)*imp_u(u(end),y(end),h,k);
        y1 = theta*exp_y(u(end),y(end),h) + (1-theta)*imp_y(u(end),y(end),h);
        
        %compute u_n+s with AB2
        vec = (h/2)*((3*prime(u(end),y(end),k)-prime(u(end-1),y(end-1),k)))+[u(end);y(end)];
        u1_ab = vec(1);
        y1_ab = vec(2);
        
        %compute E
        E = abs(c_trap/(c_ab2 - c_trap))*norm([u1;y1] - [u1_ab;y1_ab], 2);
        
        if E > h*TOL
            h = h/2; disp('changed smaller')
            %remesh
            if length(t) > 1
                t(end-1) = t(end) - h;
                u(end-1) = 0.5*(u(end-1) + u(end));
                y(end-1) = 0.5*(y(end-1) + y(end));
            end
        end
    end
    
    %accept next time step
    
    t(end+1) = t(end) + h;
    u(end+1) = u1;
    y(end+1) = y1;
    h_list(end+1) = h; %tracking step size
    
   if E <= (h/10)*TOL
       h = 2*h; disp('changed larger')
       % remesh
       if length(t) > 3
            t(end-1) = t(end) - h;
            
            %find closest time point with which we can interpolate u_n-1
            for i = 2:length(t)-1  
                if t(end-i) < t(end-1)
                    t_ref = t(end-i);
                    u_ref = u(end-i);
                    y_ref = y(end-i);
                    dist = abs(t(end) - t_ref);
                    break
                end 
            end
            weight = (t(end-1) - t_ref)/dist; %w in range [0,1]
            
            u(end-1) = (1-weight)*(u_ref) + weight*(u(end));
            y(end-1) = (1-weight)*(y_ref) + weight*(y(end));
       end
       
   end
    
end

step_size = zeros(1,length(t)-1);
for i = 2:length(t)
step_size(i) = t(i) - t(i-1);
end

figure(1)
plot(t,u)
xlabel('time')
ylabel('Population')

figure(2)
plot(log(1:length(t)), log(h_list))
xlabel('log(iterate)')
ylabel('log(step size)')
disp(length(t))

