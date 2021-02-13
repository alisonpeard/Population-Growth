%2-step Adams Bashford

k = 1;
tf = 10;
h = 0.01;
N = tf/h; 
t = 0:h:tf;

%for Crank-Nicholson, order 2 accuracy:
theta = 0.5;

%pre-allocate arrays
u = zeros(length(t),1); %population size
y = zeros(length(t),1);

u_cond = 0.1:0.1:1;
u_max = zeros(length(u_cond),1);
y_crit = zeros(length(u_cond),1);

for j = 1:length(u_cond)
    u0 = u_cond(j);
    y0 = 0;

%assign initial values
u(1) = u0; y(1) = y0;

%midpoint for first time step (same order of accuracy)
    u(2) = theta*exp_u(u(1),y(1),h,k) + (1-theta)*imp_u(u(1),y(1),h,k);
    y(2) = theta*exp_y(u(1),y(1),h) + (1-theta)*imp_y(u(2),y(1),h);

for i = 2:N
    vec = (h/2)*((3*prime(u(i),y(i),k)-prime(u(i-1),y(i-1),k)))+[u(i);y(i)];
    u(i+1) = vec(1);
    y(i+1) = vec(2);
end

hold on
plot(t,u)

[umax, maxt] = findpeaks(u);
ycrit = y(maxt);
if isempty(umax)
    umax = u0;
    ycrit = 0;
end

u_max(j) = umax;
y_crit(j) = ycrit;

end

%calculate predicted u_max and y_crit
u_num = 1 + k*log(k./(1+k-u_cond));
y_num = -k*log(k./(1+k-u_cond));


figure(2)
plot(u_cond, u_max, u_cond, u_num)
xlabel('u0')
ylabel('u_{max}')
legend('Numerical Maximum', 'Predicted Maximum')

figure(3)
plot(u_cond, y_crit, u_cond, y_num)
xlabel('u0')
ylabel('Critical toxin accumulation')
legend('Numerical y_{crit}', 'Predicted y_{crit}')
