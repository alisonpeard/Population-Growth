% Theta method with Newton
% f(y,u)
clear;clc;

figure(1)
Tfinal = 20;
N = 100;                           % # timesteps
h = Tfinal/N;                      % timestep size
k0 = 1; kstep=1; kend=5;

theta = 1/2;
% theta = 1 => implicit Euler
% theta = 0 => explicit Euler

u0 = 0.5;
u0_vec = [0; u0];
umax_vec = zeros(2,kend-k0);
umax_theor = zeros(2,kend-k0);

subplot(2,2,1);
i = 1;
for kappa = k0:kstep:kend

    Tn = 0;
    U = zeros(2,N);
    steps = (0:N).*0.2;

    U(:,1) = u0_vec;

    for n=1:N
        un = U(:,n);
        U(:,n+1) = theta_method(un, h, kappa, theta);
    end
    
    klog = kappa*log(kappa/(1+kappa-u0));
    umax_theor(:,i) = [1+klog; -klog];
    
    umax_vec(1,i) = max(U(2,:));
    ix = find(U(2,:)==umax_vec(1,i));
    umax_vec(2,i) = U(1,ix);
    
    i = i+1;

    plot(steps,U(2,:));
    hold on;
end

title('Theta method with Newton, U0 = ' + string(u0))
K = k0:kstep:kend;
legendStrings = "kappa = " + string(K);
legend(legendStrings)

hold off;

%%

subplot(2,2,3);
semilogy(K(2:end),umax_vec(1,2:end),K(2:end),umax_theor(1,2:end));
title("U_{max} numerical & theoretical, U_0 = " + string(u0))
legend("numerical","theoretical");
ylabel("log scale");
xlabel("kappa")

subplot(2,2,4);
semilogy(K(2:end),umax_vec(2,2:end),K(2:end),umax_theor(2,2:end));
title("Y_{crit} numerical & theoretical, U_0 = " + string(u0))
legend("numerical","theoretical");
ylabel("log scale");
xlabel("kappa")

hold off;


%% U0 = 1.5
subplot(2,2,2);
u0 = 1.5;
u0_vec = [u0; 0];

i = 1;
for kappa = k0:kstep:kend

    Tn = 0;
    U = zeros(2,N);
    steps = (0:N);

    U(:,1) = u0_vec;

    for n=1:N
        un = U(:,n);
        U(:,n+1) = theta_method(un, h, kappa, theta);
    end
    
    i = i+1;

    plot(steps,U(1,:));
    hold on;
end

title('Theta method with Newton, U0 = ' + string(u0))
K = k0:kstep:kend;
legendStrings = "kappa = " + string(K);
legend(legendStrings)

saveas(gcf,'newton_plots.jpeg')

hold off;

%%
figure(2)
subplot(2,1,1)
semilogy(K(2:end),abs(umax_vec(1,2:end)-umax_theor(1,2:end)));
title("absolute difference of U_{max} numerical & theoretical, U_0 = " + string(u0))
xlabel("kappa")

subplot(2,1,2)
semilogy(K(2:end),abs(umax_vec(2,2:end)-umax_theor(2,2:end)));
title("absolute difference of y_{crit}, U_0 = " + string(u0))
xlabel("kappa")

hold off;