% Theta method with Newton
% did all this as f(y,u) not f(u,y)

figure(1)
Tfinal = 20;
N = 100; % # timesteps
h = Tfinal/N; % timestep size
k0 = 1;kstep=1;kend=5;
subplot(2,2,1);

theta = 1/2;
% theta = 1 => implicit Euler
% theta = 0 => explicit Euler

U0 = 0.5;
U0_vec = [U0; 0];
Umax_vec = zeros(2,kend-k0);
Umax_theor = zeros(2,kend-k0);

i = 1;
for kappa = k0:kstep:kend

    Tn = 0;
    U = zeros(2,N);
    steps = (0:N);

    U(:,1) = U0_vec;

    for n=1:N
        Un = U(:,n);
        U(:,n+1) = mynewton(Un, h, kappa, theta);
    end
    
    klog = kappa*log(kappa/(1+kappa-U0));
    Umax_theor(:,i) = [1+klog; -klog];
    
    Umax_vec(1,i) = max(U(1,:));
    ix = find(U(1,:)==Umax_vec(1,i));
    Umax_vec(2,i) = U(2,ix);
    
    i = i+1;

    plot(steps,U(1,:));
    hold on;
end

title('Theta method with Newton, U0 = ' + string(U0))
K = k0:kstep:kend;
legendStrings = "kappa = " + string(K);
legend(legendStrings)

%%

subplot(2,2,3);
semilogy(K(2:end),Umax_vec(1,2:end),K(2:end),Umax_theor(1,2:end));
title("U_{max} numerical & theoretical, U_0 = " + string(U0))
legend("numerical","theoretical");
xlabel("kappa")

subplot(2,2,4);
semilogy(K(2:end),Umax_vec(2,2:end),K(2:end),Umax_theor(2,2:end));
title("Y_{crit} numerical & theoretical, U_0 = " + string(U0))
legend("numerical","theoretical");
xlabel("kappa")



%% U0 = 1.5
subplot(2,2,2);
U0 = 1.5;
U0_vec = [U0; 0];

i = 1;
for kappa = k0:kstep:kend

    Tn = 0;
    U = zeros(2,N);
    steps = (0:N);

    U(:,1) = U0_vec;

    for n=1:N
        Un = U(:,n);
        U(:,n+1) = mynewton(Un, h, kappa, theta);
    end
    
    i = i+1;

    plot(steps,U(1,:));
    hold on;
end

title('Theta method with Newton, U0 = ' + string(U0))
K = k0:kstep:kend;
legendStrings = "kappa = " + string(K);
legend(legendStrings)

saveas(gcf,'newton_plots.jpeg')

hold off;

%%
figure(2)
subplot(2,1,1)
semilogy(K(2:end),abs(Umax_vec(1,2:end)-Umax_theor(1,2:end)));
title("absolute difference of U_{max} numerical & theoretical, U_0 = " + string(U0))
xlabel("kappa")

subplot(2,1,2)
semilogy(K(2:end),abs(Umax_vec(2,2:end)-Umax_theor(2,2:end)));
title("absolute difference of y_{crit}, U_0 = " + string(U0))
xlabel("kappa")
