% Theta method with Newton
Tfinal = 10;
N = 100; % number of timesteps to take
h = Tfinal/N; % timestep
k0 = 0;kstep=1;kend=5;
subplot(1,2,1);

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
        U(:,n+1) = mynewton(Un,h, kappa);
    end
    
    klog = kappa*log(kappa/(1+kappa-U0));
    Umax_theor(:,i) = [1+klog; -klog];
    Umax_vec(1,i) = max(U(1,:));
    Umax_vec(2,i) = max(U(2,:));
    i = i+1;

    plot(steps,U(1,:));
    hold on;
end

title('Theta method with Newton, U0 = ' + string(U0))
K = k0:kstep:kend;
legendStrings = "kappa = " + string(K);
legend(legendStrings)

figure(2)
plot(K(2:end),Umax_vec(1,2:end),K(2:end),Umax_theor(1,2:end));
legend(["numerical max", "analytic max"])
title("max of U")

%% U0 = 1.5
subplot(1,2,2);
U0 = 1.5;
U0_vec = [U0; 0];

for kappa = k0:kstep:kend

    Tn = 0;
    U = zeros(2,N);
    steps = (0:N);

    U(:,1) = U0_vec;

    for n=1:N
        Un = U(:,n);
        U(:,n+1) = mynewton(Un,h, kappa);
    end

    plot(steps,U(1,:));
    hold on;
end

title('Theta method with Newton, U0 = ' + string(U0))
K = k0:kstep:kend;
legendStrings = "kappa = " + string(K);
legend(legendStrings)

saveas(gcf,'newton_ali.jpeg')
