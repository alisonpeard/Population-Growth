% AB2 method
% seems much less effective at umax preds?

subplot(2,2,1);

Tfinal = 20;
N = 100; % # timesteps
h = Tfinal/N; % timestep size
k0 = 1;kstep=1;kend=5;

U0 = 0.5;
U0_vec = [0; 0.5];

Umax_vec = zeros(2,kend-k0);
Umax_theor = zeros(2,kend-k0);

i=1;
for kappa = k0:kstep:kend

    Tn = 0;
    U = zeros(2,N);
    U(:,1) = U0_vec;
    U(:,2) = cn2(U0_vec,h,kappa);
    steps = (0:N);

    for n=2:N
        Un = U(:,n);
        Un_1 = U(:,n-1);
        U(:,n+1) = ab2(Un,Un_1, h, kappa);
    end
    
    klog = kappa*log(kappa/(1+kappa-U0));
    Umax_theor(:,i) = [1+klog; -klog];
    
    Umax_vec(1,i) = max(U(2,:));
    ix = find(U(2,:)==Umax_vec(1,i));
    Umax_vec(2,i) = U(1,ix);
    
    i = i+1;

    plot(steps,U(2,:));
    hold on;
end

title('Adams-bashforth method, U0 = ' + string(U0))
K = k0:kstep:kend;
legendStrings = "kappa = " + string(K);
legend(legendStrings)

%%
subplot(2,2,3);
plot(K(2:end),Umax_vec(1,2:end),K(2:end),Umax_theor(1,2:end));
legend(["numerical max", "theoretical max"])
ylabel("log scale");
title("max of U comparison, U0 = " + string(U0))

subplot(2,2,4);
plot(K(2:end),Umax_vec(2,2:end),K(2:end),Umax_theor(2,2:end));
legend(["numerical max", "theoretical max"])
ylabel("log scale");
title("max of Y comparison, U0 = " + string(U0))
