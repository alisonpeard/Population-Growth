% Explicit method
Tfinal = 10;
N = 100; % number of timesteps to take
h = Tfinal/N; % timestep
k0 = 0;kstep=1;kend=5;

subplot(1,2,1);

U0 = 0.5;
for k = k0:kstep:kend

    Tn = 0;
    U = zeros(1,N);
    steps = (0:N);

    U(1) = U0;

    for n=1:N
        Un = U(n);
        U(n+1) = Un + h * Un*(1-Un-Tn)/k;
        Tn = Tn + h * (Un + U(n+1))/2;
    end

    plot(steps,U);
    hold on;
end

title('Explicit solution with composite trapezoidal, U0 = ' + string(U0))
K = k0:kstep:kend;
legendStrings = "k = " + string(K);
legend(legendStrings)

subplot(1,2,2);
U0 = 1.5;
for k = k0:kstep:kend

    Tn = 0;
    U = zeros(1,N);
    steps = (0:N);

    U(1) = U0;

    for n=1:N
        Un = U(n);
        U(n+1) = Un + h * Un*(1-Un-Tn)/k;
        Tn = Tn + h * (Un + U(n+1))/2;
    end

    plot(steps,U);
    hold on;
end

title('Explicit solution with composite trapezoidal, U0 = ' + string(U0))
K = k0:kstep:kend;
legendStrings = "k = " + string(K);
legend(legendStrings)

saveas(gcf,'explicit_ali.jpeg')
