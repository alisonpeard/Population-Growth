%Population Growth - Sampling parameters (implicit euler method)

%u0: Initial populations of u0 < 1 experience an
%initial increase in population density before decaying to extinction.
%Smaller populations take longer to reach their peak, and longer to decay
%to extinction. Populations that start at u0 > 1 begin decaying
%immediately. 

[t1, u1] = implicit_popgrowth(0.1,10,1.5,0,1);
[t2, u2] = implicit_popgrowth(0.1,10,1,0,1);
[t3, u3] = implicit_popgrowth(0.1,10,0.8,0,1);
[t4, u4] = implicit_popgrowth(0.1,10,0.5,0,1);
[t5, u5] = implicit_popgrowth(0.1,10,0.1,0,1);
[t6, u6] = implicit_popgrowth(0.1,10,0.01,0,1);

figure(1)
plot(t1, u1, t2, u2, t3, u3, t4, u4, t5, u5, t6, u6)
title('Effect of initial population density, k = 1')
xlabel('Time'), ylabel('Population Density')
legend('u0 = 1.5', 'u0 = 1', 'u0 = 0.8', 'u0 = 0.5', 'u0 = 0.1', 'u0 = 0.01')

%effect of k: increases in the value of k decrease the peak population (if
%u0 < 1) and decrease the rate at which the population goes extinct. 

[t1, u1] = implicit_popgrowth(0.1,10,0.5,0,0.01);
[t2, u2] = implicit_popgrowth(0.1,10,0.5,0,0.1);
[t3, u3] = implicit_popgrowth(0.1,10,0.5,0,1);
[t4, u4] = implicit_popgrowth(0.1,10,0.5,0,10);
[t5, u5] = implicit_popgrowth(0.1,10,0.5,0,100);
[t6, u6] = implicit_popgrowth(0.1,10,0.5,0,1000);

figure(2)
plot(t1, u1, t2, u2, t3, u3, t4, u4, t5, u5, t6, u6)
title('Effect of parameter k, u0 = 0.5')
xlabel('Time'), ylabel('Population Density')
legend('k = 0.01', 'k = 0.1', 'k = 1', 'k = 10', 'k = 100', 'k = 1000')


