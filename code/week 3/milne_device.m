% Milne Device Rough Work
% need to change loop to be in h0 not horig
% Ab: predictor, Tr: corrector

clear; clc;

ctr = - 1/12;
cab = 5/12;

% initial parameters
u0 = 0.5;
Tfinal = 20;               % initial number of timesteps
horig = 0.1;          % timestep size
k0 = 1;kstep=1;kend=5;

% checking maximum
Umax_vec = zeros(2,kend-k0);
Umax_theor = zeros(2,kend-k0);

% algorithm
TOL = 1e-3;
U0_vec = [0; u0];

j = 1;
for kappa = k0:kstep:kend
    
    h0 = horig;
    hvec = [h0; h0];
    Un_1 = U0_vec;
    Un = theta_method(U0_vec,h0,kappa,1/2);
    U = [Un_1 Un];

    n = 2;
    t(1) = 0;
    t(2) = horig;
    
    while t(n)<Tfinal
        
        h0 = hvec(n);
        [Utr, Uab] = integrate(Un, Un_1, h0, kappa);
        E = abs(ctr/(cab-ctr)) * norm(Utr - Uab,2);

        if E > (1/10)*h0*TOL && E<=h0*TOL
        % accept step, keep timestep
            disp('E in correct range, E = ' + string(E) + ", h0*TOL = " + string(h0*TOL))
            hvec(n+1) = h0;
            Un_1 = Un;
            Un = Uab;
            U = [U Un];
            t(n+1) = t(n)+h0;
            n = n+1;
            

        elseif E < (1/10)*h0*TOL && E<=h0*TOL
        % accept step, double timestep for next time
            disp('E too small, E = ' + string(E) + ", h0*TOL = " + string(h0*TOL))
            h0 = 2*h0;
            hvec(n+1) = h0;
            hsum = 0; i = 0;
            
            % reverse interpolation
            Un_1 = 2*Un - Uab;
            Un = Uab;
            
            U = [U Un];
            t(n+1) = t(n)+h0;
            n = n+1;

        elseif E>h0*TOL
        % reject step, halve timestep and try again
            disp('E too large, E = ' + string(E) + ", h0*TOL = " + string(h0*TOL))
            hvec(n) = (1/2)*h0;
            t(n) = t(n) - hvec(n);
            Un_1=(1/2)*(Un_1 + Un);

        end
    end
    
    % calculate max and critical point
    klog = kappa*log(kappa/(1+kappa-u0));
    Umax_theor(:,j) = [1+klog; -klog];
    Umax_vec(1,j) = max(U(2,:));
    ix = find(U(2,:)==Umax_vec(1,j));
    Umax_vec(2,j) = U(1,ix);
    j = j+1;

    % plot
    M = length(hvec);
    disp('number of steps: ' + string(M));
    steps = zeros(size(hvec));
    steps(1) = 0;
    for m = 2:M
        steps(m) = steps(m-1) + hvec(m);
    end
    
    subplot(2,2,1);
    plot(steps,U(2,:))
    hold on;
    
    subplot(2,2,2)
    plot(steps, hvec,'-*')
    hold on;
    
end

subplot(2,2,1);
title('Milne device, U0 = ' + string(u0))
K = k0:kstep:kend;
legendStrings = "kappa = " + string(K);
legend(legendStrings)

subplot(2,2,2);
title('stepsize')
K = k0:kstep:kend;
legendStrings = "kappa = " + string(K);
legend(legendStrings)

hold off;

%%
subplot(2,2,3);
plot(K(2:end),Umax_vec(1,2:end),K(2:end),Umax_theor(1,2:end));
legend(["numerical max", "theoretical max"])
ylabel("log scale");
title("max of U comparison, U0 = " + string(u0))
hold off;

subplot(2,2,4);
plot(K(2:end),Umax_vec(2,2:end),K(2:end),Umax_theor(2,2:end));
legend(["numerical max", "theoretical max"])
ylabel("log scale");
title("max of Y comparison, U0 = " + string(u0))

hold off;

%%

function [Utr, Uab] = integrate(Un, Un_1, h, kappa)
    Utr = theta_method(Un, h, kappa, 1/2);
    Uab = ab2(Un_1, Un, h, kappa);
end