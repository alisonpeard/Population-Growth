% Milne Device Rough Work
% Ab: predictor, Tr: corrector
ctr = - 1/12;
cab = 5/12;

% initial parameters
u0 = 0.5;
Tfinal = 10;
N = 100; % # timesteps
h0 = Tfinal/N; % timestep size
hvec = [h0 ; h0];

kappa = 1;
TOL = 1e-3;

% algorithm
Tn = 0;
U0_vec = [0;u0];

Un_1 = U0_vec;
Un = cn2(U0_vec,h0,kappa);
U = [Un_1 Un];

n = 2;
t = 0;
while t<Tfinal
    
    h0 = hvec(n);
    [Utr, Uab] = integrate(Un, Un_1, h0, kappa);
    E = abs(ctr/(cab-ctr)) * norm(Utr - Uab,2);
    
    if E > (1/10)*h0*TOL && E<=h0*TOL % accept step, keep timestep
        hvec = [hvec; h0];
        Un_1 = Un;
        Un = Uab;
        U = [U Un];
        n = n+1;
        t = t+h0;

    elseif E < (1/10)*h0*TOL && E<=h0*TOL% accept step, double timestep for next time
        h0 = 2*h0;
        hvec = [hvec; h0];
        hsum = 0; i = 0;
        while hsum<h0
            hsum = hsum + hvec(n-i); i = i+1;
        end
        Un_1 = U(:,n-i);
        Un = Uab;
        U = [U Un];
        n = n+1;
        
    elseif E>h0*TOL
        hvec(n) = (1/2)*h0;
        Un_1=(1/2)*(Un_1 + Un);
        
    end


end

steps = zeros(1,length(hvec));
steps(1) = 0;
for n = 2:N
    steps(n) = steps(n-1) + hvec(n);
end

plot(steps,U(2,:))

function [Utr, Uab] = integrate(Un, Un_1, h, kappa)
    Utr = cn2(Un, h, kappa);
    Uab = ab2(Un,Un_1, h, kappa);
end