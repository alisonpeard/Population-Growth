function [Utr, Uab] = integrate(Un, Un_1, h, kappa,theta)
    Utr = theta_method(Un, h, kappa,theta);
    Uab = ab2(Un,Un_1, h, kappa);
end