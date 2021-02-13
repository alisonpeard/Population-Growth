function [Utr, Uab] = integrate(Un, Un_1, h, kappa)
    Utr = cn2(Un, h, kappa);
    Uab = ab2(Un,Un_1, h, kappa);
end