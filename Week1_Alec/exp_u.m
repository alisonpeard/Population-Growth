function y = exp_u(u0,T0,h,k)
y = u0+((h/k)*u0*(1-u0-T0));
end