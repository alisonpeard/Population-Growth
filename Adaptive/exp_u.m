function y = exp_u(u0,y0,h,k)
y = u0+((h/k)*u0*(1-u0-y0));
end