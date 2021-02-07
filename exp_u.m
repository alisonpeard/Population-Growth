function u = exp_u(u0, y0, h, k)
u = u0+((h/k)*u0*(1-u0-y0));
end