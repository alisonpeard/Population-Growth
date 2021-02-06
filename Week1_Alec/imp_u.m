function y = imp_u(u0,T0,h,k)
a = (h/k)+((h^2)/(2*k));
b = 1+((h/k)*T0)+((h^2)*u0/(2*k))-(h/k);
c = -u0;

y = (-b+sqrt(b^2-(4*a*c)))/(2*a);


end