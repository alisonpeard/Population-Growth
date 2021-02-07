function u = imp_u(u0, y0, h, k)
a = ((h^2)/k)+(h/k);
b = 1+((h/k)*y0)-(h/k);
c = -u0;

u = (-b+sqrt(b^2-(4*a*c)))/(2*a);

end