function dt = rhs(v,k)

u = v(1);
y = v(2);

du = (1/k)*u*(1-u-y);
dy = u;

dt = [du;dy];

end