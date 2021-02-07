function f = prime(u,y,k)

du = (1/k)*u*(1-u-y);
dy = u;
f = [du;dy];
end
