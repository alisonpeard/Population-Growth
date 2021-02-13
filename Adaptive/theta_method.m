function v1 = theta_method(v0,h,k,theta)
u0 = v0(1);
y0 = v0(2);
u1 = theta*exp_u(u0,y0,h,k) + (1-theta)*imp_u(u0,y0,h,k);
y1 = theta*exp_y(u0,y0,h) + (1-theta)*imp_y(u1,y0,h);

v1 = [u1,y1];