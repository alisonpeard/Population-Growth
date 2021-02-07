clear all;
theta = 1/2;
k = 0.2;
u0vec = [0.1:0.05:1];

Tmax = 10;

hvec = [1e-6,1e-4,1e-3,1e-2,1e-1];

Umat = zeros(1,length(hvec));
Ymat = zeros(1,length(hvec));

for index = 1:length(hvec)
disp(index);
h = hvec(index);
N = ceil(Tmax/h);
time = 0:h:Tmax;
uvec = zeros(N,1);
yvec = zeros(N,1);


for j = 1:length(u0vec)
    u0 = u0vec(j);
    y0 = 0;
    uvec(1) = u0;
    yvec(1) = y0;
for i = 2:N+1
   u1 = theta*exp_u(u0,y0,h,k) + (1-theta)*imp_u(u0,y0,h,k);
   y1 = theta*exp_y(u0,y0,h) + (1-theta)*imp_y(u1,y0,h);
    
   u0 = u1;
   y0 = y1;
   uvec(i) = u0;
   yvec(i) = y0;
end
    
    %plot(time,uvec)
    [umax,ind] = findpeaks(uvec);
    yc = yvec(ind);
    if isempty(umax)
        umax = 1;
        yc = 0;
    end
    umaxvec(j) = umax;
    ycvec(j) = yc;
    hold all
end

[yvec,uvec] = preds(k,u0vec);

unorm = norm(uvec-umaxvec);
ynorm = norm(yvec-ycvec);

Umat(index) = unorm;
Ymat(index) = ynorm;
end
clf;
figure(1)
plot(u0vec,umaxvec,u0vec,uvec)
xlabel('u0')
ylabel('u_{max}')
legend('predicted max', 'numerical max')

figure(2)
plot(u0vec,ycvec,u0vec,yvec)
xlabel('u0')
ylabel('y_c')
legend('predicted y_c', 'numerical y_c')