function C = rk(a,b,c,d,e,f,dt,dx)
% c's at n
k1 = [a;b];
% c's at n-1
k11 = [c;d];

% derivative at n
f1 = F(a,b,c,d,dx);
% derivative at n-1
f11 = F(c,d,e,f,dx);

k2 = k1 + dt*f1;

k22 = k11 + dt*f11;

f2 = F(k2(1), k2(2), k22(1) , k22(2), dx);

C = k1 + 0.5*dt*(f1 + f2);


end
