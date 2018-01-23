format long;
clear all;
clc;
% parameters
a=1.0;
T = 0.1;
J=100;      
dx = 2.0/J; 
X = -1:dx:1;
h= 0.01;
x = -1:h:1;
dt = dx/10;
N = T/dt;
c0 = zeros(J,N);
c1 = zeros(J,N);
c2 = zeros(J,N);
%  initial projection 
for j=1:J
   x1 = (X(j) + X(j+1) )/2 - 0.5*dx*sqrt(3/5);
   x2 = (X(j) + X(j+1) )/2 ;
   x3 = (X(j) + X(j+1) )/2 + 0.5*dx*sqrt(3/5);
 
    
    c0(j,1) = 0.5*(5.0*u0(x1)/9.0 + 8.0*u0(x2)/9.0 +  5.0*u0(x3)/9.0);

              
    c1(j,1) = 1.5*(-5.0*u0(x1)*(sqrt(3/5))/9.0 + 5.0*u0(x3)*sqrt(3/5)/9.0);
    
    c2(j,1) = 2.5*( 5.0*u0(x1)*0.4/9.0 - 8.0*0.5*u0(x2)/9.0 + 5.0*u0(x3)*0.4/9.0);
    
end

% forward euler method
for n = 1:N-1
    c0(1,n+1) = c0(1,n) + a*dt*(   -c0(1,n) - c1(1,n) - c2(1,n) + c0(J ,n) + c1(J,n) + c2(J,n))/dx;
    c1(1,n+1) = c1(1,n) + 3*a*dt*( -c0(1,n) - c1(1,n) - c2(1,n) - c0(J ,n) - c1(J,n) - c2(J,n) + 2*c0(1,n))/dx;
    c2(1,n+1) = c2(1,n) + 5*a*dt*( -c0(1,n) - c1(1,n) - c2(1,n) + c0(J ,n) + c1(J,n) + c2(J,n) + 2*c1(1,n))/dx;
    for j=2:J
    c0(j,n+1) = c0(j,n) + a*dt*( -c0(j,n) - c1(j,n) - c2(j,n) + c0(j-1 ,n) + c1(j-1,n) + c2(j-1,n))/dx;
    c1(j,n+1) = c1(j,n) + 3*a*dt*( -c0(j,n) - c1(j,n) - c2(j,n) - c0(j-1 ,n) - c1(j-1,n) - c2(j-1,n) + 2*c0(j,n))/dx;
    c2(j,n+1) = c2(j,n) + 5*a*dt*( -c0(j,n) - c1(j,n) - c2(j,n) + c0(j-1 ,n) + c1(j-1,n) + c2(j-1,n) + 2*c1(j,n))/dx;
    end
end
% error
err(1) = abs(c0(1,N) - c1(1,N) - Exact(-1,0.1));

for j=2:J+1
     
err(j) = abs(c0(j-1,N) + c1(j-1,N) - Exact(X(j),0.1));

end
e = min(err);
display(e)

% plot solution
plot(x, u0(x));
hold on
for j=1:J
   
   l = X(j);
   r = X(j+1);
   y = l:h:r;
   var = (2*y - l -r)/dx;
   z = c0(j,N) + c1(j,N)*var + 0.5*c2(j,N)*(3*var.*var -1);
   plot(y, z);
end
hold off