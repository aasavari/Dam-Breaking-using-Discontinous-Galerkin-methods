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
dt = dx/100;
N = T/dt;
c0 = zeros(J,N);
c1 = zeros(J,N);
c2 = zeros(J,N);
c3 = zeros(J,N);
%  initial projection

for j=1:J
   p1 =  - sqrt(3.0/7.0 + 2.0*sqrt(1.2)/7.0);
    p2 =  - sqrt(3.0/7.0 - 2.0*sqrt(1.2)/7.0);
    p3 =    sqrt(3.0/7.0 - 2.0*sqrt(1.2)/7.0);
    p4 =    sqrt(3.0/7.0 + 2.0*sqrt(1.2)/7.0);
   x1 = (X(j) + X(j+1) )/2 - 0.5*dx*sqrt(3.0/7.0 + 2.0*sqrt(1.2)/7.0);
   x2 = (X(j) + X(j+1) )/2 - 0.5*dx*sqrt(3.0/7.0 - 2.0*sqrt(1.2)/7.0);
   x3 = (X(j) + X(j+1) )/2 + 0.5*dx*sqrt(3.0/7.0 - 2.0*sqrt(1.2)/7.0);
   x4 = (X(j) + X(j+1) )/2 + 0.5*dx*sqrt(3.0/7.0 + 2.0*sqrt(1.2)/7.0); 
   b1 = 0.5 - sqrt(30)/36.0;
   b2 = 0.5 + sqrt(30)/36.0;
   b3 = 0.5 + sqrt(30)/36.0;
   b4 = 0.5 - sqrt(30)/36.0;
 c0(j,1) = 0.5*( b1*u0(x1) + b2*u0(x2) +  b3*u0(x3) + b4*u0(x4) );
 c1(j,1) = 1.5*( b1*u0(x1)*p1 + b2*u0(x2)*p2 +  b3*u0(x3)*p3 + b4*u0(x4)*p4);
 c2(j,1) = 2.5*( b1*u0(x1)*(3.0*p1*p1 - 1.0)*0.5 + ...
                   b2*u0(x2)*(3.0*p2*p2 - 1.0)*0.5 + ...
                   b3*u0(x3)*(3.0*p3*p3 - 1.0)*0.5 + ...
                   b4*u0(x4)*(3.0*p4*p4 - 1.0)*0.5);
 c3(j,1) = 3.5*( b1*u0(x1)*( 3.0*p1 - 5.0*p1*p1*p1)/2.0 +  ...
                   b2*u0(x2)*(3.0*p2 - 5.0*p2*p2*p2)/2.0 ...
                +  b3*u0(x3)*(3.0*p3 - 5.0*p3*p3*p3)/2.0 + ...
                   b4*u0(x4)*(3.0*p4 - 5.0*p4*p4*p4)/2.0 );

end

% forward euler method
for n = 1:N-1
    
    c0(1,n+1) = c0(1,n) + a*dt*( -c0(1,n) - c1(1,n) - c2(1,n) -c3(1,n)  ...
                               + c0(J ,n) + c1(J,n) + c2(J,n) + c3(J,n))/dx;                        
    c1(1,n+1) = c1(1,n) + 3.0*a*dt*( -c0(1,n) - c1(1,n) - c2(1,n) -c3(1,n) ... 
        - c0(J ,n) - c1(J,n) - c2(J,n) -c3(J,n) + 2*c0(1,n))/dx;
   
    
    c2(1,n+1) = c2(1,n) + 5.0*a*dt*( -c0(1,n) - c1(1,n) - c2(1,n) - c3(1,n)  ...
        + c0(J ,n) + c1(J,n) + c2(J,n) + c3(J,n) + 2*c1(1,n))/dx;
    
    
    c3(1,n+1) = c3(1,n) + 7.0*a*dt*( - c0(1,n) - c1(1,n) - c2(1,n) - c3(1,n) ...
                                   - c0(J,n) - c1(J,n) - c2(J,n) - c3(J,n) ... 
                                   +2*c0(1,n) + 2*c2(1,n))/dx;
    for j=2:J
    c0(j,n+1) = c0(j,n) + a*dt*( -c0(j,n) - c1(j,n) - c2(j,n) -c3(j,n) +  ...
        c0(j-1 ,n) + c1(j-1,n) + c2(j-1,n) + c3(j-1,n))/dx;
    
    c1(j,n+1) = c1(j,n) + 3.0*a*dt*( -c0(j,n) - c1(j,n) - c2(j,n) - c3(j,n) ...
        - c0(j-1 ,n) - c1(j-1,n) - c2(j-1,n) - c3(j-1,n) + 2*c0(j,n))/dx;
    
    c2(j,n+1) = c2(j,n) + 5.0*a*dt*( -c0(j,n) - c1(j,n) - c2(j,n) -c3(j,n)  ....
        + c0(j-1 ,n) + c1(j-1,n) + c2(j-1,n) + c3(j-1,n) + 2*c1(j,n))/dx;
    
    c3(j,n+1) = c3(j,n) + 7.0*a*dt*( - c0(j,n) - c1(j,n) - c2(j,n) - c3(j,n) ...
                                   - c0(j-1,n) - c1(j-1,n) - c2(j-1,n) - c3(j-1,n) ... 
                                   + 2*c0(j,n) + 2*c2(j,n))/dx;
    end
end

err = zeros(J+1,1);


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
   
   z = c0(j,N) + c1(j,N)*var + c2(j,N)*(1.5*var.*var - 0.5) + c3(j,N)*(1.5*var - 2.5*var.*var.*var);
   plot(y, z);
end
hold off





