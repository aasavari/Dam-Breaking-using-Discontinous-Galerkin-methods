format long;
clear all;
clc;
% parameters
T = 2.0;
J=100;      
dx = 2.0/J; 
X = -1:dx:1;
h= 0.01;
x = -1:h:1;
dt = dx/100;
N = T/dt;
c0 = zeros(J,N);
c1 = zeros(J,N);
% initial condition
 

%  initial projection 
for j=1:J
   x1 = (X(j) + X(j+1) )/2 - 0.5*dx/sqrt(3);
   x2 = (X(j) + X(j+1) )/2 + 0.5*dx/sqrt(3);
  c0(j,1) = 0.5*( u0(x1) + u0(x2) );
  c1(j,1) = 1.5*(-u0(x1)/(sqrt(3)) + u0(x2)/sqrt(3));

end
%forward euler for burgers
for n=1:N-1
      p = F_B((c0(1,n) + c1(1,n)) , c0(2,n)-c1(2,n) );
      q = F_B((c0(J,n) + c1(J,n)) , c0(1,n)-c1(1,n) );
        
      c0(1,n+1) = c0(1,n) - dt*(p+q)/dx ;
      c1(1,n+1) = c1(1,n) + 3*dt*(-p-q + c0(1,n)*c0(1,n) + c1(1,n)*c1(1,n)/3)/dx ;
    for j=2:J-1
        
        p = F_B( (c0(j,n) + c1(j,n)) , (c0(j+1,n)-c1(j+1,n)) );
        q = F_B((c0(j-1,n) + c1(j-1,n)) , (c0(j,n)-c1(j,n)) );
        
        c0(j,n+1) = c0(j,n) - dt*(p-q)/dx ;
        
        c1(j,n+1) = c1(j,n) + 3*dt*(-p-q + c0(j,n)*c0(j,n) + c1(j,n)*c1(j,n)/3)/dx ;
        
    end
    
     p = F_B((c0(J,n) + c1(J,n)) , c0(1,n)-c1(1,n) );
       q = F_B((c0(J-1,n) + c1(J-1,n)) , c0(J,n)-c1(J,n) );
        
       c0(J,n+1) = c0(J,n) - dt*(p+q)/dx ;
       c1(J,n+1) = c1(J,n) + 3*dt*(-p-q + c0(J,n)*c0(J,n) + c1(J,n)*c1(J,n)/3)/dx ;
end
% limiters on c1(j,n)
for n=1:N
    r =  c1(1,n);
    p = (c0(2,n)-c0(1,n));
    q = (c0(1,n)-c0(J,n));
    c1(1,n) = minmod( r, p , q);
    
    for j=2:J-1
        r = c1(j,n);
        p= c0(j+1,n)-c0(j,n);
        q= c0(j,n)-c0(j-1,n);
        c1(j,n) = minmod( r , p , q ) ;
    end
    r = c1(J,n);
    p = c0(1,n)-c0(J,n);
    q = c0(J,n)-c0(J-1,n);
    c1(J,n) = minmod( r , p , q );
end
err = zeros(J+1,1);


err(1) = abs(c0(1,N) - c1(1,N) - Exact(-1));

for j=2:J+1
     
err(j) = abs(c0(j-1,N) + c1(j-1,N) - Exact(X(j)));

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
   z = c0(j,N) + c1(j,N)*(2*y - l -r)/dx;
   plot(y, z);
end
hold off

