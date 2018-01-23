format long;
clear all;
clc;
% parameters
a=1.0;
T = 0.1;
J=50;      
dx = 2.0/J; 
X = -1:dx:1;
h= 0.01;
x = -1:h:1;
dt = dx/10;
N = T/dt;
c0 = zeros(J,N);
c1 = zeros(J,N);
%  initial projection 
for j=1:J
   x1 = (X(j) + X(j+1) )/2 - 0.5*dx/sqrt(3);
   x2 = (X(j) + X(j+1) )/2 + 0.5*dx/sqrt(3);
 
    
    c0(j,1) = 0.5*( u0(x1) + u0(x2) );

              
    c1(j,1) = 1.5*(-u0(x1)/(sqrt(3)) + u0(x2)/sqrt(3));

end

% forward euler method
% for n = 2:N
%     c0(1,n) = c0(1,n-1) - a*dt*(c0(1,n-1) + c1(1,n-1) - c0(J , n-1) - c1(J,n-1))/dx;
%     c1(1,n) = c1(1,n-1) + a*3*dt*(c0(1,n-1) - c1(1,n-1) - c0(J, n-1) - c1(J,n-1))/dx;
%     
%     for j=2:J
%     c0(j,n) = c0(j,n-1) - a*dt*(c0(j,n-1) + c1(j,n-1) - c0(j-1 , n-1)  -c1(j-1,n-1))/dx;
%     c1(j,n) = c1(j,n-1) + a*3*dt*(c0(j,n-1) - c1(j,n-1) - c0(j-1 , n-1) - c1(j-1,n-1))/dx;
%     end
% end



% RK2 
% for n=1:N-1
% 
% A = rk( c0(1,n) , c1(1,n) , c0(J,n) , c1(J,n) , c0(J-1,n) , c1(J-1,n) , dt, dx );
% c0(1,n+1) = A(1);
% c1(1,n+1) = A(2);
% 
% A = rk( c0(2,n) , c1(2,n) , c0(1,n) , c1(1,n) , c0(J,n) , c1(J,n), dt , dx );
% c0(2,n+1) = A(1);
% c1(2,n+1) = A(2);
%    
% for j=3:J
% A = rk( c0(j,n) , c1(j,n) , c0(j-1,n) , c1(j-1,n) , c0(j-2,n) , c1(j-2,n), dt , dx );
% c0(j,n+1) = A(1);
% c1(j,n+1) = A(2);
% end  
% end

% RK4
for n=1:N-1
 A = rk2(c0(1,n),c1(1,n),c0(J,n),c1(J,n),dt,dx); 
 c0(1,n+1) = A(1);
 c1(1,n+1) = A(2);  
 
 
    
for j=2:J
       
 A=rk2(c0(j,n),c1(j,n),c0(j-1,n),c1(j-1,n),dt,dx); 
 c0(j,n+1) = A(1);
 c1(j,n+1) = A(2); 
 
end
    
end
% volume integral
% V = zeros(N,1);
% for n=1:N
%     for j=1:J
%         V(n) = V(n) + dx*c0(j,n);
%     end
% end
% 
% INT = 0;
% for n=1:N-1
%     
%    INT = INT + dt*0.5*(V(n) + V(n+1));  
%    
% end
 V = Integrator(dt , dx, c0 , N , J);
 display(V)
 
% limiters on c1(j,n)
for n=1:N
    r =  c1(1,n);
    p = (c0(2,n)-c0(1,n));
    q = (c0(1,n)-c0(J,n));
    c1(1,n) = minmod( r, p , q );
    
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


err(1) = abs(c0(1,N) - c1(1,N) - Exact(-1,0.1));

for j=2:J+1
     
err(j) = abs(c0(j-1,N) + c1(j-1,N) - Exact(X(j),0.1));

end
e = min(err);
display(e)

Y = zeros(2*J+1,1);
for j=1:J
    Y(j) = Exact(x(j),0.1);

end

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




