clc; 
clear all; 
 

% parameters

TF =0.1;
J=50;      
dx = 2.0/J; 
X = -1:dx:1;
h= 0.01;
x11 = -1:h:1;
dt = dx/10;
N = TF/dt;
c0 = zeros(N,J,3);
c1 = zeros(N,J,3);


% initial conditions and boundary conditions

for j=1:J
   x1 = (X(j) + X(j+1) )/2 - 0.5*dx/sqrt(3);
   x2 = (X(j) + X(j+1) )/2 + 0.5*dx/sqrt(3);
 
    
    c0(1,j,:) = 0.5*( u0_dm(x1) + u0_dm(x2) );

              
    c1(1,j,:) = 1.5*(-u0_dm(x1)/(sqrt(3)) + u0_dm(x2)/sqrt(3));

end

for n=1:N
    
    c0(n,1,:) = [2 ; 0 ; 0];
    c1(n,1,:) = [0 ; 0 ; 0];
    
    c0(n,J,:) = [1 ; 0 ; 0];
    c1(n,J,:) = [0 ; 0 ; 0];
    
end

% project in time 
for n=1:N-1
    for j=2:J-1
 
    A = c0(n,j+1,:) - c1(n,j+1,:);
    B = c0(n,j  ,:) + c1(n,j,:);
    f1 = DER(A(1) , A(2), A(3));
    f2 = DER(B(1) , B(2 ), B(3)) ; 
   lambda1 = abs((0.5*B(2)/B(1) + 0.5*A(2)/A(1))) + sqrt(0.5*A(1) + 0.5*B(1));
   T1 = 0.5*(f1 + f2) - 0.5*lambda1*(A-B); 
   C = c0(n,j,:) - c1(n,j,:);
    D = c0(n,j-1,:) + c1(n,j-1 , :);
    f3 = DER(C(1) , C(2), C(3));
    f4 = DER(D(1) , D(2 ), D(3)) ;
    lambda2 = abs((0.5*D(2)/D(1) + 0.5*C(2)/C(1))) + sqrt((0.5*C(1) + 0.5*D(1)));
    T2 =  0.5*(f3 + f4 ) - 0.5*lambda2*(C-D); 
    c0(n+1,j,:) = c0(n,j,:) + 0.5*dt*(-T1 + T2)/dx;
   l = X(j);
    r = X(j+1); 
    q1 = c0(n,j,:) + ( 0.5*(l+ r) - 0.5*dx/sqrt(3))*c1(n,j,:);
    q2 =  c0(n,j,:) + ( 0.5*(l+ r) + 0.5*dx/sqrt(3))*c1(n,j,:);
   fq1 = DER(q1(1), q1(2) , q1(3));
    fq2 = DER(q2(1), q2(2) , q2(3));
    T3 = 0.5*(fq1 + fq2);
    c1(n+1,j,:) = c1(n,j,:) + 1.5*dt*( -T1 - T2 + T3 )/dx ;
       
    end
end
for n=1:N
    for j=1:J
        if(c0(n,j,1) > 2)
            c0(n,j,1) = 2;
        end
        if(c0(n,j,1) < 1)
            c0(n,j,1) = 1;
        end
    end
end

% limiters
for n=1:N
    
    for j=2:J-1
        r = c1(n,j,:);
        p= c0(n,j+1,:)-c0(n,j,:);
        q= c0(n,j,:)-c0(n,j-1,:);
        
         a = minmod(r(1) , p(1) , q(1) );
         b = minmod(p(2) , r(2 ) , p(2));
         c = minmod(r(3) , p(3 ) , q(3));
        
        c1(n,j,1) = a;
        c1(n,j,2) = b;
        c1(n,j,3) = c;     
    end
  
end

% plot solutions
Arr = zeros(1,201);
for j=1:201
    A=u0_dm(x11(j));
    Arr(j) = A(1);
end


hold on
plot(x11 , Arr);
for j=1:J
   l = X(j);
   r = X(j+1);
   y = l:h:r;
   var =  (2*y - l -r)/dx;
   z = c0(N,j,1) + c1(N,j,1)*var;
   plot(y, z);
end
hold off
% function to define derivative vector
function [ a , b , c] = DER( p , q , r ) 

a = q;
b = q*q/p + 0.5*p*p;
c = q*r/p;

end