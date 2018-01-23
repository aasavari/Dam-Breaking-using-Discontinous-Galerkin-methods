format long;
J = 63;              K = 63;
dx = 2.0/60.0;       dy = 2.0/60.0;
c = 0.8;             x= -1:dx:1;      y= -1:dx:1;
u0 = zeros(J,K,3);

for j=1:J
    for k= 1:K
           if( j>=17 && j<=47 && k>=17 && k<=47 )
                u0(j,k,:) = [2 0 0];
           else 
                u0(j,k,:) = [1 0 0];
           end
    end
end
plot3( x,y , u0(:,:,1));
