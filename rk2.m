function C = rk2(a,b,c,d,dt,dx)
      
      k1    = [a;b];       
      k11   = [c;d];   
      f1   = F(a,b,c,d,dx);
      k2   = k1    + 0.5*dt*f1   ;
      f2   = F(k2(1)  , k2(2)  , k11(1)   , k11(2)   , dx ) ;
      k3   = k1   + 0.5*dt*f2   ;
      f3  = F( k3(1) , k3(2) , k11(1) , k11(2)  , dx);
      k4 = k1   + dt*f3 ;
      f4 = F(k4(1), k4(2), k11(1), k11(2) , dx);
      
      C = k1 + dt*(f1/6.0 + f2/3.0 + f3/3.0 + f4/6.0);
     
end