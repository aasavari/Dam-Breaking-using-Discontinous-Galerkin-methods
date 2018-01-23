function V = Integrator(dt , dx , c0 , N , J)
V = 0;
V1 = zeros(N,1);

for n=1:N
    for j=1:J
        V1(n) = V1(n) + dx*c0(j,n);
    end
end

% trapezoidal rule

for n=1:N-1
    
   V = V + dt*0.5*(V1(n) + V1(n+1));  
   
end


end