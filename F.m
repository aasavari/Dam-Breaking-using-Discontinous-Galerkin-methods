function f = F(a,b,c,d,dx)
M1 = [-1 -1; 3 -3];
M2 = [1 1; -3 -3];

f = M1*[a;b]/dx + M2*[c;d]/dx;

end
