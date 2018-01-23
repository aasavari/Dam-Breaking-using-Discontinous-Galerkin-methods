x= -1:0.01:1;
Y = zeros(201,1);
for j=1:201
    Y(j) = Exact(x(j));
end
plot(x, Y );