function A = u0_dm(x)
if ( x>=-1 && x<0 )
    A = [ 2; 0;0];
end

if(x>=0 && x<=1)
        A = [1;0;0];
end
end