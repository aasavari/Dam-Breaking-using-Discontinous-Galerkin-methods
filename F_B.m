function c = F_B(a,b)
 if( a*b >= 0 && a>=0)
    c = a.*a/2; 
 end
 if( a*b >= 0 && a<=0)
    c = b.*b/2; 
 end
 if( a*b <= 0 && (a+b)>=0 )
    c = a.*a/2; 
 end
 if( a*b <= 0 && (a+b)<0 )
    c = b.*b/2; 
 end
end