function A = Matrix( epi , k , h )
a = i*epi/k;
beta = (epi/h)^2;
d = 20/h - 1 ;

X1 = zeros(d,d);
for s = 2:d-1
    X1(s,s) = a - 4*beta;
    X1(s,s+1)= beta;
    X1(s,s-1)= beta;
end
   X1(1,1) =  a - 4*beta;
   X1(1,2) = beta;
   X1(d,d) =  a - 4*beta;
   X1(d,d-1) = beta;
   
X2 = eye(d)*beta;
X2(1,d)= beta;

A = zeros(d^2,d^2);
for u = 2:d-1
    A( ((u-1)*d+1):u*d , ((u-1)*d+1):u*d ) = X1;
    A( ((u-1)*d+1):u*d , ((u-2)*d+1):(u-1)*d ) = X2; 
    A( ((u-1)*d+1):u*d , ( u*d+1):(u+1)*d ) = X2.'; 
end
A( 1:d , 1:d ) = X1;
A( 1:d , (d+1):2*d ) = X2.'; 
A( ((d-1)*d+1):d^2 , ((d-1)*d+1):d^2 ) = X1;
A( ((d-1)*d+1):d^2 ,  ((d-2)*d+1):(d-1)*d ) = X2; 
 
end

