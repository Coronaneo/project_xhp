epi = 0.3
k2 = -1.9718
gammay = 1.0

k = 0.00005
h = 0.1
dim = 20/h - 1 
s = h*(-(dim-1)/2:(dim-1)/2);

%initial value
U0 = zeros ( dim,dim);
for p = 1:dim
    for q = 1:dim
        U0(p,q) = 1/sqrt(pi*epi)*exp(-(s(p)^2+s(q)^2)/(2*epi));
    end
end

A = Matrix(epi , k ,h );

global width2 B2
width2 = zeros(1,40/k);
B2 = zeros(dim,dim);
%iterate
times = 40/k
DT = zeros(dim,dim);
M1 = U0;

delta1 = s;
delta2 = s;
delta1(1) = [];
delta2(dim)= [];
delta = abs(1/3.*(delta1.^3-delta2.^3));

for t = 1:times
   DT = abs(M1.^2)*k2 + 1/2*gammay*ones(dim,dim).*s.^2 + 1/2*ones(dim,dim).*s.^2.';
   DT = M1.*DT;
   b0 = reshape(DT.',1,dim^2).';
   b1 = bicgstab(A,b0);
   
  % compute condensate 
   B2 = reshape(b1,dim,dim).';
   M1 = B2;
   M2 = B2;
   M3 = B2;
   M2(1,:) = [];
   M3(dim,:) = [];
   MM = (abs(M2.^2) + abs(M3.^2))/2; 
   m = sum(MM,2).';
   width2(t) = 1/dim*sqrt( sum(delta.*m )) ;
end

B2
width2

