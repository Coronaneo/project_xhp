k = 0.00005;
h = 1/51.2;
epi = 1.0;
dim = 20/h;
DATA = zeros(dim,dim);
s = h*(-dim/2:(dim/2-1));
 
for p = 1:dim
    for q = 1:dim
        DATA(p,q) = 1/sqrt(pi*epi)*exp(-(s(p)^2+s(q)^2)/(2*epi));
    end
end

global DT 
DT = zeros(dim,dim);
for t_ = 1 : 40/k
   DT = CNSP(DATA,t_,h);
   DATA = DT;
end