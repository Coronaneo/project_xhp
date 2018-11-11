function B = CNSP(A,k,h)
  epi = 1.0;
  %k = 0.00005;
  %h = 1/51.2;
  gammay = 1.0;
  k2 = -2.0;
  M = floor(10/h) ;
  u = pi/10.*(-M:(M-1));
  s = h*(-M:(M-1));
  for p = 1:(2*M)
      for q = 1:(2*M)
          B(p,q) = k*(A(p,q)) + k*i*epi/2*( DXX(p,A(:,q).',h)+ DXX(q,A(p,:),h ))-  ...
                   k*i/(2*epi)*(A(p,q))*( (s(p))^2 + (gammay*s(q))^2  )-  ...
                   k*i*k2/(epi)*((abs(A(p,q)))^2)*(A(p,q)) ;
      end
  end
end

