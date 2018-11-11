function B = CNSP(A,k,h,k2,epi)
  gammay = 1.0;
  M = floor(10/h) ;
  u = pi/10.*(-M:(M-1));
  s = h*(-M:(M-1));
  for p = 1:(2*M)
      for q = 1:(2*M)
          B(p,q) = k*(A(p,q)) + k*i*epi/2*( DXX(p,A(:,q).',h, k2)+ DXX(q,A(p,:),h ,k2))-  ...
                   k*i/(2*epi)*(A(p,q))*( (s(p))^2 + (gammay*s(q))^2  )-  ...
                   k*i*k2/(epi)*((abs(A(p,q)))^2)*(A(p,q)) ;
      end
  end
end

