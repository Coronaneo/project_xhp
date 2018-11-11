
function  f1 = DXX( n , phi  , h ,k2) 
  gammay = 1.0;
  M = floor(10/h) ;
  u = pi/10.*(-M:(M-1));
  s = h*(-M:(M-1));
  phi_hat = zeros(2*M,1);
  for v = 1:(2*M)
     phi_hat(v) = sum(phi.*exp(-i*u(v).*(s+10)));
  end
  f1 = -sum(u.^2.*(phi_hat.').*exp(i*u.*(h*n)));
end