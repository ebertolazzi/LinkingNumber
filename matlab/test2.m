function test

  [x,y,z] = ellipse(1000) ;
  plot3(x,y,z,'-b','LineWidth',3) ;
  P = [ x'; y'; z'; ] ;

  [x,y,z] = pnts8knots(1000) ;
  hold on ;
  plot3(x,y,z,'-r','LineWidth',3) ;
  Q = [ x'; y'; z'; ] ; 

  [L,E] = lk( P, Q ) ;
  
  fprintf('L = %d, err = %g\n', L, E ) ;

end

function [x,y,z] = ellipse(npts)
  theta = [2*pi*linspace(0,1,npts)]';
  x     = 1*cos(3*theta);
  y     = 4*sin(3*theta) ;
  z     = zeros(size(theta)) ;
end

function [x,y,z] = pnts8knots(npts)
  t = [2*pi*linspace(0,1,npts)]';
  x = (2+cos(2*t)).*cos(3*t);
  y = (2+cos(2*t)).*sin(3*t);
  z = 2*sin(4*t);
end
