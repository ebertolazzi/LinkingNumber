function test

  [x,y,z] = circle(50) ;
  plot3(x,y,z,'-b','LineWidth',3) ;
  P = [ x'; y'; z' ] ;
  
  xx = z ;
  yy = x+1 ;
  zz = y ;

  hold on ;
  plot3(xx,yy,zz,'-r','LineWidth',3) ;
  Q = [ xx'; yy'; zz' ] ;
  
  [L,E] = lk( P, Q ) ;
  
  fprintf('L = %d, err = %g\n', L(1,2), E(1,2) ) ;

end

function [x,y,z] = circle(npts)
  theta = [2*pi*linspace(0,1,npts)]';
  x     = 2.5*cos(3*theta) ;
  y     = 2.5*sin(3*theta) ;
  z     = zeros(size(theta)) ;
end
