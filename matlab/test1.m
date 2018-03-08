function test

  [x,y,z] = circle(50) ;
  plot3(x,y,z,'-b','LineWidth',3) ;
  P = [ x'; y'; z' ] ;

  [x,y,z] = pnts5foil(50) ;
  hold on ;
  plot3(x,y,z,'-r','LineWidth',3) ;
  Q = [ x'; y'; z' ] ;
  
  [L,E] = lk( P, Q ) ;
  
  fprintf('L = %d, err = %g\n', L(1,2), E(1,2)) ;

end

function [x,y,z] = circle(npts)
  theta = [2*pi*linspace(0,1,npts)]';
  x     = 2.5*cos(3*theta) ;
  y     = 2.5*sin(3*theta) ;
  z     = zeros(size(theta)) ;
end

function [x,y,z] = pnts5foil(npts)
  theta = [2*pi*linspace(0,1,npts)]';
  x     = (7/3)*sin(2*theta)-(2/3)*sin(3*theta) ;
  y     = (7/3)*cos(2*theta)+(2/3)*cos(3*theta) ;
  z     = 2*sin(5*theta) ;
end
