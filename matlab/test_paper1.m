function test

  hold off ;
  [x1,y1,z1] = circle(0.2,100) ;
  plot3(x1,y1,z1,'-.b','LineWidth',3) ;
  P1 = [ x1'; y1'; z1' ] ;
  hold on ;  
  
  [x2,y2,z2] = circle(1.2,100) ;
  plot3(x2,y2,z2,'-b','LineWidth',3) ;
  P2 = [ x2'; y2'; z2' ] ;

  [x3,y3,z3] = circle(3.5,100) ;
  plot3(x3,y3,z3,'-.b','LineWidth',3) ;
  P3 = [ x3'; y3'; z3' ] ;

  [x,y,z] = trefoil(1000) ;
  plot3(x,y,z,'-r','LineWidth',3) ;
  Q = [ x'; y'; z' ] ;
  
  [L,E] = lk( P1, Q ) ;
  fprintf('L1 = %d, err = %g\n', L(1,2), E(1,2) ) ;

  [L,E] = lk( P2, Q ) ;
  fprintf('L2 = %d, err = %g\n', L(1,2), E(1,2) ) ;

  [L,E] = lk( P3, Q ) ;
  fprintf('L3 = %d, err = %g\n', L(1,2), E(1,2) ) ;
  
  axis equal;

end

function [x,y,z] = trefoil(npts)
  t = [2*pi*linspace(0,1,npts)]';
  x = sin(t)+2*sin(2*t) ;
  y = cos(t)-2*cos(2*t) ;
  z = -sin(3*t) ;
end


function [x,y,z] = circle(R,npts)
  u = [2*pi*linspace(0,1,npts)]';
  x = R*cos(u) ;
  y = R*sin(u) ;
  z = zeros(size(u)) ;
end
