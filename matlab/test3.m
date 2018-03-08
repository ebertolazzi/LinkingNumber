function test

  hold off ;

  [x,y,z] = random(500) ;
  plot3(x,y,z,'-b','LineWidth',2) ;
  P = [ x'; y'; z'; ] ;
  
  hold on ;

  [x,y,z] = random(500) ;
  plot3(x,y,z,'-r','LineWidth',2) ;
  Q = [ x'; y'; z'; ] ; 
  [L,E] = lk( P, Q );
  
  axis equal ;
  
  fprintf('L = %d, err = %g\n', L(1,2), E(1,2) ) ;
end

function [x,y,z] = random(npts)
  x = rand(npts,1);
  y = rand(npts,1);
  z = 0.1*rand(npts,1);
end
