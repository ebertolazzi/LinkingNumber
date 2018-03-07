function test

  [x,y,z] = random(5000) ;
  plot3(x,y,z,'-b','LineWidth',2) ;
  P = [ x'; y'; z'; ] ;

  [x,y,z] = random(5000) ;
  plot3(x,y,z,'-r','LineWidth',2) ;
  Q = [ x'; y'; z'; ] ; 
  [L,E] = lk( P, Q );
  
  fprintf('L = %d, err = %g\n', L, E ) ;
end

function [x,y,z] = random(npts)
  x = rand(npts,1);
  y = rand(npts,1);
  z = 0.001*rand(npts,1);
end
