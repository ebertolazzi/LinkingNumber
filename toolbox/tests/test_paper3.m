
function test

  clear all ;
  
  ns = [5, 10, 50, 150 ] ;
    
  for kk=1:4

    subplot(2,2,kk) ;

    [x,y,z] = spira();
    hold off ;
    plot3([x,x(1)],[y,y(1)],[z,z(1)],'-r.','LineWidth',3);
    hold on

    [x1,y1,z1] = spiral(ns(kk),32*ns(kk));
    plot3([x1,x1(1)],[y1,y1(1)],[z1,z1(1)],'-b.');

    P = [ x ; y ; z ] ;
    Q = [ x1; y1; z1 ] ;

    tic
    [L,E] = lk( P, Q );
    toc

    fprintf('L = %d, err = %g\n', L, E ) ;
    
    [W,E] = writhe( Q );
    fprintf('W = %g, err = %g\n', W, E ) ;

    title(sprintf('LK = %d',L));

    axis equal;
  end
  N = 100000 ;
  [x1,y1,z1] = spiral(N,32*N);
  Q = [ x1; y1; z1 ] ;

  tic
  [L,E] = lk( P, Q );
  toc

  fprintf('L = %d, err = %g\n', L, E ) ;
  
end

function [x,y,z] = spiral(N,npts)
  t = [0:N/npts:N];
  x = cos(2*pi*t);
  y = sin(2*pi*t);
  z = t/t(end);
  
  x = x + 0.8*rand(size(x));
  y = y + 0.8*rand(size(y));
  z = z + 0.4*rand(size(z))/N;

  x = [x 2 2];
  y = [y 0 0];
  z = [z 1 0];
end

function [x,y,z] = spira()
  x = [0 0 -1.1 -1.1];
  y = [0 0 0 0];
  z = [-0.1 1.1 1.1 -0.1];
end



