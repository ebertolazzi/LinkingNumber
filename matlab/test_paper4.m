function test
  close all ;
  
  data = [] ;

  subplot(2,2,1) ;
  [x,y,z] = trefoil(100,0) ;
  hold off ;
  plot3(x,y,z,'-r','LineWidth',1) ;
  axis equal;
  hold on ;

  subplot(2,2,2) ;
  for k=1:1000
    [x,y,z] = trefoil(100,0.5) ;
    plot3(x,y,z,'-r','LineWidth',1) ;
    hold on ;
    Q = [ x'; y'; z' ] ;
  
    [W,E] = Writhe( Q ) ;
    fprintf('W = %g, err = %g\n', W, E ) ;
    
    data = [ data, W ] ;
  end
  axis equal;

%  P = [ 0.0, 1.0, 0.0 ; ...
%        0.0, 0.0,-1.0 ; ...
%        0.0, 0.0, 1.0 ; ...
%        1.0, 0.0, 0.0 ; ...
%        -1.0, -1.0, 0.0 ] ;
%
%  [W,E] = lk( P.' ) ;
%  fprintf('W = %g, err = %g\n', W, E ) ;

  subplot(2,2,[3,4]) ;
  hold off ;
  histfit(data) ;
  hold on ;

end

function [x,y,z] = trefoil(npts,pert)
  t = 2*pi*linspace(0,1,npts).';
  x = sin(t)+2*sin(2*t)+pert*rand(size(t)) ;
  y = cos(t)-2*cos(2*t)+pert*rand(size(t)) ;
  z = -sin(3*t)+pert*rand(size(t)) ;
end
