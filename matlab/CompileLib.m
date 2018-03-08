clc;
clear functions;

NAMES = { 'lk', 'Writhe' } ;

LIBS = [ ...
  '-I../src ' ...
  '../src/linking_number.cc ' ...
] ;

disp('---------------------------------------------------------');
for k=1:length(NAMES)
  N=NAMES{k} ;
  fprintf(1,'Compiling: %s\n',N) ;

  CMD = [ 'while mislocked(''' N ''') ; munlock(''' N ''') ; end;'] ;
  disp(CMD);

  CMD = ['mex -output ',N,' -largeArrayDims ../src_mex/mex_',N,'.cc ', LIBS] ;
  if isunix
    if ismac
      CMD = [CMD, ' -lstdc++ CXXFLAGS="\$CXXFLAGS -Wall -O2 -g0"'] ;
    else
      CMD = [CMD, ' -lstdc++ CXXFLAGS="\$CXXFLAGS -Wall -O2 -g0"'] ;
    end
  elseif ispc
  end
  disp(CMD);
  eval(CMD);
end
disp('----------------------- DONE ----------------------------');
