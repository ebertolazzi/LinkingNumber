clc;
clear functions;

old_dir = cd(fileparts(which(mfilename)));

NAMES = { 'lk', 'Writhe' };

lst_cc = dir('src/*.cc');

name_list   = {lst_cc.name};
folder_list = {lst_cc.folder};

MROOT = matlabroot;
DEFS  = ' ';

CMDBASE = [ 'mex -c -largeArrayDims -Isrc -Isrc_lib -Isrc/Utils ', DEFS];
if isunix
  CMDBASE = [CMDBASE, 'CXXFLAGS="\$CXXFLAGS -Wall -O2 -g" '];
elseif ispc
  CMDBASE = [CMDBASE, 'COMPFLAGS="\$COMPFLAGS -O2" '];
end

LIB_OBJS = '';
for k=1:length(name_list)
  [~,dname,~] = fileparts(folder_list{k});
  [~,bname,~] = fileparts(name_list{k});
  NAME        = [dname, '/', name_list{k} ];
  if isunix
    LIB_OBJS = [ LIB_OBJS, bname, '.o ' ];
  elseif ispc
    LIB_OBJS = [ LIB_OBJS, bname, '.obj ' ];
  end
  CMD = [CMDBASE ' -c ' NAME];
  disp('---------------------------------------------------------');
  disp(CMD);
  eval(CMD);
end

for k=1:length(NAMES)
  N=NAMES{k};
  disp('---------------------------------------------------------');
  fprintf(1,'Compiling: %s\n',N);

  CMD = [ 'while mislocked(''' N '''); munlock(''' N '''); end;'];
  eval(CMD);

  CMD = [ 'mex ', DEFS, ' -Isrc -Isrc_lib -output bin/', N ];
  CMD = [ CMD, ' -largeArrayDims src_mex/mex_', N ];
  CMD = [ CMD, '.cc ', LIB_OBJS ];

  if ismac
    CMD = [CMD, ' CXXFLAGS="\$CXXFLAGS -Wall -O2 -g"'];
  elseif isunix
    % Workaround for MATLAB 2020 that force dynamic link with old libstdc++
    % solution: link with static libstdc++
    % ARCH  = computer('arch');
    % PATH1 = [MROOT, '/bin/', ARCH];
    % PATH2 = [MROOT, '/extern/bin/', ARCH];
    CMD = [ CMD, ...
      ' CXXFLAGS="\$CXXFLAGS -Wall -O2 -g"' ...
      ' LDFLAGS="\$LDFLAGS -static-libgcc -static-libstdc++"' ...
      ' LINKLIBS="-ldl -L\$MATLABROOT/bin/\$ARCH -L\$MATLABROOT/extern/bin/\$ARCH -lMatlabDataArray -lmx -lmex -lmat -lm "' ...
    ];
  elseif ispc
    CMD = [CMD, 'COMPFLAGS="\$COMPFLAGS -O2" '];
  end

  disp(CMD);
  eval(CMD);
end

if isunix
  delete *.o
else
  delete *.obj
end

cd(old_dir);

disp('----------------------- DONE ----------------------------');
