% This is a small script to compile the necessary MEX files.

% Here is the opportunity to set some additional flags here that may be
% passed to the compiler. These flags tell the GCC compiler to use the ISO
% C99 standard, and to optimize the code as much as possible. Depending
% on the compiler you use to build the MEX shared library files, you may
% want to change these variables, or set them to the empty string ('').
cflags  = '-std=gnu99 -O3 -Os';
ldflags = '-s -O3 -Os';

% These are the files containing the main functions implemented in C. Note
% that not all these files are needed to compile each of the MEX files.
Rsrcdir   = '../R/varbvs/src/'
corefiles = {'C/doublevectormatlab.c '
	     'C/singlematrixmatlab.c ' 
	     [ Rsrcdir 'vectorops.c ' ]
	     [ Rsrcdir 'sigmoid.c '   ]
	     [ Rsrcdir 'varbvs.c '    ]
	     [ Rsrcdir 'varbvsbin.c ' ]};

% These are the commands to build the build the MEX shared library files.
options = sprintf(['-O -largeArrayDims -IC -I%s ' ...
		   'COPTIMFLAGS="%s" LDOPTIMFLAGS="%s" '],...
		   Rsrcdir,cflags,ldflags);
eval(['mex ',options,'C/var1matlab.c ',corefiles{1:3}]);
eval(['mex ',options,'C/diagsqmatlab.c ',corefiles{1:3}]);
eval(['mex ',options,'C/diagsqtmatlab.c ',corefiles{1:3}]);
eval(['mex ',options,'C/varbvsupdatematlab.c ',corefiles{1:5}]);
eval(['mex ',options,'C/varbvsupdatematlab_general.c ',corefiles{1:5}]);
eval(['mex ',options,'C/varbvsupdate_fast.c ',corefiles{1:5}]);
eval(['mex ',options,'C/varbvsbinupdatematlab.c ',corefiles{[1:4 6]}]);
fprintf('Compilation of MEX files is complete.\n');
