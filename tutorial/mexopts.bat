@echo off

set MATLAB=%MATLAB%
set MW_TARGET_ARCH=win64
set PATH=C:\PROGRA~1\mingw-w64\x86_64-4.9.0-posix-seh-rt_v3-rev2\mingw64\bin;%PATH%

set COMPILER=x86_64-w64-mingw32-g++
set COMPFLAGS=-c -m64 -Wall -std=c++11 -DMATLAB_MEX_FILE
set OPTIMFLAGS=-DNDEBUG -O3 -fexpensive-optimizations
set DEBUGFLAGS=-g
set NAME_OBJECT=-o

set LINKER=x86_64-w64-mingw32-g++
set LINKFLAGS=-shared -L"%MATLAB%\extern\lib\win64\microsoft" -L"%MATLAB%\bin\win64"
set LINKFLAGSPOST=-lmx -lmex -lmat
set LINKOPTIMFLAGS=-O3
set LINKDEBUGFLAGS=-g
set LINK_FILE=
set LINK_LIB=
set NAME_OUTPUT=-o "%OUTDIR%%MEX_NAME%%MEX_EXT%"