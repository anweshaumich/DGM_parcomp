setenv('MW_MINGW64_LOC','C:\TDM-GCC-64')
mex -setup
mex -setup C++
mex -outdir potts potts/sm.c
