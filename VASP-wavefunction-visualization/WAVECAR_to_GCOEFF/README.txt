WaveTrans.f was compiled with WAVECARin.f as a subroutine to create WaveTrans.exe. If you look at WaveTrans.f, you will notice there are some defined
variables near the beginning: 'nbdim' is how many bands there are in the calculation; 'npdim' is an upper bound of the largest plane-wave basis
used by VASP; 'nwdim' is the number of k-points used; 'nsdim' is the number of spin used by the calculation. (All of these can be found by the
script by parsing the WAVECAR.)
**If any of the values are too small, you will get an error. It's normally easy to either figure out the correct values from other VASP output 
files directly, or just randomly increasing them until the script runs. Note that the Fortran scripts have to be recompiled after each edit
with 'gfortran WaveTrans.f WAVECARin.f' (this compiles it with the subroutine).

The Matlab script read_wavecar.m was just a script Fahd and I built as we were understanding how to parse the binary WAVECAR file. It might be
useful as we try to do more advanced postprocessing and move away from CMU's Fortran scripts, so I left it in here. None of my python scripts
use it.