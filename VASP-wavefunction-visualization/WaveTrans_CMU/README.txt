http://www.andrew.cmu.edu/user/feenstra/wavetrans/
(Read through 1-page site before trying to use)

WaveTrans: Real-space wavefunctions from VASP WAVECAR file

R. M. Feenstra and M. Widom
Department of Physics, Carnegie Mellon University, Pittsburgh, PA 15213

Usage

The programs read in the binary file WAVECAR that is output from VASP.
For the FORTRAN 77 version (WaveTrans.f), the dimension of the various 
arrays in the program must be sufficiently large; if they are too small
then an error message results and the user must then increase the array 
size (using the PARAMETER statements in the source code) and recompile 
the program. For the FORTRAN 90 version (WaveTrans.f90), the arrays sizes 
are set automatically.

Yousif's Notes:
Wave.Trans.f90 actually works, just because it gets the dimension
right most of the time for the WAVECAR.
To compile, run 'gfortran file.f90' in a Cygwin interface that's
been installed with the following three packages from category
"Devel": gcc-fortran, gcc-g++, make. Then make sure gfortran
is available from the Cygwin command line.
This will create executables for you, but for them to execute
from your Windows interface you need to make sure C:\cygwin64\bin
is in your Windows Path (for access to cyggfortran-3.dll and 
cygwin1.dll). Then just call the executables that the Fortran files
compile into with the WAVECAR file in your directory.