The Windows Batch Files are just for the convenience of more quickly running the python scripts (don't have to type as much with them) and being able to call them globally.



CHRG_PLOT.PY SCRIPT-----------------
'chrg_plot.py' is a python script for plotting the 2D charge density maps (like in Koiller paper).
Arguments taken: band# k-point# xdim ydim zdim
where it must be that the wave function's whose charge density is to be plotted is stored in an *.npy array file in the form
'WaveFcn.band#.k-point#.xdim.ydim.zdim.npy' which is the output from my wave_calc.py scripts (Which calculates wave function arrays of arbitrary 
dimension from GCOEFF.txt file).
So if 'WaveFcn.132.1.100.100.1.npy' is in my directory with chrg_plot.py, I would run 'python chrg_plot.py 133 1 100 100 1'.

Default is to plot in x-y plane. To change this, go to lines 60 and 62 and alter 'image = pchrg[: ,: ,0]' to whatever cut you want out of pchrg. For example, 
pchrg[0, :, :] would give y-z plane charge density.
If you look at the script, you will find there are 5 or 6 other argument types possible. These probaably won't be useful to you, but feel free to ask me
if something catches interest and it's not clear how it works.


PLOT_WAVE.PY SCRIPT-----------------
This script was made to create nice wave function figures. It can load *.mat and *.npy arrays.
Argument: 'WaveFcn.band#.k-point#.xdim.ydim.zdim.npy' or *.mat file holding wave function.
I.e., if 'WaveFcn.132.1.100.100.1.npy' is in my directory, run 'python plot_wave.py WaveFcn.132.1.100.100.1.npy' to plot the wave function.



WAVE_CALC.PY SCRIPT-----------------
This script calculates wave functions from the GCOEFF.txt file (which is output from WaveTrans.exe).
Arguments: band# k-point# xdim ydim zdim 'save'
***If 'save' is not added the wave function array will not be saved.
Note this script is a non-threaded version of the script threaded_wave_calc.py. It might be better to just use the threaded one.
Argument inputs are the same.
