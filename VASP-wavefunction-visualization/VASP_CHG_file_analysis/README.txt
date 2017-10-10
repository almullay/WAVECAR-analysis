The wave_toCHG.py file takes as input a python array that is holding wave function data. 
Normally this array would be titled 'WaveFcn.band#.kpt#.xdim.ydim.zdim.npy' as output
from either wave_calc.py or threaded_wave_calc.py (python scripts that calculate wave
function at certain band and k-point and x,y,z partition from GCOEFF.txt file, which is 
the output of a WaveTrans executable from CMU [holds G-vector indices and corresponding
coefficients]).
The output of wave_toCHG.py is the mod squared of the wave function (i.e., the charge 
density) in the form of a CHG / PARCHG file from VASP. The reason for doing this is that it is
viewable with Vesta, VMD, etc. (Though you have to rename it to *.CHGCAR to do this, even 
though it is actually in the form of the CHG file.)

read_CHG.py is essentially the (partial) inverse function of wave_toCHG.py; it takes as 
input a CHG or PARCHG file and outputs a python array holding the charge density values
at each spatial location specific by VASP in the file. This is useful for looking at 
arbitrary 2D cuts of charge density (using the chrg_plot.py script to view them).
Use format is: 'python read_CHG.py --PARCHG or CHG file name--' 
I.e., if I have a PARCHG file at the 7th band and 3rd kpoint, it would be called 
PARCHG.0007.0003 by VASP, and so I'd copy it to my directory and write
'python read_CHG.py PARCHG.0007.0003'.