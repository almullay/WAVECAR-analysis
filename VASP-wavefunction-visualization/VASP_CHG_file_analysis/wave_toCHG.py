import numpy as np
from decimal import Decimal
import sys

def getPchrg(phi):
	pchrg = np.real(np.multiply(phi,np.conj(phi)))
	return pchrg

def dec(x):
    s = "{0:.5f}".format(x)
    return s
    
phi = np.load(sys.argv[1])
chrg = getPchrg(phi)
xc,yc,zc = [len(chrg[:,0,0]),len(chrg[0,:,0]),len(chrg[0,0:])]
chrg=chrg/xc*yc*zc
all = []
for k in range(0,zc):
    for j in range(0,yc):
        for i in range(0,xc):
            all.append(dec(round(chrg[i,j,k],5)))
print(all[0])

with open("Wave.PARCHG",'w') as f:
	'''THIS HEADER WILL CHANGE DEPENDING ON THE SPECIFIC SYSTEM, BUT CAN JUST COPY AND PASTE TO REPLACE AFTER SCRIPT RUNS'''
    # Write lines with f.write("string")
    f.write('Si fcc\n   5.4310205052139802\n')     
    f.write('   0.500000    0.500000    0.000000\n')
    f.write('   0.000000    0.500000    0.500000\n')
    f.write('   0.500000    0.000000    0.500000\n')
    f.write('  Si\n')
    f.write('     2\n')
    f.write('Direct\n')
    f.write('  0.000000  0.000000  0.000000\n')
    f.write('  0.250000  0.250000  0.250000\n\n')
    f.write('  32   32   32\n')
    lines = []
    Nlines = (xc*yc*zc)//10 + 1# + 12
    for i in range(0,Nlines):
        lines.append([])
    for k in range(0,zc):
        for j in range(0,yc):
            for i in range(0,xc):
                ind = k*yc*xc + j*xc + i
                row = ind // 10
                col = ind % 10
                lines[row].append(all[ind]+'      ')
    for line in lines:
        line = ''.join(line)
        print(line)
        f.write(line+'\n')