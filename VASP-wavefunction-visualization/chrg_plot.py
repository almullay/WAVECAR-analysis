import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import sys
import scipy.io

def getPchrg(phi):
	pchrg = np.real(np.multiply(phi,np.conj(phi)))
	return pchrg

try:
	n,m,l= int(sys.argv[3]),int(sys.argv[4]),int(sys.argv[5])
	pchrg = np.zeros((n,m,1))
	image = np.zeros((n,m))
except:
	'''ONLY FOR PLOTTING SUM OF CHARGE DENSITIES'''
	pchrg = np.zeros((100,100,100))
	image = np.zeros((100,100))
	pass
	
if 'sum' in sys.argv[1]:
	'''PLOTTING SUM OF CHARGE DENSITIES'''
	#for i in range(36,42):
	for i in [535,544,553,562,571,580]:
		phi = np.load('WaveFcn.5.%i.300.300.1.npy'%i)
		#phi = np.load('CHGCAR.rho.17.00%i.npy'% i)
		pchrg=getPchrg(phi) + pchrg
elif 'chgcar' in sys.argv[1]:
	'''PLOTTING PARCHG/CHGCAR DATA'''
	str = sys.argv[1]
	phi = np.load('CHGCAR.rho.%s.npy'%str[-7::])
	print('Plotting data from CHGCAR.rho.%s.npy'%str[-7::]+'...')
	pchrg =getPchrg(phi)
elif '.mat' in sys.argv[1]:
	'''PLOTTING FROM MATLAB'''
	str = sys.argv[1]
	phi = scipy.io.loadmat('C:/Users/almul/OneDrive/Documents/MATLAB/'+str)
	phi = phi['phi']
	print('Plotting data from '+str+'...')
	pchrg =getPchrg(phi)
elif 'diff' in sys.argv[1]:
	phi = np.load(sys.argv[1])
	pchrg=getPchrg(phi)
else:
	'''PLOTTING SPECIFIC CHARGE DENSITY'''
	phi = np.load('WaveFcn.%s.%s.%s.%s.%s.npy'% tuple(sys.argv[1:6]))
	print('Plotting data from WaveFcn.%s.%s.%s.%s.%s.npy'% tuple(sys.argv[1:6]))
	pchrg=getPchrg(phi)
	### Fix axes to put into true x,y,z
	### WONT WORK FOR (CURRENT)'SUM' RUN BC PARTITION=(300,300,1)
	#pchrg = np.transpose(pchrg,axes=[1,2,0]) 
	
try:
	if '011' in sys.argv:
		l = len(phi[:,0,0])
		for k in range(l):
			image[:,k] = pchrg[:,k,k]
		print('Plotting along 011 axis...')
	else:
		image = pchrg[:,:,0]
except:
	image = pchrg[:,:,0]
'''MULTIPLY IMAGE BY 4 TO ENLARGE PATTERNS'''
#image = np.concatenate((image,image),axis=1)
#image = np.concatenate((image,image),axis=0)
np.save('image'+sys.argv[1],image)
cax=plt.imshow(image, extent=(0, 100, 0, 100),
           interpolation='nearest', cmap=cm.spectral)
plt.xticks((),())
plt.yticks((),())
plt.colorbar(cax)
plt.show()
### 1D pchrg plot to compare to PARCHG files from vasp
#plt.plot(np.arange(len(image[0,:])),image[0,:])
#plt.ylim(ymin=0)
#plt.show()