import numpy as np
from decimal import Decimal
import sys
import matplotlib.pyplot as plt

chg_file = sys.argv[1]
print('\nReading data from '+chg_file+'...')
with open(chg_file,'r+') as f:
	lines = f.readlines()
	
	line = lines[1]
	a = 5.43102050521398 #Change to whatever the cell length is. 5.431 was for Si bulk supercell.
	#print(a)
	A = np.zeros([3,3])
	for i in range(2,5):
		line = lines[i]
		A[:,i-2] = [Decimal(line[4:13]) , Decimal(line[16:25]), Decimal(line[28:37])]
	A = a*np.transpose(np.matrix(A))
	line = lines[118] #This index must be whichever line holds the x,y,z FFT dimensions in the CHG file.
	xc,yc,zc = [Decimal(line[2:6]) , Decimal(line[7:11]), Decimal(line[12:16])]
	xc,yc,zc = [int(xc),int(yc),int(zc)]
	vals = np.zeros([xc*yc*zc,1])
	ct=0
	for i in range(118+1,len(lines)): #Same as above mentioned index, plus one.
		line = lines[i]
		for j in range(0,10):
			if ct==xc*yc*zc:
				break
			ln_ind = 12*j + 1
			#print(i,ln_ind)
			vals[(i-12)*10 + j] = Decimal(line[ln_ind:ln_ind + 11])
			ct+=1
print(vals,np.shape(vals),[xc,yc,zc])
rhov = np.zeros([xc,yc,zc])
#"According to the VASP website, the NGX in the inner loop and NGZ is in the outer loop and NGY is in between."
for k in range(0,zc):
	for j in range(0,yc):
		for i in range(0,xc):
			rhov[i,j,k] = vals[k*yc*zc + j*zc + i]
print(np.sum(rhov)/(xc*yc*zc))
rhov = rhov/np.linalg.det(A) #det is triple product, doesn't matter if A is transposed
print(np.linalg.det(A))
#Maybe rho should be squared? Jingsong says no. So do vasp forums.
##rhov = np.multiply(rhov,rhov)
'''
fig = plt.figure()
dom = (a/xc)*np.arange(0,xc)

i=0
for s in [[rhov[:,0,0],311],[rhov[0,:,0],312],[rhov[0,0,:],313]]:
	ax = fig.add_subplot(s[1])
	ax.plot(dom,s[0])
	ax.set_xlim(0,a)
	if i==0:
		i+=1
		ax.set_title('Si Charge Density')
	ax.set_ylabel('Charge Density')
ax.set_xlabel('Real Space')
#plt.show()
#scipy.io.savemat('C:/Users/almul_000/Documents/MATLAB/Matlab/Research/ORNL/rhov.mat',mdict={'rhov': rhov})
'''
np.save('CHGCAR.chrg_density.%s'%chg_file[13:16],rhov)
'''
 ### 4D PLOTTING ATTEMPTS
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm, colors


# here are the x,y and respective z values
X, Y = np.meshgrid(dom, dom)

# here we create the surface plot, but pass V through a colormap
# to create a different color for each patch
norm = colors.Normalize()
chrg=rhov
chrg = norm(chrg)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.view_init(45,60)
#for i in range(0,T):	
ax.plot_surface(X, Y, chrg[:,:,0])
plt.show()'''