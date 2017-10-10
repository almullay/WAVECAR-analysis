import numpy as np
from decimal import Decimal
import sys
import matplotlib.pyplot as plt
from datetime import datetime
def find(number,Kp,lines):
	''' 
	USE: To find k-point, set number==1 and lines==f.readlines() of open file f. Take first element of output.
		 To find index of first line for band nu, set number==nu. Take second element of output.
	'''
	# Create kpt array even if not finding k-point.
	kpt = np.zeros([1,3])
	if number > 9:
		strno = str(number)
	else:
		strno =' '+str(number)
	c = int(0)
	# Start from line 3 to avoid the three other integers that could be 1 (spin, no. of kpts, no. bands)
	for i in range(3,len(lines)):
		line = lines[i]
		if '          %s' % strno in line:
			c += 1
			if c == Kp:
				IND = i
				print(IND)
				if number==1:
					# Get wavevector/k-point
					line = lines[IND-1]
					kpt[0,:] = [Decimal(line[1:13]) , Decimal(line[18:30]), Decimal(line[35:47])]
				break
			elif c == len(lines)-1:
				print("K-point out of range; not included in file.")
				sys.exit()
	return np.asmatrix(kpt), IND

def nwaves(IND,lines):
	# Get number of waves
	line = lines[IND]
	Nwaves = int(Decimal(line[18:24]))
	return Nwaves
	
def bases(lines):
	A = np.zeros([3,3])
	B = np.zeros([3,3])
	for i in range(6,9):
		line = lines[i-3]
		A[i-6,:] = [Decimal(line[2:13]) , Decimal(line[18:30]), Decimal(line[34:47])]
		line = lines[i]
		B[i-6,:] = [Decimal(line[2:13]) , Decimal(line[18:30]), Decimal(line[34:47])]
	A,B = np.asmatrix(A),np.asmatrix(B)
	return A, B

def getGCoeff(IND,Nwaves,lines):
	G = np.zeros([Nwaves, 3])
	Coeff = np.zeros([Nwaves,1], dtype=complex)
	for i in range(Nwaves):
		line = lines[i+IND+2]
		G[i] = [Decimal(line[4:6]), Decimal(line[10:12]), Decimal(line[16:18])]
		Coeff[i] = complex(Decimal(line[23:36]) , Decimal(line[40:53]))
	G=np.asmatrix(G)
	return G, Coeff

def makeWave(A,Nwaves,a0,partition):
	Tx,Ty,Tz = partition
	phi = np.zeros([Tx,Ty,Tz],dtype=complex)
	PW = np.zeros([Nwaves,1], dtype=complex)
	dr= a0/partition
	# For later:
	domx,domy,domz = dr[0]*np.arange(Tx),dr[1]*np.arange(Ty),dr[2]*np.arange(Tz)
	dr = np.asmatrix(dr)
	[A,G,Coeff] = A
	for i in range(0,Tx):
		for j in range(0,Ty):
				for k in range(0,Tz):
					r = np.transpose(np.multiply(dr,[i,j,k]))
					# Attempt for fcc basis... 
					#r = r*A
					PW = np.exp(1j*np.sum(np.matmul(G,r),axis=1,dtype=complex))
					PW = np.multiply(Coeff,PW)
					phi[i,j,k] = np.sum(PW,dtype=complex)
	return phi,domx,domy,domz

def parse_data(nu,Kp):
	with open("GCOEFF.txt",'r+') as f:
		lines = f.readlines()
		A, B = bases(lines)
		kpt=find(1,Kp,lines); kpt=kpt[0];
		print(kpt)
		IND = find(nu,Kp,lines); IND = IND[1];
		Nwaves = nwaves(IND,lines)
		print("Reciprocal cell: ", B,"Size of plane-wave basis: ",Nwaves)
		G,Coeff = getGCoeff(IND,Nwaves,lines)
		G = G+kpt
		return Nwaves,A,B,G,Coeff
		
def sort(G,Coeff):
	#KEEP OFF G = np.asmatrix(G)*np.asmatrix(B) 
	### Sort by coefficient size
	Cmag = np.absolute(Coeff)
	GC = np.array(np.concatenate((G,Cmag),axis=1))
	### When using argsort to sort whole array by a single row or column, you cannot use a matrix type.
	GCsort = GC[GC[:,3].argsort()]
	print('Sorted: \n',GCsort[0:2],'\n...\n',GCsort[-7::])
	coeff_mag = np.sum(np.multiply(Cmag,Cmag))
	print("\nSum of coeff's squared: ",coeff_mag)
	sys.exit()
	
def getPchrg(phi):
	pchrg = np.real(np.multiply(phi,np.conj(phi)))
	return pchrg

def plot(domain,A,upperbound,type,*title):	
	if type==0:
		print("\nPlotting charge density...")
		plt.plot(domain,A)
		plt.gca().set_ylim(bottom=0)
		plt.ylabel(r'$|\phi (x,0,0)|^2$')
	elif type==1:
		print("\nPlotting wave function...")
		plt.plot(domain,np.real(A),'b-')
		plt.plot(domain,np.imag(A),'r-')
		plt.ylabel(r'$\phi (x,0,0)$')
	else:
		print("Plot input error, not 0 or 1.")
	plt.xlim(0,upperbound)
	plt.xlabel('Length in Angstroms')
	plt.show()
	
def checkArgs(argv):
	try:
		if 'save'==argv[6]:
			want_save=True
		else:
			want_save=False
	except:
		want_save=False
	try:
		want_plot =('plot'==sys.argv[7])
	except:
		want_plot=False
	try:
		title = argv[8]
	except:
		title=None
	return title,want_plot,want_save
	
def main():
    # 
	#   Order of command line call arguments: 
	#   For plotting or storing phi -> band no. (int), k-point no. (int), no. points 
	#	    where no. points is three integers, i.e., 500 1 1; 32 32 32; etc.
	#	For sorting coefficient magnitudes - > band no. (int), k-point no. (int), * (doesn't matter what goes here), 'sort'
	#		Optional argument -> 'save'/any other string is 'nosave', want plot? ('plot'/anything else), plot title (string)
	#
	#
	startTime = datetime.now()
	if 'sort' in sys.argv:
		nu = int(sys.argv[1])
		Kp = int(sys.argv[2])
		Nwaves,A,B,G,Coeff = parse_data(nu,Kp)
		sort(G,Coeff)
		# Will exit system after.
	elif len(sys.argv) >= 6:
		nu = int(sys.argv[1])
		Kp = int(sys.argv[2])
		partition = [int(s) for s in sys.argv[3:6]]
		title, want_plot, want_save = checkArgs(sys.argv)
	else:
		sys.exit('\nInsufficient parameters detected. Please review script comments and use correct command line inputs.')
	
	Nwaves,A,B,G,Coeff = parse_data(nu,Kp)
	G = G*B
	a0 = 2*np.pi/B[0,0]
	phi,domx,domy,domz = makeWave([A,G,Coeff],Nwaves,a0,partition)
	pchrg = getPchrg(phi)
	if want_save:
		np.save('WaveFcn.%i.%i.%i.%i.%i'%(nu,Kp,partition[0],partition[1],
		partition[2]),phi)
	if want_plot:
		plot(domx,phi[:,0,0],a0,1,title)
	print datetime.now() - startTime
if __name__ == "__main__":
    main()