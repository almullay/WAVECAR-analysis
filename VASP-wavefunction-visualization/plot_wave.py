from pylab import *
import sys
import scipy.io
import matplotlib.pyplot as plt

file = sys.argv[1]
if '.mat' in file:
	phi = scipy.io.loadmat(file)
	phi = phi['phi']
else:
	phi = np.load(file)
print(np.shape(phi))

'''
image = np.zeros((len(phi[:,0,0]),2))
for i in range(len(phi[:,0,0])):
	image[i] = phi[0,i,i]
phi=image
'''
phi = phi[:,0,0]
#last_pt = phi[0]; phi = np.append(phi,last_pt)

print(np.shape(phi))
Rphi, Imphi = np.real(phi), np.imag(phi)
chrg = np.real(np.multiply(phi,np.conj(phi)))
print('chg max: ', max(chrg), 'phi max: ', max(Rphi),max(Imphi))
n_pts = phi.shape[0]

length = 31.4
inc = length/n_pts
x = np.arange(0, length, inc)
fig = figure()
ax = fig.add_subplot(111)

plt.rcParams['axes.linewidth'] = 100

ax.plot(x, Rphi, color='blue', linewidth=3, alpha=0.3,linestyle='-', label='linestyle="_"')
ax.plot(x, Imphi, color='red', linewidth=3, alpha=0.3, linestyle='-', label='lines tyle="_"')
ax.plot(x, chrg, color='green', linewidth=3, alpha=1, linestyle='-', label='lines tyle="_"')
ax.plot(x, np.zeros(n_pts), color='black', linewidth=1, linestyle='--', label='linestyle="--"')

min_value = 0
max_value = length

font_size = 24
ax.set_xlim(min_value, max_value)
xticks([0,max_value],['Si Core','Si Core'],fontsize=font_size, fontweight='bold')
yticks([],[]) # REMOVES Y AXIS MARKERS
ax.spines['top'].set_visible(False)
ax.spines['right'].set_linewidth(3)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_linewidth(3)
legend = legend([r"$Re(\phi)$",r"$Im(\phi)$",r"$|\phi|^2$"], loc=4);
frame = legend.get_frame()
frame.set_facecolor('1.0')
frame.set_edgecolor('1.0')

title(file[0:-3]+'Plot',fontsize=font_size, fontweight='bold')
ylabel(r'$\bf{\phi(x,0,0)}$',fontsize=font_size+6, fontweight='bold')

plotname="Plot."+file[0:-3]
if 'save' in sys.argv:
	savefig(plotname)
show()
