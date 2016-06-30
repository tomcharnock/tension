import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Ellipse
import scipy.stats as ss
import scipy.integrate as si

def read_shit(filename):
	means = np.zeros(5)
	with open('data/'+filename+'.margestats') as f:
		for line in f:
			l = [j for j in line.split(' ') if j != '']
			if l[0] == 'omegabh2':
				means[0] = float(l[1])
			if l[0] == 'omegach2':
				means[1] = float(l[1])
			if l[0] == 'theta':
				means[2] = float(l[1])
			if l[0] == 'logA':
				means[3] = float(l[1])
			if l[0] == 'ns':
				means[4] = float(l[1])
	covmat = np.zeros([5,5])
	with open('data/'+filename+'.covmat') as f:
		i = 0
		for line in f:
			if i is 0:
				covnum = np.zeros(5)
				l = [j for j in line.split(' ') if j != '']
				for j in xrange(len(l)):
					if l[j] == 'omegabh2':
						covnum[0] = float(j)
					if l[j] == 'omegach2':
						covnum[1] = float(j)
					if l[j] == 'theta':
						covnum[2] = float(j)
					if l[j] == 'logA':
						covnum[3] = float(j)
					if l[j] == 'ns':
						covnum[4] = float(j)
			covnum = [int(j) for j in covnum]
			if i in covnum:
				indices = []
				l = [float(j.replace('\n','')) for j in line.split(' ') if j != '']
				for j in covnum:
					indices.append(l[j-1])
				if i is covnum[0]:
					covmat[0,:] = np.array(indices)
				if i is covnum[1]:
					covmat[1,:] = np.array(indices)
				if i is covnum[2]:
					covmat[2,:] = np.array(indices)
				if i is covnum[3]:
					covmat[3,:] = np.array(indices)
				if i is covnum[4]:
					covmat[4,:] = np.array(indices)
			i += 1
	return means,covmat
		
def plot_cov_ellipse(cov, pos, lims, labels, nstd=2, ax=None, **kwargs):

        def eigsorted(cov):
        	vals, vecs = np.linalg.eigh(cov)
        	order = vals.argsort()[::-1]
        	return vals[order], vecs[:,order]

	if ax is None:
		makeaxis = True
		ax = []
		gs = gridspec.GridSpec(len(pos)-1,len(pos)-1)
		gs.update(hspace=0,wspace=0,left=.15,right=.96)
	else:
		makeaxis = False

	k = 0
	for i in xrange(len(pos)):
		for j in xrange(i,len(pos)):
			if i is not j:
				if makeaxis:
					ax.append(plt.subplot(gs[i,j-1]))
					ax[k].set_ylabel(labels[i],rotation=0,labelpad=20)
					plt.xticks(np.arange(lims[j][0],lims[j][1]+(lims[j][1]-lims[j][0]),(lims[j][1]-lims[j][0])/2))
					plt.yticks(np.arange(lims[i][0],lims[i][1]+(lims[i][1]-lims[i][0]),(lims[i][1]-lims[i][0])/2))
					ax[k].set_xlim(lims[j])
					ax[k].set_ylim(lims[i])
					if i > 0:
						yticks = ax[k].yaxis.get_major_ticks()
						yticks[-2].label1.set_visible(False)
					if i < len(pos)-2:
						xticks = ax[k].xaxis.get_major_ticks()
						xticks[-2].label1.set_visible(False)
					if i == len(pos)-2:
						ax[k].set_xlabel(labels[j])
					if i is not j-1:
						ax[k].set_xticklabels([])
						ax[k].set_yticklabels([])
						ax[k].set_xlabel('')
						ax[k].set_ylabel('')
				mean = np.array([pos[j],pos[i]])
				covmat = np.array([[cov[j][j],cov[j][i]],[cov[j][i],cov[i][i]]])
    				vals, vecs = eigsorted(covmat)
    				theta = np.degrees(np.arctan2(*vecs[:,0][::-1]))
    				width, height = 2 * nstd * np.sqrt(vals)
    				ellip = Ellipse(xy=mean, width=width, height=height, angle=theta, **kwargs)
    				ax[k].add_artist(ellip)
				k += 1
	return ax

def surprise_update(P1_mean,P1_cov,P2_mean,P2_cov):
	return 0.5*(np.dot((P1_mean-P2_mean),np.dot(np.linalg.inv(P1_cov),(P1_mean-P2_mean)))-np.trace(np.identity(len(P1_mean))-np.dot(P2_cov,np.linalg.inv(P1_cov))))

def surprise_independent(P1_mean,P1_cov,P2_mean,P2_cov):
	return 0.5*(np.dot((P1_mean-P2_mean),np.dot(np.linalg.inv(P1_cov),(P1_mean-P2_mean)))-np.trace(np.identity(len(P1_mean))+np.dot(P2_cov,np.linalg.inv(P1_cov))))

P1_mean,P1_cov = read_shit('CMB/CMB')
P2_mean,P2_cov = read_shit('Strong/Strong_L')
lims = []
labels = ['$\Omega_bh^2$','$\Omega_ch^2$','$\Theta_{MC}$','$\log A$','$n_s$']
domainsize = 1
for i in xrange(len(P1_mean)):
	lower = min(P1_mean[i]-5*np.sqrt(P1_cov[i][i]),P2_mean[i]-5*np.sqrt(P2_cov[i][i]))
	upper = max(P1_mean[i]+5*np.sqrt(P1_cov[i][i]),P2_mean[i]+5*np.sqrt(P2_cov[i][i]))
	domainsize = domainsize*(upper-lower)
	lims.append([lower,upper])

P1 = ss.multivariate_normal(P1_mean,P1_cov)
P2 = ss.multivariate_normal(P2_mean,P2_cov)

print 'Surprise D1 compares D2 = ',surprise_independent(P1_mean,P1_cov,P2_mean,P2_cov)
print 'Surprise D2 compares D1 = ',surprise_independent(P2_mean,P2_cov,P1_mean,P1_cov)

nmc = len(P1_mean)*5e6
means = (np.random.uniform(size=int(nmc)).reshape((int(nmc/len(P1_mean)),len(P1_mean))))
a = 0
b = 0
l = 0
m = 0

old_settings = np.geterr()

for i in xrange(len(means)):
	for j in xrange(len(means[0])):
		means[i][j] = means[i][j]*(lims[j][1]-lims[j][0])+lims[j][0]
	
	np.seterr(**old_settings)
	P1_temp = P1.pdf(means[i])
	P2_temp = P2.pdf(means[i])
	np.seterr(all='raise')
	try:
		a_temp = P2_temp*np.log(P2_temp/P1_temp)
		l_temp = 0
	except:
		a_temp = 0.0
		l_temp = 1
		pass

	a += a_temp
	l += l_temp

	try:
		b_temp = P1_temp*np.log(P1_temp/P2_temp)
		m_temp = 0
		pass
	except:
		b_temp = 0.0
		m_temp = 1
		pass
	
	b += b_temp
	m += m_temp
		
print l,m
if l != nmc/len(P1_mean):
	print 'D (D1 compares to D2) = ', a*domainsize/(nmc/len(P1_mean)-l)
else:
	print '(D1 compares to D2) all values are zero'
if m != nmc/len(P1_mean):
	print 'D (D2 compares to D1) = ', b*domainsize/(nmc/len(P1_mean)-m)
else:
	print '(D2 compares to D1) all values are zero'

axis = plot_cov_ellipse(P1_cov,P1_mean,lims,labels,nstd=2,alpha=0.5,color='blue')
plot_cov_ellipse(P2_cov,P2_mean,lims,labels,nstd=2,ax=axis,alpha=0.5,color='red')
plt.show()
