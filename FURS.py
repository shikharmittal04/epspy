#Foregrounds due to Unresolved Radio Sources (FURS)
#Created by Shikhar Mittal

import healpy as hp
import numpy as np
import random
import matplotlib.pyplot as plt
import transformcl as tcl
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import colormaps
from mpi4py import MPI
import os
import sys

def dndS(S):
	'''
	This is the distribution of flux density
	Return value is in number of sources per unit solid angle per unit flux density
	S is in units of Jy (jansky)
	'''    
	return S**-2.5*((A1*S**a1+B1*S**b1)**-1+(A2*S**a2+B2*S**b2)**-1)

def C(chi,A=7.8e-3,gam = 0.821):
	'''
	This is the correlation function from Rana & Bagla (2019).
	chi should be in radians.
	'''
	return A*(chi*180/np.pi)**(-gam)

def Tb_nu(Tb_o,beta,nu):
	'''
	Given the brightness temperature Tb_o (at reference frequency nu_o), spectral indices beta and frequency nu
	return the brightness temperature in K at requested frequency nu.
	Return value is one number.
	nu should be in Hz.
	Tb_o and beta should be of same dimensions if both are arrays of size more than 1.
	'''	
	return np.sum(Tb_o*(nu/nu_o)**-(beta))

#-------------------------------------------------------------------------------------
#Some fixed numbers ...
kB = 1.38e-23
cE = 2.998e8

#The following number sgo into dn/dS
a1,b1,a2,b2 = -0.854, 0.37, -0.856, 1.47
A1, B1 = 1.65e-4, 1.14e-4
A2A1, B2B1 = 0.24, 1.8e7
A2 = A2A1*A1
B2 = B2B1*B1

nu_o = 150e6	#Reference frequency in Hz
beta_o = 2.681	#Mean spectral index for extragalactic point sources
sigma = 0.5		#Spread in the alpha values
#-------------------------------------------------------------------------------------
'''
The following 3 numbers are required from the user.
They will be probably enter at the starting of the pipeline.
'''
path='/user1/shikhar/point_sources/'		#Path where you would like to save and load from, the Tb's and beta's.
k=7			#Number of pixels in units of log_2(Npix).
nu=50e6		#frequency (in Hz) at which you want to compute the brightness temperature map
#-------------------------------------------------------------------------------------

Tb_o_save_name = path+'Tb_o.npy'
beta_save_name = path+'beta.npy'

Nside=2**k
Npix = hp.nside2npix(Nside)

if os.path.isfile(Tb_o_save_name) and os.path.isfile(beta_save_name):
	'''
	The brightness temperatures were computed already, load them for a faster computation.  
	'''
	print("Loading pre-computed Tb_o's and beta's ...\n")
	Tb_o = np.load(Tb_o_save_name,allow_pickle=True)
	beta = np.load(beta_save_name,allow_pickle=True)
	
	N_nu = np.size(nu)
	Tb_nu_save_name = path+'Tb_nu.npy'
	if N_nu==1:
		print('Now computing the Tb at frequency {:.2f} MHz ...'.format(nu/1e6))
		Tb_nu_final = np.zeros(Npix)
		for j in range(Npix):
			Tb_nu_final[j] = Tb_nu(Tb_o[j],beta[j],nu)
		
		np.save(Tb_nu_save_name,Tb_nu_final)
		print('Done.\n File saved as',Tb_nu_save_name)

		#Since it's only one frequency, I am making a plot as well.		
		plt.rc('text', usetex=True)
		plt.rc('font', family='serif')
		hp.mollview(Tb_nu_final,title=None,unit=r'$T_{\mathrm{b}}(\nu)\,$(K)',cmap=colormaps['coolwarm'],norm='log')
		hp.graticule()
		fig_path = path+'Tb_nu.pdf'
		plt.savefig(fig_path, bbox_inches='tight')
		print('\nTb map saved as',fig_path)

	else:
		print('Now computing the Tb at multiple frequencies ...')
		Tb_nu_final = np.zeros((N_nu,Npix))	
		for i in range(N_nu):
			for j in range(Npix):
				Tb_nu_final[i,j] = Tb_nu(Tb_o[j],beta[j],nu[i])
		
		np.save(Tb_nu_save_name,Tb_nu_final)
		print('Done.\n File saved as',Tb_nu_save_name)
		print('It is an array of shape',np.shape(Tb_nu_final))
else:
	'''
	Code is being run for the first time. So we will save the brightness temperature contributed by each source.
	'''
	comm = MPI.COMM_WORLD
	cpu_ind = comm.Get_rank()
	Ncpu = comm.Get_size()
	
	if(Ncpu==1):
		print('Error: you are generating brightness temperatures for the first time.')
		print("Run this code as, say, 'mpirun -n 4 python3 %s', where 4 specifies the number of CPUs." %(sys.argv[0]))
		sys.exit()

	#-------------------------------------------------------------------------------------
	#Find the number density distribution on the master CPU and share it with all CPUs.

	low, upp = -2, 0	
	S_space = np.logspace(low,upp,1000)
	dndS_space = dndS(S_space)
	Ns_per_sr = np.trapz(dndS_space,S_space)
	Ns = 4*np.pi*Ns_per_sr
	Omega_pix = hp.nside2pixarea(Nside) #Solid angle per pixel
	nbar = Ns/Npix


	if cpu_ind==0:
		'''
		Find the number density distribution on the master CPU.
		'''
		print('\nTotal number of sources = {:d}'.format(int(Ns)))
		print('Total number of pixels, Npix =',Npix)
		print('Average number of sources per pixel = {:.2f}'.format(nbar))

		print('\nNow finding the clustered number density distribution ...')
		#corrtocl requires us to sample the angle at some specific points. Obtained by 'tcl.theta'
		th = tcl.theta(1000)
		cor = C(th)
		Cl_clus = tcl.corrtocl(cor)

		#Now calculating the clustered map fluctuation...
		del_clus = hp.synfast(Cl_clus,Nside)

		#and the corresponding number density given the fluctuation...
		n_clus = nbar*(1+del_clus)
		print('Done.\nAverage overdensity for the clustered sky (should be close to 0) =',np.mean(del_clus),'\n')

		n_clus_save_name = path+'n_clus.npy'
		np.save(n_clus_save_name,n_clus)
		print('The clustered number density has been saved into file:',n_clus_save_name)

		print("\nNow computing Tb_o's and beta's ...")
	else:
		n_clus = None

	n_clus = comm.bcast(n_clus, root=0)	#Now all CPUs have the same number density distribution.
	#-------------------------------------------------------------------------------------	
	'''
	Tb_o is an array of arrays of unequal lengths, i.e.,
	all of Tb_o[0], Tb_o[1], ..., Tb_o[Npix] are arrays of different lengths.
	The length of each Tb_o[j] tells us the number of sources, say N_j, on the jth pixel and
	Tb_o[j][0], Tb_o[j][1], ..., Tb_o[j][N_j] are the temperatures (at ref. frequency) due to 0th, 1st, ... and (N_j)th source on the jth pixel.
	'''	
	Tb_o = np.zeros(Npix, dtype=object)
	beta = np.zeros(Npix, dtype=object)
	for j in range(Npix):
		if (cpu_ind == int(j/int(Npix/Ncpu))%Ncpu):
			N = int(n_clus[j])	#no.of sources on jth pixel
			So_j = np.array(random.choices(S_space,weights=dndS_space/Ns_per_sr,k=N))	#select N flux densities for jth pixel
			beta_j = np.random.normal(loc=beta_o,scale=sigma,size=N)		#select N spectral indices for jth pixel
			
			Tb_o_j = 1e-26*So_j*cE**2/(2*kB*nu_o**2*Omega_pix)
			Tb_o[j] = Tb_o_j
			beta[j] = beta_j
		
	if cpu_ind!=0:
		'''
		I am a worker CPU. Sending my Tb's and beta's to master CPU.
		'''
		comm.send(Tb_o, dest=0, tag=13)
		comm.send(beta, dest=0, tag=29)
	else:
		'''
		I am the master CPU. Receiving all Tb's and beta's.
		I will save the Tb's and beta's as numpy arrays (in format '.npy').
		'''
		print('Done.\n')
		for i in range(1,Ncpu):
			Tb_o = Tb_o + comm.recv(source=i, tag=13)
			beta = beta + comm.recv(source=i, tag=29)
			
		
		Tb_o_save_name = path+'Tb_o.npy'
		beta_save_name = path+'beta.npy'
		np.save(Tb_o_save_name,Tb_o)
		np.save(beta_save_name,beta)
		print('The brightness temperatures have been saved into file:',Tb_o_save_name)
		print('The spectral indices have been saved into file:',beta_save_name)


