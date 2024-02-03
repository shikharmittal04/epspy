'''
The brightness temperatures were computed already. Based on those precomputed values we compute maps at a general frequency.  
'''

import healpy as hp
import numpy as np
from mpi4py import MPI

comm = MPI.COMM_WORLD
cpu_ind = comm.Get_rank()
Ncpu = comm.Get_size()


nu_o = 150e6
path='/home/hpcmitt1/rds/hpc-work/point_sources_data/'		#Path where you would like to save and load from, the Tb's and beta's.
k=7			#Number of pixels in units of log_2(Npix).
nu = np.arange(50,200)		#frequency (in Hz) at which you want to compute the brightness temperature map

Nside=2**k
Npix = hp.nside2npix(Nside)

def Tb_nu(Tb_o,beta,nu):
	'''
	Given the brightness temperature Tb_o (at reference frequency nu_o), spectral indices beta and frequency nu
	return the brightness temperature in K at requested frequency nu.
	Return value is one number.
	nu should be in Hz.
	Tb_o and beta should be of same dimensions if both are arrays of size more than 1.
	'''	
	return np.sum(Tb_o*(nu/nu_o)**-(beta))


Tb_o_save_name = path+'Tb_o.npy'
beta_save_name = path+'beta.npy'

slctd_Tb_o = np.load(Tb_o_save_name,allow_pickle=True)[cpu_ind*ppc:(cpu_ind+1)*ppc]
slctd_beta = np.load(beta_save_name,allow_pickle=True)[cpu_ind*ppc:(cpu_ind+1)*ppc]

if cpu_ind==0: print("Starting computation ...\n")
N_nu = np.size(nu)
Tb_nu_final = np.zeros((N_nu,Npix))	

ppc = int(Npix/Ncpu)	#pixels per cpu
for j in np.arange(cpu_ind*ppc,(cpu_ind+1)*ppc):
	for i in range(N_nu):
		Tb_nu_final[i,j] = Tb_nu(slctd_Tb_o[j-cpu_ind*ppc],slctd_beta[j-cpu_ind*ppc],nu[i])

del slctd_Tb_o
del slctd_beta

#An additional short loop is required if Npix/Ncpu is not an integer. We do the remaining pixels on rank 0.
if cpu_ind==0:
	slctd_Tb_o = np.load(Tb_o_save_name,allow_pickle=True)[Ncpu*ppc:Npix]
	slctd_beta = np.load(beta_save_name,allow_pickle=True)[Ncpu*ppc:Npix]
	for j in np.arange(Ncpu*ppc,Npix):
		for i in range(N_nu):
			Tb_nu_final[i,j] = Tb_nu(slctd_Tb_o[j-Ncpu*ppc],slctd_beta[j-Ncpu*ppc],nu[i])


#-------------------------------------------------------------------------------------
#Now all CPUs have done their jobs of calculating the Tb's and beta's. 
if cpu_ind!=0:
	'''
	I am a worker CPU. Sending my Tb's to master CPU.
	'''
	comm.send(Tb_nu_final, dest=0, tag=11)
else:
	'''
	I am the master CPU. Receiving all Tb's.
	'''
	print("Done. Receiving Tb's from all CPUs ...\n")
	for i in range(1,Ncpu):
		Tb_nu_final = Tb_nu_final + comm.recv(source=i, tag=11)

	Tb_nu_save_name = path+'Tb_nu.npy'
	np.save(Tb_nu_save_name,Tb_nu_final)
	
	print('Done.\n File saved as',Tb_nu_save_name)
	print('It is an array of shape',np.shape(Tb_nu_final))




