'''
Foregrounds due to Unresolved Radio Sources (FURS)
Created by Shikhar Mittal

This code will save the power law index and brightness temperature contributed by each source. The outputs are going to be very large.
'''

import healpy as hp
import numpy as np
import random
import transformcl as tcl
from mpi4py import MPI
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

#-------------------------------------------------------------------------------------
#Some fixed numbers ...
kB = 1.38e-23
cE = 2.998e8

#The following numbers go into dn/dS
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
path='/home/hpcmitt1/rds/hpc-work/point-sources-data/'		#Path where you would like to save and load from, the Tb's and beta's.
k=7			#Number of pixels in units of log_2(Npix).
nu=150e6		#frequency (in Hz) at which you want to compute the brightness temperature map
#-------------------------------------------------------------------------------------

Nside=2**k
Npix = hp.nside2npix(Nside)


comm = MPI.COMM_WORLD
cpu_ind = comm.Get_rank()
Ncpu = comm.Get_size()

if(Ncpu==1):
    print('Error: you are generating brightness temperatures for the first time.')
    print("Run this code as, say, 'mpirun -n 4 python3 %s', where 4 specifies the number of CPUs." %(sys.argv[0]))
    sys.exit()

#-------------------------------------------------------------------------------------
#Find the number density distribution on the master CPU and share it with all CPUs.

low, upp = -6, -1	
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
Tb_o_individual is an array of arrays of unequal lengths, i.e.,
all of Tb_o_individual[0], Tb_o_individual[1], ..., Tb_o_individual[Npix] are arrays of different lengths.
The length of Tb_o[j] tells us the number of sources, say N_j, on the jth pixel and
Tb_o_individual[j][0], Tb_o_individual[j][1], ..., Tb_o_individual[j][N_j] are the temperatures (at ref. frequency) due to 0th, 1st,...(N_j)th source on the jth pixel.
'''	
ppc = int(Npix/Ncpu)	#pixels per cpu

Tb_o_local = np.zeros(ppc)
Tb_o_local_individual = np.zeros(ppc, dtype=object)
beta_local = np.zeros(ppc, dtype=object)

for j in np.arange(cpu_ind*ppc,(cpu_ind+1)*ppc):
    N = int(n_clus[j])	#no.of sources on jth pixel
    So_j = np.array(random.choices(S_space,weights=dndS_space/Ns_per_sr,k=N))	#select N flux densities for jth pixel, in Jy
    beta_local[j-cpu_ind*ppc] = np.random.normal(loc=beta_o,scale=sigma,size=N)		#select N spectral indices for jth pixel
		
    Tb_o_local_individual[j-cpu_ind*ppc] = 1e-26*So_j*cE**2/(2*kB*nu_o**2*Omega_pix)
    Tb_o_local[j-cpu_ind*ppc] = np.sum(Tb_o_local_individual[j-cpu_ind*ppc])

#An additional short loop is required if Npix/Ncpu is not an integer. Do the remaining pixels ('rm_pix') on rank 0.
rm_pix=Npix%Ncpu
if cpu_ind==0 and rm_pix!=0:
    Tb_o_remain_individual = np.zeros(rm_pix, dtype=object)
    Tb_o_remain = np.zeros(rm_pix)
    beta_remain = np.zeros(rm_pix, dtype=object)
    for j in np.arange(Ncpu*ppc,Npix):
        N = int(n_clus[j])  #no.of sources on jth pixel
        So_j = np.array(random.choices(S_space,weights=dndS_space/Ns_per_sr,k=N))   #select N flux densities for jth pixel
        beta_remain[j-Ncpu*ppc] = np.random.normal(loc=beta_o,scale=sigma,size=N)   #select N spectral indices for jth pixel
        Tb_o_remain_individual[j-Ncpu*ppc] = 1e-26*So_j*cE**2/(2*kB*nu_o**2*Omega_pix)
        Tb_o_remain[j-Ncpu*ppc] = np.sum(Tb_o_remain_individual[j-Ncpu*ppc])

#-------------------------------------------------------------------------------------
#Now all CPUs have done their jobs of calculating the Tb's and beta's. 
if cpu_ind!=0:
    '''
    I am a worker CPU. Sending my Tb's and beta's to master CPU.
    '''
    comm.send(Tb_o_local, dest=0, tag=11)
    comm.send(Tb_o_local_individual, dest=0, tag=13)
    comm.send(beta_local, dest=0, tag=29)
else:
    '''
    I am the master CPU. Receiving all Tb's and beta's.
    I will save the Tb's and beta's as numpy arrays (in format '.npy').
    '''
    print('Done.\n')
    Tb_o = Tb_o_local
    Tb_o_individual = Tb_o_local_individual
    beta = beta_local
    for i in range(1,Ncpu):
        Tb_o = np.concatenate((Tb_o,comm.recv(source=i, tag=11)))
        Tb_o_individual = np.concatenate((Tb_o_individual,comm.recv(source=i, tag=13)))
        beta = np.concatenate((beta,comm.recv(source=i, tag=29)))
		
    if rm_pix!=0:
        Tb_o = np.concatenate((Tb_o,Tb_o_remain))
        Tb_o_individual = np.concatenate((Tb_o_individual,Tb_o_remain_individual))
        beta = np.concatenate((beta,beta_remain))

    Tb_o_individual_save_name = path+'Tb_o_individual.npy'
    Tb_o_save_name = path+'Tb_o.npy'
    beta_save_name = path+'beta.npy'

    np.save(Tb_o_individual_save_name,Tb_o_individual)
    np.save(Tb_o_save_name,Tb_o)
    np.save(beta_save_name,beta)
	
    print('The brightness temperatures for each source individually have been saved into file:',Tb_o_individual_save_name)
    print('The pixel wise brightness temperatures have been saved into file:',Tb_o_save_name)
    print('The spectral indices for each source individually have been saved into file:',beta_save_name)


