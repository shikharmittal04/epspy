import healpy as hp
import numpy as np

nu_o = 150e6
path='/shikhar/point_sources/'		#Path where you would like to save and load from, the Tb's and beta's.
k=7			#Number of pixels in units of log_2(Npix).
nu=150e6		#frequency (in Hz) at which you want to compute the brightness temperature map

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

'''
The brightness temperatures were computed already, load them for a faster computation.  
'''


Tb_o_individual_save_name = path+'Tb_o_individual.npy'
Tb_o_save_name = path+'Tb_o.npy'
beta_save_name = path+'beta.npy'
	
print("Loading pre-computed Tb_o's and beta's ...\n")
Tb_o_individual = np.load(Tb_o_individual_save_name,allow_pickle=True)
Tb_o = np.load(Tb_o_save_name)
beta = np.load(beta_save_name,allow_pickle=True)

N_nu = np.size(nu)
Tb_nu_save_name = path+'Tb_nu.npy'
if N_nu==1:
	#If there is only one frequency at which you want to calculate Tb... 
	print('Now computing the Tb at frequency {:.2f} MHz ...'.format(nu/1e6))
	Tb_nu_final = np.zeros(Npix)
	for j in range(Npix):
		Tb_nu_final[j] = Tb_nu(Tb_o_individual[j],beta[j],nu)
	
	np.save(Tb_nu_save_name,Tb_nu_final)
	print('Done.\n File saved as',Tb_nu_save_name)

else:
	#If you want to compute Tb at multiple frequencies ...
	print('Now computing the Tb at multiple frequencies ...')
	Tb_nu_final = np.zeros((N_nu,Npix))	
	for i in range(N_nu):
		for j in range(Npix):
			Tb_nu_final[i,j] = Tb_nu(Tb_o[j],beta[j],nu[i])
	
	np.save(Tb_nu_save_name,Tb_nu_final)
	print('Done.\n File saved as',Tb_nu_save_name)
	print('It is an array of shape',np.shape(Tb_nu_final))
