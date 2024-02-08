import healpy as hp
import numpy as np
from scipy.interpolate import CubicSpline
import random
import transformcl as tcl
import sys
from mpi4py import MPI

#Some fixed numbers ...
kB = 1.38e-23   #Boltzmann constant in J/K units
cE = 2.998e8    #Speed of light in m/s
Tcmb_o = 2.725  #CMB temperature today in K

class extragalactic():
    def __init__(self, log2Npix=6, low=-6,upp=-1, nu_o=150e6, beta_o=2.681,sigma_beta=0.5, A=7.8e-3,gam=0.821, path=''):
        self.nu_o = nu_o        #Reference frequency in Hz

        self.beta_o = beta_o    #Mean spectral index for extragalactic point sources
        self.sigma_beta = sigma_beta   #Spread in the beta values

        self.A = A      #Amplitude of the power-law 2-point angular correlation function
        self.gam = gam  #-exponent of the power-law 2-point angular correlation function

        self.low = low  #log_10(S_min), where S_max is in Jy
        self.upp = upp  #log_10(S_max)
        
        self.path = path        #Path where you would like to save and load from, the Tb's and beta's
                            
        self.log2Nside = log2Nside    #Number of divisions in units of log_2
        global Nside, Npix
        Nside=2**log2Nside
        Npix = hp.nside2npix(Nside) #Actual number of pixels
    #End of function __init__()
	
	def num_sources():
		S_space = np.logspace(self.low,self.upp,1000)
        dndS_space = dndS(S_space)

		Ns_per_sr = np.trapz(dndS_space,S_space)
        Ns = 4*np.pi*Ns_per_sr
		print('Number of sources =',Ns)
		return Ns

    def ref_freq(self):
        '''
        Tb_o_individual is an array of arrays of unequal lengths, i.e.,
        all of Tb_o_individual[0], Tb_o_individual[1], ..., Tb_o_individual[Npix] are arrays of different lengths.
        The length of Tb_o[j] tells us the number of sources, say N_j, on the jth pixel and
        Tb_o_individual[j][0], Tb_o_individual[j][1], ..., Tb_o_individual[j][N_j] are the temperatures (at ref. frequency) due to 0th, 1st,...(N_j)th source on the jth pixel.
        '''
        def dndS(S):
            '''
            This is the distribution of flux density.
            Return value is in number of sources per unit solid angle per unit flux density.
            S is in units of Jy (jansky).
            I have taken the functional form and the numbers from Gervasi et al (2008) ApJ.
            '''
            a1,b1,a2,b2 = -0.854, 0.37, -0.856, 1.47
            A1, B1 = 1.65e-4, 1.14e-4
            A2A1, B2B1 = 0.24, 1.8e7
            A2 = A2A1*A1
            B2 = B2B1*B1
            return S**-2.5*((A1*S**a1+B1*S**b1)**-1+(A2*S**a2+B2*S**b2)**-1)

        def C(chi):
            '''
            This is the correlation function from Rana & Bagla (2019).
            chi should be in radians.
            '''
            return self.A*(chi*180/np.pi)**(-self.gam)
        
        #-------------------------------------------------------------------------------------
        comm = MPI.COMM_WORLD
        cpu_ind = comm.Get_rank()
        Ncpu = comm.Get_size()

        if Ncpu==1: print("Best to run this code in parallel. Run as 'mpirun -np 4 python3 %s', where 4 specifies the number of CPUs." %(sys.argv[0]))
            
        #-------------------------------------------------------------------------------------
        #Find the number density distribution on the master CPU and share it with all CPUs.
        
        S_space = np.logspace(self.low,self.upp,1000)
        dndS_space = dndS(S_space)
        
        Omega_pix = hp.nside2pixarea(Nside) #Solid angle per pixel

        if cpu_ind==0:
            '''
            Find the number density distribution on the master CPU.
            '''
			Ns_per_sr = np.trapz(dndS_space,S_space)
        	Ns = 4*np.pi*Ns_per_sr
			nbar = Ns/Npix

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
            print('Done.\nAverage overdensity for the clustered sky (should be close to 0) = {:.3f}'.format(np.mean(del_clus)))

            n_clus_save_name = self.path+'n_clus.npy'
            np.save(n_clus_save_name,n_clus)
            print('The clustered number density has been saved into file:',n_clus_save_name)

            print("\nNow computing Tb_o's and beta's ...")
        else:
            n_clus = None

        n_clus = comm.bcast(n_clus, root=0) #Now all CPUs have the same number density distribution.
        #-------------------------------------------------------------------------------------
        #Visit each pixel on the sky, for each source on that pixel, assign fluxes and spectral indices to them.
        
        ppc = int(Npix/Ncpu)	#pixels per cpu

        Tb_o_local = np.zeros(ppc)
        Tb_o_local_individual = np.zeros(ppc, dtype=object)
        beta_local = np.zeros(ppc, dtype=object)

        for j in np.arange(cpu_ind*ppc,(cpu_ind+1)*ppc):
            N = int(n_clus[j])	#no.of sources on jth pixel
            So_j = np.array(random.choices(S_space,weights=dndS_space/Ns_per_sr,k=N))   #select N flux densities for jth pixel, in Jy
            beta_local[j-cpu_ind*ppc] = np.random.normal(loc=self.beta_o,scale=self.sigma_beta,size=N) #select N spectral indices for jth pixel
                        
            Tb_o_local_individual[j-cpu_ind*ppc] = 1e-26*So_j*cE**2/(2*kB*self.nu_o**2*Omega_pix)
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
                beta_remain[j-Ncpu*ppc] = np.random.normal(loc=self.beta_o,scale=self.sigma_beta,size=N)   #select N spectral indices for jth pixel
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

            Tb_o_individual_save_name = self.path+'Tb_o_individual.npy'
            Tb_o_save_name = self.path+'Tb_o.npy'
            beta_save_name = self.path+'beta.npy'

            np.save(Tb_o_individual_save_name,Tb_o_individual)
            np.save(Tb_o_save_name,Tb_o)
            np.save(beta_save_name,beta)
                
            print('The brightness temperatures for each source individually have been saved into file:',Tb_o_individual_save_name)
            print('The pixel wise brightness temperatures have been saved into file:',Tb_o_save_name)
            print('The spectral indices for each source individually have been saved into file:',beta_save_name)
        return None
    #End of function ref_freq()


    def gen_freq(self, nu=1e6*np.arange(50,200)):
        '''
        If you are running this function you must have run ref_freq().    
        This function computes the map(s) at general frequency(ies) based on the precomputed values from ref_freq().
        nu is the frequency (in Hz) at which you want to evaluate the brightness temperature map.
        nu can be one number or an array.
        '''

        comm = MPI.COMM_WORLD
        cpu_ind = comm.Get_rank()
        Ncpu = comm.Get_size()

        #making nu_eval global so that it can be used in the function plotter
        global nu_glob
        nu_glob = nu
        
        def Tb_nu(Tb_ref,beta,nu):
            '''
            Given the brightness temperature Tb_ref (at reference frequency nu_o), spectral indices beta and frequency nu
            return the brightness temperature in K at requested frequency nu.
            Return value is one number.
            nu should be in Hz. (Best to give one nu at a time.)
            Tb_o and beta should be of same dimensions if both are arrays of size more than 1.
            '''	
            return np.sum(Tb_ref*(nu/self.nu_o)**-(beta))


        Tb_o_individual_save_name = self.path+'Tb_o_individual.npy'
        beta_save_name = self.path+'beta.npy'

        if cpu_ind==0: print("\nStarting computation ...\n")
        N_nu = np.size(nu)
        Tb_nu_final = np.zeros((Npix,N_nu))	

        ppc = int(Npix/Ncpu)    #pixels per cpu
        slctd_Tb_o = np.load(Tb_o_individual_save_name,allow_pickle=True)[cpu_ind*ppc:(cpu_ind+1)*ppc]
        slctd_beta = np.load(beta_save_name,allow_pickle=True)[cpu_ind*ppc:(cpu_ind+1)*ppc]

        for j in np.arange(cpu_ind*ppc,(cpu_ind+1)*ppc):
            for i in range(N_nu):
                Tb_nu_final[j,i] = Tb_nu(slctd_Tb_o[j-cpu_ind*ppc],slctd_beta[j-cpu_ind*ppc],nu[i])

        del slctd_Tb_o
        del slctd_beta

        #An additional short loop is required if Npix/Ncpu is not an integer. We do the remaining pixels on rank 0.
        if cpu_ind==0 and Npix%Ncpu!=0:
            slctd_Tb_o = np.load(Tb_o_individual_save_name,allow_pickle=True)[Ncpu*ppc:Npix]
            slctd_beta = np.load(beta_save_name,allow_pickle=True)[Ncpu*ppc:Npix]
            for j in np.arange(Ncpu*ppc,Npix):
                for i in range(N_nu):
                    Tb_nu_final[j,i] = Tb_nu(slctd_Tb_o[j-Ncpu*ppc],slctd_beta[j-Ncpu*ppc],nu[i])


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

            Tb_nu_save_name = self.path+'Tb_nu.npy'
            np.save(Tb_nu_save_name,Tb_nu_final)
            
            print('Done.\nFile saved as',Tb_nu_save_name)
            print('It is an array of shape',np.shape(Tb_nu_final))

        return Tb_nu_final
    #End of function gen_freq()

    def visual(self, nu=None, skymap=False, spectrum=True, xlog=False,ylog=True):
        '''
        Use this function to plotting a sky map at a given freqeuncy ('skymap') and/or
        the global extragalactic foregrounds as a function of frequency ('spectrum').
        'nu' is required only for making the sky map. It should be one number in Hz.
        By default we only plot the spectrum and not the skymap.
        'xlog' and 'ylog' are the boolean values deciding the scale of x and y axis, respectively.
        '''
        Tb_o = np.load(self.path+'Tb_o.npy')
        Tb_nu = np.load(self.path+'Tb_nu.npy')
        Tb_o_mean = np.mean(Tb_o)
        
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')
        if skymap:    
            if np.size(nu)==1:
                if nu==None:
                    print("No frequency given with 'skymap=True'. Creating sky map at the reference frequency ...")
                    nu=self.nu_o
                    Tb_plot = Tb_o
                else:
                    ind = np.where(nu_glob==nu)
                    if ind==None:
                        print('Given frequency unavailable in gen_freq(). Interpolating ...')
                        if nu<np.min(nu_glob):
                            print('Warning! Given frequency outside the range. Using the lowest available frequency; {:.2f} MHz ...'.format(np.min(nu_glob)/1e6))
                            Tb_plot = Tb_nu[:,0]
                        elif nu>np.max(nu_glob):
                            print('Warning! Given frequency outside the range. Using the highest available frequency; {:.2f} MHz ...'.format(np.max(nu_glob)/1e6))
                            Tb_plot = Tb_nu[:,0]
                        else:
                            print("Creating sky map at {:.2f} ...".format(nu/1e6))
                            spl = CubicSpline(nu_glob, Tb_nu)
                            Tb_plot = spl(nu)
                    else:
                        print("Creating sky map at {:.2f} ...".format(nu/1e6))
                        Tb_plot = Tb_nu[:,ind]

                print('\nGenerating the sky map at frequency = {:.2f} MHz ...'.format(nu/1e6))
            else:
                print("Warning! Multiple frequencies given with 'skymap=True'. Plotting only at the reference frequency ...")
                Tb_plot = Tb_o

            hp.mollview(Tb_plot,title=None,unit=r'$T_{\mathrm{b}}^{\mathrm{eg}}\,$(K)',cmap=colormaps['coolwarm'],min=0.05,max=200,norm='log')
            hp.graticule()
            fig_path = self.path+'Tb_nu_map_'+str(int(nu/1e6))+'-MHz.pdf'
            plt.savefig(fig_path, bbox_inches='tight')
            print('Done. Tb map saved as',fig_path)

        if spectrum:
            nu=nu_glob
                    
            Tb_mean = Tb_o_mean*(nu/self.nu_o)**-self.beta_o
            Tb_glob = np.mean(Tb_nu,axis=0)
            
            print('\nCreating Tb vs nu plot ...')
            left=0.12
            fs=22
            fig,ax=plt.subplots(figsize=(8, 7.9))
            fig.subplots_adjust(left=left, bottom=0.06, right=1-left, top=0.94)
            
            ax.axhline(y=Tcmb_o,color='k',ls='--',lw=1.5, label='CMB')
            ax.plot(nu,Tb_mean,color='r',lw=1.5,ls=':',label=r'$\beta= $%.2f'%beta_o)
            ax.plot(nu,Tb_glob,color='b',lw=1.5,label='Extragalactic')

            if xlog:
                ax.set_xscale('log')
            if ylog:
                ax.set_yscale('log')
            
            ax.set_xlabel(r'$\nu\,$(Hz)',fontsize=fs)
            ax.set_ylabel(r'$T_{\mathrm{b}}^{\mathrm{eg}}\,$(K)',fontsize=fs)

            ax.minorticks_on()
            ax.yaxis.set_ticks_position('both')
            ax.tick_params(axis='both', which='major', labelsize=fs)
            ax.legend(fontsize=18,frameon=False)

            secax = ax.secondary_xaxis('top', functions=(nu2z,z2nu))
            secax.set_xlabel(r'$z$',fontsize=fs, labelpad=12)
            secax.tick_params(which='major', labelsize=fs)

            #plt.xlim([0.1,1e2])
            #plt.ylim([1,3e4])
            ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
            fig_path = self.path+'Tb_vs_nu.pdf'
            plt.savefig(fig_path)
            print('Done. Tb vs frequency saved as',fig_path)
    #End of function plotter()
#End of class run()

