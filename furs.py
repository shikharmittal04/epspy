import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colormaps
from scipy.interpolate import CubicSpline
from scipy.special import legendre
import random
import transformcl as tcl
import healpy as hp
from mpi4py import MPI
from mpi4py.util import pkl5
import sys
import os

#Some fixed numbers ...
kB = 1.38e-23   #Boltzmann constant in J/K
cE = 2.998e8    #Speed of light in m/s
Tcmb_o = 2.725  #CMB temperature today in K
nu21 = 1420     #21-cm frequency in MHz

np.seterr(all='ignore')

def print_banner():
    banner = """\033[94m
    ███████╗ ██╗   ██╗ ██████╗  ███████╗
    ██╔════╝ ██║   ██║ ██╔══██╗ ██╔════╝
    █████╗   ██║   ██║ ██████╔╝ ███████╗
    ██╔══╝   ██║   ██║ ██╔══██╗ ╚════██║
    ██║      ╚██████╔╝ ██║  ██║ ███████║
    ╚═╝       ╚═════╝  ╚═╝  ╚═╝ ╚══════╝
    \033[00m"""                                
    print(banner)
    return None

#The following 2 functions are required for adding a secondary x-axis at the top for figures.
def nu2z(nu):
    return nu21/nu-1

def z2nu(z):
    return nu21/(1+z)


class extragalactic():
    def __init__(self, log2Nside=6, logSmin=-2,logSmax=-1, nu_o=150e6, beta_o=2.681,sigma_beta=0.5, amp=7.8e-3,gam=0.821, path=''):
        self.nu_o = nu_o        #Reference frequency in Hz

        self.beta_o = beta_o    #Mean spectral index for extragalactic point sources
        self.sigma_beta = sigma_beta   #Spread in the beta values

        self.amp = amp  #Amplitude of the power-law 2-point angular correlation function (2PACF)
        self.gam = gam  #-exponent of the power-law 2-point angular correlation function

        self.logSmin = logSmin  #log_10(S_min), where S_min is in Jy
        self.logSmax = logSmax  #log_10(S_max)
        
        self.path = path        #Path where you would like to save and load from, the Tb's and beta's
                            
        self.log2Nside = log2Nside    #Number of divisions in units of log_2
    #End of function __init__()
    
    def print_input(self):
    	print("\n\033[93mnu_o =",self.nu_o)
    	print("beta_o =",self.beta_o)
    	print("sigma_beta =",self.sigma_beta)
    	print("amp =",self.amp)
    	print("gam =",self.gam)
    	print("log2Nside =",self.log2Nside)
    	print("logSmax =",self.logSmax)
    	print("logSmin =",self.logSmin)
    	print("path =",self.path,"\033[00m\n")

    	return None
    	
    def dndS(self, S):
        '''
        This is the distribution of flux density.
        I have taken the functional form and the numbers from Gervasi et al (2008) ApJ.
        Input S is in units of Jy (jansky). Can be 1 value or an array.
        Output is in number of sources per unit solid angle per unit flux density. 1 value or an array depending on input.
        '''
        a1,b1,a2,b2 = -0.854, 0.37, -0.856, 1.47
        A1, B1 = 1.65e-4, 1.14e-4
        A2A1, B2B1 = 0.24, 1.8e7
        A2 = A2A1*A1
        B2 = B2B1*B1
        return S**-2.5*((A1*S**a1+B1*S**b1)**-1+(A2*S**a2+B2*S**b2)**-1)

    def acf(self, chi):
        '''
        This is the popular form of the 2PACF; a power law.
        The default values for amplitude and index are from Rana & Bagla (2019).
        Input chi should be in radians. One number or an array.
        Output is pure number or an array accordingly as chi is a number or an array. 
        '''
        return self.amp*(chi*180/np.pi)**(-self.gam)

    '''
    #Uncomment this function for method 2 for converting 2PACF, i.e., C(\chi), to APS, i.e., C_\ell.
    def acf2Cl(ell):
        def P_ell(ell,chi):
            poly = legendre(ell)
            return poly(np.cos(chi))

        return 2*np.pi*scint.quad(lambda chi: self.acf(chi)*P_ell(ell,chi)*np.sin(chi),0.0,np.pi)[0]
    '''
    
    def num_sources(self):
        '''
        This function gives the total number of unresolved point sources on the full sky.
        This is specifically for the flux density distribution defined in dndS(),
        and the minimum and maximum S values are set during the initialisation of the class object.
        No input is required.
        Output is a pure number.
        '''
        S_space = np.logspace(self.logSmin,self.logSmax,1000)
        dndS_space = self.dndS(S_space)

        Ns_per_sr = np.trapz(dndS_space,S_space)
        Ns = 4*np.pi*Ns_per_sr
        print('\nTotal number of sources = {:d}'.format(round(Ns)))
        return Ns

    def num_den(self):
        '''
        This function calculates the number density function n_clus for the 2PACF defined in acf().
        No input is required.
        Output is in units of number per pixel. It will be an array of length Npix.
        '''
        
        #Method 1:
        '''
        This methods requires transformcl package. Methods 2 does the same thing manually.
        tcl.corrtocl requires us to sample the angle at some specific points; obtained by 'tcl.theta'
        I have chosen 200 randomly. (Using a large like 1000 may result in del_clus<-1.)
        '''
        th = tcl.theta(200)
        cor = self.acf(th)
        Cl_clus = tcl.corrtocl(cor)
        

        '''
        #Method 2:
        #Here we manually compute the C_\ell's without using transformcl. This is much slower but more accurate. 
        Cl_clus = np.zeros(50)
        for i in range(50):
            Cl_clus[i] = self.acf2cl(i)
        '''
        Nside= 2**self.log2Nside
        Npix = hp.nside2npix(Nside) #number of pixels

        #Now calculating the clustered map fluctuation...
        del_clus = hp.synfast(Cl_clus,Nside)
        
        where_n_neg = np.where(del_clus<-1.0)[0]
        Npix_n_neg = 100*np.size(where_n_neg)/Npix
        if Npix_n_neg!=0:
            print('\n\033[31mError! Your choice of 2PACF parameters is NOT valid.')
            print('{:.2f}% pixels have negative number of sources!'.format(Npix_n_neg))
            print('Terminating ...\033[00m\n')
            sys.exit()
        
        print('Done.\nAverage overdensity for the clustered sky (should be ~ 0) = {:.3f}'.format(np.mean(del_clus)))

        Ns = self.num_sources()
        nbar = Ns/Npix
        print('Total number of pixels, Npix =',Npix)
        print('Average number of sources per pixel = {:.2f}'.format(nbar))

        #and the corresponding clustered number density function given the fluctuation...
        n_clus = nbar*(1+del_clus)
        
        where_n_less_than_1 = np.where(np.round(n_clus)<1.0)
        Npix_less_than_1 = 100*np.size(where_n_less_than_1)/Npix
        if Npix_less_than_1>50:
            print('\n\033[91m{:.2f}% pixels have no sources!'.format(Npix_less_than_1))
            print('This can happen either because there are very few sources in your chosen flux density range or your resolution is too high.')
            print('Recommendation: either increase (`logSmax`-`logSmin`) or decrease `log2Nside`.\n\033[00m')
        
        return n_clus
    #End of function num_den()
    
    def ref_freq(self):
        '''
        Tb_o_individual is an array of arrays of unequal lengths, i.e.,
        all of Tb_o_individual[0], Tb_o_individual[1], ..., Tb_o_individual[Npix] are arrays of different lengths.
        The length of Tb_o[j] tells us the number of sources, say N_j, on the jth pixel and
        Tb_o_individual[j][0], Tb_o_individual[j][1], ..., Tb_o_individual[j][N_j] are the temperatures (at ref. frequency) due to 0th, 1st,...(N_j)th source on the jth pixel.
        '''
        #-------------------------------------------------------------------------------------
        comm = pkl5.Intracomm(MPI.COMM_WORLD)
        cpu_ind = comm.Get_rank()
        Ncpu = comm.Get_size()
        
        if cpu_ind==0: print_banner()
        if Ncpu==1: print("\033[91mBetter to parallelise. Eg. 'mpirun -np 4 python3 %s', where 4 specifies the number of tasks.\033[00m" %(sys.argv[0]))
            
        #-------------------------------------------------------------------------------------
        Nside= 2**self.log2Nside
        Npix = hp.nside2npix(Nside) #number of pixels

        #Find the number density distribution on the master CPU and share it with all CPUs.
        if cpu_ind==0:
            '''
            Find the number density distribution on the master CPU.
            '''
            if os.path.isdir(self.path)==False:
                print('The requested directory does not exist. Creating one ...')
                os.mkdir(self.path)

            print('\n\033[94mRunning ref_freq() ...\033[00m\n')          
            print('Finding the clustered number density distribution ...')
            n_clus = self.num_den()
            
            
            
            n_clus_save_name = self.path+'n_clus.npy'
            np.save(n_clus_save_name,n_clus)
            print('\033[32mThe clustered number density has been saved into file:\n',n_clus_save_name,'\033[00m')

            print("\nAssigning flux density and power-law index ...")
        else:
            n_clus = None

        n_clus = comm.bcast(n_clus, root=0) #Now all CPUs have the same number density distribution.
        #-------------------------------------------------------------------------------------
        #For each pixel on the sky and for each source on that pixel, assign flux and spectral index.
        
        S_space = np.logspace(self.logSmin,self.logSmax,1000)
        dndS_space = self.dndS(S_space)
        Ns_per_sr = np.trapz(dndS_space,S_space)

        ppc = int(Npix/Ncpu)    #pixels per cpu
        Omega_pix = hp.nside2pixarea(Nside) #Solid angle per pixel
        
        Tb_o_local = np.zeros(ppc)
        Tb_o_local_individual = np.zeros(ppc, dtype=object)
        beta_local = np.zeros(ppc, dtype=object)

        for j in np.arange(cpu_ind*ppc,(cpu_ind+1)*ppc):
            N = round(n_clus[j])    #no.of sources on jth pixel
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
                N = round(n_clus[j])  #no.of sources on jth pixel
                So_j = np.array(random.choices(S_space,weights=dndS_space/Ns_per_sr,k=N))   #select N flux densities for jth pixel
                beta_remain[j-Ncpu*ppc] = np.random.normal(loc=self.beta_o,scale=self.sigma_beta,size=N)   #select N spectral indices for jth pixel
                Tb_o_remain_individual[j-Ncpu*ppc] = 1e-26*So_j*cE**2/(2*kB*self.nu_o**2*Omega_pix)
                Tb_o_remain[j-Ncpu*ppc] = np.sum(Tb_o_remain_individual[j-Ncpu*ppc])

        #-------------------------------------------------------------------------------------
        comm.Barrier()
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
            I will save the Tb's and beta's as hickle objects (in format '.hkl').
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
            Tb_o_save_name = self.path+'Tb_o_map.npy'
            beta_save_name = self.path+'beta.npy'
            
           
            np.save(Tb_o_individual_save_name,Tb_o_individual)
            #hkl.dump(Tb_o_individual, Tb_o_individual_save_name, mode='w')
            np.save(Tb_o_save_name,Tb_o)
            np.save(beta_save_name,beta)
            #hkl.dump(beta, beta_save_name, mode='w')
                
            print('\033[32mThe brightness temperature (at reference frequency) for each source saved into:\n',Tb_o_individual_save_name,'\033[00m\n')
            print('\033[32mThe pixel wise brightness temperature (at reference frequency) saved into:\n',Tb_o_save_name,'\033[00m\n')
            print('\033[32mThe spectral index for each source saved into:\n',beta_save_name,'\033[00m')
            print('\n\033[94m================ End of function ref_freq(). ================\033[00m\n')
            
            mempertask = 2e-6*os.path.getsize(Tb_o_individual_save_name)
            if mempertask > 2000:
                print("\033[96mRecommendation for '--mem-per-task' to run gen_freq() {:d} MB\n\033[00m".format(round(mempertask)))
        comm.Barrier()
        return None
    #End of function ref_freq()


    def gen_freq(self, nu=1e6*np.arange(50,200)):
        '''
        If you are running this function you must have run ref_freq().    
        This function computes the map(s) at general frequency(ies) based on the precomputed values from ref_freq().
        nu is the frequency (in Hz) at which you want to evaluate the brightness temperature map.
        nu can be one number or an array.
        '''
        #-------------------------------------------------------------------------------------
        comm = MPI.COMM_WORLD
        cpu_ind = comm.Get_rank()
        Ncpu = comm.Get_size()
        #-------------------------------------------------------------------------------------
        Nside= 2**self.log2Nside
        Npix = hp.nside2npix(Nside) #number of pixels

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

        if cpu_ind==0:
            print_banner()
            print('\n\033[94mRunning gen_freq() ...\033[00m\n')
            print("Beginning scaling extragalactic maps to general frequency ...")
        N_nu = np.size(nu)
        Tb_nu_final = np.zeros((Npix,N_nu),dtype='float64')    

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
        comm.Barrier()
        #Now all CPUs have done their job of calculating the Tb's and beta's. 
        if cpu_ind!=0:
            '''
            I am a worker CPU. Sending my Tb to master CPU.
            '''
            comm.Send([Tb_nu_final,MPI.FLOAT], dest=0, tag=11)
        else:
            '''
            I am the master CPU. Receiving all Tb's.
            '''
            for i in range(1,Ncpu):
                receive_local = np.empty((Npix, N_nu),dtype='float64')
                comm.Recv([receive_local,MPI.FLOAT],source=i, tag=11)
                Tb_nu_final = Tb_nu_final + receive_local

            Tb_nu_save_name = self.path+'Tb_nu_map.npy'
            np.save(Tb_nu_save_name,Tb_nu_final)
            
            Tb_nu_glob = np.mean(Tb_nu_final,axis=0)
            tbnuglob = self.path+'Tb_nu_glob.npy'
            np.save(tbnuglob,Tb_nu_glob)

            print('Done.\n\033[32mFile saved as',Tb_nu_save_name,'\033[00m')
            print('It is an array of shape',np.shape(Tb_nu_final))
            print('\n\033[94m================ End of function gen_freq(). ================\033[00m\n')

            nu_save_name = self.path+'nu_glob.npy'
            np.save(nu_save_name,nu)
        comm.Barrier()
        return None    
    #End of function gen_freq()

    def visual(self, nu_skymap=None, t_skymap=False, spectrum=True, n_skymap=False, xlog=False,ylog=True, fig_ext = 'pdf'):
        '''
        Use this function for creating a sky map at a given freqeuncy ('skymap') and/or
        the global extragalactic foregrounds as a function of frequency ('spectrum').
        'nu_skymap' is required only for making the sky map. It should be one number in Hz.
        By default we only plot the spectrum and not the skymap.
        'xlog' and 'ylog' are the boolean values deciding the scale of x and y axis, respectively.
        '''
        #-------------------------------------------------------------------------------------
        comm = MPI.COMM_WORLD
        cpu_ind = comm.Get_rank()
        Ncpu = comm.Get_size()
        #-------------------------------------------------------------------------------------
        
        if cpu_ind==0:
            print('\n\033[94mRunning visual() ...\033[00m\n')
            if Ncpu>1:
                print("\033[91m'visual' does not require parallelisation.\033[00m")
                print("\033[91mYou can run as 'python3 %s'.\033[00m\n" %(sys.argv[0]))
            nu = np.load(self.path+'nu_glob.npy')
            Tb_o_map = np.load(self.path+'Tb_o_map.npy')
            Tb_o_glob = np.mean(Tb_o_map)
            
            Tb_nu_map = np.load(self.path+'Tb_nu_map.npy')
            Tb_nu_glob = np.load(self.path+'Tb_nu_glob.npy')
            
            if np.size(nu_skymap)==1 and nu_skymap is not None: t_skymap = True
            
            plt.rc('text', usetex=True)
            plt.rc('font', family='serif')
            
            if t_skymap:    
                if np.size(nu_skymap)==1:
                    if nu_skymap==None:
                        print("No frequency given with 't_skymap=True'. Creating sky map at the reference frequency ...")
                        nu_skymap=self.nu_o
                        Tb_plot = Tb_o_map
                    else:
                        ind = np.where(nu==nu_skymap)[0]
                        if np.size(ind)==0:
                            if nu_skymap<np.min(nu):
                                print('\033[91mGiven frequency outside the range.\033[00m')
                                print('Using the lowest available frequency; {:.2f} MHz ...'.format(np.min(nu)/1e6))
                                Tb_plot = Tb_nu_map[:,0]
                                nu_skymap = np.min(nu)
                            elif nu_skymap>np.max(nu):
                                print('\033[91mGiven frequency outside the range.\033[00m')
                                print('Using the highest available frequency; {:.2f} MHz ...'.format(np.max(nu)/1e6))
                                Tb_plot = Tb_nu_map[:,-1]
                                nu_skymap = np.max(nu)
                            else:
                                print('Given frequency unavailable in gen_freq(). Interpolating ...')
                                print("Creating sky map at {:.2f} MHz...".format(nu_skymap/1e6))
                                spl = CubicSpline(nu, Tb_nu_map, axis=1)
                                Tb_plot = spl(nu_skymap)
                        else:
                            print("Creating sky map at {:.2f} MHz...".format(nu_skymap/1e6))
                            Tb_plot = Tb_nu_map[:,ind[0]]

                else:
                    print("\033[91mMultiple values given for 'nu_skymap' with 't_skymap=True'. Plotting only at the reference frequency ...\033[00m")
                    Tb_plot = Tb_o_map

                hp.mollview(Tb_plot,title=None,unit=r'$T_{\mathrm{b}}^{\mathrm{eg}}\,$(K)',cmap=colormaps['coolwarm'],norm='log') #,min=0.05,max=200
                hp.graticule()
                
                fig_path = self.path+'Tb_nu_map_'+str(int(nu_skymap/1e6))+'-MHz.'+fig_ext
                plt.savefig(fig_path, bbox_inches='tight')

                print('Done.\n\033[32mTb map saved as:\n',fig_path,'\033[00m')
                
            if n_skymap:
                print('Creating number density map ...')
                n_clus = np.load(self.path+'n_clus.npy')
                nmax = int(np.max(n_clus))+1
                nmin = int(np.min(n_clus))
                
                hp.mollview(n_clus,title=None,unit='$n$',min=nmin,max=nmax)
                hp.graticule()
                
                fig_path = self.path+'n_clus.' + fig_ext
                plt.savefig(fig_path, bbox_inches='tight')
                print('Done.\n\033[32mnumber density map saved as:\n',fig_path,'\033[00m')
            if spectrum:
                Tb_mean = Tb_o_glob*(nu/self.nu_o)**-self.beta_o
                
                print('\nCreating Tb vs nu plot ...')
                left=0.12
                fs=22
                fig,ax=plt.subplots(figsize=(8, 7.9))
                fig.subplots_adjust(left=left, bottom=0.06, right=1-left, top=0.94)
                
                ax.axhline(y=Tcmb_o,color='k',ls='--',lw=1.5, label='CMB')
                ax.plot(nu/1e6,Tb_mean,color='r',lw=1.5,ls=':',label=r'$\beta= $ %.2f'%self.beta_o)
                ax.plot(nu/1e6,Tb_nu_glob,color='b',lw=1.5,label='Extragalactic')

                if xlog:
                    ax.set_xscale('log')
                if ylog:
                    ax.set_yscale('log')
                
                ax.set_xlabel(r'$\nu\,$(MHz)',fontsize=fs)
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
                fig_path = self.path+'Tb_vs_nu.' + fig_ext
                plt.savefig(fig_path)
                
                print('Done.\n\033[32mTb vs frequency saved as:\n',fig_path,'\n\033[00m')
            
            print('\n\033[94m================ End of function visual(). ================\033[00m\n')
    #End of function visual()
#End of class extragalactic()

