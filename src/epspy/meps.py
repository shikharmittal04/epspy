#Copyright 2024, Shikhar Mittal.
'''
Application programming interface

This is the module for extragalactic point sources (hence the name 'meps'). It contains public functions :func:`save_eps`, :func:`load_eps` and a class :class:`eps`. 
'''
import numpy as np
from numpy.polynomial import polynomial as pol
import pickle
import matplotlib.pyplot as plt
from matplotlib import colormaps
from scipy.interpolate import CubicSpline
from scipy.special import legendre
import scipy.optimize as op
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

def _print_banner():
    banner = """\033[94m
    ███████╗██████╗ ███████╗
    ██╔════╝██╔══██╗██╔════╝
    █████╗  ██████╔╝███████╗
    ██╔══╝  ██╔═══╝ ╚════██║
    ███████╗██║     ███████║
    ╚══════╝╚═╝     ╚══════╝
    Copyright 2024, Shikhar Mittal.
    \033[00m"""                                
    print(banner)
    return None

#--------------------------------------------------------------------------------------------
#The following 2 functions are required for adding a secondary x-axis at the top for figures.
def _nu2z(nu):
    return nu21/nu-1

def _z2nu(z):
    return nu21/(1+z)

#--------------------------------------------------------------------------------------------
#The following 2 functions will be useful if you want to save and load your class object.
def save_eps(obj, filename):
    '''Saves the class object :class:`eps`.
    
    Save the class object :class:`eps` for later use. It will save the object in the path where you have all the other outputs
    from this package.
    
    Parameters
    ~~~~~~~~~~

    obj : class
        This should be the class object you want to save.
        
    filename : str
        Give a filename to your object. It will be saved in the ``obj.path`` directory.
    
    '''
    try:
        comm = MPI.COMM_WORLD
        cpu_ind = comm.Get_rank()
        Ncpu = comm.Get_size()
    except:
        cpu_ind=0
    if cpu_ind==0:
        if filename[-4:]!='.pkl': filename=filename+'.pkl'
        fullpath = obj.path+filename
        with open(fullpath, 'wb') as outp:  # Overwrites any existing file.
            pickle.dump(obj, outp, pickle.HIGHEST_PROTOCOL)
        print('\033[32m',filename,'saved as',fullpath,'\033[00m\n')
    return None
    
def load_eps(filename):
    '''To load the class object :class:`eps`.
    
    Parameters
    ~~~~~~~~~~

    filename : str
        This should be the name of the file you gave in :func:`save_eps()` for saving class object :class:`eps`.
        Important: provide the full path for ``filename`` with the extension ``.pkl``.
        
    Returns
    ~~~~~~~

    class object    
    '''
    try:
        comm = MPI.COMM_WORLD
        cpu_ind = comm.Get_rank()
        Ncpu = comm.Get_size()
    except:
        cpu_ind=0
    if cpu_ind==0:
        with open(filename, 'rb') as inp:
            epsobj = pickle.load(inp)
        print('Loaded the eps class object.\n')
    return epsobj
#--------------------------------------------------------------------------------------------


class eps():
    '''
    This is the class for extragalactic point sources.
    
    Parameters
    ~~~~~~~~~~

    nu_o : float, optional
        Reference frequency, :math:`\\nu_0`, in Hz
    
    beta_o : float, optional
        Mean spectral index for extragalactic point sources; :math:`\\beta_0`
    
    sigma_beta : float, optional
        Standard deviation, :math:`\\sigma_{\\beta}`, of the Gaussian distribution of :math:`\\beta`

    amp : float, optional
        Amplitude, :math:`A`, of the power-law 2-point angular correlation function (2PACF)
    
    gam : float, optional
        Negative exponent, :math:`\\gamma`, of the power-law 2PACF
    
    logSmin : float, optional
        :math:`\\log_{10}(S_{\\mathrm{min}})`, where :math:`S_{\\mathrm{min}}` is in Jy
    
    logSmax : float, optional
        :math:`\\log_{10}(S_{\\mathrm{max}})`, where :math:`S_{\\mathrm{max}}` is in Jy
    
    dndS_form : int, optional
        Choose the functional form for :math:`\\mathrm{d}n/\\mathrm{d}S`. Available options -> 0 (default),1 or 2. For more details see the API for :func:`dndS`.
    
    log2Nside : int, optional
        Number of divisions of the base HEALPix pixel in units of :math:`\\log_2`. It sets the sky pixelisation. The number of pixels is :math:`N_{\\mathrm{pix}} = 12\\times 2^{2k}`, where :math:`k=` ``log2Nside``
    
    path : str, optional
        Path where you would like to put all your outputs
            
    lbl : str, optional
        Append an extra string to all the output file names
    
    Derived attributes
    ~~~~~~~~~~~~~~~~~~
    
    Nside : int
        The :code:`nside` required in healpy package
        
    Npix : int
        Number of pixels; :math:`N_{\\mathrm{pix}}`
    
    Nps : float
        Total number of extragalactic point sources on the sky in the :math:`S` range :math:`S_{\\mathrm{min}}` and :math:`S_{\\mathrm{max}}`; :math:`N_{\\mathrm{ps}}`
        
    Methods
    ~~~~~~~
    
    '''
    def __init__(self, beta_o=2.681,sigma_beta=0.5, logSmin=-2,logSmax=-1,dndS_form=0, log2Nside=6, nu_o=150e6, amp=7.8e-3,gam=0.821, path='',lbl=''):
        self.nu_o = nu_o

        self.beta_o = beta_o    
        self.sigma_beta = sigma_beta

        self.amp = amp  
        self.gam = gam  

        self.logSmin = logSmin
        self.logSmax = logSmax
        self.dndS_form = dndS_form
        
        self.log2Nside = log2Nside
        
        self.path = path
        self.lbl = lbl
        
        self.Nside= 2**self.log2Nside
        self.Npix = hp.nside2npix(self.Nside)
        
        
        S_space = np.logspace(self.logSmin,self.logSmax,1000)
        dndS_space = self.dndS(S_space)

        Nps_per_sr = np.trapz(dndS_space,S_space)
        self.Nps = 4*np.pi*Nps_per_sr
        mempertask = 1.1*16e-6*self.Nps
        
        try:
            comm = MPI.COMM_WORLD
            cpu_ind = comm.Get_rank()
            Ncpu = comm.Get_size()
        except:
            cpu_ind=0
        if cpu_ind==0:
            if mempertask > 2000:
                print("\033[96mRecommendation for '--mem-per-task' {:d} MB\n\033[00m".format(round(mempertask)))
        
    #End of function __init__()
    
    def print_input(self):
        '''Prints the input parameters you gave.'''
        
        print("\n\033[93mnu_o =",self.nu_o)
        print("beta_o =",self.beta_o)
        print("sigma_beta =",self.sigma_beta)
        print("amp =",self.amp)
        print("gam =",self.gam)
        print("log2Nside =",self.log2Nside)
        print("logSmax =",self.logSmax)
        print("logSmin =",self.logSmin)
        print("dndS_form =",self.dndS_form)
        print("path =",self.path)
        print("lbl =",self.lbl,"\033[00m\n")
        
        return None
        
    def dndS(self, S):
        ''':math:`\\mathrm{d}n/\\mathrm{d}S=\\mathrm{d}n/\\mathrm{d}S(S)`
        
        Distribution of flux density, S.
        
        The default choice (0) is by `Gervasi et al (2008) <https://iopscience.iop.org/article/10.1086/588628/meta>`__. It is approximately a power law. See paper for the exact form.
        
        Form 1 is by `Mandal et al. (2021) <https://www.aanda.org/articles/aa/full_html/2021/04/aa39998-20/aa39998-20.html>`__. It is a :math:`7^{\\mathrm{th}}` order log-log polynomial fit to the Euclidean number count.
        
        Form 2 is by `Intema et al. (2017) <https://www.aanda.org/articles/aa/full_html/2017/02/aa28536-16/aa28536-16.html>`__. It is a :math:`5^{\\mathrm{th}}` order log-log polynomial fit to the Euclidean number count.
        
        Parameters
        ^^^^^^^^^^
        
        S : float
            Flux density in units of Jy (jansky). Can be 1 value or a :mod:`numpy` array.
        
        Returns
        ^^^^^^^
        
        float
            Number of point sources per unit solid angle per unit flux density :math:`(\\mathrm{sr}^{-1}\\mathrm{Jy}^{-1})`. 1 value or an array depending on input.
        '''

        if self.dndS_form==1:
            P = pol.Polynomial((1.655, -0.115,0.2272,0.51788,-0.449661,0.160265,-0.028541,0.002041))
            return S**-2.5*10**(P(np.log10(S*1e3)))        

        elif self.dndS_form==2:
            P = pol.Polynomial((3.5142, 0.3738, -0.3138, -0.0717, 0.0213, 0.0097))
            return S**-2.5*10**P(np.log10(S))
        else:
            if self.dndS_form!=0:
                print("\033[31mInvalid option! Using sum-of-2-double-inverse-power-law form ... \033[00m")
            
            a1,b1,a2,b2 = -0.854, 0.37, -0.856, 1.47
            A1, B1 = 1.65e-4, 1.14e-4
            A2A1, B2B1 = 0.24, 1.8e7
            A2 = A2A1*A1
            B2 = B2B1*B1
            return S**-2.5*((A1*S**a1+B1*S**b1)**-1+(A2*S**a2+B2*S**b2)**-1)

    def acf(self, chi):
        ''':math:`C(\\chi)=A\\chi^{-\\gamma}`
        
        2-point angular correlation function (2PACF). The amplitude :math:`A` and negative power-law index :math:`\\gamma` are set when you initialise the class object. If you want to simulate an isotropic sky, set ``amp=0``. The default values for amplitude and index are from `Rana & Bagla (2019) <https://academic.oup.com/mnras/article/485/4/5891/5420431>`__.
                
        Parameters
        ^^^^^^^^^^

        chi : float
            Angle at which you want to get the 2PACF, should be in radians. One number or an array.
        
        Returns
        ^^^^^^^

        float
            A pure number or an array accordingly as ``chi`` is a number or an array. 
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
    
    def num_den(self):
        ''':math:`n_{\\mathrm{ps}}=n_{\\mathrm{ps}}(\\hat{n})`
        
        The number density function. It is dependent on the choice of 2PACF and :math:`\\mathrm{d}n/\\mathrm{d}S`.
        
        Returns
        ^^^^^^^
        
        float
            An array of length :math:`N_{\\mathrm{pix}}` whose elements are the number of point sources on the corresponding pixel. The array will also be saved as ``n_ps.npy`` in the path you gave during initialisation.
        '''
        
        nbar = self.Nps/self.Npix
        print('\nTotal number of point sources, Nps = {:d}'.format(round(self.Nps)))
        print('Total number of pixels, Npix =',self.Npix)
        print('Average number of sources per pixel, n_bar = {:.2f}'.format(nbar))
        
        if self.amp==0:
            print('Simulating an isotropic sky ...')
            n_ps = np.random.poisson(lam=nbar,size=self.Npix)
        else:
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
            
            print('Simulating a clustered sky ...')
            del_clus = hp.synfast(Cl_clus,self.Nside)
            
            where_n_neg = np.where(del_clus<-1.0)[0]
            Npix_n_neg = 100*np.size(where_n_neg)/self.Npix
            if Npix_n_neg!=0:
                print('\n\033[31mError! Your choice of 2PACF parameters is NOT valid.')
                print('{:.2f}% pixels have negative number of sources!'.format(Npix_n_neg))
                print('Terminating ...\033[00m\n')
                sys.exit()
            
            print('Done.\nAverage overdensity for the clustered sky (should be ~ 0) = {:.3f}'.format(np.mean(del_clus)))

            #and the corresponding clustered number density function given the fluctuation...
            n_ps = nbar*(1+del_clus)
        
        where_n_less_than_1 = np.where(np.round(n_ps)<1.0)
        Npix_less_than_1 = 100*np.size(where_n_less_than_1)/self.Npix
        if Npix_less_than_1>50:
            print('\n\033[91m{:.2f}% pixels have no sources!'.format(Npix_less_than_1))
            print('This can happen either because there are very few sources in your chosen flux density range or your resolution is too high.')
            print('Recommendation: either increase (`logSmax`-`logSmin`) or decrease `log2Nside`.\n\033[00m')
        
        n_ps_save_name = self.path+'n_ps'+self.lbl
        np.save(n_ps_save_name,n_ps)
        print('\033[32mThe number density file has been saved as:\n',n_ps_save_name,'\033[00m')
            
        return n_ps
    #End of function num_den()
    
    def ref_freq(self):
        '''Generates the brightness temperature and spectral indices at reference frequency. 
        
        3 files will be generated ``Tb_o_individual.npy``, ``Tb_o_map.npy`` and ``beta.npy``
        
        To understand the structure of these output files see :ref:`ref-freq`.
        
        '''
        #-------------------------------------------------------------------------------------
        comm = pkl5.Intracomm(MPI.COMM_WORLD)
        cpu_ind = comm.Get_rank()
        Ncpu = comm.Get_size()
        
        if cpu_ind==0:
            _print_banner()
            self.print_input()
        
        if Ncpu==1: print("\033[91mBetter to parallelise. Eg. 'mpirun -np 4 python3 %s', where 4 specifies the number of tasks.\033[00m" %(sys.argv[0]))
            
        #-------------------------------------------------------------------------------------
        #Find the number density distribution on the master CPU and share it with all CPUs.
        if cpu_ind==0:
            #Find the number density distribution on the master CPU.
            if os.path.isdir(self.path)==False:
                print('The requested directory does not exist. Creating one ...')
                os.mkdir(self.path)

            print('\n\033[94mRunning eps.ref_freq() ...\033[00m\n')          
            print('Finding the number density distribution ...')
            n_ps = self.num_den()

            print("\nAssigning flux density and power-law index ...")
        else:
            n_ps = None

        n_ps = comm.bcast(n_ps, root=0) #Now all CPUs have the same number density distribution.
        #-------------------------------------------------------------------------------------
        #For each pixel on the sky and for each source on that pixel, assign flux and spectral index.
        
        S_space = np.logspace(self.logSmin,self.logSmax,1000)
        dndS_space = self.dndS(S_space)
        Nps_per_sr = np.trapz(dndS_space,S_space)

        ppc = int(self.Npix/Ncpu)    #pixels per cpu
        Omega_pix = hp.nside2pixarea(self.Nside) #Solid angle per pixel
        
        Tb_o_local = np.zeros(ppc)
        Tb_o_local_individual = np.zeros(ppc, dtype=object)
        beta_local = np.zeros(ppc, dtype=object)

        for j in np.arange(cpu_ind*ppc,(cpu_ind+1)*ppc):
            N = round(n_ps[j])    #no.of sources on jth pixel
            So_j = np.array(random.choices(S_space,weights=dndS_space/Nps_per_sr,k=N))   #select N flux densities for jth pixel, in Jy
            beta_local[j-cpu_ind*ppc] = np.random.normal(loc=self.beta_o,scale=self.sigma_beta,size=N) #select N spectral indices for jth pixel
                        
            Tb_o_local_individual[j-cpu_ind*ppc] = 1e-26*So_j*cE**2/(2*kB*self.nu_o**2*Omega_pix)
            Tb_o_local[j-cpu_ind*ppc] = np.sum(Tb_o_local_individual[j-cpu_ind*ppc])

        #An additional short loop is required if Npix/Ncpu is not an integer. Do the remaining pixels ('rm_pix') on rank 0.
        rm_pix=self.Npix%Ncpu
        if cpu_ind==0 and rm_pix!=0:
            Tb_o_remain_individual = np.zeros(rm_pix, dtype=object)
            Tb_o_remain = np.zeros(rm_pix)
            beta_remain = np.zeros(rm_pix, dtype=object)
            for j in np.arange(Ncpu*ppc,self.Npix):
                N = round(n_ps[j])  #no.of sources on jth pixel
                So_j = np.array(random.choices(S_space,weights=dndS_space/Nps_per_sr,k=N))   #select N flux densities for jth pixel
                beta_remain[j-Ncpu*ppc] = np.random.normal(loc=self.beta_o,scale=self.sigma_beta,size=N)   #select N spectral indices for jth pixel
                Tb_o_remain_individual[j-Ncpu*ppc] = 1e-26*So_j*cE**2/(2*kB*self.nu_o**2*Omega_pix)
                Tb_o_remain[j-Ncpu*ppc] = np.sum(Tb_o_remain_individual[j-Ncpu*ppc])

        #-------------------------------------------------------------------------------------
        comm.Barrier()
        #Now all CPUs have done their jobs of calculating the Tb's and beta's. 
        if cpu_ind!=0:
            #I am a worker CPU. Sending my Tb's and beta's to master CPU.
            comm.send(Tb_o_local, dest=0, tag=11)
            comm.send(Tb_o_local_individual, dest=0, tag=13)
            comm.send(beta_local, dest=0, tag=29)
        else:
            #I am the master CPU. Receiving all Tb's and beta's.
            #I will save the Tb's and beta's as .npy objects.
            
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

            Tb_o_individual_save_name = self.path+'Tb_o_individual'+self.lbl
            Tb_o_save_name = self.path+'Tb_o_map'+self.lbl
            beta_save_name = self.path+'beta'+self.lbl
            
           
            np.save(Tb_o_individual_save_name,Tb_o_individual)
            np.save(Tb_o_save_name,Tb_o)
            np.save(beta_save_name,beta)
                
            print('\033[32mThe brightness temperature (at reference frequency) for each source saved into:\n',Tb_o_individual_save_name,'\033[00m\n')
            print('\033[32mThe pixel wise brightness temperature (at reference frequency) saved into:\n',Tb_o_save_name,'\033[00m\n')
            print('\033[32mThe spectral index for each source saved into:\n',beta_save_name,'\033[00m')
            print('\n\033[94m================ End of function eps.ref_freq() ================\033[00m\n')
    
        comm.Barrier()
        return None
    #End of function ref_freq()


    def gen_freq(self, nu=1e6*np.arange(50,201)):
        '''Scale the brightness temperature at reference frequency to a general frequency.
        
        If you are running this function you must have run :func:`ref_freq`.   
         
        This function computes the map(s) at general frequency(ies) based on the precomputed values from :func:`ref_freq`.
        
        Parameters
        ^^^^^^^^^^

        nu : float
            Frequency (in Hz) at which you want to evaluate the foregrounds due to extragalactic point sources. Can be one number or an array.
            (Default = `1e6*np.arange(50,201)`)
            
        3 files will be generated namely, ``Tb_nu_glob.npy``, ``Tb_nu_glob.npy``, and ``nu_glob.npy``.
        
        To understand the structure of these output files see :ref:`gen-freq`.
        '''
        #-------------------------------------------------------------------------------------
        comm = MPI.COMM_WORLD
        cpu_ind = comm.Get_rank()
        Ncpu = comm.Get_size()
        #-------------------------------------------------------------------------------------

        def Tb_nu(Tb_ref,beta,nu):
            '''
            Given the brightness temperature Tb_ref (at reference frequency nu_o), spectral indices beta and frequency nu
            return the brightness temperature in K at requested frequency nu.
            Return value is one number.
            nu should be in Hz. (Best to give one nu at a time.)
            Tb_o and beta should be of same dimensions if both are arrays of size more than 1.
            '''    
            return np.sum(Tb_ref*(nu/self.nu_o)**-(beta))


        Tb_o_individual_save_name = self.path+'Tb_o_individual'+self.lbl+'.npy'
        beta_save_name = self.path+'beta'+self.lbl+'.npy'

        if cpu_ind==0:
            _print_banner()
            print('\n\033[94mRunning eps.gen_freq() ...\033[00m\n')
            print("Beginning scaling extragalactic maps to general frequency ...")
        N_nu = np.size(nu)
        Tb_nu_final = np.zeros((self.Npix,N_nu),dtype='float64')    

        ppc = int(self.Npix/Ncpu)    #pixels per cpu
        slctd_Tb_o = np.load(Tb_o_individual_save_name,allow_pickle=True)[cpu_ind*ppc:(cpu_ind+1)*ppc]
        slctd_beta = np.load(beta_save_name,allow_pickle=True)[cpu_ind*ppc:(cpu_ind+1)*ppc]

        for j in np.arange(cpu_ind*ppc,(cpu_ind+1)*ppc):
            for i in range(N_nu):
                Tb_nu_final[j,i] = Tb_nu(slctd_Tb_o[j-cpu_ind*ppc],slctd_beta[j-cpu_ind*ppc],nu[i])

        del slctd_Tb_o
        del slctd_beta

        #An additional short loop is required if Npix/Ncpu is not an integer. We do the remaining pixels on rank 0.
        if cpu_ind==0 and self.Npix%Ncpu!=0:
            slctd_Tb_o = np.load(Tb_o_individual_save_name,allow_pickle=True)[Ncpu*ppc:self.Npix]
            slctd_beta = np.load(beta_save_name,allow_pickle=True)[Ncpu*ppc:self.Npix]

            for j in np.arange(Ncpu*ppc,self.Npix):
                for i in range(N_nu):
                    Tb_nu_final[j,i] = Tb_nu(slctd_Tb_o[j-Ncpu*ppc],slctd_beta[j-Ncpu*ppc],nu[i])

        #-------------------------------------------------------------------------------------
        comm.Barrier()
        #Now all CPUs have done their job of calculating the Tb's and beta's. 
        if cpu_ind!=0:
            #I am a worker CPU. Sending my Tb to master CPU.
            comm.Send([Tb_nu_final,MPI.FLOAT], dest=0, tag=11)
        else:
            #I am the master CPU. Receiving all Tb's.
            for i in range(1,Ncpu):
                receive_local = np.empty((self.Npix, N_nu),dtype='float64')
                comm.Recv([receive_local,MPI.FLOAT],source=i, tag=11)
                Tb_nu_final = Tb_nu_final + receive_local

            Tb_nu_save_name = self.path+'Tb_nu_map'+self.lbl
            np.save(Tb_nu_save_name,Tb_nu_final)
            
            Tb_nu_glob = np.mean(Tb_nu_final,axis=0)
            tbnuglob = self.path+'Tb_nu_glob'+self.lbl
            np.save(tbnuglob,Tb_nu_glob)

            print('Done.\n\033[32mFile saved as',Tb_nu_save_name,'\033[00m')
            print('It is an array of shape',np.shape(Tb_nu_final))
            print('\n\033[94m================ End of function eps.gen_freq() ================\033[00m\n')

            nu_save_name = self.path+'nu_glob.npy'
            np.save(nu_save_name,nu)
        comm.Barrier()
        return None    
    #End of function gen_freq()
    
    def couple2D(self, bd=None):
        '''Couple the foregrounds generated by :func:`gen_freq` to the antenna beam directivity.
        
        ``Tb_nu_map.npy`` generated by :func:`gen_freq` does not account for chromatic distortions in the antenna. To use this function one must provide an array, which should be in the shape of :math:`N_{\\mathrm{pix}} \\times N_{\\nu}`. An output file will be generated called ``T_ant.npy``.
        
        Also, the best-fitting parameters (along with :math:`1\\sigma` uncertainty) :math:`T_{\\mathrm{f}}, \\beta_{\\mathrm{f}}` and :math:`\\Delta\\beta_{\\mathrm{f}}` based on a simple least-squares fitting of power-law-with-a-running-spectral-index function to the antenna temperature data will be printed.

        Parameters
        ^^^^^^^^^^

        bd : str 
            Full path to the beam directivity array. By default we call this array as 'D' and assume it to be in the path where you have all other outputs.
        '''
        #-------------------------------------------------------------------------------------
        try:
            comm = MPI.COMM_WORLD
            cpu_ind = comm.Get_rank()
            Ncpu = comm.Get_size()
        except:
            cpu_ind=0
            Ncpu=1
        #-------------------------------------------------------------------------------------
        def fitter(nu,T_f,beta_f,Dbeta_f):
            return T_f*(nu/self.nu_o)**(-beta_f+Dbeta_f*np.log(nu/self.nu_o))


        
        if cpu_ind==0:
            print('\n\033[94mRunning eps.couple2D() ...\033[00m\n')
            
            if bd==None:
                D_file_path = self.path+'D.npy'
            else:
                D_file_path = bd

            D = np.load(D_file_path)
            
            nu_save_name = self.path+'nu_glob.npy'
            nu = np.load(nu_save_name)
            
            Tb_nu_save_name = self.path+'Tb_nu_map'+self.lbl+'.npy'
            Tb_nu_map = np.load(Tb_nu_save_name)
            
            if np.shape(Tb_nu_map)!=np.shape(D):
                print('\033[31mError:')
                print('Directivity pattern array should be of shape (Npix, Nnu),')
                print('where Npix =',self.Npix, 'and Nnu =',len(nu),'.\n')
                print('Current given D has shape',np.shape(D))
                print('Terminating ...\033[00m\n')
                sys.exit()
            
            print('Performing integral of foregrounds weighted by directivity over the sky ...')
            T_ant = np.mean(Tb_nu_map*D,axis=0)
            
            Tb_o_map = np.load(self.path+'Tb_o_map'+self.lbl+'.npy')
            Tb_o_glob = np.mean(Tb_o_map)
            popt,pcov = op.curve_fit(fitter,nu,T_ant,p0=[Tb_o_glob,self.beta_o,0.0])
            
            print('\nT_f, beta_f, Dbeta_f = {:.3f} {:.3f} {:.3f}'.format(popt[0],popt[1],popt[2]))
            print('1sigma uncertainty = ',np.sqrt(np.diag(pcov)),'\n')
            
            T_ant_save_name = self.path+'T_ant'+self.lbl
            np.save(T_ant_save_name,T_ant)
            print('Done.\n\033[32mFile saved as\n',T_ant_save_name,'\033[00m')
            print('It is an array of shape',np.shape(T_ant))
            print('\n\033[94m================ End of function eps.couple2D() ================\033[00m\n')
        return None
    #End of function couple2D    
        
    def visual(self, t_skymap=False, nu_skymap=None, aps=False, n_skymap=False, dndS_plot = False, spectrum=True, antenna = False, xlog=False,ylog=True, fig_ext = 'pdf'):
        '''Plotting function.
        
        This function can produce several figures such as:-
        
        * flux density distribution function
        
        * number density map
        
        * angular power spectrum
        
        * map of foregrounds due to extragalactic point sources (FEPS)
        
        * sky-averaged FEPS as function of frequency
        
        * antenna temperature as function of frequency
        
        Parameters
        ^^^^^^^^^^
        
        t_skymap : bool, optional
            Want to plot the FEPS map (a Mollweide projection plot)? (Default = ``False``).
            
        nu_skymap : float, optional
            Frequency in Hz at which you want to construct the FEPS map. Relevant only when you give ``t_skymap = True``. Currently, only 1 value supported.
            (Default = ``nu_o``)
            
        aps : bool, optional
            Want to plot the angular power spectrum? (Default = ``False``)
            
        n_skymap : bool, optional
            Want to plot the number density map (a Mollweide projection plot)? (Default = ``False``).
            
        dndS : bool, optional
            Want to plot the flux density distribution? (Default = ``False``). The form of :math:`\\mathrm{d}n/\\mathrm{d}S` is set during initialisation.
        
        spectrum : bool, optional
            Want to plot the sky-averaged FEPS? (Default = ``True``).
        
        antenna : bool, optional
            Add the antenna temperature? (Default = ``False``). You should have run :func:`couple2D` to use this.
             
        xlog : bool, optional
            Set the x-axis scale of spectrum plot in log? (Default = ``False``)
            
        ylog : bool, optional
            Set the y-axis scale of spectrum plot in log? (Default = ``True``)
        
        fig_ext : str, optional
            What should be the format of the figure files? Common choices include png, pdf or jpg. (Default = ``pdf``)
        '''
        #-------------------------------------------------------------------------------------
        try:
            comm = MPI.COMM_WORLD
            cpu_ind = comm.Get_rank()
            Ncpu = comm.Get_size()
        except:
            cpu_ind=0
            Ncpu=1
        #-------------------------------------------------------------------------------------
        fs=22
        if cpu_ind==0:
            print('\n\033[94mRunning eps.visual() ...\033[00m\n')
            if Ncpu>1:
                print("\033[91m'visual' does not require parallelisation.\033[00m")
                print("\033[91mYou can run as 'python3 %s'.\033[00m\n" %(sys.argv[0]))
            nu = np.load(self.path+'nu_glob.npy')
            Tb_o_map = np.load(self.path+'Tb_o_map'+self.lbl+'.npy')
            Tb_o_glob = np.mean(Tb_o_map)
            
            Tb_nu_map = np.load(self.path+'Tb_nu_map'+self.lbl+'.npy')
            Tb_nu_glob = np.load(self.path+'Tb_nu_glob'+self.lbl+'.npy')
            
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

                hp.mollview(Tb_plot,title=None,unit=r'$T_{\mathrm{ps}}\,$(K)',cmap=colormaps['coolwarm'],norm = 'log') #,min=0.05,max=200
                hp.graticule()
                
                fig_path = self.path+'Tb_map_'+str(int(nu_skymap/1e6))+'-MHz.'+fig_ext
                plt.savefig(fig_path, bbox_inches='tight')

                print('Done.\n\033[32mTb map saved as:\n',fig_path,'\033[00m\n')
                plt.close()
                
            if aps:
                n_ps = np.load(self.path+'n_ps'+self.lbl+'.npy')
                nbar = self.Nps/self.Npix
                
                fig,ax=plt.subplots(figsize=(8.3,7.5),dpi=300)
                fig.subplots_adjust(left=0.12, bottom=0.07, right=0.88, top=0.97)
                
                if self.amp!=0:
                    del_clus = (n_ps - nbar)/nbar
                    Cl_rec = hp.anafast(del_clus)
                    th = tcl.theta(200)
                    cor = self.acf(th)
                    Cl_clus = tcl.corrtocl(cor)
                    n_poisson = np.random.poisson(lam=nbar,size=self.Npix)
                    ax.loglog(range(200),Cl_clus,'b--',label='Input (clustered)')
                    ax.loglog(range(self.Nside*3),Cl_rec,'limegreen',label='Recovered (clustered)')
                else:
                    n_poisson=n_ps
                
                del_poisson = (n_poisson - nbar)/nbar
                Cl_poisson = hp.anafast(del_poisson)
                ax.loglog(range(self.Nside*3),Cl_poisson,'r:',label='Poisson (isotropic)')
                
                ax.set_xlabel(r'$\ell$',fontsize=fs)
                ax.set_ylabel(r'$C_{\ell}$',fontsize=fs)
                
                ax.yaxis.set_ticks_position('both')
                ax.xaxis.set_ticks_position('both')
                ax.tick_params(axis='both', which='major', length=5, width=1, labelsize=fs,direction='in',pad=8)
                ax.tick_params(axis='both', which='minor', length=3, width=1,direction='in')
                ax.legend(fontsize=18,frameon=False)
                ax.minorticks_on()
                ax.set_xlim([1,200])
                ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
                
                fig_path = self.path+'APS.'+fig_ext
                plt.savefig(fig_path)
                print('Done.\n\033[32mAPSs saved as:\n',fig_path,'\033[00m\n')
                
                plt.close()
                
            if n_skymap:
                print('Creating number density map ...')
                
                n_ps = np.load(self.path+'n_ps'+self.lbl+'.npy')
                nmax = int(np.max(n_ps))+1
                nmin = int(np.min(n_ps))
                
                hp.mollview(n_ps,title=None,unit=r'$n_{\mathrm{ps}}$',min=nmin,max=nmax)
                hp.graticule()
                
                fig_path = self.path+'n_ps.' + fig_ext
                plt.savefig(fig_path, bbox_inches='tight')
                print('Done.\n\033[32mNumber density map saved as:\n',fig_path,'\033[00m\n')
                plt.close()
            
            if dndS_plot:
                print('Creating flux density distribution function ...')
                fig,ax=plt.subplots(figsize=(8.3,7.5),dpi=300)
                fig.subplots_adjust(left=0.12, bottom=0.07, right=0.88, top=0.97)
                
                S_space = np.logspace(self.logSmin,self.logSmax,1000)
                dndS_space = self.dndS(S_space)
                ax.loglog(S_space,dndS_space,'b')

                ax.set_xlabel(r'$S\,$(Jy)',fontsize=fs)
                ax.set_ylabel(r'$\mathrm{d}n/\mathrm{d}S\,(\mathrm{Jy}^{-1}\mathrm{sr}^{-1})$',fontsize=fs)
                ax.tick_params(axis='both', which='major', length=5, width=1, labelsize=fs,direction='in',pad=8)
                ax.tick_params(axis='both', which='minor', length=3, width=1, direction='in')
                ax.minorticks_on()
                ax.set_xlim([10**self.logSmin,10**self.logSmax])
                ax.yaxis.set_ticks_position('both')
                ax.xaxis.set_ticks_position('both')
                ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
                
                fig_path = self.path+'dndS.' + fig_ext
                plt.savefig(fig_path)
                print('Done.\n\033[32mFlux density distribution saved as:\n',fig_path,'\033[00m\n')
                plt.close()
            
            if spectrum:                
                print('\nCreating T vs nu plot ...')
                left=0.12
                fig,ax=plt.subplots(figsize=(8, 7.9))
                fig.subplots_adjust(left=left, bottom=0.06, right=1-left, top=0.94)
                
                if antenna==False:
                    ax.plot(nu/1e6,Tb_nu_glob,color='b',lw=1.5)
                    ax.set_ylabel(r'$\langle T_{\mathrm{ps}}\rangle\,$(K)',fontsize=fs)
                else:
                    T_ant = np.load(self.path+'T_ant'+self.lbl+'.npy')
                    ax.plot(nu/1e6,Tb_nu_glob,color='b',lw=1.5,label=r'$\langle T_{\mathrm{ps}}\rangle$')
                    ax.plot(nu/1e6,T_ant,color='limegreen',lw=1.5,label=r'$T_{\mathrm{A,ps}}$')
                    ax.set_ylabel(r'$T\,$(K)',fontsize=fs)
                    ax.legend(fontsize=18,frameon=False)
                    print('Added the antenna temperature to the figure.')
                
                if xlog:
                    ax.set_xscale('log')
                if ylog:
                    ax.set_yscale('log')
                
                ax.set_xlabel(r'$\nu\,$(MHz)',fontsize=fs)
                ax.minorticks_on()
                ax.yaxis.set_ticks_position('both')
                ax.tick_params(axis='both', which='major', length=5, width=1, labelsize=fs,direction='in')
                ax.tick_params(axis='both', which='minor', length=3, width=1,direction='in')

                secax = ax.secondary_xaxis('top', functions=(_nu2z,_z2nu))
                secax.set_xlabel(r'$z$',fontsize=fs, labelpad=12)
                secax.tick_params(axis='both', which='major', length=5, width=1, labelsize=fs,direction='in')
                secax.tick_params(axis='both', which='minor', length=3, width=1,direction='in')
                secax.minorticks_on()

                #plt.xlim([0.1,1e2])
                #plt.ylim([1,3e4])
                ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
                fig_path = self.path+'T_vs_nu.' + fig_ext
                plt.savefig(fig_path)
                plt.close()
                print('Done.\n\033[32mT vs frequency saved as:\n',fig_path,'\n\033[00m')
            
            print('\n\033[94m================ End of function eps.visual() ================\033[00m\n')
    #End of function visual()
#End of class eps()

