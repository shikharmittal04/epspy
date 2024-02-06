import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colormaps

nu21=1420e6 #21-cm frequency in Hz

def nu2z(nu):
    return nu21/nu-1

def z2nu(z):
    return nu21/(1+z)

np.seterr(all='ignore')


def plotter(Tb_nu,Tb_o=None, nu=None, path='/home/hpcmitt1/rds/hpc-work/point-sources-data/',beta_o=2.681,Tcmb_o=2.725):
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    if np.size(nu)==1:
        #No or only one frequency given, that mean you want to plot only the Mollweide plot at a single frequency.

        if nu!=None: print('\nGenerating the sky map at frequency = {:.2f} MHz ...'.format(nu/1e6))
        hp.mollview(Tb_nu,title=None,unit=r'$T_{\mathrm{b}}^{\mathrm{eg}}\,$(K)',cmap=colormaps['coolwarm'],min=0.05,max=200,norm='log')
        hp.graticule()
        fig_path = path+'Tb_nu_map_'+str(int(nu/1e6))+'-MHz.pdf'
        plt.savefig(fig_path, bbox_inches='tight')
        print('Done. Tb map saved as',fig_path)
    else:
        print('\nCreating Tb vs nu plot ...')
        left=0.12
        fs=22
        fig,ax=plt.subplots(figsize=(8, 7.9))
        fig.subplots_adjust(left=left, bottom=0.06, right=1-left, top=0.94)
        
        try:
            Tb_o_mean = np.mean(Tb_o)
        except:
            print('Error! Suppy the Tb at reference frequency.')
            sys.exit()

        Tb_mean = Tb_o_mean*(nu/nu_o)**-beta_o
        Tb_glob = np.mean(Tb_nu,axis=1)
        
        ax.axhline(y=Tcmb_o,color='k',ls='--',lw=1.5, label='CMB')
        ax.loglog(nu,Tb_mean,color='r',lw=1.5,ls=':',label=r'$\beta= $%.2f'%beta_o)
        ax.loglog(nu,Tb_glob,color='b',lw=1.5,label='Extragalactic')

        ax.set_xlabel(r'$\nu\,$(Hz)',fontsize=fs)
        ax.set_ylabel(r'$T_{\mathrm{b}}^{\mathrm{eg}}\,$(K)',fontsize=fs)

        ax.minorticks_on()
        ax.yaxis.set_ticks_position('both')
        ax.tick_params(axis='both', which='major', labelsize=fs)
        ax.legend(fontsize=18,loc=1,frameon=False)

        secax = ax.secondary_xaxis('top', functions=(nu2z,z2nu))
        secax.set_xlabel(r'$z$',fontsize=fs, labelpad=12)
        secax.tick_params(which='major', labelsize=fs)

        #plt.xlim([0.1,1e2])
        #plt.ylim([1,3e4])
        ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
        fig_path = path+'Tb_vs_nu.pdf'
        plt.savefig(fig_path)
        print('Done. Tb vs frequency saved as',fig_path)
        
nu_o=150e6
nu = 1e6*np.arange(50,200)
path='/home/hpcmitt1/rds/hpc-work/point-sources-data/'

Tb_o = np.load(path+'Tb_o.npy')
Tb_nu = np.load(path+'Tb_nu.npy')

plotter(Tb_nu=Tb_nu,Tb_o=Tb_o,nu=nu)

plotter(Tb_nu=Tb_nu[149,:],nu=199e6)
