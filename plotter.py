import healpy as hp
import numpy as np
import matplotlib.pyplot as plt

def nu2z(nu):
    return 1420/nu-1

def z2nu(z):
    return 1420/(1+z)

np.seterr(all='ignore')


def plotter(nu,Tb_nu,path='/home/hpcmitt1/rds/hpc-work/point-sources-data/',beta0=2.681,To=2.725):
    Nnu=np.size(nu)
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    if Nnu==1:
        #Since it's only one frequency, I am making a plot as well.		
        hp.mollview(Tb_nu,title=None,unit=r'$T_{\mathrm{b}}(\nu)\,$(K)',cmap=colormaps['coolwarm'],norm='log')
        hp.graticule()
        fig_path = path+'Tb_nu_map.pdf'
        plt.savefig(fig_path, bbox_inches='tight')
        print('\nTb map saved as',fig_path)
    else:
        fig,ax=plt.subplots(figsize=(11, 11))
        fig.subplots_adjust(left=0.12, bottom=0.1, right=0.88, top=0.9)
        Tb_mean = Tb_o*(nu/nu_o)**-beta_o
        Tb_glob = np.mean(Tb_nu,axis=1)

        ax.axhline(y=To,color='k',ls='--',lw=1.5, label='CMB')
        ax.loglog(nu,Tb_mean,color='r',lw=1.5,ls=':')
        ax.loglog(nu,Tb_glob,color='b',lw=1.5,label='Extragalactic')

        ax.set_xlabel(r'$\nu\,$(MHz)',fontsize=24)
        ax.set_ylabel(r'$T_{\mathrm{b}}\,$(K)',fontsize=24)

        ax.minorticks_on()
        ax.yaxis.set_ticks_position('both')
        ax.tick_params(axis='both', which='major', labelsize=22)
        ax.legend(fontsize=18,loc=1,frameon=False)

        secax = ax.secondary_xaxis('top', functions=(nu2z,z2nu))
        secax.set_xlabel(r'$z$',fontsize=24, labelpad=12)
        secax.tick_params(which='major', labelsize=22)

        #plt.xlim([0.1,1e2])
        #plt.ylim([1,3e4])
        ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
        fig_path = path+'Tb_vs_nu.pdf'
        plt.savefig(fig_path)
        print('\nTb vs frequency saved as',fig_path)
        plt.savefig(fig_path)

nu=np.arange(50,200)
path='/home/hpcmitt1/rds/hpc-work/point-sources-data/Tb_nu.npy'
Tb_nu = np.load(path)

plotter(nu,Tb_nu)
