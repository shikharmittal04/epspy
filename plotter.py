import healpy as hp

#Since it's only one frequency, I am making a plot as well.		
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
hp.mollview(Tb_nu_final,title=None,unit=r'$T_{\mathrm{b}}(\nu)\,$(K)',cmap=colormaps['coolwarm'],norm='log')
hp.graticule()
fig_path = path+'Tb_nu.pdf'
plt.savefig(fig_path, bbox_inches='tight')
print('\nTb map saved as',fig_path)
