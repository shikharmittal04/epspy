Use this code to generate the extragalactic Foregrounds due to Unresolved Radio Sources (FURS).

If you use this code please cite Mittal et al (2024).

## Why do you need this code?
A cosmological global 21-cm signal hides under foregrounds due to galactic and extragalactic emissions. These foregrounds can easily be 4 to 5 orders of magnitude higher than the signal of interest. For a reliable inference it is important to accurately model these foregrounds. While we have a reasonable understanding of galactic emission (typically fit as log-log polynomial), we do not understand the extragalactic contributions. Based on previous work, this code models the foregrounds due to unresolved extragalactic radio sources.

## Installation and requirements
This package can be installed as 

```
pip install furs
```

It is recommended to work on a Python version > 3.8. Additional packages required are `healpy` and `transformcl`.

## Quick start
The code is run in two main steps:
* Assign the unresolved sources flux densities (at a chosen reference frequency) and spectral indices.
* Then generate the sky maps at desired frequencies of observation.

The following code captures the main functionalities of this package.
```
from furs import furs
obj = furs.extragalactic()

obj.ref_freq()

obj.gen_freq()
```
Save the above code as (say) `eg_script.py` and run it as

`python eg_script.py`

The default values have been chosen such that the code can be run on a PC. Since modern PCs have at least 4 cores, for a better performance one could also run the code as

`mpirun -np 4 python eg_script.py`

However, in general and for more realistic flux density ranges and high resolution maps, it is recommended to run the code on HPCs.

## Detailed explanation
#### Initialisation
`furs.extragalactic()` initialises the class object with default settings. There are total 9 available optional arguments as follows:

1. `nu_o`
    * reference frequency, $\nu_0$ (Hz)
    * type *float*
    * default **150e6**
2. `beta_o`
    * mean spectral index, $\beta_0$
    * type *float*
    * default **2.681**
3. `sigma_beta`
    * spread in the Gaussian distribution of spectral index, $\sigma_\beta$
    * type *float*
    * default **0.5**
4. `amp`
    * amplitude of the power-law 2-point angular correlation function (2PACF), $A$
    * type *float*
    * default **7.8e-3**
5. `gam`
    * negative of the exponent of the 2PACF, $\gamma$
    * type *float*
    * default **0.821**
6. `logSmin`
    * $\log_{10}(S_{\mathrm{min}})$, where $S_{\mathrm{min}}$ is in Jansky (Jy)
    * type *float*
    * default **-2.0**
7. `logSmax`
    * $\log_{10}(S_{\mathrm{max}})$, where $S_{\mathrm{max}}$ is in Jansky (Jy)
    * type *float*
    * default **-1.0**
8. `path`
    * path where you would like to save all output files
    * type *string*
    * default **''**
9. `log2Nside`
    * Number of divisions of each side of the pixel for `HEALPix` maps in units of log_2
    * type *int*
    * default **6**

Thus, if you want to chose a different set of parameters, you must initialise the object as

`obj = furs.extragalactic(log2Nside=6, logSmin=-2,logSmax=-1, nu_o=150e6, beta_o=2.681,sigma_beta=0.5, amp=7.8e-3,gam=0.821, path='')`

(Replace the above values by values of your choice.)

#### Reference frequency
The function `ref_freq` does 3 tasks:-

* Calcuates the total number of unresolved sources corresponding to your specified `logSmin` and `logSmin`.
* Creates a 'clustered' sky density of unresolved radio sources, fluctuation for which follows the 2PACF whose parameters are set by `amp` and `gam`.
* Finally, it visits each pixel on the sky and assigns each source a flux density chosen from a flux distribution function, $\mathrm{d}n/\mathrm{d}S$ (see eq. in Mittal et al 2024) and a spectral index which is normally distributed. The normal distribution is set by `beta_o` and `sigma_beta`.

Sky pixelisation is set by `log2Nside`. The number of pixels is $N_{\mathrm{pix}} = 12 \times 2^{2k}$, where $k=$`log2Nside`.
 
The function does not return anything, but produces 4 output files, namely `n_clus.npy`, `Tb_o_individual.npy`, `Tb_o_map.npy`, and `beta.npy` in the path specified by `path` during initialisation. The files are described below.

1. `n_clus.npy` is a 1D array which stores number density of unresolved radio sources as number per pixel. `n_clus[i]` gives the number of sources on the $i^{\mathrm{th}}$ pixel, where $i=0,1,\ldots,N_{\mathrm{pix}}-1$. Note that in general `n_clus[i]` will not be a natural number; we simulate for a rounded-off value.

2. Both `Tb_o_individual.npy` and `beta.npy` are array of arrays of unequal sizes and share equal amount of memory. Typically, these files will be huge (for default settings they will of size ~ 17 MB each). Each of `Tb_o_individual[0]`, `Tb_o_individual[1]`, ..., is an array and they are total $N_{\mathrm{pix}}$ in number corresponding to $N_{\mathrm{pix}}$ pixels. The last array is `Tb_o_individual[Npix-1]`. The length of array `Tb_o_individual[i]` is equal to the number of sources on the $i^{\mathrm{th}}$ pixel, which is `round(n_clus[i])`. The values itself are the brightness temperature contributed by each source in kelvin at reference frequency.

3. The structure of `beta` is same as `Tb_o_individual`. The values itself are the spectral indices assigned to each source.

4. `Tb_o_map` is an array similar in structure to `n_clus`. It is the pixel wise brightness temperature contributed by the extragalactic radio sources at the reference frequency. Thus, `Tb_o_map[i] = numpy.sum(Tb_o_individual[i])`. 

#### General frequency
The next important task is performed by the function `gen_freq`. It scales the brightness temperature at reference frequency for each source according to a power law to a desired range of frequencies. The desired frequencies should be supplied (in Hz) as a numpy array to this function. For example
```
obj.gen_freq(nu = 1e6*numpy.arange(50,200))
```
The default value is as given in the above command. This function does not return anything but produces 3 files namely `Tb_nu_map.npy`, `Tb_nu_glob.npy`, and `nu_glob.npy` in the path specified by `path` during initialisation. The files are described below.

1. `Tb_nu_map` is a 2D array of shape $N_{\mathrm{pix}}\times N_{\nu}$, so that `Tb_nu_map[i,j]` gives the brightness temperature on the $i^{\mathrm{th}}$ pixel at `nu[j]` frequency. $N_{\nu}$ is the number of frequencies you gave in the argument of `gen_freq()`.

2. `Tb_nu_glob` is derived directly from `Tb_nu_map`. It is the sky average of the map at each frequency and is thus a 1D array. It is calculated as `Tb_nu_glob = numpy.mean(Tb_nu_map,axis=0)`.

3. `nu_glob.npy` is simply the frequency array you gave else it is the default value.

Note that this function loads `Tb_o_individual.npy` and `beta.npy`. These files can easily be 10s of GB in size for 'realistic' `logSmin` and `logSmax`. Common personal computers have ~ 4 GB RAM. It is thus recommended to run this code on supercomputers. For job submission scipt users are requested to specify `#SBATCH --mem-per-cpu=[size in MB]`, where a recommendation for `size in MB` will be printed by `ref_freq()` function.


#### Visualisation
The final part of the code is to visualise the results. Main data for inspection is in the file `Tb_nu_map.npy`. Each of `Tb_nu_map[:,j]` is an array in the standard ring ordered `HEALPix` format and is thus ready for visualisation as a Mollweide projection. You may also be interested in inspecting the global spectrum of extragalactic emission, i.e, temperature as a function of frequency. This is simply the data in the file `Tb_nu_glob.npy` generated by `gen_freq()`.

You may use the function `visual()` for both the above purposes. This function is again a method of class object extragalactic and is thus called as

```
obj = furs.extragalactic()
obj.visual()
```

The following optional arguments are available for this function:-
1. `nu_skymap`
	* the frequency at which you want to produce a Mollweide projection of extragalactic foregrounds
	* type *float*
	* default `nu_o`
2. `t_skymap`
	* Create a sky map of extragalactic foregrounds?
	* type *bool*
	* default `False`
3. `spectrum`
	* Create the foreground spectrum?
	* type *bool*
	* default `True`
4.	`n_skymap`
	* Create a sky map of number density of unresolved radio sources?
	* type *bool*
	* default `False`
5. `xlog`
	* Set x-axis in log scale?
	* type *bool*
	* default `False`
6. `ylog`
	* Set y-axis in log scale?
	* type *bool*
	* default `True`
7.  `fig_ext`
	* Choose your format of figure file; popular choices include `pdf`, `jpeg`, `png`
	* type *string*
	* default `pdf`

This function will produce figures in the path specficied during initialisation. 

## Other functions

There are 4 additional useful methods of the class `extragalactic`. These are:-

1. `acf(chi)`
	* returns the 2PACF, $C(\chi)$
	* requires one argument, the angle $\chi$ in radians; can be a number or an array
	* output is a dimensionless quantity
2. `dndS(S)`
	* returns flux distribution, $\mathrm{d}n/\mathrm{d}S$
	* requires one argument, the flux density $S$ in Jy; can be a number or an array
	* output is in units of number per unit flux density per unit solid angle, i.e. $\mathrm{Jy^{-1}sr^{-1}}$
3. `num_den()`
	* returns the clustered number density as number per pixel, $n_{\mathrm{clus}}$
	* no arguments required
	* output is an array of length $N_{\mathrm{pix}}$
4. `num_sources()`
	* returns the total number of unresolved extragalactic radio sources for the full sky, $N_{\mathrm{s}}$
	* no arguments required
	* output is a pure number
5. `print_input()`
	* If you want to print the all raw parameter values you gave, you may use this function to print them
	* no arguments required
	* no return value

Example usage: to find the number of sources between $10^{-6}$ and $10^{-1}\mathrm{Jy}$ do
```
obj = furs.extragalactic(logSmin=-6,logSmax=-1)
Ns = obj.num_sources()
```

## General remarks

Users do not have to run `ref_freq()` everytime. If they want to use the same data for source distribution (`n_clus.npy`), flux density (`Tb_o_individual.npy`) and spectral index (`beta.npy`) assignments at reference frequency to generate spectrum and sky maps for a different frequency range, then run only `gen_freq()` for a new choice of `nu`.
