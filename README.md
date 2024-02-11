Use this code to generate the extragalactic Foregrounds due to Unresolved Radio Sources (FURS).

If you use this code please cite Mittal et al (2024).

##Why do you need this code?
A cosmological global 21-cm signal hides under foregrounds due to galactic and extragalactic emission. These foregrounds can easily be 4 to 5 orders of magnitude higher than the signal of interest. For a reliable inference it is important to accurately model these foregrounds. While we have a reasonable understanding of galactic emission (typically log-log polynomial), we do not understand the extragalactic contributions. Based on previous work, this code models the foregrounds due to unresolved extragalactic radio sources.

##Quick start
The code is run in two main steps:
* Assign the unresolved sources flux densities (at a chosen reference frequency) and spectral indices.
* Then generate the sky maps a desired frequencies of observations.

The following code captures the main functionality of this package.
```
from furs import furs
obj = furs.extragalactic()

obj.ref_freq()

obj.gen_freq()
```
Save the above code as (say) `eg_script.py` and run it as

`mpirun -np 4 python eg_script.py`

##Detailed explanation
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
    * spread in the Gaussian distribution of spectral index
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
6. `low`
    * $\log_{10}(S_{\mathrm{min}})$, where $S_{\mathrm{min}}$ is in Jansky (Jy)
    * type *float*
    * default **-6.0**
7. `upp`
    * $\log_{10}(S_{\mathrm{max}})$, where $S_{\mathrm{max}}$ is in Jansky (Jy)
    * type *float*
    * default **-1.0**
8. `path`
    * path where you would like to save all output files
    * type *string*
    * default **''**
9. `log2Nside`
    * Number of divisions of each side of the pixel for HEALPix maps in units of log_2
    * type *int*
    * default **6**

Thus, if you want to chose a different set of parameters, you must initialise the object as

`obj = furs.extragalactic(log2Nside=6, low=-6,upp=-1, nu_o=150e6, beta_o=2.681,sigma_beta=0.5, amp=7.8e-3,gam=0.821, path='')`.
(Replace the above values by values of your choice)

The function `ref_freq` does 3 tasks

* Calcuates the total number of unresolved sources corresponding to `low` and `upp`
* Creates a 'clustered' sky density of unresolved radio sources, fluctuation for which follows the 2PACF whose parameters are set by `amp` and `gam`
* Finally, it visits each pixel on the sky and assigns each source a flux density chosen from a sum-of-two-double-inverse-power-law function (see eq.() in Mittal et al (2024)) and a spectral index which is normally distributed. The normal distribution is set by `beta_o` and `sigma_beta`.

Sky pixelisation is set by `log2Nside`. The number of pixels is $N_{\mathrm{pix}} = 12 \times 2^{2k}$, where $k=$`log2Nside`.
 



