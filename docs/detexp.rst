Detailed explanation
--------------------

Initialisation
^^^^^^^^^^^^^^
The first step to use this package is to initialise the properties of the unresolved sources. This is done using the class 
:class:`furs.furs`. If you give no arguments default settings are assumed.

There are total 11 available optional arguments. So if you want to chose a different set of parameters, you must
initialise the object as

.. code:: python

   from furs import furs
   obj = furs.furs(log2Nside=6, logSmin=-2,logSmax=-1,dndS_form=0, nu_o=150e6, beta_o=2.681,sigma_beta=0.5, amp=7.8e-3,gam=0.821, path='', lbl='')

(Replace the above values by values of your choice.) See also the API reference section.

Reference frequency
^^^^^^^^^^^^^^^^^^^

The function :func:`ref_freq` does 3 tasks:-

1. Calcuates the total number of unresolved sources corresponding to your specified ``logSmin`` and ``logSmin``. Then it creates a 'clustered' sky density of unresolved radio sources, fluctuation for which follows the 2PACF whose parameters are set by ``amp`` and ``gam``.
   
2. Next, it visits each pixel on the sky and assigns each source a flux density chosen from a flux distribution function, :math:`\mathrm{d}n/\mathrm{d}S`

3. Finally, it assigns a spectral index to each source which is normally distributed. The normal distribution is set by ``beta_o`` and ``sigma_beta``.

Sky pixelisation is set by ``log2Nside``. The number of pixels is
:math:`N_{\mathrm{pix}} = 12\times 2^{2k}`, where :math:`k=` ``log2Nside``.

The function does not return anything, but produces 4 output files,
namely ``n_clus.npy``, ``Tb_o_individual.npy``, ``Tb_o_map.npy``, and
``beta.npy`` in the path specified by ``path`` during initialisation.
The files are described below.

- ``n_clus.npy`` is a 1D array which stores number density of unresolved radio sources as number per pixel. ``n_clus[i]`` gives the number of sources on the :math:`i^{\mathrm{th}}` pixel, where :math:`i=0,1,\ldots,N_{\mathrm{pix}}-1`. Note that in general ``n_clus[i]`` will not be a natural number; we simulate for a rounded-off value.

- Both ``Tb_o_individual.npy`` and ``beta.npy`` are array of arrays of unequal sizes and share equal amount of memory. Typically, these files will be huge (for default settings they will of size ~ 17 MB each). Each of ``Tb_o_individual[0]``, ``Tb_o_individual[1]``, ..., is an array and they are total :math:`N_{\mathrm{pix}}` in number corresponding to :math:`N_{\mathrm{pix}}` pixels. The last array is ``Tb_o_individual[Npix-1]``. The length of array ``Tb_o_individual[i]`` is equal to the number of sources on the :math:`i^{\mathrm{th}}` pixel, which is ``round(n_clus[i])``. The values itself are the brightness temperature contributed by each source in kelvin at reference frequency.

- The structure of ``beta`` is same as ``Tb_o_individual``. The values itself are the spectral indices assigned to each source.

- ``Tb_o_map`` is an array similar in structure to ``n_clus``. It is the pixel wise brightness temperature contributed by the extragalactic radio sources at the reference frequency. Thus, ``Tb_o_map[i] = numpy.sum(Tb_o_individual[i])``.


.. code:: python

   from furs import furs
   
   obj = furs.furs()
   
   obj.ref_freq()

General frequency
^^^^^^^^^^^^^^^^^

The next important task is performed by the function :func:`gen_freq`. It
scales the brightness temperature at reference frequency for each source
according to a power law to a desired range of frequencies. The desired
frequencies should be supplied (in Hz) as a :mod:`numpy` array to this
function. For example

.. code:: python

   from furs import furs
   
   obj = furs.furs()
   
   obj.ref_freq()

   obj.gen_freq(nu = 1e6*numpy.arange(50,201))

The default value is as given in the above command. This function does
not return anything but produces 3 files namely ``Tb_nu_map.npy``,
``Tb_nu_glob.npy``, and ``nu_glob.npy`` in the path specified by
``path`` during initialisation. The files are described below.

1. ``Tb_nu_map`` is a 2D array of shape :math:`N_{\mathrm{pix}}\times
   N_{\nu}`, so that ``Tb_nu_map[i,k]`` gives the brightness temperature
   on the :math:`i^{\mathrm{th}}` pixel at ``nu[k]`` frequency. :math:`N_{\nu}` is
   the number of frequencies you gave in the argument of ``gen_freq()``.

2. ``Tb_nu_glob`` is derived directly from ``Tb_nu_map``. It is the sky
   average of the map at each frequency and is thus a 1D array. It is
   calculated as ``Tb_nu_glob = numpy.mean(Tb_nu_map,axis=0)``.

3. ``nu_glob.npy`` is simply the frequency array you gave else it is the
   default value.

Note that this function loads ``Tb_o_individual.npy`` and ``beta.npy``.
These files can easily be 10s of GB in size for 'realistic' ``logSmin``
and ``logSmax``. Common personal computers have ~ 4 GB RAM. It is thus
recommended to run this code on supercomputers. For job submission scipt
users are requested to specify ``#SBATCH --mem-per-cpu=[size in MB]``,
where a recommendation for ``size in MB`` will be printed by
:func:`ref_freq` function.

Chromatic distortions
^^^^^^^^^^^^^^^^^^^^^

``Tb_nu_map`` and hence ``Tb_nu_glob`` so generated do NOT account for
chromatic distortions. They are simply the model outputs for foregrounds
due to unresolved radio sources. However, in reality because of the
chromatic nature of the antenna beam the actual foregrounds spectrum
registered will be different. You can use the function
:func:`chromatisize()` to account for the chromaticity.

Since this is experiment specific you will need to provide an external
data file: the beam directivity pattern, :math:`D`. This should be a 2D array
of shape :math:`N_{\mathrm{pix}}\times N_{\nu}`, such that ``D[i,k]`` should
give the beam directivity at :math:`i^{\mathrm{th}}` pixel at ``nu[k]'' frequency.
The frequencies at which you generate your data :math:`D` should be the same
as the frequencies you gave in ``gen_freq()``. (In case you forgot,
:func:`gen_freq` will have saved the frequency array in your ``obj.path``
path.) Put this array :math:`D` in your ``obj.path`` path by the name of
``D.npy``.

Only after running :func:`ref_freq` and :func:`gen_freq`, run :func:`chromatisize`
as

.. code:: python

   from furs import furs

   obj = furs.furs()

   obj.ref_freq()

   obj.gen_freq()
   
   #If you have already ran ref_freq and gen_freq previously then comment
   #obj.ref_freq() and obj.gen_freq(). 
   obj.chromatisize()

No input argument is required. The return value is ``None``. This
function will generate a file called ``T_ant.npy`` in your path. This
will be a 1D array with length of number of frequencies.

Visualisation
^^^^^^^^^^^^^

The final part of the code is to visualise the results. Main data for
inspection is in the file ``Tb_nu_map.npy``. Each of ``Tb_nu_map[:,k]``
is an array in the standard ring ordered ``HEALPix`` format and is thus
ready for visualisation as a Mollweide projection. You may also be
interested in inspecting the global spectrum of extragalactic emission,
i.e, temperature as a function of frequency. This is simply the data in
the file ``Tb_nu_glob.npy`` generated by :func:`gen_freq`.

You may use the function :func:`visual` for both the above purposes. It is
possible to make several other additional figures by simply setting the
optional arguments to ``True`` (see below). This function is again a
method of class object :class:`furs.furs` and is thus called as

.. code:: python
   
   obj = furs.furs()

   obj.ref_freq()

   obj.gen_freq()

   obj.chromatisize()

   #comment out obj.ref_freq(), obj.gen_freq(), obj.chromatisize() if you have already run them.
   obj.visual()

For all the available options for this function see the API reference section. This function will produce figures in the path specficied during initialisation.
