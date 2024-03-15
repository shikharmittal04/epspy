Other useful functions
----------------------

The main functionality of the package was already discussed in previous sections. Here we will learn about 5 additional useful methods of the class :class:`furs.furs`. These are 

1. The 2-point angular correlation function (2PACF, :math:`C=C(\chi)`), :func:`acf`

2. The flux density distribution function (:math:`\mathrm{d}n/\mathrm{d}S`), :func:`dndS`. As this function is a method of class :class:`furs.furs`, the choice of the form of :math:`\mathrm{d}n/\mathrm{d}S` will have been set when you initialise your class :class:`furs.furs` object as

.. code:: python
   
   from furpy import furs
   
   obj = furs.furs()

3. Total number of sources on the full sky ( :math:`N_{\mathrm{s}}` ), :func:`num_sources()`. This function has direct correspondence with the form of :math:`\mathrm{d}n/\mathrm{d}S` and the choice of :math:`S_{\mathrm{min}}` and :math:`S_{\mathrm{max}}`.  

4. Number of source per pixel ( :math:`n_{\mathrm{cl}}` ), :func:`num_den()`. As in the :func:`num_sources()`, the total number density returned will correspond to the chosen form of :math:`\mathrm{d}n/\mathrm{d}S` the choice of :math:`S_{\mathrm{min}}` and :math:`S_{\mathrm{max}}`.
   
5. :func:`print_input()` can be useful if you want see what parameters you are currently running with.

Note that these functions are computationally cheap and thus can be run directly in a jupyter notebook or interactively on the terminal.

**Examples**

Suppose you want to find the number of sources between :math:`S=10^{-6}` and :math:`10^{-1}\mathrm{Jy}` for a :math:`7^{\mathrm{th}}` order log-log polynomial fit to Euclidean number count then do

.. code::

   >>> from furspy import furs
   >>> obj = furs.furs(dndS_form=1, logSmin=-6,logSmax=-1)
   >>> Ns = obj.num_sources()
   >>> Ns
   249825819.67727068

Let us calculate the number density.

.. code::

   >>> ncl = obj.num_den()
   Total number of sources, Ns = 249825820
   Total number of pixels, Npix = 49152
   Average number of sources per pixel, n_bar = 5082.72
   Done.
   Average overdensity for the clustered sky (should be ~ 0) = 0.007
   The clustered number density has been saved into file:
    n_clus 

(Note that since this is a random process you might get different numbers.)

Let us also check if the sum of the elements of array ``ncl`` is ``Ns``

.. code::

   >>> ncl.sum()
   251605169.7207427

There some discrepancy but do not be worried as this is just a numerical artefact. In fact the fractional error is already printed and is about 0.7%, which is sufficiently small.

