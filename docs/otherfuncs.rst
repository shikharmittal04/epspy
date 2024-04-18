More functions
--------------

Functions of class :class:`eps.eps`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The main functionality of the package was already discussed in previous sections. Here we will learn about 5 additional useful methods of the class :class:`eps.eps`. These are 

1. The 2-point angular correlation function (2PACF, :math:`C=C(\chi)`), :func:`acf`

2. The flux density distribution function (:math:`\mathrm{d}n/\mathrm{d}S`), :func:`dndS`. As this function is a method of class :class:`eps.eps`, the choice of the form of :math:`\mathrm{d}n/\mathrm{d}S` will have been set when you initialise your class :class:`eps.eps` object as

.. code:: python
   
   from furpy import eps
   
   obj = eps.eps()

3. Number of source per pixel ( :math:`n_{\mathrm{ps}}` ), :func:`num_den()`. This function has direct correspondence with the form of :math:`\mathrm{d}n/\mathrm{d}S` and the choice of :math:`S_{\mathrm{min}}` and :math:`S_{\mathrm{max}}`.
   
4. :func:`print_input()` can be useful if you want see what parameters you are currently running with.

Note that these functions are computationally cheap and thus can be run directly in a jupyter notebook or interactively on the terminal.

**Examples**

Suppose you want to find the number of sources between :math:`S=10^{-6}` and :math:`10^{-1}\mathrm{Jy}` for a :math:`7^{\mathrm{th}}` order log-log polynomial fit to Euclidean number count then do

.. code::

   >>> from furspy import eps
   >>> obj = eps.eps(dndS_form=1, logSmin=-6,logSmax=-1)
   >>> Nps = round(obj.Nps)
   >>> Nps
   249825820

Let us calculate the number density.

.. code::

   >>> nps = obj.num_den()
   Total number of sources, Nps = 249825820
   Total number of pixels, Npix = 49152
   Average number of sources per pixel, n_bar = 5082.72
   Done.
   Average overdensity for the clustered sky (should be ~ 0) = 0.007
   The clustered number density has been saved into file:
    n_ps 

(Note that since this is a random process you might get different numbers.)

Let us also check if the sum of the elements of array ``nps`` is ``Nps``

.. code::

   >>> ncl.sum()
   251605169.7207427

There some discrepancy but do not be worried as this is just a numerical artefact. In fact the fractional error is already printed and is about 0.7%, which is sufficiently small.


Functions of module :mod:`eps`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In case you forgot what data set you generated with what parameter specifications, you can always save your class object using the function :func:`save_eps` in the directory where all other outputs are saved and load back using :func:`load_eps`. (Both functions are part of module ``eps.py``.)

Thus, after initialising your class object (i.e. ``obj = eps.eps([YOUR SPECIFICATIONS])``), you can add to your script ``eps.save_eps(obj,'myobj')``.

**Examples**

.. code:: python
   
   from furspy import eps
   
   obj = eps.eps()
   eps.save_eps(obj,'myobj')

Now check if your there is a file called ``myobj.pkl`` in ``obj.path`` directory. 

When you came back next time you can load you class object as

.. code:: python
   
   from furspy import eps
   obj=eps.load_eps('/give/full/path/to/myobj.pkl')

Remember to give the full path to the ``myobj`` with the extension ``.pkl``. 

You may now check that indeed the specifications are correctly loaded by printing them using function :func:`print_input()`.

.. code:: python
   
   from furspy import eps
   obj=eps.load_eps('/give/full/path/to/myobj.pkl')
   obj.print_input()


There is also an argument ``lbl``, which you can use to put an extra label to you output files. For example,

.. code:: python
   
   from furspy import eps
   
   obj = eps.eps(lbl='_mylabel')

Now all files names will have `_mylabel` appended to them. For example, when you run :func:`num_den`, the output file name will be called ``n_ps_mylabel``.
