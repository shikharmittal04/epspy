General remarks
---------------

Users do not have to run :func:`ref_freq()` everytime. If they want to use the same data for source distribution (``n_clus.npy``), flux density (``Tb_o_individual.npy``) and spectral index (``beta.npy``) assignments at reference frequency to generate spectrum and sky maps for a different frequency range, then run only :func:`gen_freq()` for a new choice of ``nu``.

Similarly, if you have already run :func:`gen_freq()` and are happy with the specifications of the model then you can directly jump to the :func:`visual()` function.

In case you forgot what data set you generated with what parameter specifications, you can always save your class object using the function :func:`save_furs` in the directory where all other outputs are saved and load back using :func:`load_furs`. (Both functions are part of module ``furs.py``.)

Thus, after initialising your class object (i.e. ``obj = furs.furs([YOUR SPECIFICATIONS])``), you can add to your script ``furs.save_furs(obj,'myobj')``.

**Examples**

.. code:: python
   
   from furspy import furs
   
   obj = furs.furs()
   furs.save_furs(obj,'myobj')

So when you came back next time you can load it as

.. code:: python
   
   from furspy import furs
   obj=furs.load_furs('/give/full/path/to/myobj.pkl')

Remember to give the full path to the ``myobj`` with the extension ``.pkl``. 

You may now check that indeed the specfications are correctly loaded by printing them using function :func:`print_input()`.

.. code:: python
   
   from furspy import furs
   obj=furs.load_furs('/give/full/path/to/myobj.pkl')
   obj.print_input()


There is also an argument ``lbl``, which you can use to put an extra label to you output files. For example,

.. code:: python
   
   from furspy import furs
   
   obj = furs.furs(lbl='mylabel')

