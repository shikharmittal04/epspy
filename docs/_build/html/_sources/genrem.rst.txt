General remarks
---------------

Users do not have to run ``ref_freq()`` everytime. If they want to use
the same data for source distribution (``n_clus.npy``), flux density
(``Tb_o_individual.npy``) and spectral index (``beta.npy``) assignments
at reference frequency to generate spectrum and sky maps for a different
frequency range, then run only ``gen_freq()`` for a new choice of
``nu``.

Similarly, if you have already run ``gen_freq()`` and are happy with the
specifications of the model then you can directly jump to the
``visual()`` function.

In case you forgot what data set you generated with what specifications,
you can always save your class object using the function
``save_furs(class_object,'file_name.pkl')`` in the directory where all
other outputs are saved and load back using ``load_furs``. (Both
functions are part of module ``furs.py``.)

Thus, after initialising your class object (i.e.
``obj = furs([YOUR SPECIFICATIONS])``), you can add to your script

::

   furs.save_furs(obj,'myobj.pkl')

So when you came back next time you can load it as

::

   obj=furs.load_furs('/give/full/path/to/myobj.pkl')

You can check that indeed the specfications are correctly loaded by
printing them via command ``obj.print_input()``.
