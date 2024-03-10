Basics
======

Overview
--------

:Name: Foregrounds due to Unresolved Radio Sources
:Author: Shikhar Mittal
:Homepage: https://github.com/shikharmittal04/furs

Why do you need this code?
--------------------------

Use this code to generate the Foregrounds due to Unresolved Radio Sources (FURS).

A cosmological global 21-cm signal hides under foregrounds due to galactic and extragalactic emissions. These foregrounds can easily be 4 to 5 orders of magnitude higher than the signal of interest. For a reliable inference it is important to accurately model these foregrounds. While we have a reasonable understanding of galactic emission (typically fit as log-log polynomial), we do not understand the extragalactic contributions. Based on existing models, this code models the foregrounds due to unresolved extragalactic radio sources.

Read more about it in the paper `Mittal et al (2024) <https://arxiv.org/abs/2311.03447>`_.

Installation and requirements
-----------------------------

This package can be installed as

.. code:: bash

   pip install furs

It is recommended to work on a Python version > 3.8. Packages required are 

- `numpy <https://pypi.org/project/numpy/>`_
- `scipy <https://pypi.org/project/scipy/>`_
- `matplotlib <https://pypi.org/project/matplotlib/>`_
- `mpi4py <https://pypi.org/project/mpi4py/>`_
- `healpy <https://pypi.org/project/healpy/>`_
- `transformcl <https://pypi.org/project/transformcl/>`_


Quick start
-----------

The code is run in two main steps:

-  Assign the unresolved sources flux densities (at a chosen reference frequency) and spectral indices.
-  Then generate the sky maps at desired frequencies of observation.

The following code captures the main functionalities of this package.

.. code:: python

   from furs import furs

   #Step-1 initialise the object with default settings
   obj = furs.furs()

   #Step-2 generate the data at the reference frequency
   obj.ref_freq()

   #Step-3 generate the sky maps at multiple frequencies as well as their sky average
   obj.gen_freq()

   #Step-4 finally, generate a sky averaged spectrum vs frequency figure
   obj.visual()


Save the above code as (say) ``eg_script.py`` and run it as

.. code:: bash

    python eg_script.py

Running the code will generate several files. The terminal messages will
guide you to these output files. The most important of all files of your
interest will be ``Tb_nu_map.npy``. However, you may never have to deal
with them yourself. To visualise your outputs use the function
:func:`visual`. See :ref:`api` for the available features for :func:`visual`.

The default values have been chosen such that the above script can be
run on a PC. Since modern PCs have at least 4 cores, for a better
performance one could also run the code as

.. code:: bash

    mpirun -np 4 python eg_script.py

However, in general and for more realistic flux density ranges and high
resolution maps, it is recommended to run the code on HPCs.

License and citation
--------------------
The software is free to use on the MIT open source license. If you use the software then please consider citing `Mittal et al (2024) <https://arxiv.org/abs/2311.03447>`_.
