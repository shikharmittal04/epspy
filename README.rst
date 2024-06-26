Overview
--------

:Name: Extragalactic Point Sources
:Author: `Shikhar Mittal <https://sites.google.com/view/shikharmittal/home>`_
:Paper: `Mittal et al. (2024) <https://arxiv.org/abs/2406.17031>`_

Why do you need this code?
--------------------------

Use this code to generate the foregrounds due to extragalactic point sources.

A cosmological global 21-cm signal hides under foregrounds due to galactic and extragalactic emissions. These foregrounds can easily be 4 to 5 orders of magnitude higher than the signal of interest. For a reliable inference it is important to accurately model these foregrounds. While we have a reasonable understanding of galactic emission (typically fit as log-log polynomial), we do not understand the extragalactic contributions. This package models the foregrounds due to extragalactic radio sources.

Read more about it in the paper `Mittal et al (2024) <https://arxiv.org/abs/2406.17031>`_.

Installation and requirements
-----------------------------

This package can be installed as

.. code:: bash

   pip install epspy

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

-  Assign the point sources flux densities (at a chosen reference frequency) and spectral indices.
-  Then generate the sky maps at desired frequencies of observation.

The following code captures the main functionalities of this package.

.. code:: python

   from epspy import meps

   #Step-1 initialise the object with default settings
   obj = meps.eps()

   #Step-2 generate the data at the reference frequency
   obj.ref_freq()

   #Step-3 generate the sky maps at multiple frequencies as well as their sky average
   obj.gen_freq()

   #Step-4 finally, generate a sky averaged spectrum vs frequency figure
   obj.visual()


Save the above code as (say) ``eg_script.py`` and run it as

.. code:: bash

    python eg_script.py

Running the above will generate several files. The terminal messages will guide you to these output files. The most important of all files of your interest will be ``Tb_nu_map.npy``. To visualise your outputs use the function ``visual()``. Refer to the documentation for more details. To learn about the physics of this package see our `paper <https://arxiv.org/abs/2406.17031>`_.

The default values in this package have been chosen such that the users can run their scripts on a PC. Since modern PCs have at least 4 cores, for a better performance, one could also run ``eg_script.py`` parallely as 

.. code:: bash

    mpirun -np 4 python eg_script.py

In general, and for more realistic flux density ranges and high resolution maps, it is recommended to use the code on high performance clusters.

Documentation
-------------

For more details on the working of the package and understanding the output files refer to the documentation. 

License and citation
--------------------

The software is free to use on the MIT open source license. If you use the software then please consider citing `Mittal et al (2024) <https://arxiv.org/abs/2406.17031>`_.

Contact
-------

In case of any confusion or suggestions for improvement please do not hesitate to contact me.

