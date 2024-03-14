Frequently Asked Questions
--------------------------

1. When I run ``visual()``, I am getting the following error 

``RuntimeError: Failed to process string with tex because latex could not be found``

How to rectify this?


Solution: This is probably because you do not have a proper LaTeX installation. Run the following commands on your terminal

.. code:: bash

   sudo apt install texlive texlive-latex-extra texlive-fonts-recommended dvipng
   pip install latex

Solution borrowed from `Stack Overflow -failed-to-process-string <https://stackoverflow.com/questions/58121461/runtimeerror-failed-to-process-string-with-tex-because-latex-could-not-be-found>`_.


2. I am getting the following error

``ERROR: Could not build wheels for mpi4py, which is required to install pyproject.toml-based projects``

when I install :mod:`mpi4py`.

Solution: Run the following commands on your terminal 

.. code:: bash

   sudo apt update
   sudo apt-get install libopenmpi-dev 
   
Solution borrowed from `Stack Overflow -could-not-build-wheels <https://stackoverflow.com/questions/74427664/error-could-not-build-wheels-for-mpi4py-which-is-required-to-install-pyproject>`_.   


3. Will this package run on windows?

No, because it uses :mod:`healpy` and since :mod:`healpy` is not (yet) supported on windows, this package cannot be used on windows. However, there is still a workaround without having to dual boot your PC with ubuntu. You have use the Windows Subsystem for Linux (WSL). See the official `ubuntu <https://ubuntu.com/desktop/wsl>`_ page and `this <https://learn.microsoft.com/en-us/windows/wsl/install>`_ page.


