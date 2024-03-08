Other functions
---------------

There are 5 additional useful methods of the class ``furs``. These are:-

1. ``acf(chi)``

   -  returns the 2PACF, $C(\chi)$
   -  requires one argument, the angle $\chi$ in radians; can be a
      number or an array
   -  output is a dimensionless quantity

2. ``dndS(S)``

   -  returns flux distribution, $\mathrm{d}n/\mathrm{d}S$. The
      functional form will be according to your choice for ``dndS_form``
      you gave during initialisation. Default is 0.
   -  requires one argument, the flux density $S$ in Jy; can be a number
      or an array
   -  output is in units of number per unit flux density per unit solid
      angle, i.e. $\mathrm{Jy^{-1}sr^{-1}}$

3. ``num_den()``

   -  returns the clustered number density as number per pixel,
      $n_{\mathrm{clus}}$
   -  no arguments required
   -  output is an array of length $N_{\mathrm{pix}}$

4. ``num_sources()``

   -  returns the total number of unresolved extragalactic radio sources
      for the full sky, $N_{\mathrm{s}}$
   -  no arguments required
   -  output is a pure number

5. ``print_input()``

   -  If you want to print the all raw parameter values you gave, you
      may use this function to print them
   -  no arguments required
   -  no return value

Example usage: to find the number of sources between $10^{-6}$ and
$10^{-1}\mathrm{Jy}$ do

::

   obj = furs(logSmin=-6,logSmax=-1)
   Ns = obj.num_sources()

