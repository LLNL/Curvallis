Some experimental fitters have been added to Curvallis, but with limited success in achieving a reasonable fit.  These 
are outlined here for those who might want to try them.

# GammaPoly and GammaPolyV
Reference [Curvallis enhancement request 5](https://github.com/LLNL/Curvallis/issues/5)

## GammaPoly
The GammaPoly fitter tries to fit the Gr&#252;neisen parameter (gamma) as a function of density (rho) using separate 
polynomials (both of the same degree) for high and low pressure regimes.  To determine whether any particular density is 
high or low pressure, it uses a reference density rho0, which must be set with parameter `rho0_guess`.  Low 
pressure values are where rho <= rho0 and high pressure values are where rho >= rho0.  Note that if rho0 is included in 
the data set, the point will be used by both fitters.  The plotted value for rho0 will be that of the low pressure fit.

Like the polynomial fitter (`fit_type poly[1-12]`), Curvallis recognizes any of`fit_type gammapoly[1-12]` as 
a valid parameter value, and will take the integer at the end of the `fit_type` as the appropriate degree of polynomial 
to use (e.g., `fit_type gammapoly5` will result in a GammaPoly fit of degree 5). The GammaPoly fitter takes all data 
in the low pressure regime and fits a polynomial of the specified degree using x = rho0/rho as the independent variable.  
It then separately takes all data in the high pressure regime and fits a polynomial of the same degree using x = rho/rho0 
as the independent variable.  Note that although the console output for this fitter shows polynomials in terms of x as 
defined above, the plot displayed by Curvallis has density on the x-axis, not x.

In an (unsuccessful) attempt to reduce the discontinuity between the fit for high and low pressure regimes, the value
of Curvallis's `overlap` parameter can be used to include some points from the other regime when creating the fit.

Since GammaPoly already uses the reciprocal values in a specific way, the `v_axis` parameter should not be used with 
this fitter, although the program does not explicitly prevent it.  Rather, one should use the GammaPolyV fitter 
described below.

## GammaPolyV
This fitter behaves very similarly to the GammaPoly fitter described above, except it fits gamma as a function of unit 
volume (v), where v = 1/rho.  The value of the `rho0_guess` parameter for this fitter should be in terms of unit volume 
rather than in density; this value will henceforth be referred to as v0.  Thus high pressure values are 
those where v <= v0, and low pressure values are those where v >= v0. The polynomial fit for high pressure is in terms 
of x = v0/v, while x = v/v0 is used as the independent variable for the low pressure polynomial fit.  The plotted value 
will still be that of the low pressure fit, and the the plot will have volume on the x-axis rather than x as defined in
this paragraph.

Similar to GammaPoly, set the `fit_type` parameter to `gammapolyv[1-12]`, and both polynomials will use the degree
indicated by the `fit_type` name  (e.g. `fit_type gammapolyv3` will result in a GammaPolyV fit of degree 3).  

Example configuration and data files for this fitter are provided in this directory (v_gamma.ini, v_gamma.dat).  The 
`rho0_guess` parameter in this config file is set to a value near the unit volume of the peak gamma value in the data 
file.  Note that because unit volumes in this data set extend over 13 degrees of magnitude, it is recommended to use 
Curvallis's interactive commands to view the resulting plot with a logarithmic scale on the x-axis.

# B-spline fitting
The Jupyter notebook "B-spline Fitting with Sparse Data" gives an overview of B-splines and describes the basic 
functionality exposed in the SciPy library for Python.  The notebook can either be viewed in GitHub or run directly.  

To run the notebook, first install the following:

* [Python 3](https://www.python.org/downloads/)
* Python libraries
  * [numpy](https://numpy.org/install/)
  * [scipy](https://scipy.org/install.html)
  * [matplotlib](https://matplotlib.org/users/installing.html)
* [Jupyter](https://jupyter.org/install)

Copy the notebook .ipynb file to your computer, then navigate to the directory where the file is stored.
Open a command window and run command `jupyter notebook`, which will start a local Jupyter server and will open the notebook 
directory in a browser.  Click on the notebook to open it.  Click the Run button to run cells one at a time 
(recommended), or the fast forward button to run all cells at once. 

To close the notebook, use options on the File menu.  Then go to the command window and stop the server process using
Ctrl+C, or use the command `jupyter notebook stop` in a separate command window.