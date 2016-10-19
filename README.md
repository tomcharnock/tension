# Calculate the tension between cosmological parameters obtained from two datasets.

Parameters and default parameter values in the `params.py` file in `params/` are
```
parameters['method'] = ['combine'] 
parameters['chain_dir'] = 'chains/'
parameters['CMB'] = 'CMB'
parameters['CMB_chains'] = 6
parameters['LSS'] = 'Strong_L'
parameters['LSS_chains'] = 6
parameters['params'] = ['omegabh2', 'omegach2', 'theta', 'logA', 'ns']
parameters['bins'] = 40
parameters['plot_dir'] = None
parameters['save_dir'] = None
parameters['load'] = None
parameters['priors'] = None
parameters['smoothing'] = 0
parameters['interpolate_bins'] = None
parameters['integration_bounds'] = None
parameters['plot_interpolate'] = None
parameters['sampling_method'] = None
parameters['sampling_parameter'] = None
parameters['sampling_constraints'] = None
```
Possible methods are
```
parameters['method'] = ['combine', 'surprise', 'bhattacharyya', 'overlap_coefficient', 'ibi', 'marshall', 'verde', 'difference_vector']
```
Note that `'difference_vector'` cannot be run in the list with others.
Plot are not made unless specifed. A `plots/` directory is present and can be used by setting
```
parameters['plot_dir'] = 'plots/'
```
The arrays are also not saved unless specifed. To save all arrays a `saves/` directory is present and can be used by setting
```
parameters['save_dir'] = 'saves/`
```
If loading the arrays when using `'difference_vector'` then
```
parameters['load'] = ['directory_where_difference_vector_is_saved/name_of_difference_vector']
```
If loading the arrays when using anything other than `'difference_vector'` or `'verde'` then
```
parameters['load'] = ['directory_where_CMB_array_is_saved/name_of_CMB_array', 'directory_where_LSS_array_is_saved/name_of_LSS_array']
```
When loading the arrays using `'verde'` then
```
parameters['load'] = ['directory_where_CMB_array_is_saved/name_of_CMB_array', 'directory_where_LSS_array_is_saved/name_of_LSS_array', 'directory_where_shifted_array_is_saved/name_of_shifted_array']
```
If using `'marshall'` or `'verde'` the priors will be uniform unless the correct dimension array containing the prior is loaded by setting
```
parameters['priors'] = 'directory_where_priors_are_saved/name_of_priors'
```
For Gaussian smoothing or arrays (advised) then a smoothing of 1 should be used
```
parameters['smoothing'] = 1
```
To set optimum bins for the interpolated grid then the `'interpolate_bins'` parameter should be unset. This gives better results, but can use up HUGE amounts of memory - might be better to just go with smoothing.
```
parameters['interpolate_bins'] = None
```
When using `'ibi'` then the integration bounds need to be set, for 3 sigma bounds use
```
parameters['integration_bounds'] = 0.997
```
To use interpolation when plotting (for whatever reason) the number of bins can be set (to 10 for this example) using
```
parameters['plot_interpolate'] = 10
```
To impliment importance sampling on the chains before processing the sampling method, parameter and constraints need to be set. For Planck2013 Gaussian priors on theta and ns this could be
```
parameters['sampling_method'] = ['gaussian', 'gaussian']
parameters['sampling_parameter'] = ['theta', 'ns']
parameters['sampling_constraints'] = [[1.04086, 0.00029], [0.9666, 0.0040]]
```
For Gaussian priors the `'sampling_constraints'` parameter must be a list with the mean in the first position and the standard deviation in the second position. Flat priors can also be applied using
```
parameters['sampling_method'] = ['uniform']
parameters['sampling_parameter'] = ['omegabh2']
parameters['sampling_constraints'] = [[0.0224, 0.0225]]
```
where the first index of the list is the lower bound and the second position in the list is the upper bound. Both Gaussian and uniform priors can be applied to different parameters at the same time.
##To run
To run the code (with parameters saved in the premade parameter file, `params/params.py` then use
```
python measures.py
```
To define new parameter file, say `params/new_params.py` the file must contain a dictionary called `parameters`. An example could be
```
#params/new_params.py
parameters = {}

parameters['method'] = ['difference_vector'] 
parameters['chain_dir'] = 'chains/'
parameters['CMB'] = 'CMB'
parameters['CMB_chains'] = 2
parameters['LSS'] = 'Weak_L'
parameters['LSS_chains'] = 2
parameters['params'] = ['omegabh2', 'omegach2', 'theta', 'logA', 'ns']
parameters['bins'] = 20
parameters['plot_dir'] = 'plots/'
parameters['save_dir'] = 'saves/'
parameters['smoothing'] = 1
```
This is then run with
```
python measures.py -file new_params
```
##utils.py
This file contains all the functions to get probability distributions and mess with arrays etc. It also contains the default parameters and the class which contains them. For some reason I decided it should be in alphabetical order...

##measures.py
This file contains all the calculation of the methods and whether to plot. It basically calls a load of stuff from `utils.py`.

##plot.py
This file contains the plotting functions and one function for working out the isocountours which is also used in `'ibi'`.

##params/params.py
The default input file.

##plots/
Directory for plots

##saves/
Directory for saving arrays
