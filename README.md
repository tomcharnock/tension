# Calculate the tension between cosmological parameters obtained from two datasets.

Parameters and default parameter values in the `params.py` file in `params/` are
```
parameters['method'] = None 
parameters['chain_dir'] = 'chains/'
parameters['CMB'] = 'CMB'
parameters['CMB_chains'] = 6
parameters['LSS'] = 'Strong_L'
parameters['LSS_chains'] = 6
parameters['params'] = ['omegabh2', 'omegach2', 'theta', 'logA', 'ns']
parameters['bins'] = 40
parameters['plot_dir'] = None
parameters['save_dir'] = None
parameters['load_file'] = None
parameters['load_priors'] = None
parameters['filename'] = None
parameters['smoothing'] = 0.5
parameters['interpolate_bins'] = None
parameters['integration_bounds'] = None
parameters['sampling_method'] = None
parameters['sampling_parameter'] = None
parameters['sampling_constraints'] = None
```
Possible methods are
```
parameters['method'] = ['surprise', 'bhattacharyya', 'overlap_coefficient', 'ibi', 'marshall', 'verde', 'difference_vector']
```
Note that `'difference_vector'` cannot be run in the list with others.
Plot are not made unless specifed. A `plots/` directory is present and can be used by setting
```
parameters['plot_dir'] = 'plots'
```
To give a filename to any of the saves or plots this can be set using
```
parameters['filename'] = 'test'
```
The arrays are also not saved unless specifed. To save all arrays a `saves/` directory is present and can be used by setting
```
parameters['save_dir'] = 'saves`
```
If loading the arrays when using `'difference_vector'` then if saved as `saves/test_Strong_L_difference_vector_hist.npy` using `parameters['filename'] = 'test'` and `parameters['LSS'] = 'Strong_L'` then use
```
parameters['load_file'] = ['test_Strong_L']
```
If loading the arrays when using anything other than `'difference_vector'` then use
```
parameters['load_file'] = ['name_of_CMB_array', 'name_of_LSS_array']
```
Note that `['save_dir']` must be set when loading parameters.
If using `'marshall'` or `'verde'` the priors will be uniform unless the correct dimension array containing the prior is loaded by setting
```
parameters['load_priors'] = 'name_of_priors'
```
For Gaussian smoothing or arrays (advised) then a smoothing can be applied using
```
parameters['smoothing'] = 1
```
To set optimum bins for the interpolated grid then the `'interpolate_bins'` parameter should be unset. This gives better results, but can use up HUGE amounts of memory - might be better to just go with smoothing.
```
parameters['interpolate_bins'] = 50
```
When using `'ibi'` then the integration bounds need to be set, for 3 sigma bounds use
```
parameters['integration_bounds'] = 0.997
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
where the first index of the list is the lower bound and the second position in the list is the upper bound. Both Gaussian and uniform priors can be applied to different parameters at the same time. Priors can also be placed on parameters not included in the analysis, but which are contained in the chains. For example, placing the Planck2016+lowE optical depth constraints to the CMB chain only can be acheived by using
```
parameters['sampling_method'] = ['gaussian']
parameters['sampling_parameter'] = ['tau']
parameters['sampling_constraints'] = [[0.058, 0.012]]
```
`tau` is not in the LSS chains and so no importance sampling of `tau` occurs. 
##To run
To run the code (with parameters saved in the premade parameter file, `params/params.py` then use
```
python tension.py
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
python tensions.py -file new_params
```
This can also be imported as a module, for example
```
import tension
parameters = {}
parameters['params'] = ['omegabh2', 'omegach2']
parameters['bins'] = 50
parameters['interpolate_bins'] = 50
parameters['filename'] = 'surprise_only'
parameters['plot_dir'] = 'plots'
parameters['plot_show'] = True
parameters['save_dir'] = 'saves'
parameters['smoothing'] = 1
parameters['sampling_method'] = ['gaussian']
parameters['sampling_parameter'] = ['tau']
parameters['sampling_constraints'] = [[0.058, 0.012]]

T = tension.tension(parameters)
T.method = ['surprise']
T.combine()
T.D1_int, T.S1 = T.get_surprise(T.CMB_hist, T.LSS_hist)
T.D2_int, T.S2 = T.get_surprise(T.LSS_hist, T.CMB_hist)
filename = T.get_filename(T.LSS, string = 'surprise_S1_' + str(T.S1) + '_S2_' + str(T.S2))
T.plot_surprise(filename)
```

##tension.py
Module file.

##params/params.py
The default input file.

##plots/
Directory for plots.

##saves/
Directory for saving arrays.
