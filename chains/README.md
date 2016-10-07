#Folder for the chains
Should contain chains (from `COSMOMC`) named along the lines of `name_of_chain_i.txt` where `name_of_chain = CMB` for example and `_i = _1` or `= _2` for the number of chains.
These are passed to the code as
```
parameters['CMB'] = 'name_of_CMB_chain'
parameters['CMB_chains'] = number_of_CMB_chains
parameters['LSS'] = 'name_of_LSS_chain'
parameters['LSS_chains'] = number_of_LSS_chains
```
Parameter name files should also be contained here and called `name_of_chain.paramnames`.

