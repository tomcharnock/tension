# Calculate the tension between cosmological parameters obtained from two datasets.

## Int P1(theta)P2(theta) dtheta
To calculate the integral of the probability distribution P1 multiplied by P2 you can use the code in the tension folder. This includes a python version and a Fortran90 code.

###Python 
*uses chains.*
```
python tension/python/likelihoods.py
```

###Fortran
*uses multivariate gaussian.*

This example compares `data/CMB/CMB` results with `data/Strong/Strong_L` results and outputs into `output/Strong_L.txt`.
```
cd tension/fortran/
ifort -openmp tension.f90 -o tension
cd ../../
./tension/fortran/tension CMB/CMB Strong/Strong_L Strong_L
```


##Amara surprise
Calculates the surprise and amount of information gained when comparing constraints from independent sources.

###Python
*uses multivariate gaussian.*

This also plots the 2D contours for each of the parameters.
```
python surprise/python/surprise.py
```

###Fortran
*uses multivariate gaussian.*

This example compares `data/CMB/CMB` results with `data/Strong/Strong_L` results and vice versa and outputs into `output/Strong_L.txt`.
```
cd surprise/fortran/
ifort -openmp surprise.f90 -o surprise
cd ../..
./surprise/fortran/surprise CMB/CMB Strong/Strong_L Strong_L
```

##Difference vector
Calculates the difference between CMB and LSS and then calculates the integral of the contour above the pdf value at the origin.

###Python
*uses chains directly.*

This also plots the 2D contours for each of the parameters. Needs chains in `data/.../chains/`.

This example finds the difference vector of the `data/CMB/chains/CMB_*.txt` data to the `data/Strong/chains/Strong_L_*.txt` data using 10 bins in the histogram.
```
python difference_vector/difference_vector.py 10 CMB CMB Strong Strong_L
```





