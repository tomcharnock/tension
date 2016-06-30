Calculate the tension between cosmological parameters obtained from two datasets.

Int P1P2
To calculate the integral of the probability distribution P1 multiplied by P2 you can use the code in the tension folder. This includes a python version and a Fortran90 code.

Python (uses chains)
python tension/python/likelihoods.py

Fortran (uses multivariate gaussian)
cd tension/fortran/
ifort -openmp tension.f90 -o tension
cd ../../
./tension/fortran/tension CMB/CMB Strong/Strong_L Strong_L
* This example compares data/CMB/CMB results with data/Strong/Strong_L results and outputs into output/Strong_L.txt


Amara surprise
Calculates the surprise and amount of information gained when comparing constraints from independent sources.

Python (uses multivariate gaussian)
python surprise/python/surprise.py
* This also plots the 2D contours for each of the parameters

Fortran (uses multivariate gaussian)
cd surprise/fortran/
ifort -openmp surprise.f90 -o surprise
cd ../..
./surprise/fortran/surprise CMB/CMB Strong/Strong_L Strong_L
* This example compares data/CMB/CMB results with data/Strong/Strong_L results and vice versa and outputs into output/Strong_L.txt







