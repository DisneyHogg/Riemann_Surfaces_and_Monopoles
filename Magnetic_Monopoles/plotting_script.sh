#!/bin/bash

# Example usage of monopole_plotting.py. The program find the Higgs field corresponding to some Nahm data
# on a specified spatial extent, using a discretised grid. 

# -T is the name of a nahmdata.py file that contains a create_nahmdata function. This determines what Nahm data
# is being integrated

# -p is an integer giving the number of parallel processes to run the program with

# -o is a string which gives a prefix for the name of the output files

# -x, -y, and -z determine the spatial extent on which to find the Higgs field, e.g if -x xmax
# is taken the x spatial extent is -xmax to xmax.

# -s is the grid size in the spatial extent, uniform in each direction

# Any remaining parameters are passed to the nahmdata.py function, in this case we pass k, alpha, and sgn. 

python monopole_plotting.py -T V4_nahmdata.py -p 1 -o V4_scattering -x 5 -y 3 -z 5 -s 0.77 -2.0 1.0

# Output is two files, numpy arrays of the same size storing the spatial extent (x_array) and the Higgs
# field value at each corresponding spatial point (higgs_array).
