# Symmetries of Riemann Surfaces and Magnetic Monopoles 
This repository contains code giving calculations from my PhD thesis. The files are split according to which of the two main chapters the code is relevant for.
One section of the thesis is based upon [*Bring's Curve: Old and New*](https://arxiv.org/abs/2208.13692) by H. W. Braden and Linden Disney-Hogg, and the code relevant for that 
section is available [here](https://github.com/DisneyHogg/Brings\_Curve). 

## Requirements
I shall not give a complete list of system requirements, but as a rough guideline:
* the bulk of the code was written in [Sage](https://www.sagemath.org/) version 9.4,
* the code training a classifier on theta characteristic orbit data used [Python](https://www.python.org/) and the additional packages [numpy](https://numpy.org/), [pandas](https://pandas.pydata.org/) and [scikit-learn](https://scikit-learn.org/stable/index.html),
* additional python packages were required for the plotting of monopole energy density isosurfaces; the corresponding packages may be deduced by reading the imported functions,
* (optional) and further functionality can be added with a licensed copy of [Maple](https://www.maplesoft.com/) which [can be ran via Sage](https://doc.sagemath.org/html/en/reference/interfaces/sage/interfaces/maple.html),

## Acknowledgements
Not all the code contained in this repository is entirely written by me, there are two execptions.
* I had no part in writing `grupos_que_actuan.sage` and `polyB.sage`; these were programmed by [Antonino Behn and Anita Rojas](https://sites.google.com/a/u.uchile.cl/mat-ciencias-prof-anita-rojas/home/proyectos) for the papers [*Adapted hyperbolic polygons and symplectic representations for group actions on Riemann surfaces*](https://doi.org/10.1016/j.jpaa.2012.06.030) by Antonio Behn, Rubí E. Rodríguez and Anita M. Rojas, and [*A SAGE package for equisymmetric stratification in Mg*](https://doi.org/10.1080/10586458.2020.1763872) by A. Behn, A. M. Rojas, M. Tello-Carrera. I have included them here for ease of access, as I use them to compute data on orbits of theta characteristics on which a classifier is trained.
* The files `monopole_plotting.py`, `V4_nahmdata.py` and `minimal_plotting_from_file.ipynb` were developed from an initial source provided to me by [Paul Sutcliffe](https://www.maths.dur.ac.uk/users/p.m.sutcliffe/index.html).
