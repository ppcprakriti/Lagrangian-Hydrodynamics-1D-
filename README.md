# Lagrangian-Hydrodynamics-1D-
A modular Python (some modules are cythonized for speed) code to evolve Lagrangian gas shells in presence of a halo gravitational potential, radiative cooling, and a prescribed model for thermal heating. 
In CL5_M14_n_3_2, there is a script setup.py which compiles the cythonized modules as following: python setup.py build_ext --inplace. This creates the .so files from .pyx files and from .so, any other python script can import functions as usual. There are 3 .pyx scripts altogether: (a)subcycle.pyx, (b) main_step.pyx, (c)find_min_tme.pyx.
The paths to Lagrangian_1D need to be changed in mypath.py
Things to change: CL5_M14_n_3_2/parameters.py contains some initial parameters (e.g., halo mass). CL5_M14_n_3_2/CL_new.py has cooling=True which can be turned to False to shut off radiative cooling and heating. 
The .py scripts in Lagrangian_1D/ are related to modeling the gravitational profile for any halo mass and can be used even if one makes a copy of CL5_M14_n_3_2/ to CLfresh (say) and changes the halo mass in CLfresh/parameters.py in order to run the simulation of a smaller or larger halo. 
