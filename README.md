# Lagrangian-Hydrodynamics-1D-
* A modular Python (some modules are cythonized for speed) code to evolve Lagrangian gas shells in presence of a halo gravitational potential, radiative cooling, and a prescribed model for thermal heating. 
* In CL5_M14_n_3_2, there is a script setup.py which compiles the cythonized modules as following: python setup.py build_ext --inplace. This creates the .so files from .pyx files and from .so, any other python script can import functions as usual. There are 3 .pyx scripts altogether: (a)subcycle.pyx [the cython functions for cooling, heating, black hole growth], (b) main_step.pyx [the main hydro step], (c)find_min_tme.pyx [finding the minimum dt for hydro step].
* The paths to Lagrangian_1D need to be changed in mypath.py
* Things that can be changed: CL5_M14_n_3_2/parameters.py contains some initial parameters (e.g., halo mass - M_halo, number of shells- N_r). CL5_M14_n_3_2/CL_new.py has cooling=True which can be turned to False to shut off radiative cooling and heating. Note CL5_M14_n_3_2/CL_new.py also has restart=False. The restart number can be added in parameters.py. Data dumping is in CLnew.py and user may modify it in a more convenient format (I dumped data by file number but e.g., one can dump data if 0.95*t_dump <= t <= 1.05*t_dump, etc). 
* The .py scripts in Lagrangian_1D/ are related to modeling the gravitational profile for any halo mass and can be used even if one makes a copy of CL5_M14_n_3_2/ to CLfresh (say) and changes the halo mass in CLfresh/parameters.py in order to run the simulation of a smaller or larger halo. In latter case, mypath.py needs to be edited in the new CLfresh directory. Lagrangian_1D/ also has constants.py for any global constants, new ones can be added. This is imported in almost all the scripts. 
* Finally, of course, run the entire combination by python CLnew.py from within the directory of a specific case (CL5_M14_n_3_2 or CLfresh or whatever).
