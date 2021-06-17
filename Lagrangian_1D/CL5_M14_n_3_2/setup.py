from setuptools import Extension, setup
from Cython.Build import cythonize
import numpy

extensions1 = [
       Extension("subcycle", 
                sources=["subcycle.pyx"], 
                libraries=["m"],
                define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")]
                )
]
setup(
    ext_modules=cythonize(extensions1),
    include_dirs=[numpy.get_include()]
)


extensions2 = [
       Extension("main_step",
                sources=["main_step.pyx"],
                libraries=["m"],
                define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")]
                )
]



setup(
    ext_modules=cythonize(extensions2),
    include_dirs=[numpy.get_include()]
)



extensions3 = [
       Extension("find_min_tme",
                sources=["find_min_tme.pyx"],
                libraries=["m"],
                define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")]
                )
]

setup(
    ext_modules=cythonize(extensions3),
    include_dirs=[numpy.get_include()]
)
