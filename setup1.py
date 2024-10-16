from distutils.core import setup 
from Cython.Build import cythonize
import numpy

setup(ext_modules = cythonize('Cable_eq_current_cyt.pyx'),include_dirs=[numpy.get_include()])