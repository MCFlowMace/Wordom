#! /usr/bin/env python

# System imports
from distutils.core import *
from distutils      import sysconfig

# Third-party modules - we depend on numpy
import numpy
# in order to check whether lapack are present ...
import numpy.distutils.system_info as sysinfo

# Obtain the numpy include directory. This works across numpy versions.
try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

# wordom extension module
if len(sysinfo.get_info('lapack')) == 0:
  _wordom = Extension("_wordom",
                     ["wordom.i","fileio.c","tools.c","qcprot.c", "xdrfile.c", "xdrfile_xtc.c"],
                     )
else:
  _wordom = Extension("_wordom",
                     ["wordom.i","fileio.c","tools.c","qcprot.c", "xdrfile.c", "xdrfile_xtc.c"],
                     include_dirs = [numpy_include],
                     extra_compile_args = ["-D LAPACK"],
                     libraries = [ 'lapack', 'blas' ]
                     )

# NumyTypemapTests setup
setup(  name        = "wordom",
        description = "wordom is a molecular structure and data manipulation program/library",
        author      = "Michele Seeber & colleagues",
        url         = "http://wordom.sf.net",
        author_email= "mseeber@gmail.com",
        license     = "GPL",
        version     = "0.23",
        ext_modules = [_wordom],
        py_modules  = ['wordom']
        )

