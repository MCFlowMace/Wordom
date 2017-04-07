#! /usr/bin/env python

# System imports
from distutils.core import *
from distutils      import sysconfig

# wordom extension module
_wordom = Extension("_wordom",
                   ["wordom.i","fileio.c","tools.c","qcprot.c", "xdrfile.c", "xdrfile_xtc.c"],
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

