from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension
setup(name="metrics", ext_modules=cythonize('metrics.pyx'),)

#ext_modules=[Extension("demo",sources=["demo.pyx"],libraries=["m"] # Unix-like specific)]

#setup( name = "Demos",ext_modules = cythonize(ext_modules))

