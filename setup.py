from setuptools import setup
from distutils.core import setup
from distutils.extension import Extension

def readme():
        with open('README.rst') as f:
                    return f.read()
# try:
#     from Cython.Distutils import build_ext
#     from Cython.Build import cythonize
# except ImportError:
#     use_cython = False
# else:
#     use_cython = True
#
# cmdclass = {}
# # ext_modules = [cythonize('cython/metrics.pyx') ]
# ext_modules = []
#
# if use_cython:
#     ext_modules += [Extension("metrics", ["metrics.pyx"]), ]
#     cmdclass.update({'build_ext': build_ext})
# else:
#     ext_modules += [Extension("metrics", ["metrics.c"]), ]
setup(
    name='correlcalc',
    version='0.97rc1',
    description='Two-point correlation function (2pCF) calculation',
    long_description=readme(),
    url='http://github.com/rohinkumar/correlcalc',
    author='Rohin Kumar Y',
    author_email='yrohinkumar@gmail.com',
    license='MIT',
    packages=['correlcalc'],
    install_requires=['numpy', 'scipy', 'astropy', 'cython', 'tqdm', 'matplotlib', 'pymangle', 'sklearn'],
    dependency_links=['https://github.com/esheldon/healpix_util/master/tarball#egg=package-0.1'],
    # cmdclass = cmdclass,
    # ext_modules=ext_modules,
    zip_safe=False)
