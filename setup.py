import setuptools
from distutils.core import Extension


with open("README.md", "r") as fh:
    long_description = fh.read()

pyblitzdg = Extension('pyblitzdg',
                        sources=['src/pyblitzdg/pyblitzdg.cpp'],
                        include_dirs=['/usr/local/include', 'include'],
                        library_dirs=['/usr/local/lib/boost', '/usr/local/lib', 'lib'],
                        runtime_library_dirs=['/usr/local/lib/boost'],
                        libraries=['boost_python3', 'blitzdg', 'blitz'])

setuptools.setup(
     name='pyblitzdg',  
     version='0.1.3',
     scripts=[] ,
     author="Waterloo Quantitative Consulting Group",
     author_email="dsteinmo@wqcg.ca",
     description="Discontinuous Galerkin Finite Element Library and Solvers",
     long_description=long_description,
     long_description_content_type="text/markdown",
     url="https://github.com/WQCG/blitzdg",
     platforms=['manylinux2010'],
     ext_modules=[pyblitzdg],
     packages = [],
     package_data = {},
     classifiers=[
         "Programming Language :: Python :: 3",
         "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
         "Operating System :: POSIX :: Linux",
     ],
 )
