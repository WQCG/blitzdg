import setuptools
from distutils.core import Extension


with open("README.md", "r") as fh:
    long_description = fh.read()



pyblitzdg_shared = Extension('pyblitzdg.pyblitzdg',
                        sources=['src/pyblitzdg/pyblitzdg.cpp'],
                        include_dirs=['/usr/local/include', 'include'],
                        library_dirs=['/usr/local/lib/boost', 'lib'],
                        runtime_library_dirs=['/usr/local/lib/boost'],
                        libraries=['boost_python3', 'lblitzdg'])

setuptools.setup(
     name='pyblitzdg',  
     version='0.1.0',
     scripts=[] ,
     author="Waterloo Quantitative Consulting Group",
     author_email="dsteinmo@wqcg.ca",
     description="Discontinuous Galerkin Finite Element Library and Solvers",
     ext_modules=[pyblitzdg_shared],
     long_description=long_description,
   long_description_content_type="text/markdown",
     url="https://github.com/WQCG/blitzdg",
     packages=setuptools.find_packages(),
     classifiers=[
         "Programming Language :: Python :: 3",
         "License :: GPLv3",
         "Operating System :: Linux/OSX",
     ],
 )
