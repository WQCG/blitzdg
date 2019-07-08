import setuptools
from distutils.core import Extension


with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
     name='pyblitzdg',  
     version='0.1.2',
     scripts=[] ,
     author="Waterloo Quantitative Consulting Group",
     author_email="dsteinmo@wqcg.ca",
     description="Discontinuous Galerkin Finite Element Library and Solvers",
     long_description=long_description,
     long_description_content_type="text/markdown",
     url="https://github.com/WQCG/blitzdg",
     platforms=['linux_x86_64'],
     packages=[
         'pyblitzdg'
     ],
     package_data= {
        'pyblitzdg': ['pyblitzdg.so', 'libblitzdg.so', 'libblitz.so.0']
    },
     classifiers=[
         "Programming Language :: Python :: 3",
         "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
         "Operating System :: POSIX :: Linux",
     ],
 )
