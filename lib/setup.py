import setuptools
with open("../README.md", "r") as fh:
    long_description = fh.read()


setuptools.setup(
     name='pyblitzdg',  
     version='0.1.0',
     scripts=[] ,
     author="Waterloo Quantitative Consulting Group",
     author_email="dsteinmo@wqcg.ca",
     description="Discontinuous Galerkin Finite Element Library and Solvers",
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
