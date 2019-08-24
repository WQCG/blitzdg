import setuptools
from distutils.core import Extension


with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
     name='pyblitzdg',  
     version='0.1.5',
     scripts=[] ,
     author="Waterloo Quantitative Consulting Group",
     author_email="dsteinmo@wqcg.ca",
     description="Discontinuous Galerkin Finite Element Library and Solvers",
     long_description=long_description,
     long_description_content_type="text/markdown",
     url="https://github.com/WQCG/blitzdg",
     platforms=['linux_x86_64', 'win64'],
     ext_modules=[Extension('pyblitzdg', ['src/pyblitzdg/pyblitzdg.cpp'], 
        include_dirs=[ 'C:/ProgramData/Anaconda3/pkgs/py-boost-1.67.0-py37h8300f20_4/Library/include', 'C:/ProgramData/Anaconda3/pkgs/libboost-1.67.0-hfd51bdf_4/Library/include/', 'include'],
        library_dirs=["lib", "C:/dev/blitzdg/pyblitzdg/Win64"],
        libraries=['blitzdg', 'vtkIOXML-7.1', 'vtkCommonCore-7.1', 'vtkCommonExecutionModel-7.1', 'vtkCommonDataModel-7.1', 'vtkCommonMisc-7.1', 'vtkCommonSystem-7.1', 'vtkCommonTransforms-7.1', 'vtkexpat-7.1', 'vtkIOCore-7.1', 'vtkIOGeometry-7.1', 'vtkIOXML-7.1', 'vtkIOXMLParser-7.1', 'vtksys-7.1', 'vtkzlib-7.1'],
        define_macros=[("VTKCOMMONCORE_STATIC_DEFINE", None),
            ("VTKCOMMONEXECUTIONMODEL_STATIC_DEFINE", None),
            ("VTKIOGEOMETRY_STATIC_DEFINE", None),
            ("VTKCOMMONDATAMODEL_STATIC_DEFINE", None),
            ("VTKIOXML_STATIC_DEFINE", None),
            ("VTKRENDERINGCORE_STATIC_DEFINE", None),
            ("PY_MAJOR_VERSION", "3"),
            ("PY_MINOR_VERSION", "7")
            ]
        )],
     classifiers=[
         "Programming Language :: Python :: 3",
         "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
         "Operating System :: POSIX :: Linux",
     ],
 )
