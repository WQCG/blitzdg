import setuptools
from distutils.core import Extension


with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
     name='pyblitzdg',  
     version='0.3.0',
     scripts=[] ,
     author="Waterloo Quantitative Consulting Group",
     author_email="dsteinmo@wqcg.ca",
     description="Discontinuous Galerkin Finite Element Library and Solvers",
     long_description=long_description,
     long_description_content_type="text/markdown",
     url="https://github.com/WQCG/blitzdg",
     platforms=['manylinux2010', 'win64', 'macosx'],
     ext_modules=[Extension('pyblitzdg', ['src/pyblitzdg/pyblitzdg.cpp'], 
        include_dirs=[ "include", "include/igloo", "/usr/include/python3.7m/", "/usr/local/include/", "/usr/local/Cellar/python/3.7.4_1/Frameworks/Python.framework/Versions/3.7/include/python3.7m/",
            "C:\\ProgramData\\Anaconda3\\pkgs\\python-3.7.3-h8c8aaf0_0\\include\\",
            "C:\\ProgramData\\Anaconda3\\pkgs\\py-boost-1.67.0-py37h8300f20_4\\Library\\include\\",
            "C:/ProgramData/Anaconda3/Library/include\\",
            "C:/dev/blitzdg/include/"
        ],
        library_dirs=[ "/root/project/", "/root/project/pyblitzdg", "/usr/lib64", "/usr/local/lib64", "/usr/lib/x86_64-linux-gnu/", "/usr/local/lib", "/usr/local/Cellar/python/3.7.4_1/Frameworks/Python.framework/Versions/3.7/lib/",
            "C:/ProgramData/Anaconda3/pkgs/py-boost-1.67.0-py37h8300f20_4/Library/lib/",
            "C:\\dev\\blitzdg\\Win64\\",
            "C:\\dev\\blitzdg\\lib\\",
            "/miniconda3/pkgs/python-3.7.3-h0371630_0/lib/",
            "/usr/local/Cellar/python/3.7.4_1/Frameworks/Python.framework/Versions/3.7/lib/",
            "/usr/local/Cellar/gcc/8.2.0/lib/gcc/8/"
        ],
        libraries=['blitzdg', 'vtkIOXML-7.1', 'vtkCommonCore-7.1', 'vtkCommonExecutionModel-7.1', 'vtkCommonDataModel-7.1', 'vtkCommonMisc-7.1', 'vtkCommonSystem-7.1', 'vtkCommonTransforms-7.1', 'vtkexpat-7.1', 'vtkIOCore-7.1', 'vtkIOGeometry-7.1', 'vtkIOXML-7.1', 'vtkIOXMLParser-7.1', 'vtksys-7.1', 'vtkzlib-7.1', 'cholmod', 'cxsparse'],
        # libraries=['blitzdg', 'vtkIOXML-7.1', 'vtkCommonCore-7.1', 'vtkCommonExecutionModel-7.1', 'vtkCommonDataModel-7.1', 'vtkCommonMisc-7.1', 'vtkCommonSystem-7.1', 'vtkCommonTransforms-7.1', 'vtkexpat-7.1', 'vtkIOCore-7.1', 'vtkIOGeometry-7.1', 'vtkIOXML-7.1', 'vtkIOXMLParser-7.1', 'vtksys-7.1', 'vtkzlib-7.1', 'gfortran.5'],
        define_macros=[("VTKCOMMONCORE_STATIC_DEFINE", None),
            ("VTKCOMMONEXECUTIONMODEL_STATIC_DEFINE", None),
            ("VTKIOGEOMETRY_STATIC_DEFINE", None),
            ("VTKCOMMONDATAMODEL_STATIC_DEFINE", None),
            ("VTKIOXML_STATIC_DEFINE", None),
            ("VTKRENDERINGCORE_STATIC_DEFINE", None),
            ("PY_MAJOR_VERSION", "3"),
            ("PY_MINOR_VERSION", "7")
            ],
        language = 'c++14',
        extra_compile_args = ["-std=c++14"],
	    extra_link_args = ["-Wl,-rpath,/usr/local/lib"]
        )],
     classifiers=[
         "Programming Language :: Python :: 3",
         "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
         "Operating System :: POSIX :: Linux",
     ],
 )
