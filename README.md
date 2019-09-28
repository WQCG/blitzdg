# blitzdg

blitzdg is an open-source library offering discontinuous Galerkin (dg) solvers for common partial differential equations systems using blitz++ for array and tensor manipulations in a C++ environment or NumPy as a Python 3 library.

<img alt="shallow water wave example" src="https://raw.githubusercontent.com/WQCG/blitzdg/master/example/sw_coarsebox.gif" />

Shallow Water Wave Example ```blitzdg``` output

[![Build Status](https://travis-ci.org/dsteinmo/blitzdg.svg?branch=master)](https://travis-ci.org/WQCG/blitzdg) [![Coverage Status](https://coveralls.io/repos/github/WQCG/blitzdg/badge.svg?branch=master)](https://coveralls.io/github/WQCG/blitzdg?branch=master)

Support blitzdg: <a href="https://www.paypal.com/cgi-bin/webscr?cmd=_s-xclick&hosted_button_id=3RM7VGV28NEHU&source=url" ><img alt="Donate to blitzdg development" src="https://dantheman827.github.io/images/donate-button.svg" width="75" /></a>
<img alt="" border="0" src="https://www.paypal.com/en_CA/i/scr/pixel.gif" width="1" height="1" />


## Building From Source

Build and development support has broadened from linux only to Mac OSX and Windows systems. Tested with GNU make (written to be cross-platform) and `g++` on linux/MinGW64, `clang++` on Mac OSX Sierra, and MSVC on Windows. Our build system depends on the cross-platform `cmake` tooling for Makefile generation.

1. `git clone https://github.com/dsteinmo/blitzdg.git`
2. `cd blitzdg && ./pull-deps.sh`
3. `cmake . && make advec1d && ./bin/advec1d` (This binary is a 1D advection equation solver.).
4. Run unit tests with `make test`.

### Running with Docker

You can also run the build and tests inside a docker (linux) container. The container is based on an ubuntu 18.04 image.

1. `git clone https://github.com/dsteinmo/blitzdg.git && cd blitzdg`
2. `docker build -t blitzdg .`
3. `docker run -t blitzdg`

## Dependencies

* `cmake`
* `blitz++`
* `SuiteSparse (umfpack, cxsparse)`
* `LAPACK`
* `metis`
* `boost`
* `igloo` for BDD-style testing.
* `vtk` for visualization in Paraview.
* `boost-python` for python bindings.

Dependency installation is outlined in `pull-deps.sh` (tested on Ubuntu and Mac OSX).

### Windows

Our windows distribution was recently switched from MinGW to MSVC in order to achieve better support for graphics APIs, so support is lacking at the moment. Instructions will be made available here soon.

## Contributing

We accept pull requests from public forks, and we use pull requests as the primary delivery mechanism of any new code within the base repository.

If you add code, please write tests using the igloo testing framework that is included as a project dependency. Your code additions will be subject to peer review and will be run through our Travis-CI continuous integration process.

Interested developers should consult the [Guidelines for Contributing](https://github.com/WQCG/blitzdg/blob/master/CONTRIBUTING.md "Contributing Markdown") before getting started.

## Maintainer

* [Derek Steinmoeller](https://github.com/dsteinmo)

## Lead Developers

* [Derek Steinmoeller](https://github.com/dsteinmo)
* [Killian Miller](https://github.com/k7miller)

## Documentation

We actively maintain an interactive set of docs using Doxygen for end-user consumption.

The documentation is available on github pages at [https://wqcg.github.io/blitzdg](https://wqcg.github.io/blitzdg "blitzdg Documentation") and is kept synchronized with the master branch via automation.

## License

This project is licensed under the [GNU Public License Version 3](https://www.gnu.org/licenses/gpl-3.0.en.html "GPLv3 License").

Our license choice is driven by the desire to keep this project and any of its derivative works open-source for public consumption by developers, mathematicians, scientists, engineers, and anyone else who might be interested in this project.

## Contact

Any questions regarding the project may be addressed via email to the [project maintainer](mailto:dsteinmo@wqcg.ca).
