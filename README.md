# blitzdg

blitzdg is an open-source project aiming to implement parallel discontinuous Galerkin (dg) solvers for common partial differential equations systems using blitz++ for array and tensor manipulations and MPI for distributed parallelism.

[![Build Status](https://travis-ci.org/dsteinmo/blitzdg.svg?branch=master)](https://travis-ci.org/WQCG/blitzdg)  [![Coverage Status](https://coveralls.io/repos/github/WQCG/blitzdg/badge.svg?branch=master)](https://coveralls.io/github/WQCG/blitzdg?branch=master)

## Getting Started

Build and development support has broadened from linux only to Mac OSX and Windows systems. Tested with GNU make (written to be cross-platform) and `g++` on linux/MinGW64, `c++` on Mac OSX Sierra.

1. `git clone https://github.com/dsteinmo/blitzdg.git`
2. `cd blitzdg && ./pull-deps.sh`
3. `make && ./bin/advec1d` (This binary is a 1D advection equation solver.).
4. Run unit tests with `make test`.

### Running with Docker

You can also run the build and tests inside a docker (linux) container. The container is based on an ubuntu 16.04 image.

1. `git clone https://github.com/dsteinmo/blitzdg.git && cd blitzdg`
2. `docker build -t blitzdg .`
3. `docker run -t blitzdg`

## Dependencies

So far:

* `blitz++`
* `SuiteSparse (umfpack, cxsparse)`
* `LAPACK`
* `metis`
* `boost`
* `igloo` for BDD-style testing.

Dependency installation is outlined in `pull-deps.sh` (tested on Ubuntu and Mac OSX).

### Windows Dependencies

The Windows build requires [MinGW/MinGW64](http://www.mingw.org/wiki/Getting_Started "MinGW Installation Instructions"). We are testing the windows build in AppVeyor using the `mingw64-x86_64-7.3.0-posix-seh-rt_v5-rev0` distribution of MinGW. There is currently not a Visual Studio/MSVC build.

Run `.\pull-deps.ps1` in Powershell (4+ or Powershell Core (`pwsh`)).

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