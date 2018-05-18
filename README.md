# blitzdg

blitzdg is an open-source project aiming to implement parallel discontinuous Galerkin (dg) solvers for common partial differential equations systems using blitz++ for array and tensor manipulations and MPI for distributed parallelism.

[![Build Status](https://travis-ci.org/dsteinmo/blitzdg.svg?branch=master)](https://travis-ci.org/WQCG/blitzdg)  [![Coverage Status](https://coveralls.io/repos/github/WQCG/blitzdg/badge.svg?branch=master)](https://coveralls.io/github/WQCG/blitzdg?branch=master) [![Windows Build Status](https://ci.appveyor.com/api/projects/status/pmx725yhsrnq3thu?svg=true)](https://ci.appveyor.com/project/WQCG/blitzdg)

## Getting Started

Build and development support has broadened from linux only to Mac OSX and Windows systems. Tested with GNU make (written to be cross-platform) and g++ on linux/MinGW64/Mac OSX Sierra.

1. `git clone https://github.com/dsteinmo/blitzdg.git`
2. `cd blitzdg && ./pull-deps.sh`
3. `make && ./bin/blitzdg` (The binary currently doesn't do much).
4. Run unit tests with `make test`.

### Running with Docker

You can also run the build and tests inside a docker (linux) container.

1. `git clone https://github.com/dsteinmo/blitzdg.git && cd blitzdg`
2. `docker build -t blitzdg .`
3. `docker run -t blitzdg`

## Dependencies

So far:

* `blitz++`
* `SuiteSparse (umfpack)`
* `LAPACK`
* `metis`
* `boost`
* `igloo` for BDD-style testing.

Dependency installation is outlined in `pull-deps.sh` (tested on Ubuntu) and `pull-deps-mac.sh`.

### Windows Dependencies

The Windows build requires [MinGW/MinGW64](http://www.mingw.org/wiki/Getting_Started "MinGW Installation Instructions") and is currently tested in the AppVeyor CI process. There is currently not a Visual Studio/MSVC build.

Run `.\pull-deps.ps1` in Powershell (4+ or Powershell Core (`pwsh`)).

## Contributing

We accept pull requests.

If you add code, please write tests using the igloo testing framework that is included as a project dependency.

Interested developers should consult the [guidelines for contributing](https://github.com/WQCG/blitzdg/blob/master/CONTRIBUTING.md "Contributing Markdown").

## Maintainer

* [Derek Steinmoeller](https://github.com/dsteinmo)

## Lead Developers

* [Derek Steinmoeller](https://github.com/dsteinmo)
* [Killian Miller](https://github.com/k7miller)

## Documentation

[Here](https://wqcg.github.io/blitzdg "blitzdg Documentation")

## License

[GNU Public License Version 3](https://www.gnu.org/licenses/gpl-3.0.en.html "GPLv3 License")
