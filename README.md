# blitzdg

blitzdg is an open-source project aiming to implement parallel discontinuous Galerkin (dg) solvers for common partial differential equations systems using blitz++ for array and tensor manipulations and MPI for distributed parallelism.

[![Build Status](https://travis-ci.org/dsteinmo/blitzdg.svg?branch=master)](https://travis-ci.org/dsteinmo/blitzdg)  [![Coverage Status](https://coveralls.io/repos/github/dsteinmo/blitzdg/badge.svg)](https://coveralls.io/github/dsteinmo/blitzdg)

[![Windows Build Status](https://ci.appveyor.com/api/projects/status/pmx725yhsrnq3thu?svg=true)](https://ci.appveyor.com/project/dsteinmo/blitzdg)

## Running

Currently only supporting running/development on linux systems, primarily ubuntu. Tested with GNU g++ compiler and GNU make.

1. `git clone https://github.com/dsteinmo/blitzdg.git`
2. `cd blitzdg && ./pull-deps.sh`
3. `make && ./bin/blitzdg` (The binary currently doesn't do much).
4. Run unit tests with `make test`.

### Running with Docker

You can also run the build and tests inside a docker container.

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

Dependency installation is outlined in `pull-deps.sh` (tested on Ubuntu) and `pull-deps-mac.sh` (test on Mac OSX Sierra).

## Contributing

We accept pull requests. 

If you add code, please write tests using the igloo testing framework that is included as a project dependency.

## Maintainer

* [Derek Steinmoeller](https://github.com/dsteinmo)

## Documentation

[Here](https://dsteinmo.github.io/blitzdg "blitzdg Documentation")

## License

MIT
