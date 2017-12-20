# blitzdg

blitzdg is an open-source project aiming to implement parallel discontinuous Galerkin (dg) solvers for common partial differential equations systems using blitz++ for array and tensor manipulations and MPI for distributed parallelism.

## Running

Currently only supporting running/development on linux systems, primarily ubuntu. Tested with GNU g++ compiler and GNU make.

1. `git clone https://github.com/dsteinmo/blitzdg.git`
2. `cd blitzdg && ./pull-deps.sh`
3.  `make && ./bin/blitzdg` (The binary currently doesn't do much).
4. Run unit tests with `make tests && ./bin/test`.

## Dependencies

So far: 

* `blitz++`
* `SuiteSparse (umfpack)`
* `metis`
* `libgmsh`
* `igloo` for BDD-style testing.

Dependency installation is outlined in `pull-deps.sh` (tested on Ubuntu).

## Contributing

We accept pull requests. 

If you add code, please write tests using the igloo testing framework that is included as a project dependency.

## Maintainer

* [Derek Steinmoeller](https://github.com/dsteinmo)

## License

MIT