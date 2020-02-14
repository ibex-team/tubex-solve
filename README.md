# [tubex-solve](http://simon-rohou.fr/research/tubex-lib) [![Build Status](https://travis-ci.org/ibex-team/tubex-solve.svg)](https://travis-ci.org/ibex-team/tubex-solve)


### Installation
--------------------------------------

The [`Tubex library`](https://github.com/SimonRohou/tubex-lib) have to be installed first. See the [installation notice](http://simon-rohou.fr/research/tubex-lib/02_installation/index.html).

Then the library `tubex-solve` can be compiled using the following commands:
```bash
git clone https://github.com/ibex-team/tubex-solve
cd tubex-solve
mkdir -p make && cd make
make
```

The problems are compiled together with the `Solver` class.
One can try the first problem with:
```bash
VIBes-viewer &
./problems/01_picard/01_picard
```