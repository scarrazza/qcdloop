[![Build Status](https://travis-ci.org/scarrazza/qcdloop.svg?branch=master)](https://travis-ci.org/scarrazza/qcdloop)

![alt text](https://raw.githubusercontent.com/scarrazza/qcdloop/master/extra/logo.png "Logo QCDLoop")

QCDLoop: an object-oriented one-loop scalar Feynman integrals framework

# General information

Homepage with library description: https://qcdloop.web.cern.ch

If you use this code in your publication, please cite
[arXiv:0712.1851](http://arxiv.org/abs/0712.1851) and
[arXiv:1605.03181](http://arxiv.org/abs/1605.03181).

## Download

You can obtain QCDLoop releases directly from the github repository:

https://github.com/scarrazza/qcdloop/releases

For the last development version you can clone the master code:

```Shell
git clone https://github.com/scarrazza/qcdloop.git
```

For the latest tag:

```Shell
git tag -l
git checkout tags/tag_name
```

## Installation

### From source

Checkout the code and compile the code using the
following procedure:

```Shell
mkdir build
cd build
cmake ..
make && make install
```

By the default, if prefix is not set the program is installed in
`/usr/local`. If you define a custom prefix, use the `-DCMAKE_INSTALL_PREFIX` option and
remember to export `qcdloop/lib` to the `LD_LIBRARY_PATH`. QCDLoop requires a compiler with
C++11 and `quadmath` features (e.g. `gcc >= 5`).

Other qcdloop cmake options are:
- `ENABLE_EXAMPLES`, build examples in C++, default `OFF`.
- `ENABLE_FORTRAN_WRAPPER`, include fortran wrapper in the library, default `ON`.

The fortran wrapper follows the previous syntax in `qcdloop`, see details in table 2 of https://arxiv.org/pdf/1605.03181.pdf.

### From conda-forge

QCDLoop is packaged and distributed on [conda-forge](https://github.com/conda-forge/qcdloop-feedstock/) for the following platforms:

[![Conda Version](https://img.shields.io/conda/vn/conda-forge/qcdloop.svg)](https://anaconda.org/conda-forge/qcdloop)
[![Conda Platforms](https://img.shields.io/conda/pn/conda-forge/qcdloop.svg)](https://anaconda.org/conda-forge/qcdloop)

To install and add QCDLoop to a project with [`pixi`](https://pixi.sh/), from the project directory run

```
pixi add qcdloop
```

and to install into a particular conda environment with [`conda`](https://docs.conda.io/projects/conda/), in the activated environment run

```
conda install --channel conda-forge qcdloop
```

## Contact Information

Maintainer: Stefano Carrazza
