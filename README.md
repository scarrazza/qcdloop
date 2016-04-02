![alt text](https://github.com/scarrazza/qcdloop/raw/master/resources/logo.png "Logo QCDLoop")

QCDLoop: an object-oriented one-loop scalar Feynman integrals framework

# General information

Homepage with library description: http://cern.ch/qcdloop

If you use this code in your publication, please cite
[arXiv:0712.1851](http://arxiv.org/abs/0712.1851) and
[arXiv:16xx.xxxxx](http://arxiv.org/abs/16xx.xxxxx).
 
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

Checkout the code and compile the code using the
following procedure:

```Shell
cd qcdloop
./configure --prefix=/where/install/qcdloop #(optional)
make && make install
```

By the default, if prefix is not set the program is installed in
/usr/local. If you define a custom prefix, remember to export
qcdloop/lib to the LD_LIBRARY_PATH. QCDLoop requires a compiler with
C++11 and quadmath features (e.g. gcc >= 4.7)

## Contact Information

Maintainer: Stefano Carrazza