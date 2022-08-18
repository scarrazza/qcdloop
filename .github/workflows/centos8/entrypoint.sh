#!/bin/sh -l
set -x
uname -a 
cat /etc/issue
yum -y install  epel-release dnf-plugins-core
dnf config-manager --set-enabled powertools
yum -y install  gcc gcc-c++ gcc-gfortran make which cmake cmake-data cmake-filesystem 

cmake -S . -B BUILD -DCMAKE_INSTALL_PREFIX=$(pwd)/INSTALL
cmake --build BUILD
cmake --install BUILD

out=$?
echo ::set-output name=out::$out
