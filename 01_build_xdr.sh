#!/bin/bash
FC=ifort 
CXX=icc 


tar -xzvf xdrfile-1.1.4.tar.gz
cd xdrfile-1.1.4
CDIR=`pwd`

./configure FC=${FC} CXX=${CXX} --prefix=${CDIR}/../xdr

make
make install

cd ../xdr

##
