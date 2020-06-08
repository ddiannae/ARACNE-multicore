#!/bin/bash

cd src
tar xzvf ARACNE2.tgz
cd ARACNE
make clean
make
cd ..
ln -s ARACNE ../bin