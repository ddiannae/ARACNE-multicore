#!/bin/bash

cd src
tar xzvf ARACNE.src.tar.gz
cd ARACNE
make clean
make
cd ..
ln -s ../bin ARACNE 