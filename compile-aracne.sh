#!/bin/bash

cd src
tar xzvf ARACNE.src.tar.gz
cd ARACNE
make clean
make
cd ..
ln -s ../bin ARACNE

[[ -f "../bin/usage.txt"]] && echo "Finished successfully." && exit 0
echo "Something was wrong"
exit 15
