Bootstrap: docker
From: continuumio/miniconda3

%post
    apt-get update
    apt-get install -y build-essential git

    /opt/conda/bin/conda install -y numpy pandas
    
    git clone https://github.com/ddiannae/ARACNE-multicore
    cd ARACNE-multicore/src
    tar xf ARACNE.src.tar.gz
    cd ARACNE
    make clean
    make
    cd ../../
    python3 setup.py install

%runscript 
    python $@

%environment
    export ARACNEHOME=/ARACNE-multicore/src/ARACNE

