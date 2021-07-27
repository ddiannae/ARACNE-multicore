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

%runscript 
	 if [ $# -ne 3 ]; then
        echo "Please provide the name of your matrix, a p-value. and the number of processors"
        exit 1
    fi
    awk '{print $1}' $1 > gene.list
	mkdir output
	python /ARACNE-multicore/parallel/aracne-par.py $1 gene.list $2 $3 output

%environment
    export ARACNEHOME=/ARACNE-multicore/src/ARACNE

