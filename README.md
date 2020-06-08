# ARACNE-multicore


- ``ARACNE2.tgz`` -> Codigo fuente para ARACNE 2.
- ``parallel`` -> scripts para correr aracne sobre una plataforma multicore
- ``work`` -> Carpeta para ejecutar los trabajos

## Prerequisitos

- python 3
- pandas
- multiprocessing 
- functools
- os, sys

## Ejecución

Copy  matrix gene expression into directory ``launch``:

``bash run.sh norm-Stage1.tsv &> salida &``

## Compilación de ARACNE

```bash
bash compile-aracne.sh
```

## Install bioconda

```bash
bash install-miniconda.sh
```
