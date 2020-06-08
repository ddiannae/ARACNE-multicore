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

Mover a la carpeta ``work``:

``bash run.sh norm-Stage1.tsv &> salida &``

Como parámetro sólo se necesita la matriz de expresión.


## Compilación de ARACNE

```bash
bash compile-aracne.sh
```


