# ARACNE-multicore


- ``ARACNE2.tgz`` -> Codigo fuente para ARACNE 2.
- ``parallel`` -> scripts para correr aracne sobre una plataforma multicore
- ``work`` -> Carpeta para ejecutar los trabajos

## Prerequisitos

python 2.7
import pandas as pd
from multiprocessing import Pool
from functools import partial
import os, sys

## Ejecución

Mover a la carpeta ``work``:

``bash run.sh norm-Stage1.tsv &> salida &``

Como parámetro sólo se necesita la matriz de expresión.


## Compilación de ARACNE

```bash
tar xzvf ARACNE2.tgz
cd ARACNE
make clean
make
```


