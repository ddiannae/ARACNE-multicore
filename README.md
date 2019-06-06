# ARACNE-multicore


- ``ARACNE2.tgz`` - Codigo fuente para ARACNE 2.
- parallel - scripts para correr aracne sobre una plataforma multicore
- work - Carpeta para ejecutar los trabajos


## Ejecución

Mover a la carpeta ``work``:

``bash run.sh norm-Stage1.tsv &> salida &``

como parametro solo se necesita la matriz de expresion.


## Compilación de ARACNE

```bash
tar xzvf ARACNE2.tgz
cd ARACNE
make clean
make
```
