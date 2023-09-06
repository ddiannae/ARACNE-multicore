# -*- coding: utf-8 -*-
#!/usr/bin/env python

from multiprocessing import Pool
from functools import partial
import pandas as pd
import numpy as np
import os
import time

class Aracne:

    def __init__(self, matrix_file, genes_file=None):

        self.matrix_file = matrix_file

        if not genes_file:
            self.genes = []
            with open(matrix_file) as file_in:
                for line in file_in:
                    self.genes.append(line.split("\t")[0])
            # First line has headers
            self.genes.pop(0)
        else:    
            with open(genes_file, "r") as file:
                self.genes = [line.strip() for line in file]

        self.genes.sort()
            
    def run(self, outdir, processes=None, pval=1, aracnehome=None):  

        if not os.path.exists(outdir) or not os.path.isdir(outdir):
            raise IOError(f"Output directory does not exist: {outdir}")
        
        if not os.environ.get('ARACNEHOME') and not aracnehome:
            raise AttributeError(f"Path for aracne not provided. Set ARACNEHOME env variable.")

        print(f"Running aracne computations for {len(self.genes)} features" )
        p = Pool(processes)
        fparam = partial(self.run_gene,  pval=pval, outdir=outdir)
        result = p.map_async(fparam, self.genes)
        p.close()
        p.join()
        print(f"\nAracne computations done")

    def run_gene(self, gene, outdir, pval=1):

        if not os.path.exists(outdir) or not os.path.isdir(outdir):
            raise IOError(f"Output directory does not exist: {outdir}")
        
        print(".", end="")
        cmd = (f"{os.environ.get('ARACNEHOME')}/aracne2 -H"
            f" {os.environ.get('ARACNEHOME')} -i {self.matrix_file}"
            f" -p {pval} -h {gene}" 
            f" -o {outdir}/{gene.replace('/', '-') if '/' in gene else gene}.adj >"
            f" {outdir}/{gene.replace('/', '-') if '/' in gene else gene}.log")
        os.system(cmd)
    
    ### Build full matrix 
    def build_triu_missing_genes(self, outdir, outfile):
        
        if not os.path.exists(outdir) or not os.path.isdir(outdir):
            raise IOError(f"Output directory does not exist: {outdir}")
        
        file_names = [outdir+"/"+fn for fn in os.listdir(outdir) if fn.endswith(".adj")]

        mis = []

        print("Reading files")
        start_time = time.time()

        for gene in self.genes:
            fname = outdir+"/"+gene+".adj"
            if not os.path.exists(fname):
                raise IOError(f"Missing aracne file: {gene}.adj")

            with open(fname, "r") as fh:
                lines = fh.readlines()                
                if lines:
                    last_line = lines[-1]
                    if not last_line.startswith(">"):
                        line = last_line.split("\t")
                        mis.append(pd.Series([float(line[i+1]) for i in range(1, len(line), 2)] + [1], 
                                  index=[line[i] for i in range(1, len(line), 2)] + [gene], name=gene))
                        
        print("--- %s seconds ---" % (time.time() - start_time))

        print("Building data frame")
        start_time = time.time()
        mi_df = pd.concat(mis, axis=1)
        mi_df.sort_index(axis=0, inplace=True)
        mi_df.sort_index(axis=1, inplace=True)
        mi_matrix = mi_df.to_numpy()
        mi_matrix = np.triu(mi_matrix)
        mi_df = pd.DataFrame(data=mi_matrix, index=self.genes, columns=self.genes)
        print("--- %s seconds ---" % (time.time() - start_time))
        mi_df.to_csv(outfile, index=False)

    def build_triu(self, outdir, outfile):
        
        if not os.path.exists(outdir) or not os.path.isdir(outdir):
            raise IOError(f"Output directory does not exist: {outdir}")
        
        file_names = [outdir+"/"+fn for fn in os.listdir(outdir) if fn.endswith(".adj")]

        mi_matrix = []

        print("Reading files and building matrix")
        start_time = time.time()

        for gene in self.genes:
            fname = outdir+"/"+gene+".adj"
            if not os.path.exists(fname):
                raise IOError(f"Missing aracne file: {gene}.adj")

            with open(fname, "r") as fh:
                lines = fh.readlines()                
                if lines:
                    last_line = lines[-1]
                    if not last_line.startswith(">"):
                        line = last_line.split("\t")
                        serie = pd.Series([float(line[i+1]) for i in range(1, len(line), 2)], 
                                  index=[line[i] for i in range(1, len(line), 2)])
                        serie = pd.concat([pd.Series([1], index=[gene]), serie])
                        serie.sort_index(inplace=True)
                        mi_matrix.append(np.array(serie.array, ndmin=2))

        print("--- %s seconds ---" % (time.time() - start_time))

        print("Building data frame")
        mi_matrix = np.concatenate(mi_matrix)
        mi_matrix = np.triu(mi_matrix)
        mi_df = pd.DataFrame(data=mi_matrix, index=self.genes, columns=self.genes)
        start_time = time.time()
        print("--- %s seconds ---" % (time.time() - start_time))
        mi_df.to_csv(outfile, index=False)
    
