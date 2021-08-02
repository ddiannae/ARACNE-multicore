# -*- coding: utf-8 -*-
#!/usr/bin/env python

from multiprocessing import Pool
from functools import partial
import pandas as pd
import numpy as np
import os

class MultiAracne:

    def __init__(self, matrix_file, genes_file=None):

        self.matrix_file = matrix_file

        if not genes_file:
            self.genes = []
            with open(matrix_file) as file_in:
                for line in file_in:
                    self.genes.append(line.split("\t")[0])
        else:    
            with open(genes_file, "r") as file:
                self.genes = [line.strip() for line in file]

    def run_aracne(self, outdir, processes=None, pval=1, aracnehome=None):  

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
            f" -o {outdir}/{gene}.adj > {outdir}/{gene}.log")
        os.system(cmd)

    def join_adj(self, outdir, outfile="output.adj"):

        if not os.path.exists(outdir) or not os.path.isdir(outdir):
            raise IOError(f"Output directory does not exist: {outdir}")

        file_names = [outdir+"/"+fn for fn in os.listdir(outdir) if fn.endswith(".adj")]
        print(f"Joining {len(file_names)}*.adj files in {outdir}") 

        wlines = []
        for fname in file_names:
            fh = open(fname, "r")
            lines = fh.read().splitlines()
            if lines:
                last_line = lines[-1]
                wlines.append(last_line)
            fh.close()
        
        with open(outfile, 'w') as f:
            for l in wlines:
                if not l.startswith(">"):
                    f.write("%s\n" % l)
    
    @staticmethod
    def adj_to_matrix(adjfile, output, genes_file=None):
        
        print("Reading file")
        with open(adjfile) as file_in:
            lines = []
            for line in file_in:
                lines.append(line.split("\t"))
        
        if not genes_file:
             genes = [l[0] for l in lines]
        else:     
            with open(genes_file) as file_in:
                genes = []
                for line in file_in:
                    genes.append(line.strip())
        
        print("Building dictionaries")
        mi_genes = [l[0] for l in lines]
        mi_dicts = [dict([(l[i], float(l[i+1])) for i in range(1, len(l), 2)]) for l in lines]
        
        print("Building all genes series")
        mi_series = [pd.Series(mid) for mid in mi_dicts]
        genes_serie = pd.Series(np.nan, index=genes)
        sum_series = [s.add(genes_serie, fill_value = 0) for s in mi_series]
        
        print("Building data frame")
        mi_df = pd.concat(sum_series, axis=1)
        mi_df.columns = mi_genes
        
        print("Adding missing columns")
        missing_columns = [c for c in mi_df.index if c not in mi_df.columns]
        mi_df[missing_columns] = np.nan
        
        print("Sorting data frame")
        mi_df.sort_index(axis=0, inplace=True)
        mi_df.sort_index(axis=1, inplace=True)
        cols = mi_df.columns
        
        print("Getting triangular matrix")
        mi_matrix = mi_df.to_numpy()
        mi_matrix = np.triu(mi_matrix)
        
        print("Saving data frame")
        mi_df = pd.DataFrame(data=mi_matrix, index=cols, columns=cols)
        mi_df.replace(0, np.nan, inplace=True)
        mi_df.to_csv(output, index=False)
        
    @staticmethod
    def adj_to_unsif(adjfile, output):
        
        print("Reading values")
        with open(adjfile) as file_in:
            mi_vals = []
            for line in file_in:
                line_vals = line.split("\t")
                for i in range(1, len(line_vals), 2):
                    st = [line_vals[0], line_vals[i]]
                    st.sort()
                    st.append(float(line_vals[i+1]))
                    mi_vals.append(st)
        
        print(len(mi_vals))
        print("Building sif")
        df = pd.DataFrame(mi_vals, columns=["source", "target", "mi"])
        df.sort_values(by = ["mi"], ascending = False, inplace = True)
        df.drop_duplicates(inplace = True)
        print(df.shape)
        print("Saving data frame")
        df.to_csv(output, sep="\t", index=False, header=False)


