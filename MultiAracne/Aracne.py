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

        # Check for existing .adj files and skip those genes
        existing_files = set()
        if os.path.exists(outdir):
            for fn in os.listdir(outdir):
                if fn.endswith(".adj"):
                    # Remove .adj extension to get gene name
                    gene_name = fn[:-4]
                    existing_files.add(gene_name)
        
        # Filter out genes that already have .adj files
        genes_to_process = [gene for gene in self.genes if gene not in existing_files]
        
        if len(existing_files) > 0:
            print(f"Found {len(existing_files)} existing .adj files, skipping those genes")
        
        if len(genes_to_process) == 0:
            print("All genes already have .adj files. Nothing to compute.")
            return
        
        print(f"Running aracne computations for {len(genes_to_process)} features" )
        p = Pool(processes)
        fparam = partial(self.run_gene,  pval=pval, outdir=outdir)
        result = p.map_async(fparam, genes_to_process)
        p.close()
        p.join()
        print(f"\nAracne computations done")

    def run_gene(self, gene, outdir, pval=1):

        if not os.path.exists(outdir) or not os.path.isdir(outdir):
            raise IOError(f"Output directory does not exist: {outdir}")
        
        print(".", end="", flush=True)
        cmd = (f"{os.environ.get('ARACNEHOME')}/aracne2 -H"
            f" {os.environ.get('ARACNEHOME')} -i {self.matrix_file}"
            f" -p {pval} -h {gene}" 
            f" -o {outdir}/{gene.replace('/', '-') if '/' in gene else gene}.adj >"
            f" {outdir}/{gene.replace('/', '-') if '/' in gene else gene}.log")
        os.system(cmd)
    
    ### Build full matrix 
    def build_triu_missing_genes(self, outdir, outfile, row_names=False, chunk_size=100):
        
        if not os.path.exists(outdir) or not os.path.isdir(outdir):
            raise IOError(f"Output directory does not exist: {outdir}")
        
        print("Building upper triangular matrix in chunks")
        start_time = time.time()
        
        # First pass: collect all unique genes that appear in the network
        all_genes_set = set()
        gene_data = {}  # Store parsed data: {gene: {target_gene: mi_value}}
        
        print(f"Reading {len(self.genes)} adjacency files...")
        for idx, gene in enumerate(self.genes):
            if (idx + 1) % 100 == 0:
                print(f"  Processed {idx + 1}/{len(self.genes)} files")
            
            fname = outdir + "/" + gene + ".adj"
            if not os.path.exists(fname):
                raise IOError(f"Missing aracne file: {gene}.adj")

            with open(fname, "r") as fh:
                lines = fh.readlines()
                if lines:
                    last_line = lines[-1]
                    if not last_line.startswith(">"):
                        line = last_line.split("\t")
                        # Parse gene interactions: odd indices are gene names, even are MI values
                        interactions = {}
                        for i in range(1, len(line), 2):
                            if i + 1 < len(line):
                                target_gene = line[i]
                                mi_value = float(line[i + 1])
                                interactions[target_gene] = mi_value
                                all_genes_set.add(target_gene)
                        
                        # Add self-interaction
                        interactions[gene] = 1.0
                        gene_data[gene] = interactions
                        all_genes_set.add(gene)
        
        # Create sorted list of all genes
        all_genes_sorted = sorted(list(all_genes_set))
        n_genes = len(all_genes_sorted)
        gene_to_idx = {gene: idx for idx, gene in enumerate(all_genes_sorted)}
        
        print(f"Total unique genes in network: {n_genes}")
        print(f"--- Reading completed in {time.time() - start_time:.2f} seconds ---")
        
        # Write output in chunks to minimize memory usage
        print("Writing triangular matrix to file...")
        write_start = time.time()
        
        with open(outfile, 'w') as out_fh:
            # Write header
            if row_names:
                out_fh.write(',')
            out_fh.write(','.join(all_genes_sorted) + '\n')
            
            # Process and write row by row
            for row_idx, row_gene in enumerate(all_genes_sorted):
                if (row_idx + 1) % 100 == 0:
                    print(f"  Writing row {row_idx + 1}/{n_genes}")
                
                # Build row with upper triangular constraint
                row_values = []
                for col_idx, col_gene in enumerate(all_genes_sorted):
                    if col_idx < row_idx:
                        # Lower triangle: set to NaN
                        row_values.append('')
                    elif row_gene in gene_data and col_gene in gene_data[row_gene]:
                        # Upper triangle or diagonal: use actual value
                        row_values.append(str(gene_data[row_gene][col_gene]))
                    elif col_gene in gene_data and row_gene in gene_data[col_gene]:
                        # Check symmetric entry
                        row_values.append(str(gene_data[col_gene][row_gene]))
                    else:
                        # No edge: set to NaN
                        row_values.append('')
                
                # Write row
                if row_names:
                    out_fh.write(row_gene + ',')
                out_fh.write(','.join(row_values) + '\n')
        
        print(f"--- Writing completed in {time.time() - write_start:.2f} seconds ---")
        print(f"--- Total time: {time.time() - start_time:.2f} seconds ---")

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
        start_time = time.time()
        mi_matrix = np.concatenate(mi_matrix)
        mi_matrix = np.triu(mi_matrix)
        mi_df = pd.DataFrame(data=mi_matrix, index=self.genes, columns=self.genes)
        mi_df.replace(0, np.nan, inplace=True)
        print("--- %s seconds ---" % (time.time() - start_time))
        mi_df.to_csv(outfile, index=False)

     ### Build full matrix 
    def build_nm_matrix(self, outdir, outfile, genes_filter):
        
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
                                  index=[line[i] for i in range(1, len(line), 2)] + [gene], name=gene).loc[genes_filter])
                        
        print("--- %s seconds ---" % (time.time() - start_time))

        print("Building data frame")
        start_time = time.time()
        mi_df = pd.concat(mis, axis=1)
        mi_df.sort_index(axis=0, inplace=True)
        mi_df.sort_index(axis=1, inplace=True)
        print("--- %s seconds ---" % (time.time() - start_time))
        mi_df.to_csv(outfile, index=True)

