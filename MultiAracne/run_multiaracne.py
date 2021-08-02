from multiaracne import MultiAracne
import sys, os

if __name__=='__main__':
    exp_matrix = sys.argv[1]
    outdir = sys.argv[2]
    outmatrix = sys.argv[3]
    procs = int(sys.argv[4])

    if not os.path.exists(outdir):
        os.mkdir(outdir)
        
    #ma = MultiAracne(exp_matrix)
    #ma.run_aracne(processes=procs, outdir=outdir, pval=1)
    #ma.join_adj(outdir, outdir+"/matrix.adj")
    MultiAracne.adj_to_matrix(outdir+"/matrix.adj", outmatrix)
