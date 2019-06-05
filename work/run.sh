# Batch run
# for i in $(ls *.tsv); do bash run.sh $i; done > salida &

RAIZ="/HD-BIO/ARACNE2/parallel"

[[ $1 == "" ]] && echo "nead tsv file with expression" && exit 15

ftsv=$1
# ftsv="norm-Control.tsv"
nom=$(echo $ftsv | cut -d. -f 1)

awk '{print $1}' $ftsv > node.list
cname=$(head -1 node.list)
echo "Column name: $cname"

SECONDS=0
python $RAIZ/aracne-par.py $ftsv node.list $cname $(nproc) &> aracne.log 
echo "ARACNe time: $(echo $SECONDS/60 | bc -l) minutes."

SECONDS=0
n=$( (cd adj; ls) | head -1 | cut -d'.' -f 2 )
echo "Parameters to join: $nom $n node.list $cname"
python $RAIZ/joinadj.py $nom $n node.list $cname
echo "join ADJ matriz time: $(echo $SECONDS/60 | bc -l) minutes."

echo "Moving adjancy matrix"
mv adj/mat.adj .

SECONDS=0
python $RAIZ/adj2sif.py > ${nom}.sif
echo "Creating SIF: $(echo $SECONDS/60 | bc -l) minutes."

SECONDS=0
sort -r -k3,3 ${nom}.sif | head -10000 > ${nom}-IM-1e5.txt
echo "Sorting: $(echo $SECONDS/60 | bc -l) minutes."

rm -rf adj log mat.adj node.list
