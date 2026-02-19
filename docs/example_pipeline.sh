# pip install . 

# conda install mafft -c bioconda 
# conda install -c bioconda iqtree
# conda install 'jclusterfunk>=0.0.25' -c bioconda -c conda-forge

raccoon seq-qc examples/mev/mev_sample.fasta -o examples/mev/seq-qc/mev_sample.seq_qc.fasta --metadata examples/mev/mev_sample.tsv --metadata-id-field accessionVersion --metadata-location-field geoLocCountry --metadata-date-field sampleCollectionDate --min-length 14000 --max-n-content 0.2

mkdir examples/mev/aln-qc/

mafft examples/mev/seq-qc/mev_sample.seq_qc.fasta > examples/mev/aln-qc/mev_sample.aln.fasta

raccoon aln-qc examples/mev/aln-qc/mev_sample.aln.fasta -d examples/mev/aln-qc/

raccoon mask examples/mev/aln-qc/mev_sample.aln.fasta --mask-file examples/mev/aln-qc/mask_sites.csv -d examples/mev/masked/
# realign if sequence removed? 

iqtree -s examples/mev/masked/mev_sample.aln.masked.fasta -m HKY -czb -blmin 0.00000001 -asr  -o 'PP_003MAAS.2||2019' -redo

jclusterfunk prune  -i "examples/mev/masked/mev_sample.aln.masked.fasta.treefile" -t 'PP_003MAAS.2||2019' -o 'examples/mev/masked/mev_sample.pruned.tree'

raccoon tree-qc --phylogeny 'examples/mev/masked/mev_sample.pruned.tree' --asr-state examples/mev/masked/mev_sample.aln.masked.fasta.state --alignment examples/mev/masked/mev_sample.aln.masked.fasta -d examples/mev/tree-qc/
