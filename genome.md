# Reference genomes

## Arabidopsis thaliana Col-0

```shell script
mkdir -p ~/data/organelles/genome/col_0
cd ~/data/organelles/genome/col_0

wget -N ftp://ftp.ensemblgenomes.org/pub/release-45/plants/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa.gz

faops order Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa.gz \
    <(for chr in {1,2,3,4,5,Mt,Pt}; do echo $chr; done) \
    genome.fa

bwa index genome.fa

```

