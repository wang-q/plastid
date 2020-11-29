# Reference genomes

## Arabidopsis thaliana Col-0

```shell script
mkdir -p ~/data/organelles/genome/col_0
cd ~/data/organelles/genome/col_0

wget -N ftp://ftp.ensemblgenomes.org/pub/release-45/plants/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa.gz

faops order Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa.gz \
    <(for chr in {1,2,3,4,5,Mt,Pt}; do echo $chr; done) \
    genome.fa

# bwa index
bwa index genome.fa

# chr.sizes
faops size genome.fa > chr.sizes

```


## Oryza sativa Japonica Group Cultivar Nipponbare

```shell script
mkdir -p ~/data/organelles/genome/nip
cd ~/data/organelles/genome/nip

wget -N ftp://ftp.ensemblgenomes.org/pub/release-45/plants/fasta/oryza_sativa/dna/Oryza_sativa.IRGSP-1.0.dna_sm.toplevel.fa.gz

faops order Oryza_sativa.IRGSP-1.0.dna_sm.toplevel.fa.gz \
    <(for chr in $(seq 1 1 12) Mt Pt; do echo $chr; done) \
    genome.fa

# bwa index
bwa index genome.fa

# chr.sizes
faops size genome.fa > chr.sizes

```

