# Reference genomes

## *Arabidopsis thaliana* Col-0

```shell script
mkdir -p ~/data/plastid/genome/col_0
cd ~/data/plastid/genome/col_0

wget -N ftp://ftp.ensemblgenomes.org/pub/release-45/plants/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa.gz

faops order Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa.gz \
    <(for chr in {1,2,3,4,5,Mt,Pt}; do echo $chr; done) \
    genome.fa

# chr.sizes
faops size genome.fa > chr.sizes

```


## *Oryza sativa* Japonica Group Cultivar Nipponbare

```shell script
mkdir -p ~/data/plastid/genome/nip
cd ~/data/plastid/genome/nip

wget -N ftp://ftp.ensemblgenomes.org/pub/release-45/plants/fasta/oryza_sativa/dna/Oryza_sativa.IRGSP-1.0.dna_sm.toplevel.fa.gz

faops order Oryza_sativa.IRGSP-1.0.dna_sm.toplevel.fa.gz \
    <(for chr in $(seq 1 1 12) Mt Pt; do echo $chr; done) \
    genome.fa

# chr.sizes
faops size genome.fa > chr.sizes

```

## *Medicago truncatula* A17

```shell script
mkdir -p ~/data/plastid/genome/a17
cd ~/data/plastid/genome/a17

for ACCESSION in "NC_003119" "NC_029641"; do
    URL=$(printf "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=%s&id=%s&retmode=text" "fasta" "${ACCESSION}");
    curl $URL -o ${ACCESSION}.fa
done

aria2c -x 4 -s 2 -c ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/219/495/GCF_000219495.3_MedtrA17_4.0/GCF_000219495.3_MedtrA17_4.0_genomic.fna.gz

TAB=$'\t'
cat <<EOF > replace.tsv
NC_016407.2${TAB}1
NC_016408.2${TAB}2
NC_016409.2${TAB}3
NC_016410.2${TAB}4
NC_016411.2${TAB}5
NC_016412.2${TAB}6
NC_016413.2${TAB}7
NC_016414.2${TAB}8
NC_003119.8${TAB}Pt
NC_029641.1${TAB}Mt
EOF

gzip -dcf GCF_000219495.3_MedtrA17_4.0_genomic.fna.gz NC_003119.fa NC_029641.fa |
    faops replace stdin replace.tsv stdout |
    faops order stdin <(for chr in $(seq 1 1 8) Mt Pt; do echo $chr; done) genome.fa

# chr.sizes
faops size genome.fa > chr.sizes

```

## *Prunus persica* PLov2-2N (a double haploid genotype of the peach cv. Lovell)

```shell script
mkdir -p ~/data/plastid/genome/lovell
cd ~/data/plastid/genome/lovell

aria2c -x 4 -s 2 -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/346/465/GCF_000346465.2_Prunus_persica_NCBIv2/GCF_000346465.2_Prunus_persica_NCBIv2_genomic.fna.gz

TAB=$'\t'
cat <<EOF > replace.tsv
NC_034009.1${TAB}G1
NC_034010.1${TAB}G2
NC_034011.1${TAB}G3
NC_034012.1${TAB}G4
NC_034013.1${TAB}G5
NC_034014.1${TAB}G6
NC_034015.1${TAB}G7
NC_034016.1${TAB}G8
NC_014697.1${TAB}Pt
EOF

gzip -dcf GCF*_genomic.fna.gz |
    faops replace -s stdin replace.tsv genome.fa

# chr.sizes
faops size genome.fa > chr.sizes

```

## *Solanum lycopersicum* Cultivar: Heinz 1706

```shell script
mkdir -p ~/data/plastid/genome/h1706
cd ~/data/plastid/genome/h1706

aria2c -x 4 -s 2 -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/188/115/GCF_000188115.4_SL3.0/GCF_000188115.4_SL3.0_genomic.fna.gz

TAB=$'\t'
cat <<EOF > replace.tsv
NC_015438.3${TAB}1
NC_015439.3${TAB}2
NC_015440.3${TAB}3
NC_015441.3${TAB}4
NC_015442.3${TAB}5
NC_015443.3${TAB}6
NC_015444.3${TAB}7
NC_015445.3${TAB}8
NC_015446.3${TAB}9
NC_015447.3${TAB}10
NC_015448.3${TAB}11
NC_015449.3${TAB}12
NC_035963.1${TAB}Mt
NC_007898.3${TAB}Pt
EOF

gzip -dcf GCF*_genomic.fna.gz |
    faops replace stdin replace.tsv stdout |
    faops order stdin <(for chr in $(seq 1 1 12) Mt Pt; do echo $chr; done) genome.fa

# chr.sizes
faops size genome.fa > chr.sizes

```

## *Glycine max* Williams 82

```shell script
mkdir -p ~/data/plastid/genome/w82
cd ~/data/plastid/genome/w82

aria2c -x 4 -s 2 -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/515/GCF_000004515.5_Glycine_max_v2.1/GCF_000004515.5_Glycine_max_v2.1_genomic.fna.gz

TAB=$'\t'
cat <<EOF > replace.tsv
NC_016088.3${TAB}1
NC_016089.3${TAB}2
NC_016090.3${TAB}3
NC_016091.3${TAB}4
NC_038241.1${TAB}5
NC_038242.1${TAB}6
NC_038243.1${TAB}7
NC_038244.1${TAB}8
NC_038245.1${TAB}9
NC_038246.1${TAB}10
NC_038247.1${TAB}11
NC_038248.1${TAB}12
NC_038249.1${TAB}13
NC_038250.1${TAB}14
NC_038251.1${TAB}15
NC_038252.1${TAB}16
NC_038253.1${TAB}17
NC_038254.1${TAB}18
NC_038255.1${TAB}19
NC_038256.1${TAB}20
NC_007942.1${TAB}Pt
NC_020455.1${TAB}Mt
EOF

gzip -dcf GCF_000004515.5_Glycine_max_v2.1_genomic.fna.gz |
    faops replace stdin replace.tsv stdout |
    faops order stdin <(for chr in $(seq 1 1 20) Mt Pt; do echo $chr; done) genome.fa

# chr.sizes
faops size genome.fa > chr.sizes

```
