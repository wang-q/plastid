# Mapping to genomes on folds

[TOC levels=1-3]: # ""

- [Mapping to genomes on folds](#mapping-to-genomes-on-folds)
    - [Reference genomes](#reference-genomes)
    - [Download fastq files from ENA](#download-fastq-files-from-ena)
    - [基本信息](#基本信息)
    - [Symlink](#symlink)
    - [Trim, cutoff and mapping](#trim-cutoff-and-mapping)
    - [Combine chromosomes](#combine-chromosomes)
    - [Merge all results](#merge-all-results)
    - [Remove intermediate files](#remove-intermediate-files)

## Tools

* https://www.ibm.com/aspera/connect/
* https://www.biostars.org/p/9528910/
* https://github.com/PRIDE-Archive/pride-inspector/blob/master/aspera/etc/asperaweb_id_dsa.openssh
* `~/.aspera/connect/bin/ascp`

```shell
curl -L https://delivery04-mul.dhe.ibm.com/sar/CMA/OSA/0cz8y/0/ibm-aspera-connect_4.2.14.855-HEAD_linux_x86_64.tar.gz |
    tar xvz

bash ibm-aspera-connect*

cat >> $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh <<EOF
-----BEGIN DSA PRIVATE KEY-----
MIIBuwIBAAKBgQDkKQHD6m4yIxgjsey6Pny46acZXERsJHy54p/BqXIyYkVOAkEp
KgvT3qTTNmykWWw4ovOP1+Di1c/2FpYcllcTphkWcS8lA7j012mUEecXavXjPPG0
i3t5vtB8xLy33kQ3e9v9/Lwh0xcRfua0d5UfFwopBIAXvJAr3B6raps8+QIVALws
yeqsx3EolCaCVXJf+61ceJppAoGAPoPtEP4yzHG2XtcxCfXab4u9zE6wPz4ePJt0
UTn3fUvnQmJT7i0KVCRr3g2H2OZMWF12y0jUq8QBuZ2so3CHee7W1VmAdbN7Fxc+
cyV9nE6zURqAaPyt2bE+rgM1pP6LQUYxgD3xKdv1ZG+kDIDEf6U3onjcKbmA6ckx
T6GavoACgYEAobapDv5p2foH+cG5K07sIFD9r0RD7uKJnlqjYAXzFc8U76wXKgu6
WXup2ac0Co+RnZp7Hsa9G+E+iJ6poI9pOR08XTdPly4yDULNST4PwlfrbSFT9FVh
zkWfpOvAUc8fkQAhZqv/PE6VhFQ8w03Z8GpqXx7b3NvBR+EfIx368KoCFEyfl0vH
Ta7g6mGwIMXrdTQQ8fZs
-----END DSA PRIVATE KEY-----
EOF

chmod 600 $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh

cd $HOME/bin
ln -s ../.aspera/connect/bin/ascp ascp

```

## Reference genomes

* Ensembl Release 114 (May 2025)
* Ensembl Plants release 61

* *Arabidopsis thaliana* Col-0

```shell
mkdir -p ~/data/plastid/genome/col0
cd ~/data/plastid/genome/col0

wget -N https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-61/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa.gz

hnsm order Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa.gz \
    <(for chr in {1,2,3,4,5,Mt,Pt}; do echo $chr; done) \
    -o genome.fa

# chr.sizes
hnsm size genome.fa -o chr.sizes

```

* *Oryza sativa* Japonica Group Cultivar Nipponbare

```shell
mkdir -p ~/data/plastid/genome/nip
cd ~/data/plastid/genome/nip

wget -N https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-61/fasta/oryza_sativa/dna/Oryza_sativa.IRGSP-1.0.dna_sm.toplevel.fa.gz

hnsm order Oryza_sativa.IRGSP-1.0.dna_sm.toplevel.fa.gz \
    <(for chr in $(seq 1 1 12) Mt Pt; do echo $chr; done) \
    -o genome.fa

# chr.sizes
hnsm size genome.fa -o chr.sizes

```

* *Medicago truncatula* A17

The Ensembl version lacks chromosome naming and does not include chloroplast and mitochondrial
genomes, which need to be constructed manually.

```shell
#mkdir -p ~/data/plastid/genome/a17
#cd ~/data/plastid/genome/a17
#
#wget -N https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-61/fasta/medicago_truncatula/dna/Medicago_truncatula.MtrunA17r5.0_ANR.dna_sm.toplevel.fa.gz
#
#hnsm order Medicago_truncatula.MtrunA17r5.0_ANR.dna_sm.toplevel.fa.gz \
#    <(for chr in $(seq 1 1 12) Mt Pt; do echo $chr; done) \
#    -o genome.fa
#
## chr.sizes
#hnsm size genome.fa -o chr.sizes

```

```shell
mkdir -p ~/data/plastid/genome/a17
cd ~/data/plastid/genome/a17

#for ACCESSION in "NC_003119" "NC_029641"; do
#    URL=$(printf "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=%s&id=%s&retmode=text" "fasta" "${ACCESSION}");
#    curl $URL -o ${ACCESSION}.fa
#done

aria2c -x 4 -s 2 -c \
    https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/473/485/GCF_003473485.1_MtrunA17r5.0-ANR/GCF_003473485.1_MtrunA17r5.0-ANR_genomic.fna.gz

TAB=$'\t'
cat <<EOF > replace.tsv
NC_053042.1${TAB}1
NC_053043.1${TAB}2
NC_053044.1${TAB}3
NC_053045.1${TAB}4
NC_053046.1${TAB}5
NC_053047.1${TAB}6
NC_053048.1${TAB}7
NC_053049.1${TAB}8
NC_003119.8${TAB}Pt
NC_029641.1${TAB}Mt
EOF

gzip -dcf GCF_003473485.1_MtrunA17r5.0-ANR_genomic.fna.gz |
    hnsm replace stdin replace.tsv |
    hnsm order stdin <(for chr in $(seq 1 1 8) Mt Pt; do echo $chr; done) -o genome.fa

# chr.sizes
hnsm size genome.fa -o chr.sizes

```

* *Solanum lycopersicum* Micro-Tom

```shell
mkdir -p ~/data/plastid/genome/microtom
cd ~/data/plastid/genome/microtom

aria2c -x 4 -s 2 -c \
    https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/036/512/215/GCF_036512215.1_SLM_r2.1/GCF_036512215.1_SLM_r2.1_genomic.fna.gz

TAB=$'\t'
cat <<EOF > replace.tsv
NC_090800.1${TAB}1
NC_090801.1${TAB}2
NC_090802.1${TAB}3
NC_090803.1${TAB}4
NC_090804.1${TAB}5
NC_090805.1${TAB}6
NC_090806.1${TAB}7
NC_090807.1${TAB}8
NC_090808.1${TAB}9
NC_090809.1${TAB}10
NC_090810.1${TAB}11
NC_090811.1${TAB}12
NC_035963.1${TAB}Mt
NC_007898.3${TAB}Pt
EOF

gzip -dcf GCF_036512215*_genomic.fna.gz |
    hnsm replace stdin replace.tsv |
    hnsm order stdin <(for chr in $(seq 1 1 12) Mt Pt; do echo $chr; done) -o genome.fa

# chr.sizes
hnsm size genome.fa -o chr.sizes

```

* *Solanum lycopersicum* Cultivar: Heinz 1706

```shell
# mkdir -p ~/data/plastid/genome/h1706
# cd ~/data/plastid/genome/h1706

# aria2c -x 4 -s 2 -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/188/115/GCF_000188115.5_SL3.1/GCF_000188115.5_SL3.1_genomic.fna.gz

# TAB=$'\t'
# cat <<EOF > replace.tsv
# NC_015438.3${TAB}1
# NC_015439.3${TAB}2
# NC_015440.3${TAB}3
# NC_015441.3${TAB}4
# NC_015442.3${TAB}5
# NC_015443.3${TAB}6
# NC_015444.3${TAB}7
# NC_015445.3${TAB}8
# NC_015446.3${TAB}9
# NC_015447.3${TAB}10
# NC_015448.3${TAB}11
# NC_015449.3${TAB}12
# NC_035963.1${TAB}Mt
# NC_007898.3${TAB}Pt
# EOF

# gzip -dcf GCF*_genomic.fna.gz |
#     faops replace stdin replace.tsv stdout |
#     faops order stdin <(for chr in $(seq 1 1 12) Mt Pt; do echo $chr; done) genome.fa

# # chr.sizes
# faops size genome.fa > chr.sizes

```

* *Prunus persica* PLov2-2N (a double haploid genotype of the peach cv. Lovell)

```shell
mkdir -p ~/data/plastid/genome/lovell
cd ~/data/plastid/genome/lovell

aria2c -x 4 -s 2 -c \
    https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/346/465/GCF_000346465.2_Prunus_persica_NCBIv2/GCF_000346465.2_Prunus_persica_NCBIv2_genomic.fna.gz

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

gzip -dcf GCF_000346465*_genomic.fna.gz |
    hnsm replace -s stdin replace.tsv -o genome.fa

# chr.sizes
hnsm size genome.fa -o chr.sizes

```

* *Glycine max* Williams 82

```shell
mkdir -p ~/data/plastid/genome/w82
cd ~/data/plastid/genome/w82

aria2c -x 4 -s 2 -c \
    https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/515/GCF_000004515.6_Glycine_max_v4.0/GCF_000004515.6_Glycine_max_v4.0_genomic.fna.gz

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

gzip -dcf GCF*_genomic.fna.gz |
    hnsm replace stdin replace.tsv |
    hnsm order stdin <(for chr in $(seq 1 1 20) Pt Mt; do echo $chr; done) -o genome.fa

# chr.sizes
hnsm size genome.fa -o chr.sizes

```

## Download fastq files from ENA

```shell
mkdir -p ~/data/plastid/ena
cd ~/data/plastid/ena

cat << EOF > source.csv
SRX202246,Atha_Col_0_1,HiSeq 2000 PE100
SRX2527206,Atha_Col_0_2,MiSeq 2000 PE300
SRR616965,Atha_Ler_0,HiSeq 2000 PE100
SRX179254,Osat_Nip,HiSeq 2000 PE100
SRX025260,Osat_Nip_2,Osat_50
SRX673852,Mtru_A17,HiSeq 2000 PE150
SRX150254,Pper_Lovell,Illumina Genome Analyzer IIx PE100
SRX698770,Slyc_H1706,Illumina HiSeq 2000 PE100
SRX7009428,Gmax_W82,HiSeq X Ten
EOF

anchr ena info | perl - -v source.csv > ena_info.yml
anchr ena prep | perl - ena_info.yml --ascp

rgr md ena_info.tsv --fmt

cat ena_info.ascp.sh |
    parallel --no-run-if-empty -j 2 "{}"

#aria2c -j 4 -x 4 -s 1 -c -i ena_info.ftp.txt
md5sum --check ena_info.md5.txt

# Bad quality
# SRR1542422

```

| name         | srx        | platform | layout | ilength | srr         |       spots | bases  |
|--------------|------------|----------|--------|---------|-------------|------------:|--------|
| Atha_Col_0_1 | SRX202246  | ILLUMINA | PAIRED | 450     | SRR616966   |  24,851,796 | 4.63G  |
| Atha_Col_0_1 | SRX202246  | ILLUMINA | PAIRED | 450     | SRR611086   |  49,891,349 | 9.29G  |
| Atha_Col_0_2 | SRX2527206 | ILLUMINA | PAIRED |         | SRR5216995  |  26,893,065 | 14.46G |
| Atha_Ler_0   | SRX202247  | ILLUMINA | PAIRED | 450     | SRR616965   |  25,436,255 | 4.74G  |
| Atha_Ler_0   | SRX202247  | ILLUMINA | PAIRED | 450     | SRR611087   |  50,791,450 | 9.46G  |
| Gmax_W82     | SRX7009428 | ILLUMINA | PAIRED |         | SRR10296600 | 162,110,355 | 45.29G |
| Mtru_A17     | SRX673852  | ILLUMINA | PAIRED | 360     | SRR1542423  |  29,663,436 | 8.34G  |
| Mtru_A17     | SRX673852  | ILLUMINA | PAIRED | 360     | SRR1542422  |  99,418,334 | 16.67G |
| Osat_Nip     | SRX179254  | ILLUMINA | PAIRED | 300     | SRR545231   |  85,251,097 | 16.04G |
| Osat_Nip     | SRX179254  | ILLUMINA | PAIRED | 300     | SRR545059   |  85,148,124 | 7.93G  |
| Osat_Nip_2   | SRX025260  | ILLUMINA | PAIRED | 463     | SRR063638   |  29,784,011 | 5.55G  |
| Pper_Lovell  | SRX150254  | ILLUMINA | PAIRED | 400     | SRR502985   | 123,590,441 | 23.25G |
| Slyc_H1706   | SRX698770  | ILLUMINA | PAIRED |         | SRR1572628  |  24,198,345 | 4.51G  |

## 基本信息

`cutoff = FOLD * DEPTH`

* FOLD 0, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64

## Symlink

```shell
mkdir -p ~/data/plastid/evaluation
cd ~/data/plastid/evaluation

SRRS=(
    'SRR616966::col0'  # Col-0
    'SRR611086::col0'
    'SRR5216995::col0'
    'SRR616965::col0'  # Ler-0
    'SRR611087::col0'
    'SRR545231::nip'    # Nipponbare
    'SRR063638::nip'
    'SRR1542423::a17'   # A17
    'SRR1572628::h1706' # Heinz 1706
)
FOLDS=(0 0.25 0.5 1 2 4 8 16 32 64)

for item in "${SRRS[@]}"; do
    SRR="${item%%::*}"
    GENOME="${item##*::}"

    for FOLD in "${FOLDS[@]}"; do
        BASE_NAME=${SRR}_${FOLD}

        mkdir -p ${BASE_NAME}/1_genome
        pushd ${BASE_NAME}/1_genome

        ln -fs ../../../genome/${GENOME}/genome.fa genome.fa
        cp ../../../genome/${GENOME}/chr.sizes chr.sizes
        popd

        mkdir -p ${BASE_NAME}/2_illumina
        pushd ${BASE_NAME}/2_illumina

        ln -fs ../../../ena/${SRR}_1.fastq.gz R1.fq.gz
        ln -fs ../../../ena/${SRR}_2.fastq.gz R2.fq.gz
        popd

    done
done

```

## Trim, cutoff and mapping

* Rsync to hpcc

```shell
rsync -avP \
    ~/data/plastid/ \
    wangq@202.119.37.251:data/plastid

rsync -avP \
    -e 'ssh -p 8804' \
    ~/data/plastid/ \
    wangq@58.213.64.36:data/plastid

# rsync -avP wangq@202.119.37.251:data/plastid/ ~/data/plastid

```

```shell
cd ~/data/plastid/evaluation

SRRS=(
    'SRR616966::41'  # Col-0
    'SRR611086::78'
    'SRR5216995::121'
    'SRR616965::42'  # Ler-0
    'SRR611087::79'
    'SRR545231::46'  # Nipponbare
    'SRR063638::15'
    'SRR1542423::23' # A17
    'SRR1572628::5'  # Heinz 1706
)
FOLDS=(0 0.25 0.5 1 2 4 8 16 32 64)

for item in "${SRRS[@]}"; do
    SRR="${item%%::*}"
    DEPTH="${item##*::}"

    for FOLD in "${FOLDS[@]}"; do
        CUTOFF=$(bc <<< "(${DEPTH} * ${FOLD}) / 1")

        echo 1>&2 "==> ${item} ${FOLD}"

        BASE_NAME=${SRR}_${FOLD}
        pushd ${BASE_NAME}

        if [ ! -f 3_bwa/join.tsv ]; then
            rm *.sh
            if [[ "${FOLD}" == "0" ]]; then
                anchr template \
                    --genome $(tsv-summarize 1_genome/chr.sizes --sum 2) \
                    --parallel 24 \
                    --xmx 80g \
                    \
                    --fastqc \
                    --insertsize \
                    --kat \
                    \
                    --trim "--dedupe" \
                    --qual "25" \
                    --len "60" \
                    --filter "adapter artifact" \
                    \
                    --bwa Q25L60

                bsub -q mpi -n 24 -J "${BASE_NAME}" "
                    bash 2_fastqc.sh
                    bash 2_insert_size.sh
                    bash 2_kat.sh
                    bash 2_trim.sh
                    bash 9_stat_reads.sh
                    bash 3_bwa.sh
                "
            else
                anchr template \
                    --genome $(tsv-summarize 1_genome/chr.sizes --sum 2) \
                    --parallel 24 \
                    --xmx 80g \
                    \
                    --trim "--dedupe --cutoff ${CUTOFF} --cutk 31" \
                    --qual "25" \
                    --len "60" \
                    --filter "adapter artifact" \
                    \
                    --bwa Q25L60

                bsub -q mpi -n 24 -J "${BASE_NAME}" "
                    bash 2_trim.sh
                    bash 9_stat_reads.sh
                    bash 3_bwa.sh
                "
            fi

        fi

        popd

    done
done

```

## Combine chromosomes

```shell
cd ~/data/plastid/evaluation

SRRS=(
    'SRR616966::5'   # Col-0
    'SRR611086::5'
    'SRR5216995::5'
    'SRR616965::5'   # Ler-0
    'SRR611087::5'
    'SRR545231::12'  # Nipponbare
    'SRR063638::12'
    'SRR1542423::8'  # A17
    'SRR1572628::12' # Heinz 1706
)
FOLDS=(0 0.25 0.5 1 2 4 8 16 32 64)

for item in "${SRRS[@]}"; do
    SRR="${item%%::*}"
    CHR_NUM="${item##*::}"

    for FOLD in "${FOLDS[@]}"; do
        BASE_NAME=${SRR}_${FOLD}

        echo 1>&2 "==> ${BASE_NAME}"

        pushd ${BASE_NAME}/3_bwa

        cat join.tsv |
            grep -v "^Mt" |
            grep -v "^Pt" |
            tsv-summarize -H --sum chrLength,covLength,bases --min min --max max |
            sed '1d' |
            perl -e '
                my $line = <>;
                chomp $line;
                my ($chrLength, $covLength, $bases, $min, $max, ) = split qq(\t), $line;
                my $covRate = sprintf qq(%.4f), $covLength / $chrLength;
                my $mean = sprintf qq(%.2f), $bases / $chrLength;
                print join qq(\t), (
                    "Nc", $chrLength, $covLength, $covRate, $bases, $mean, $min, $max,
                );
                print qq(\n);
            ' |
            (cat join.tsv | sed "2,$((CHR_NUM+1))d" && cat) \
            > combine.tsv

        popd

    done
done

```

## Merge all results

```shell
cd ~/data/plastid/evaluation

SRRS=(
    'SRR616966::Col-0'      # Col-0
    'SRR611086::Col-0'
    'SRR5216995::Col-0'
    'SRR616965::Ler-0'      # Ler-0
    'SRR611087::Ler-0'
    'SRR545231::Nipponbare' # Nipponbare
    'SRR063638::NP'
    'SRR1542423::A17'       # A17
    'SRR1572628::Heinz1706' # Heinz 1706
)
FOLDS=(0 0.25 0.5 1 2 4 8 16 32 64)

for item in "${SRRS[@]}"; do
    SRR="${item%%::*}"
    STRAIN="${item##*::}"

    for FOLD in "${FOLDS[@]}"; do
        BASE_NAME=${SRR}_${FOLD}

        pushd ${BASE_NAME}/3_bwa > /dev/null

        echo -e "Fold\tchrom\n${FOLD}\tNc\n${FOLD}\tMt\n${FOLD}\tPt" |
            tsv-join -H --filter-file combine.tsv --key-fields chrom --append-fields 2-8

        popd > /dev/null

    done |
        tsv-uniq \
        > ${SRR}_folds.tsv

    echo
    echo "Table: ${STRAIN} ${SRR} Folds"
    echo
    for PART in Nc Mt Pt; do
        cat ${SRR}_folds.tsv |
            tsv-filter -H --str-eq chrom:${PART} |
            mlr --itsv --omd cat
        echo
        echo
    done

done

```

Table: Col-0 SRR616966 Folds

| Fold | chrom | chrLength | covLength | covRate | bases      | mean  | min | max   |
|:-----|:------|:----------|:----------|:--------|:-----------|:------|:----|:------|
| 0    | Nc    | 119146348 | 118877276 | 0.9977  | 2503233613 | 21.01 | 0   | 22463 |
| 0.25 | Nc    | 119146348 | 118728768 | 0.9965  | 2487364808 | 20.88 | 0   | 22596 |
| 0.5  | Nc    | 119146348 | 78370392  | 0.6578  | 1156569560 | 9.71  | 0   | 22690 |
| 1    | Nc    | 119146348 | 13893050  | 0.1166  | 373927474  | 3.14  | 0   | 22620 |
| 2    | Nc    | 119146348 | 9240571   | 0.0776  | 316297152  | 2.65  | 0   | 22666 |
| 4    | Nc    | 119146348 | 6528463   | 0.0548  | 277049848  | 2.33  | 0   | 22679 |
| 8    | Nc    | 119146348 | 3868287   | 0.0325  | 219278076  | 1.84  | 0   | 22798 |
| 16   | Nc    | 119146348 | 2629096   | 0.0221  | 196224459  | 1.65  | 0   | 22623 |
| 32   | Nc    | 119146348 | 1661524   | 0.0139  | 177713329  | 1.49  | 0   | 22504 |
| 64   | Nc    | 119146348 | 1081156   | 0.0091  | 167328831  | 1.40  | 0   | 22098 |

| Fold | chrom | chrLength | covLength | covRate | bases    | mean   | min | max |
|:-----|:------|:----------|:----------|:--------|:---------|:-------|:----|:----|
| 0    | Mt    | 366924    | 363783    | 0.9914  | 53132850 | 144.81 | 0   | 361 |
| 0.25 | Mt    | 366924    | 363652    | 0.9911  | 53074874 | 144.65 | 0   | 360 |
| 0.5  | Mt    | 366924    | 363794    | 0.9915  | 53015767 | 144.49 | 0   | 360 |
| 1    | Mt    | 366924    | 363740    | 0.9913  | 53063245 | 144.62 | 0   | 360 |
| 2    | Mt    | 366924    | 363421    | 0.9905  | 53022871 | 144.51 | 0   | 360 |
| 4    | Mt    | 366924    | 356213    | 0.9708  | 49676263 | 135.39 | 0   | 360 |
| 8    | Mt    | 366924    | 95163     | 0.2594  | 6738481  | 18.36  | 0   | 306 |
| 16   | Mt    | 366924    | 28210     | 0.0769  | 630212   | 1.72   | 0   | 256 |
| 32   | Mt    | 366924    | 25910     | 0.0706  | 544902   | 1.49   | 0   | 250 |
| 64   | Mt    | 366924    | 17448     | 0.0476  | 373665   | 1.02   | 0   | 256 |

| Fold | chrom | chrLength | covLength | covRate | bases     | mean    | min | max  |
|:-----|:------|:----------|:----------|:--------|:----------|:--------|:----|:-----|
| 0    | Pt    | 154478    | 154478    | 1.0000  | 390524494 | 2528.03 | 4   | 3444 |
| 0.25 | Pt    | 154478    | 154478    | 1.0000  | 389956792 | 2524.35 | 4   | 3442 |
| 0.5  | Pt    | 154478    | 154478    | 1.0000  | 389862800 | 2523.74 | 4   | 3442 |
| 1    | Pt    | 154478    | 154478    | 1.0000  | 389924494 | 2524.14 | 4   | 3442 |
| 2    | Pt    | 154478    | 154478    | 1.0000  | 389867227 | 2523.77 | 4   | 3442 |
| 4    | Pt    | 154478    | 154478    | 1.0000  | 389838717 | 2523.59 | 4   | 3442 |
| 8    | Pt    | 154478    | 154478    | 1.0000  | 389852569 | 2523.68 | 4   | 3442 |
| 16   | Pt    | 154478    | 154478    | 1.0000  | 389918324 | 2524.10 | 4   | 3442 |
| 32   | Pt    | 154478    | 154478    | 1.0000  | 389773309 | 2523.16 | 2   | 3442 |
| 64   | Pt    | 154478    | 148548    | 0.9616  | 258998421 | 1676.60 | 0   | 3437 |

Table: Col-0 SRR611086 Folds

| Fold | chrom | chrLength | covLength | covRate | bases      | mean  | min | max   |
|:-----|:------|:----------|:----------|:--------|:-----------|:------|:----|:------|
| 0    | Nc    | 119146348 | 118911572 | 0.9980  | 5014032993 | 42.08 | 0   | 46514 |
| 0.25 | Nc    | 119146348 | 118742272 | 0.9966  | 4991705640 | 41.90 | 0   | 46534 |
| 0.5  | Nc    | 119146348 | 72377756  | 0.6075  | 1909924348 | 16.03 | 0   | 46478 |
| 1    | Nc    | 119146348 | 15877379  | 0.1333  | 716200662  | 6.01  | 0   | 46530 |
| 2    | Nc    | 119146348 | 11158262  | 0.0937  | 602522244  | 5.06  | 0   | 46337 |
| 4    | Nc    | 119146348 | 8351847   | 0.0701  | 526308824  | 4.42  | 0   | 46473 |
| 8    | Nc    | 119146348 | 5431558   | 0.0456  | 403060220  | 3.38  | 0   | 46467 |
| 16   | Nc    | 119146348 | 4068756   | 0.0341  | 357141776  | 3.00  | 0   | 46295 |
| 32   | Nc    | 119146348 | 3068787   | 0.0258  | 319157424  | 2.68  | 0   | 45961 |
| 64   | Nc    | 119146348 | 1992580   | 0.0167  | 297862635  | 2.50  | 0   | 45354 |

| Fold | chrom | chrLength | covLength | covRate | bases     | mean   | min | max |
|:-----|:------|:----------|:----------|:--------|:----------|:-------|:----|:----|
| 0    | Mt    | 366924    | 363907    | 0.9918  | 110294430 | 300.59 | 0   | 729 |
| 0.25 | Mt    | 366924    | 364126    | 0.9924  | 110259316 | 300.50 | 0   | 728 |
| 0.5  | Mt    | 366924    | 364079    | 0.9922  | 110195760 | 300.32 | 0   | 728 |
| 1    | Mt    | 366924    | 363961    | 0.9919  | 110124043 | 300.13 | 0   | 728 |
| 2    | Mt    | 366924    | 363810    | 0.9915  | 110177032 | 300.27 | 0   | 728 |
| 4    | Mt    | 366924    | 362531    | 0.9880  | 106881513 | 291.29 | 0   | 728 |
| 8    | Mt    | 366924    | 122927    | 0.3350  | 14628914  | 39.87  | 0   | 611 |
| 16   | Mt    | 366924    | 62492     | 0.1703  | 1428313   | 3.89   | 0   | 616 |
| 32   | Mt    | 366924    | 57643     | 0.1571  | 1206441   | 3.29   | 0   | 609 |
| 64   | Mt    | 366924    | 41074     | 0.1119  | 979520    | 2.67   | 0   | 620 |

| Fold | chrom | chrLength | covLength | covRate | bases     | mean    | min | max  |
|:-----|:------|:----------|:----------|:--------|:----------|:--------|:----|:-----|
| 0    | Pt    | 154478    | 154478    | 1.0000  | 693275576 | 4487.86 | 2   | 5904 |
| 0.25 | Pt    | 154478    | 154478    | 1.0000  | 692237254 | 4481.14 | 2   | 5898 |
| 0.5  | Pt    | 154478    | 154478    | 1.0000  | 692046597 | 4479.90 | 2   | 5898 |
| 1    | Pt    | 154478    | 154478    | 1.0000  | 692189465 | 4480.83 | 2   | 5898 |
| 2    | Pt    | 154478    | 154478    | 1.0000  | 692117486 | 4480.36 | 2   | 5898 |
| 4    | Pt    | 154478    | 154478    | 1.0000  | 692006900 | 4479.65 | 2   | 5898 |
| 8    | Pt    | 154478    | 154478    | 1.0000  | 691899452 | 4478.95 | 2   | 5898 |
| 16   | Pt    | 154478    | 154478    | 1.0000  | 692042081 | 4479.87 | 2   | 5898 |
| 32   | Pt    | 154478    | 154431    | 0.9997  | 691313793 | 4475.16 | 0   | 5898 |
| 64   | Pt    | 154478    | 151268    | 0.9792  | 412296703 | 2668.97 | 0   | 5879 |

Table: Col-0 SRR5216995 Folds

| Fold | chrom | chrLength | covLength | covRate | bases      | mean  | min | max    |
|:-----|:------|:----------|:----------|:--------|:-----------|:------|:----|:-------|
| 0    | Nc    | 119146348 | 118940483 | 0.9983  | 7036461038 | 59.06 | 0   | 105509 |
| 0.25 | Nc    | 119146348 | 118925235 | 0.9981  | 6984132507 | 58.62 | 0   | 105324 |
| 0.5  | Nc    | 119146348 | 54619245  | 0.4584  | 2290944232 | 19.23 | 0   | 105580 |
| 1    | Nc    | 119146348 | 11854942  | 0.0995  | 1045474074 | 8.77  | 0   | 104938 |
| 2    | Nc    | 119146348 | 7844827   | 0.0658  | 910697502  | 7.64  | 0   | 104715 |
| 4    | Nc    | 119146348 | 5310899   | 0.0446  | 772500832  | 6.48  | 0   | 104335 |
| 8    | Nc    | 119146348 | 3231231   | 0.0271  | 680141555  | 5.71  | 0   | 103665 |
| 16   | Nc    | 119146348 | 2249146   | 0.0189  | 629499845  | 5.28  | 0   | 102273 |
| 32   | Nc    | 119146348 | 1403503   | 0.0118  | 581422777  | 4.88  | 0   | 99633  |
| 64   | Nc    | 119146348 | 623905    | 0.0052  | 537191586  | 4.51  | 0   | 92620  |

| Fold | chrom | chrLength | covLength | covRate | bases    | mean   | min | max |
|:-----|:------|:----------|:----------|:--------|:---------|:-------|:----|:----|
| 0    | Mt    | 366924    | 362889    | 0.9890  | 94469130 | 257.46 | 0   | 617 |
| 0.25 | Mt    | 366924    | 362398    | 0.9877  | 93651000 | 255.23 | 0   | 615 |
| 0.5  | Mt    | 366924    | 362644    | 0.9883  | 93792065 | 255.62 | 0   | 615 |
| 1    | Mt    | 366924    | 362508    | 0.9880  | 93720659 | 255.42 | 0   | 615 |
| 2    | Mt    | 366924    | 362521    | 0.9880  | 93529324 | 254.90 | 0   | 615 |
| 4    | Mt    | 366924    | 156529    | 0.4266  | 19678194 | 53.63  | 0   | 587 |
| 8    | Mt    | 366924    | 23498     | 0.0640  | 457700   | 1.25   | 0   | 592 |
| 16   | Mt    | 366924    | 18081     | 0.0493  | 434297   | 1.18   | 0   | 577 |
| 32   | Mt    | 366924    | 19262     | 0.0525  | 430357   | 1.17   | 0   | 575 |
| 64   | Mt    | 366924    | 7556      | 0.0206  | 296209   | 0.81   | 0   | 582 |

| Fold | chrom | chrLength | covLength | covRate | bases     | mean    | min  | max  |
|:-----|:------|:----------|:----------|:--------|:----------|:--------|:-----|:-----|
| 0    | Pt    | 154478    | 154478    | 1.0000  | 895029108 | 5793.89 | 1368 | 8188 |
| 0.25 | Pt    | 154478    | 154478    | 1.0000  | 890591936 | 5765.17 | 1357 | 8168 |
| 0.5  | Pt    | 154478    | 154478    | 1.0000  | 889908392 | 5760.75 | 1301 | 8166 |
| 1    | Pt    | 154478    | 154478    | 1.0000  | 889556246 | 5758.47 | 1298 | 8166 |
| 2    | Pt    | 154478    | 154478    | 1.0000  | 889454944 | 5757.81 | 1330 | 8166 |
| 4    | Pt    | 154478    | 154478    | 1.0000  | 889391520 | 5757.40 | 1347 | 8166 |
| 8    | Pt    | 154478    | 154478    | 1.0000  | 889366585 | 5757.24 | 1346 | 8166 |
| 16   | Pt    | 154478    | 154478    | 1.0000  | 889340183 | 5757.07 | 1324 | 8166 |
| 32   | Pt    | 154478    | 154478    | 1.0000  | 888571719 | 5752.09 | 1213 | 8166 |
| 64   | Pt    | 154478    | 83678     | 0.5417  | 286346173 | 1853.64 | 0    | 7844 |

Table: Ler-0 SRR616965 Folds

| Fold | chrom | chrLength | covLength | covRate | bases      | mean  | min | max   |
|:-----|:------|:----------|:----------|:--------|:-----------|:------|:----|:------|
| 0    | Nc    | 119146348 | 112843627 | 0.9471  | 2393195441 | 20.09 | 0   | 23605 |
| 0.25 | Nc    | 119146348 | 112773888 | 0.9465  | 2380158986 | 19.98 | 0   | 23694 |
| 0.5  | Nc    | 119146348 | 50167841  | 0.4211  | 817795168  | 6.86  | 0   | 23632 |
| 1    | Nc    | 119146348 | 11881091  | 0.0997  | 383657385  | 3.22  | 0   | 23653 |
| 2    | Nc    | 119146348 | 8406746   | 0.0706  | 328370381  | 2.76  | 0   | 23656 |
| 4    | Nc    | 119146348 | 6090497   | 0.0511  | 290987442  | 2.44  | 0   | 23694 |
| 8    | Nc    | 119146348 | 4449987   | 0.0373  | 260503330  | 2.19  | 0   | 23550 |
| 16   | Nc    | 119146348 | 3065811   | 0.0257  | 182678949  | 1.53  | 0   | 23378 |
| 32   | Nc    | 119146348 | 2325855   | 0.0195  | 160846975  | 1.35  | 0   | 23116 |
| 64   | Nc    | 119146348 | 1876951   | 0.0158  | 152689486  | 1.28  | 0   | 22457 |

| Fold | chrom | chrLength | covLength | covRate | bases     | mean   | min | max |
|:-----|:------|:----------|:----------|:--------|:----------|:-------|:----|:----|
| 0    | Mt    | 366924    | 351472    | 0.9579  | 109519085 | 298.48 | 0   | 822 |
| 0.25 | Mt    | 366924    | 351476    | 0.9579  | 109367942 | 298.07 | 0   | 822 |
| 0.5  | Mt    | 366924    | 351297    | 0.9574  | 109310192 | 297.91 | 0   | 820 |
| 1    | Mt    | 366924    | 350967    | 0.9565  | 109281709 | 297.83 | 0   | 820 |
| 2    | Mt    | 366924    | 350815    | 0.9561  | 109308035 | 297.90 | 0   | 820 |
| 4    | Mt    | 366924    | 350948    | 0.9565  | 109197324 | 297.60 | 0   | 820 |
| 8    | Mt    | 366924    | 348567    | 0.9500  | 103864413 | 283.07 | 0   | 820 |
| 16   | Mt    | 366924    | 114147    | 0.3111  | 8703600   | 23.72  | 0   | 622 |
| 32   | Mt    | 366924    | 79795     | 0.2175  | 1378446   | 3.76   | 0   | 619 |
| 64   | Mt    | 366924    | 76285     | 0.2079  | 1374260   | 3.75   | 0   | 616 |

| Fold | chrom | chrLength | covLength | covRate | bases     | mean    | min | max  |
|:-----|:------|:----------|:----------|:--------|:----------|:--------|:----|:-----|
| 0    | Pt    | 154478    | 154478    | 1.0000  | 767923755 | 4971.09 | 244 | 6898 |
| 0.25 | Pt    | 154478    | 154478    | 1.0000  | 767928954 | 4971.12 | 244 | 6927 |
| 0.5  | Pt    | 154478    | 154478    | 1.0000  | 767914496 | 4971.03 | 244 | 6968 |
| 1    | Pt    | 154478    | 154478    | 1.0000  | 767908001 | 4970.99 | 244 | 7018 |
| 2    | Pt    | 154478    | 154478    | 1.0000  | 767817464 | 4970.40 | 244 | 6957 |
| 4    | Pt    | 154478    | 154478    | 1.0000  | 767871785 | 4970.75 | 244 | 6949 |
| 8    | Pt    | 154478    | 154478    | 1.0000  | 767821009 | 4970.42 | 243 | 6945 |
| 16   | Pt    | 154478    | 154478    | 1.0000  | 767834842 | 4970.51 | 243 | 6970 |
| 32   | Pt    | 154478    | 154478    | 1.0000  | 767808163 | 4970.34 | 243 | 6924 |
| 64   | Pt    | 154478    | 154478    | 1.0000  | 767857310 | 4970.66 | 177 | 6968 |

Table: Ler-0 SRR611087 Folds

| Fold | chrom | chrLength | covLength | covRate | bases      | mean  | min | max   |
|:-----|:------|:----------|:----------|:--------|:-----------|:------|:----|:------|
| 0    | Nc    | 119146348 | 113259476 | 0.9506  | 4673921200 | 39.23 | 0   | 39590 |
| 0.25 | Nc    | 119146348 | 113146168 | 0.9496  | 4648634221 | 39.02 | 0   | 39520 |
| 0.5  | Nc    | 119146348 | 52628352  | 0.4417  | 1533820866 | 12.87 | 0   | 39626 |
| 1    | Nc    | 119146348 | 14001234  | 0.1175  | 723111708  | 6.07  | 0   | 39631 |
| 2    | Nc    | 119146348 | 10346152  | 0.0868  | 611879662  | 5.14  | 0   | 39715 |
| 4    | Nc    | 119146348 | 7900444   | 0.0663  | 536577336  | 4.50  | 0   | 39595 |
| 8    | Nc    | 119146348 | 6113437   | 0.0513  | 477635606  | 4.01  | 0   | 39420 |
| 16   | Nc    | 119146348 | 4507030   | 0.0378  | 320100221  | 2.69  | 0   | 39164 |
| 32   | Nc    | 119146348 | 3748466   | 0.0315  | 275464674  | 2.31  | 0   | 38765 |
| 64   | Nc    | 119146348 | 3197041   | 0.0268  | 257588137  | 2.16  | 0   | 37975 |

| Fold | chrom | chrLength | covLength | covRate | bases     | mean   | min | max  |
|:-----|:------|:----------|:----------|:--------|:----------|:-------|:----|:-----|
| 0    | Mt    | 366924    | 352249    | 0.9600  | 219983375 | 599.53 | 0   | 1653 |
| 0.25 | Mt    | 366924    | 351845    | 0.9589  | 220062812 | 599.75 | 0   | 1653 |
| 0.5  | Mt    | 366924    | 351508    | 0.9580  | 219961213 | 599.47 | 0   | 1653 |
| 1    | Mt    | 366924    | 351482    | 0.9579  | 219931622 | 599.39 | 0   | 1653 |
| 2    | Mt    | 366924    | 351366    | 0.9576  | 219947740 | 599.44 | 0   | 1653 |
| 4    | Mt    | 366924    | 351373    | 0.9576  | 219703022 | 598.77 | 0   | 1653 |
| 8    | Mt    | 366924    | 351158    | 0.9570  | 216422346 | 589.83 | 0   | 1653 |
| 16   | Mt    | 366924    | 168021    | 0.4579  | 19925031  | 54.30  | 0   | 1568 |
| 32   | Mt    | 366924    | 131237    | 0.3577  | 3010889   | 8.21   | 0   | 1248 |
| 64   | Mt    | 366924    | 128565    | 0.3504  | 2769396   | 7.55   | 0   | 1262 |

| Fold | chrom | chrLength | covLength | covRate | bases      | mean    | min | max  |
|:-----|:------|:----------|:----------|:--------|:-----------|:--------|:----|:-----|
| 0    | Pt    | 154478    | 154478    | 1.0000  | 1102940670 | 7139.79 | 20  | 9621 |
| 0.25 | Pt    | 154478    | 154478    | 1.0000  | 1103004796 | 7140.21 | 20  | 9580 |
| 0.5  | Pt    | 154478    | 154478    | 1.0000  | 1102908208 | 7139.58 | 20  | 9639 |
| 1    | Pt    | 154478    | 154478    | 1.0000  | 1103005489 | 7140.21 | 20  | 9561 |
| 2    | Pt    | 154478    | 154478    | 1.0000  | 1102906662 | 7139.57 | 20  | 9614 |
| 4    | Pt    | 154478    | 154478    | 1.0000  | 1102892044 | 7139.48 | 20  | 9607 |
| 8    | Pt    | 154478    | 154478    | 1.0000  | 1102687420 | 7138.15 | 20  | 9568 |
| 16   | Pt    | 154478    | 154478    | 1.0000  | 1102868223 | 7139.32 | 20  | 9569 |
| 32   | Pt    | 154478    | 154478    | 1.0000  | 1102811593 | 7138.96 | 20  | 9656 |
| 64   | Pt    | 154478    | 154445    | 0.9998  | 1102669960 | 7138.04 | 0   | 9611 |

Table: Nipponbare SRR545231 Folds

| Fold | chrom | chrLength | covLength | covRate | bases       | mean  | min | max  |
|:-----|:------|:----------|:----------|:--------|:------------|:------|:----|:-----|
| 0    | Nc    | 373245519 | 368402032 | 0.9870  | 12983797798 | 34.79 | 0   | 8186 |
| 0.25 | Nc    | 373245519 | 366583918 | 0.9822  | 12964540235 | 34.73 | 0   | 8170 |
| 0.5  | Nc    | 373245519 | 356144812 | 0.9542  | 12665609206 | 33.93 | 0   | 8173 |
| 1    | Nc    | 373245519 | 174199288 | 0.4667  | 4466878641  | 11.97 | 0   | 8192 |
| 2    | Nc    | 373245519 | 120942207 | 0.3240  | 3151703014  | 8.44  | 0   | 8172 |
| 4    | Nc    | 373245519 | 98773643  | 0.2646  | 2556632934  | 6.85  | 0   | 8140 |
| 8    | Nc    | 373245519 | 82490186  | 0.2210  | 2129813140  | 5.71  | 0   | 8140 |
| 16   | Nc    | 373245519 | 67192387  | 0.1800  | 1740850801  | 4.66  | 0   | 7965 |
| 32   | Nc    | 373245519 | 50568041  | 0.1355  | 1272065380  | 3.41  | 0   | 7797 |
| 64   | Nc    | 373245519 | 33468069  | 0.0897  | 798576696   | 2.14  | 0   | 7212 |

| Fold | chrom | chrLength | covLength | covRate | bases    | mean   | min | max |
|:-----|:------|:----------|:----------|:--------|:---------|:-------|:----|:----|
| 0    | Mt    | 490520    | 488459    | 0.9958  | 62088079 | 126.58 | 0   | 465 |
| 0.25 | Mt    | 490520    | 488638    | 0.9962  | 62070637 | 126.54 | 0   | 530 |
| 0.5  | Mt    | 490520    | 488594    | 0.9961  | 62057759 | 126.51 | 0   | 489 |
| 1    | Mt    | 490520    | 488642    | 0.9962  | 62097034 | 126.59 | 0   | 466 |
| 2    | Mt    | 490520    | 488581    | 0.9960  | 61976607 | 126.35 | 0   | 484 |
| 4    | Mt    | 490520    | 432170    | 0.8810  | 46324406 | 94.44  | 0   | 501 |
| 8    | Mt    | 490520    | 88134     | 0.1797  | 7604591  | 15.50  | 0   | 458 |
| 16   | Mt    | 490520    | 43163     | 0.0880  | 4466508  | 9.11   | 0   | 498 |
| 32   | Mt    | 490520    | 33508     | 0.0683  | 3760523  | 7.67   | 0   | 447 |
| 64   | Mt    | 490520    | 6881      | 0.0140  | 486583   | 0.99   | 0   | 405 |

| Fold | chrom | chrLength | covLength | covRate | bases    | mean   | min | max  |
|:-----|:------|:----------|:----------|:--------|:---------|:-------|:----|:-----|
| 0    | Pt    | 134525    | 132943    | 0.9882  | 56062147 | 416.74 | 0   | 1746 |
| 0.25 | Pt    | 134525    | 132903    | 0.9879  | 56268862 | 418.28 | 0   | 1746 |
| 0.5  | Pt    | 134525    | 132771    | 0.9870  | 56165888 | 417.51 | 0   | 1746 |
| 1    | Pt    | 134525    | 132892    | 0.9879  | 56113124 | 417.12 | 0   | 1746 |
| 2    | Pt    | 134525    | 132957    | 0.9883  | 56090073 | 416.95 | 0   | 1746 |
| 4    | Pt    | 134525    | 132841    | 0.9875  | 56059193 | 416.72 | 0   | 1746 |
| 8    | Pt    | 134525    | 132867    | 0.9877  | 56121372 | 417.18 | 0   | 1746 |
| 16   | Pt    | 134525    | 132857    | 0.9876  | 55970896 | 416.06 | 0   | 1746 |
| 32   | Pt    | 134525    | 115373    | 0.8576  | 29714063 | 220.88 | 0   | 1743 |
| 64   | Pt    | 134525    | 18374     | 0.1366  | 2141712  | 15.92  | 0   | 450  |

Table: NP SRR063638 Folds

| Fold | chrom | chrLength | covLength | covRate | bases      | mean | min | max  |
|:-----|:------|:----------|:----------|:--------|:-----------|:-----|:----|:-----|
| 0    | Nc    | 373245519 | 328139822 | 0.8792  | 2902528518 | 7.78 | 0   | 3467 |
| 0.25 | Nc    | 373245519 | 326604651 | 0.8750  | 2898173424 | 7.76 | 0   | 3446 |
| 0.5  | Nc    | 373245519 | 310867716 | 0.8329  | 2792272341 | 7.48 | 0   | 3466 |
| 1    | Nc    | 373245519 | 174870691 | 0.4685  | 1350011316 | 3.62 | 0   | 3421 |
| 2    | Nc    | 373245519 | 97467718  | 0.2611  | 769361428  | 2.06 | 0   | 3452 |
| 4    | Nc    | 373245519 | 78270101  | 0.2097  | 644317110  | 1.73 | 0   | 3477 |
| 8    | Nc    | 373245519 | 64422539  | 0.1726  | 563854945  | 1.51 | 0   | 3494 |
| 16   | Nc    | 373245519 | 50283322  | 0.1347  | 476517152  | 1.28 | 0   | 3475 |
| 32   | Nc    | 373245519 | 36156052  | 0.0969  | 393014785  | 1.05 | 0   | 3462 |
| 64   | Nc    | 373245519 | 23672609  | 0.0634  | 324573449  | 0.87 | 0   | 3451 |

| Fold | chrom | chrLength | covLength | covRate | bases    | mean  | min | max |
|:-----|:------|:----------|:----------|:--------|:---------|:------|:----|:----|
| 0    | Mt    | 490520    | 488175    | 0.9952  | 46766188 | 95.34 | 0   | 665 |
| 0.25 | Mt    | 490520    | 488206    | 0.9953  | 46725421 | 95.26 | 0   | 646 |
| 0.5  | Mt    | 490520    | 488220    | 0.9953  | 46702662 | 95.21 | 0   | 695 |
| 1    | Mt    | 490520    | 488115    | 0.9951  | 46594712 | 94.99 | 0   | 645 |
| 2    | Mt    | 490520    | 488166    | 0.9952  | 46678496 | 95.16 | 0   | 690 |
| 4    | Mt    | 490520    | 487797    | 0.9944  | 46668284 | 95.14 | 0   | 647 |
| 8    | Mt    | 490520    | 474402    | 0.9671  | 44866194 | 91.47 | 0   | 669 |
| 16   | Mt    | 490520    | 247375    | 0.5043  | 17454207 | 35.58 | 0   | 650 |
| 32   | Mt    | 490520    | 82202     | 0.1676  | 5115594  | 10.43 | 0   | 640 |
| 64   | Mt    | 490520    | 61990     | 0.1264  | 4177739  | 8.52  | 0   | 673 |

| Fold | chrom | chrLength | covLength | covRate | bases    | mean   | min | max  |
|:-----|:------|:----------|:----------|:--------|:---------|:-------|:----|:-----|
| 0    | Pt    | 134525    | 133283    | 0.9908  | 68261865 | 507.43 | 0   | 2615 |
| 0.25 | Pt    | 134525    | 133119    | 0.9895  | 68409621 | 508.53 | 0   | 2617 |
| 0.5  | Pt    | 134525    | 133196    | 0.9901  | 68233413 | 507.22 | 0   | 2616 |
| 1    | Pt    | 134525    | 133116    | 0.9895  | 68022282 | 505.65 | 0   | 2616 |
| 2    | Pt    | 134525    | 133340    | 0.9912  | 68225954 | 507.16 | 0   | 2616 |
| 4    | Pt    | 134525    | 133057    | 0.9891  | 68210977 | 507.05 | 0   | 2615 |
| 8    | Pt    | 134525    | 133107    | 0.9895  | 68090615 | 506.16 | 0   | 2617 |
| 16   | Pt    | 134525    | 133226    | 0.9903  | 68244550 | 507.30 | 0   | 2616 |
| 32   | Pt    | 134525    | 133094    | 0.9894  | 68281992 | 507.58 | 0   | 2618 |
| 64   | Pt    | 134525    | 133204    | 0.9902  | 68314292 | 507.82 | 0   | 2616 |

Table: A17 SRR1542423 Folds

| Fold | chrom | chrLength | covLength | covRate | bases      | mean | min | max    |
|:-----|:------|:----------|:----------|:--------|:-----------|:-----|:----|:-------|
| 0    | Nc    | 384466993 | 348072517 | 0.9053  | 3237657237 | 8.42 | 0   | 107159 |
| 0.25 | Nc    | 384466993 | 319243960 | 0.8304  | 3077659704 | 8.01 | 0   | 92341  |
| 0.5  | Nc    | 384466993 | 197683490 | 0.5142  | 2210392114 | 5.75 | 0   | 92229  |
| 1    | Nc    | 384466993 | 76576135  | 0.1992  | 1292591729 | 3.36 | 0   | 91636  |
| 2    | Nc    | 384466993 | 51550888  | 0.1341  | 1079035637 | 2.81 | 0   | 108024 |
| 4    | Nc    | 384466993 | 39755248  | 0.1034  | 942208106  | 2.45 | 0   | 89569  |
| 8    | Nc    | 384466993 | 31023999  | 0.0807  | 848499067  | 2.21 | 0   | 108386 |
| 16   | Nc    | 384466993 | 23226239  | 0.0604  | 720018690  | 1.87 | 0   | 87693  |
| 32   | Nc    | 384466993 | 16235103  | 0.0422  | 618721163  | 1.61 | 0   | 103831 |
| 64   | Nc    | 384466993 | 11134978  | 0.0290  | 507222214  | 1.32 | 0   | 105054 |

| Fold | chrom | chrLength | covLength | covRate | bases    | mean   | min | max |
|:-----|:------|:----------|:----------|:--------|:---------|:-------|:----|:----|
| 0    | Mt    | 271618    | 271618    | 1.0000  | 49932049 | 183.83 | 1   | 726 |
| 0.25 | Mt    | 271618    | 271618    | 1.0000  | 49142557 | 180.93 | 1   | 721 |
| 0.5  | Mt    | 271618    | 271618    | 1.0000  | 49020706 | 180.48 | 1   | 721 |
| 1    | Mt    | 271618    | 271618    | 1.0000  | 48878964 | 179.95 | 1   | 720 |
| 2    | Mt    | 271618    | 271618    | 1.0000  | 48861022 | 179.89 | 1   | 718 |
| 4    | Mt    | 271618    | 271495    | 0.9995  | 48383733 | 178.13 | 0   | 717 |
| 8    | Mt    | 271618    | 253902    | 0.9348  | 38889573 | 143.18 | 0   | 708 |
| 16   | Mt    | 271618    | 95498     | 0.3516  | 9127266  | 33.60  | 0   | 706 |
| 32   | Mt    | 271618    | 21853     | 0.0805  | 578709   | 2.13   | 0   | 571 |
| 64   | Mt    | 271618    | 14959     | 0.0551  | 40418    | 0.15   | 0   | 127 |

| Fold | chrom | chrLength | covLength | covRate | bases    | mean   | min | max  |
|:-----|:------|:----------|:----------|:--------|:---------|:-------|:----|:-----|
| 0    | Pt    | 124033    | 124031    | 1.0000  | 80757838 | 651.10 | 0   | 3314 |
| 0.25 | Pt    | 124033    | 124032    | 1.0000  | 80570481 | 649.59 | 0   | 3313 |
| 0.5  | Pt    | 124033    | 124029    | 1.0000  | 80085169 | 645.68 | 0   | 3292 |
| 1    | Pt    | 124033    | 124032    | 1.0000  | 79931815 | 644.44 | 0   | 3285 |
| 2    | Pt    | 124033    | 124032    | 1.0000  | 79714333 | 642.69 | 0   | 3285 |
| 4    | Pt    | 124033    | 124024    | 0.9999  | 79818760 | 643.53 | 0   | 3281 |
| 8    | Pt    | 124033    | 123963    | 0.9994  | 79719423 | 642.73 | 0   | 3280 |
| 16   | Pt    | 124033    | 123771    | 0.9979  | 79182159 | 638.40 | 0   | 3281 |
| 32   | Pt    | 124033    | 117909    | 0.9506  | 75454137 | 608.34 | 0   | 3284 |
| 64   | Pt    | 124033    | 82647     | 0.6663  | 54141806 | 436.51 | 0   | 3278 |

Table: Heinz1706 SRR1572628 Folds

| Fold | chrom | chrLength | covLength | covRate | bases      | mean | min | max   |
|:-----|:------|:----------|:----------|:--------|:-----------|:-----|:----|:------|
| 0    | Nc    | 807224664 | 687361649 | 0.8515  | 3180562052 | 3.94 | 0   | 20886 |
| 0.25 | Nc    | 807224664 | 687367973 | 0.8515  | 3180568663 | 3.94 | 0   | 20936 |
| 0.5  | Nc    | 807224664 | 678841856 | 0.8410  | 3165122034 | 3.92 | 0   | 20924 |
| 1    | Nc    | 807224664 | 512507255 | 0.6349  | 2438099474 | 3.02 | 0   | 21022 |
| 2    | Nc    | 807224664 | 204747604 | 0.2536  | 1110973726 | 1.38 | 0   | 21151 |
| 4    | Nc    | 807224664 | 137292726 | 0.1701  | 819204674  | 1.01 | 0   | 20919 |
| 8    | Nc    | 807224664 | 111626374 | 0.1383  | 699778951  | 0.87 | 0   | 20716 |
| 16   | Nc    | 807224664 | 92796197  | 0.1150  | 605505425  | 0.75 | 0   | 20469 |
| 32   | Nc    | 807224664 | 75908434  | 0.0940  | 515050993  | 0.64 | 0   | 20284 |
| 64   | Nc    | 807224664 | 61797609  | 0.0766  | 436774399  | 0.54 | 0   | 20279 |

| Fold | chrom | chrLength | covLength | covRate | bases    | mean  | min | max |
|:-----|:------|:----------|:----------|:--------|:---------|:------|:----|:----|
| 0    | Mt    | 446257    | 446257    | 1.0000  | 43828430 | 98.21 | 2   | 828 |
| 0.25 | Mt    | 446257    | 446257    | 1.0000  | 43807272 | 98.17 | 1   | 799 |
| 0.5  | Mt    | 446257    | 446255    | 1.0000  | 43837766 | 98.23 | 0   | 858 |
| 1    | Mt    | 446257    | 446257    | 1.0000  | 43778626 | 98.10 | 2   | 846 |
| 2    | Mt    | 446257    | 446254    | 1.0000  | 43750446 | 98.04 | 0   | 902 |
| 4    | Mt    | 446257    | 446257    | 1.0000  | 43804934 | 98.16 | 1   | 875 |
| 8    | Mt    | 446257    | 446257    | 1.0000  | 43789320 | 98.13 | 3   | 846 |
| 16   | Mt    | 446257    | 436794    | 0.9788  | 39605474 | 88.75 | 0   | 861 |
| 32   | Mt    | 446257    | 209356    | 0.4691  | 13305819 | 29.82 | 0   | 840 |
| 64   | Mt    | 446257    | 31133     | 0.0698  | 1332256  | 2.99  | 0   | 863 |

| Fold | chrom | chrLength | covLength | covRate | bases     | mean   | min | max  |
|:-----|:------|:----------|:----------|:--------|:----------|:-------|:----|:-----|
| 0    | Pt    | 155461    | 155461    | 1.0000  | 132760191 | 853.98 | 3   | 2245 |
| 0.25 | Pt    | 155461    | 155461    | 1.0000  | 132746965 | 853.89 | 3   | 2245 |
| 0.5  | Pt    | 155461    | 155461    | 1.0000  | 132756549 | 853.95 | 5   | 2245 |
| 1    | Pt    | 155461    | 155461    | 1.0000  | 132674798 | 853.43 | 2   | 2243 |
| 2    | Pt    | 155461    | 155461    | 1.0000  | 132594962 | 852.91 | 4   | 2241 |
| 4    | Pt    | 155461    | 155461    | 1.0000  | 132593950 | 852.91 | 9   | 2241 |
| 8    | Pt    | 155461    | 155461    | 1.0000  | 132612435 | 853.03 | 8   | 2241 |
| 16   | Pt    | 155461    | 155461    | 1.0000  | 132611693 | 853.02 | 4   | 2241 |
| 32   | Pt    | 155461    | 155461    | 1.0000  | 132639793 | 853.20 | 5   | 2241 |
| 64   | Pt    | 155461    | 155461    | 1.0000  | 132292659 | 850.97 | 6   | 2241 |

```shell
cd ~/data/plastid/evaluation

SRRS=(
    'SRR616966::Col-0'      # Col-0
    'SRR611086::Col-0'
    'SRR5216995::Col-0'
    'SRR616965::Ler-0'      # Ler-0
    'SRR611087::Ler-0'
    'SRR545231::Nipponbare' # Nipponbare
    'SRR063638::NP'
    'SRR1542423::A17'       # A17
    'SRR1572628::Heinz1706' # Heinz 1706
)
FOLDS=(0 0.25 0.5 1 2 4 8 16 32 64)

for item in "${SRRS[@]}"; do
    SRR="${item%%::*}"
    STRAIN="${item##*::}"

    for FOLD in "${FOLDS[@]}"; do
        BASE_NAME=${SRR}_${FOLD}

        cat ${BASE_NAME}/statReads.md |
            mdtable2csv |
            mlr --icsv --otsv cat |
            grep -v "^Name" |
            grep -v "^Genome" |
            tsv-select -f 1,3 |
            (echo -e "Fold\t${FOLD}" && cat) |
            datamash transpose
    done |
        tsv-uniq \
        > ${SRR}_reads.tsv

    echo
    echo "Table: ${STRAIN} ${SRR} Reads"
    echo
    cat ${SRR}_reads.tsv |
        mlr --itsv --omd cat
    echo

done

```

Table: Col-0 SRR616966 Reads

| Fold | Illumina.R | trim.R  | Q25L60  |
|:-----|:-----------|:--------|:--------|
| 0    | 4.97G      | 4.52G   | 4.17G   |
| 0.25 | 4.97G      | 3.67G   | 3.41G   |
| 0.5  | 4.97G      | 2.13G   | 1.99G   |
| 1    | 4.97G      | 1.26G   | 1.17G   |
| 2    | 4.97G      | 1.19G   | 1.1G    |
| 4    | 4.97G      | 1.14G   | 1.06G   |
| 8    | 4.97G      | 1.02G   | 945.53M |
| 16   | 4.97G      | 981.34M | 912.91M |
| 32   | 4.97G      | 957.19M | 890.32M |
| 64   | 4.97G      | 761.94M | 710.47M |

Table: Col-0 SRR611086 Reads

| Fold | Illumina.R | trim.R | Q25L60 |
|:-----|:-----------|:-------|:-------|
| 0    | 9.98G      | 9.04G  | 8.4G   |
| 0.25 | 9.98G      | 7.24G  | 6.78G  |
| 0.5  | 9.98G      | 3.74G  | 3.51G  |
| 1    | 9.98G      | 2.42G  | 2.26G  |
| 2    | 9.98G      | 2.28G  | 2.13G  |
| 4    | 9.98G      | 2.19G  | 2.04G  |
| 8    | 9.98G      | 1.93G  | 1.81G  |
| 16   | 9.98G      | 1.86G  | 1.74G  |
| 32   | 9.98G      | 1.81G  | 1.69G  |
| 64   | 9.98G      | 1.39G  | 1.3G   |

Table: Col-0 SRR5216995 Reads

| Fold | Illumina.R | trim.R | Q25L60 |
|:-----|:-----------|:-------|:-------|
| 0    | 15.53G     | 13.42G | 12.14G |
| 0.25 | 15.53G     | 12.92G | 11.82G |
| 0.5  | 15.53G     | 5.44G  | 4.96G  |
| 1    | 15.53G     | 3.51G  | 3.15G  |
| 2    | 15.53G     | 3.3G   | 2.95G  |
| 4    | 15.53G     | 2.96G  | 2.65G  |
| 8    | 15.53G     | 2.78G  | 2.49G  |
| 16   | 15.53G     | 2.69G  | 2.41G  |
| 32   | 15.53G     | 2.59G  | 2.33G  |
| 64   | 15.53G     | 1.53G  | 1.35G  |

Table: Ler-0 SRR616965 Reads

| Fold | Illumina.R | trim.R | Q25L60 |
|:-----|:-----------|:-------|:-------|
| 0    | 5.09G      | 4.28G  | 4.06G  |
| 0.25 | 5.09G      | 4.21G  | 3.99G  |
| 0.5  | 5.09G      | 2.44G  | 2.31G  |
| 1    | 5.09G      | 1.96G  | 1.85G  |
| 2    | 5.09G      | 1.89G  | 1.79G  |
| 4    | 5.09G      | 1.85G  | 1.74G  |
| 8    | 5.09G      | 1.81G  | 1.71G  |
| 16   | 5.09G      | 1.61G  | 1.52G  |
| 32   | 5.09G      | 1.57G  | 1.48G  |
| 64   | 5.09G      | 1.56G  | 1.47G  |

Table: Ler-0 SRR611087 Reads

| Fold | Illumina.R | trim.R | Q25L60 |
|:-----|:-----------|:-------|:-------|
| 0    | 10.16G     | 8.15G  | 7.7G   |
| 0.25 | 10.16G     | 7.99G  | 7.55G  |
| 0.5  | 10.16G     | 4.44G  | 4.19G  |
| 1    | 10.16G     | 3.53G  | 3.32G  |
| 2    | 10.16G     | 3.4G   | 3.2G   |
| 4    | 10.16G     | 3.31G  | 3.11G  |
| 8    | 10.16G     | 3.24G  | 3.04G  |
| 16   | 10.16G     | 2.82G  | 2.65G  |
| 32   | 10.16G     | 2.74G  | 2.58G  |
| 64   | 10.16G     | 2.71G  | 2.55G  |

Table: Nipponbare SRR545231 Reads

| Fold | Illumina.R | trim.R | Q25L60  |
|:-----|:-----------|:-------|:--------|
| 0    | 17.22G     | 14.56G | 13.72G  |
| 0.25 | 17.22G     | 14.52G | 13.69G  |
| 0.5  | 17.22G     | 14.12G | 13.36G  |
| 1    | 17.22G     | 5.1G   | 4.81G   |
| 2    | 17.22G     | 3.66G  | 3.45G   |
| 4    | 17.22G     | 2.99G  | 2.81G   |
| 8    | 17.22G     | 2.48G  | 2.33G   |
| 16   | 17.22G     | 2.04G  | 1.92G   |
| 32   | 17.22G     | 1.49G  | 1.4G    |
| 64   | 17.22G     | 931.8M | 872.28M |

Table: NP SRR063638 Reads

| Fold | Illumina.R | trim.R  | Q25L60  |
|:-----|:-----------|:--------|:--------|
| 0    | 5.96G      | 3.73G   | 3.33G   |
| 0.25 | 5.96G      | 3.68G   | 3.28G   |
| 0.5  | 5.96G      | 3.53G   | 3.16G   |
| 1    | 5.96G      | 1.78G   | 1.59G   |
| 2    | 5.96G      | 1.1G    | 975.86M |
| 4    | 5.96G      | 946.63M | 840.97M |
| 8    | 5.96G      | 845.48M | 751.19M |
| 16   | 5.96G      | 703.88M | 626.12M |
| 32   | 5.96G      | 586.72M | 522M    |
| 64   | 5.96G      | 501.67M | 446.5M  |

Table: A17 SRR1542423 Reads

| Fold | Illumina.R | trim.R | Q25L60 |
|:-----|:-----------|:-------|:-------|
| 0    | 8.96G      | 5.61G  | 4.97G  |
| 0.25 | 8.96G      | 5.33G  | 4.74G  |
| 0.5  | 8.96G      | 4.15G  | 3.69G  |
| 1    | 8.96G      | 2.91G  | 2.58G  |
| 2    | 8.96G      | 2.6G   | 2.3G   |
| 4    | 8.96G      | 2.43G  | 2.15G  |
| 8    | 8.96G      | 2.27G  | 2.01G  |
| 16   | 8.96G      | 2.07G  | 1.83G  |
| 32   | 8.96G      | 1.9G   | 1.68G  |
| 64   | 8.96G      | 1.71G  | 1.51G  |

Table: Heinz1706 SRR1572628 Reads

| Fold | Illumina.R | trim.R  | Q25L60  |
|:-----|:-----------|:--------|:--------|
| 0    | 4.84G      | 4.13G   | 3.81G   |
| 0.25 | 4.84G      | 4.13G   | 3.81G   |
| 0.5  | 4.84G      | 4.1G    | 3.79G   |
| 1    | 4.84G      | 3.24G   | 3G      |
| 2    | 4.84G      | 1.67G   | 1.52G   |
| 4    | 4.84G      | 1.32G   | 1.2G    |
| 8    | 4.84G      | 1.17G   | 1.06G   |
| 16   | 4.84G      | 1.05G   | 945.98M |
| 32   | 4.84G      | 902.96M | 813.96M |
| 64   | 4.84G      | 790.91M | 711.65M |

## Remove intermediate files

```shell
cd ~/data/plastid/evaluation

find . -type d -name "trim" | xargs rm -fr
find . -type f -path "*3_bwa/genome.fa*" | xargs rm
find . -type f -name "*.ba[mi]" | xargs rm
find . -type f -name "*.per-base.bed.gz" | xargs rm

find . -type f -name "*.tadpole.contig.*" | xargs rm

find . -type f -name "core.*" | xargs rm
find . -type f -name "output.*" | xargs rm

```

