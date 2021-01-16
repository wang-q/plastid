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


## Reference genomes

* *Arabidopsis thaliana* Col-0

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

* *Oryza sativa* Japonica Group Cultivar Nipponbare

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

* *Medicago truncatula* A17

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

* *Prunus persica* PLov2-2N (a double haploid genotype of the peach cv. Lovell)

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

* *Solanum lycopersicum* Cultivar: Heinz 1706

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

* *Glycine max* Williams 82

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

## Download fastq files from ENA

```shell script
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

mlr --icsv --omd cat ena_info.csv

cat ena_info.ascp.sh |
    parallel --no-run-if-empty -j 1 "{}"

```

| name         | srx        | platform | layout | ilength | srr         | spots     | bases  |
|:-------------|:-----------|:---------|:-------|:--------|:------------|:----------|:-------|
| Atha_Col_0_1 | SRX202246  | ILLUMINA | PAIRED | 450     | SRR611086   | 49891349  | 9.29G  |
| Atha_Col_0_1 | SRX202246  | ILLUMINA | PAIRED | 450     | SRR616966   | 24851796  | 4.63G  |
| Atha_Col_0_2 | SRX2527206 | ILLUMINA | PAIRED |         | SRR5216995  | 26893065  | 14.46G |
| Atha_Ler_0   | SRX202247  | ILLUMINA | PAIRED | 450     | SRR611087   | 50791450  | 9.46G  |
| Atha_Ler_0   | SRX202247  | ILLUMINA | PAIRED | 450     | SRR616965   | 25436255  | 4.74G  |
| Gmax_W82     | SRX7009428 | ILLUMINA | PAIRED |         | SRR10296600 | 162110355 | 45.29G |
| Mtru_A17     | SRX673852  | ILLUMINA | PAIRED | 360     | SRR1542422  | 99418334  | 16.67G |
| Mtru_A17     | SRX673852  | ILLUMINA | PAIRED | 360     | SRR1542423  | 29663436  | 8.34G  |
| Osat_Nip     | SRX179254  | ILLUMINA | PAIRED | 300     | SRR545059   | 85148124  | 7.93G  |
| Osat_Nip     | SRX179254  | ILLUMINA | PAIRED | 300     | SRR545231   | 85251097  | 16.04G |
| Osat_Nip_2   | SRX025260  | ILLUMINA | PAIRED | 463     | SRR063638   | 29784011  | 5.55G  |
| Pper_Lovell  | SRX150254  | ILLUMINA | PAIRED | 400     | SRR502985   | 123590441 | 23.25G |
| Slyc_H1706   | SRX698770  | ILLUMINA | PAIRED |         | SRR1572628  | 24198345  | 4.51G  |


## 基本信息

`cutoff = 倍数因子 * 覆盖深度`

* 因子值 0, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64

* Install `mdtable2csv`

```shell script
mkdir -p ${HOME}/bin
curl -fsSL $(
    curl -fsSL https://api.github.com/repos/515hikaru/mdtable2csv/releases/latest |
        jq -r '.assets[] | select(.name == "mdtable2csv_linux_x86_64.tar.gz").browser_download_url'
    ) |
    tar xvz mdtable2csv
mv mdtable2csv ${HOME}/bin

```

## Symlink

```shell script
mkdir -p ~/data/plastid/evaluation
cd ~/data/plastid/evaluation

SRRS=(
    'SRR616966::col_0' # Col-0
#    'SRR611086::col_0'
#    'SRR5216995::col_0'
    'SRR616965::col_0' # Ler-0
#    'SRR611087::col_0'
#    'SRR545231::nip'   # Nipponbare
#    'SRR063638::nip'
#    'SRR1542423::a17'  # A17
#    'SRR1542422::a17'
)
FOLDS=(0 0.25 0.5 1 2 4 8 16 32 64)

for item in "${SRRS[@]}"; do
    SRR="${item%%::*}"
    STRAIN="${item##*::}"

    for FOLD in "${FOLDS[@]}"; do
        BASE_NAME=${SRR}_${FOLD}

        mkdir -p ${BASE_NAME}/1_genome
        pushd ${BASE_NAME}/1_genome

        ln -fs ../../../genome/${STRAIN}/genome.fa genome.fa
        cp ../../../genome/${STRAIN}/chr.sizes chr.sizes
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

```shell script
rsync -avP \
    ~/data/plastid/evaluation/ \
    wangq@202.119.37.251:data/plastid/evaluation

# rsync -avP wangq@202.119.37.251:data/plastid/evaluation/ ~/data/plastid/evaluation

```

```shell script
cd ~/data/plastid/evaluation

SRRS=(
    'SRR616966::41'  # Col-0
#    'SRR611086::78'
#    'SRR5216995::121'
    'SRR616965::42'  # Ler-0
#    'SRR611087::79'
#    'SRR545231::46'  # Nipponbare
#    'SRR063638::15'
#    'SRR1542423::23' # A17
#    'SRR1542422::43'
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

```shell script
cd ~/data/plastid/evaluation

SRRS=(
    'SRR616966::5' # Col-0
#    'SRR611086::col_0'
#    'SRR5216995::col_0'
    'SRR616965::col_0' # Ler-0
#    'SRR611087::col_0'
#    'SRR545231::nip'   # Nipponbare
#    'SRR063638::nip'
#    'SRR1542423::a17'  # A17
#    'SRR1542422::a17'
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

```shell script
cd ~/data/plastid/evaluation

SRRS=(
    'SRR616966::5' # Col-0
#    'SRR611086::col_0'
#    'SRR5216995::col_0'
    'SRR616965::5' # Ler-0
#    'SRR611087::col_0'
#    'SRR545231::nip'   # Nipponbare
#    'SRR063638::nip'
#    'SRR1542423::a17'  # A17
#    'SRR1542422::a17'
)
FOLDS=(0 0.25 0.5 1 2 4 8 16 32 64)

for item in "${SRRS[@]}"; do
    SRR="${item%%::*}"
    CHR_NUM="${item##*::}"

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
    echo "Table: ${SRR} Folds"
    echo
    for PART in Nc Mt Pt; do
        cat SRR616966_folds.tsv |
            tsv-filter -H --str-eq chrom:${PART} |
            mlr --itsv --omd cat
        echo
        echo
    done

done

```

Table: SRR616966 Folds

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


Table: SRR616965 Folds

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


```shell script
cd ~/data/plastid/evaluation

SRRS=(
    'SRR616966::5' # Col-0
#    'SRR611086::col_0'
#    'SRR5216995::col_0'
    'SRR616965::5' # Ler-0
#    'SRR611087::col_0'
#    'SRR545231::nip'   # Nipponbare
#    'SRR063638::nip'
#    'SRR1542423::a17'  # A17
#    'SRR1542422::a17'
)
FOLDS=(0 0.25 0.5 1 2 4 8 16 32 64)

for item in "${SRRS[@]}"; do
    SRR="${item%%::*}"
    CHR_NUM="${item##*::}"

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
    echo "Table: ${SRR} Reads"
    echo
    cat ${SRR}_reads.tsv |
        mlr --itsv --omd cat
    echo

done

```

Table: SRR616966 Reads

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


Table: SRR616965 Reads

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


## Remove intermediate files

```shell script
cd ~/data/plastid/evaluation

find . -type d -name "trim" | xargs rm -fr
find . -type f -path "*3_bwa/genome.fa*" | xargs rm
find . -type f -name "*.ba[mi]" | xargs rm
find . -type f -name "*.per-base.bed.gz" | xargs rm

find . -type f -name "*.tadpole.contig.*" | xargs rm

find . -type f -name "core.*" | xargs rm
find . -type f -name "output.*" | xargs rm

```

