# VCFs on folds


[TOC levels=1-3]: # ""

- [VCFs on folds](#vcfs-on-folds)
  - [Genomes](#genomes)
  - [Symlink](#symlink)
  - [Trim, bwa and gatk](#trim-bwa-and-gatk)
  - [Remove intermediate files](#remove-intermediate-files)
  - [VCF processing](#vcf-processing)

## Genomes

```shell script
mkdir -p ~/data/plastid/gatk/genome
cd ~/data/plastid/gatk/genome

cp ../../Atha_1001/genome/genome.fa col_0.fa
cp ../../Osat_50/genome/genome.fa nip.fa
cp ../../Mtru_384/genome/genome.fa a17.fa

```

## Symlink

```shell script
mkdir -p ~/data/plastid/gatk/
cd ~/data/plastid/gatk/

SRRS=(
    'SRR616966::col_0' # Col-0
    'SRR611086::col_0'
    'SRR5216995::col_0'
    'SRR616965::col_0' # Ler-0
    'SRR611087::col_0'
    'SRR545231::nip'   # Nipponbare
    'SRR063638::nip'
    'SRR1542423::a17'  # A17
)
FOLDS=(0 0.25 0.5 1 2 4 8 16 32 64)

for item in "${SRRS[@]}"; do
    SRR="${item%%::*}"
    STRAIN="${item##*::}"

    for FOLD in "${FOLDS[@]}"; do
        BASE_NAME=${SRR}_${FOLD}

        mkdir -p ${BASE_NAME}/1_genome
        pushd ${BASE_NAME}/1_genome

        cp ../../genome/${STRAIN}.fa genome.fa
        popd

        mkdir -p ${BASE_NAME}/2_illumina
        pushd ${BASE_NAME}/2_illumina

        ln -fs ../../../ena/${SRR}_1.fastq.gz R1.fq.gz
        ln -fs ../../../ena/${SRR}_2.fastq.gz R2.fq.gz
        popd

    done
done

```

## Trim, bwa and gatk

* Rsync to hpcc

```shell script
rsync -avP \
    ~/data/plastid/gatk/ \
    wangq@202.119.37.251:data/plastid/gatk

# rsync -avP wangq@202.119.37.251:data/plastid/gatk/ ~/data/plastid/gatk

```

```shell script
cd ~/data/plastid/gatk/

SRRS=(
    'SRR616966::41'  # Col-0
    'SRR611086::78'
    'SRR5216995::121'
    'SRR616965::42'  # Ler-0
    'SRR611087::79'
    'SRR545231::46'  # Nipponbare
    'SRR063638::15'
    'SRR1542423::23' # A17
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

        if [ ! -f 3_gatk/R.filtered.vcf ]; then
            rm *.sh
            anchr template \
                --genome 1000000 \
                --parallel 24 \
                --xmx 80g \
                \
                --trim "--dedupe --cutoff ${CUTOFF} --cutk 31" \
                --qual "25" \
                --len "60" \
                --filter "adapter artifact" \
                \
                --bwa Q25L60 \
                --gatk

            bsub -q mpi -n 24 -J "${BASE_NAME}" "
                bash 2_trim.sh
                bash 3_bwa.sh
                bash 3_gatk.sh
            "
        fi

        popd

    done
done

```

## Remove intermediate files

```shell script
cd ~/data/plastid/gatk/

find . -type d -name "trim" | xargs rm -fr
find . -type f -path "*3_bwa/genome.fa*" | xargs rm
find . -type f -name "*.ba[mi]" | xargs rm
find . -type f -name "*.per-base.bed.gz" | xargs rm

find . -type f -name "core.*" | xargs rm
find . -type f -name "output.*" | xargs rm

```

## VCF processing


```shell script
cd ~/data/plastid/gatk/

SRRS=(
    'SRR616966::col_0' # Col-0
    'SRR611086::col_0'
    'SRR5216995::col_0'
    'SRR616965::col_0' # Ler-0
    'SRR611087::col_0'
    'SRR545231::nip'   # Nipponbare
    'SRR063638::nip'
    'SRR1542423::a17'  # A17
)
FOLDS=(0 0.25 0.5 1 2 4 8 16 32 64)

for item in "${SRRS[@]}"; do
    SRR="${item%%::*}"
    STRAIN="${item##*::}"

    mkdir -p ${SRR}

    for FOLD in "${FOLDS[@]}"; do
        BASE_NAME=${SRR}_${FOLD}

        echo 1>&2 "==> ${item} ${FOLD}"

        bcftools reheader ${BASE_NAME}/3_gatk/R.filtered.vcf \
            --samples <(echo ${BASE_NAME}) |
            bcftools view \
                --apply-filters PASS --types snps --max-alleles 2 --targets Pt -Oz |
            bcftools view --include 'AF>0.01' -Oz -o ${SRR}/${BASE_NAME}.vcf.gz

        bcftools index -f ${SRR}/${BASE_NAME}.vcf.gz

    done

    bcftools merge --merge all $(
            for FOLD in "${FOLDS[@]}"; do echo "${SRR}/${SRR}_${FOLD}.vcf.gz"; done
        ) \
        > ${SRR}.vcf

done

```
